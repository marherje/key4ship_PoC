"""The ask/tell optimizer: search space in, unit-cube proposals out.

Knows nothing about how a point is evaluated, or by whom — `ask()` returns
coordinates and `observe()` takes metrics back. That is what lets the same
object drive a synthetic test function in a unit test and a night of HTCondor
jobs without changing a line.

Two policies live here and nowhere else:

* the initial design is Sobol until `n_init` points have been *proposed* (not
  completed — in an asynchronous loop waiting for completions would idle the
  cluster), after which qLogNEHVI takes over as soon as there are enough
  observations to fit a GP on;
* a failed trial is dropped, not modelled. Feeding a sentinel value to the GP
  would invent a landscape that does not exist; learning the feasible region
  properly needs a classifier, which is Phase 5 material.

The loop never dies from a numerical problem: if fitting the GP or optimizing
the acquisition raises, the optimizer logs it and falls back to Sobol.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import torch

from .acquisition import propose_qlognehvi, sobol_points
from .models import DTYPE, build_model
from .pareto import hypervolume, pareto_mask, ref_point_tensor
from .search_space import SearchSpace
from .types import ObjectiveSpec, Trial, TrialStatus, derive_seed

log = logging.getLogger(__name__)

# The driver runs on a shared login node. Left to itself torch spreads the
# acquisition optimization over every core it sees, and 2026-09-18 that was 20
# CPU-hours in an evening; one thread keeps it to one core, at nice 19, which
# is the only footprint acceptable on lxplus (CLAUDE.md). MOBO_TORCH_THREADS
# raises it on a machine that is actually ours.
torch.set_num_threads(max(1, int(os.environ.get("MOBO_TORCH_THREADS", "1"))))


@dataclass
class OptimizerConfig:
    """Knobs of the optimizer. Everything here is settable from Hydra."""

    n_init: int = 0  # 0 -> 2 * (d + 1)
    batch_size: int = 1
    seed: int = 0
    num_restarts: int = 10
    raw_samples: int = 512
    mc_samples: int = 128
    sequential: bool = True
    prune_baseline: bool = True
    # Below this many observations a GP is a random number generator with extra
    # steps; keep drawing Sobol points instead.
    min_observations: int = 3
    # Per-observation noise from the metrics' standard errors, when the
    # evaluator reports them (`<metric>_stderr`). Off for the proxy phase: a hit
    # count has no meaningful stderr from a single run.
    fixed_noise: bool = False
    # An INFEASIBLE proposal (the geometry does not build: over the
    # longitudinal budget, say) enters the GP as an observation AT THE
    # REFERENCE POINT -- the worst value still worth having on every
    # objective, i.e. zero hypervolume gain -- so the acquisition learns to
    # stay away. Without this the point was merely counted and the model,
    # which only ever saw the good side of the boundary, proposed the same
    # corner again: 146 rejections in one evening (snd_genie_gen_10k,
    # 2026-09-18). Run failures (crashes) are NOT modelled this way: they are
    # not a property of the geometry.
    penalize_infeasible: bool = True
    # Circuit breaker: after this many infeasible proposals in a row the next
    # one is Sobol, whatever the model says. Each rejected ask still costs a
    # full acquisition optimization, and the penalty above needs a refit to
    # take effect.
    max_consecutive_infeasible: int = 3

    @classmethod
    def from_config(cls, cfg: Any) -> OptimizerConfig:
        fields = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in dict(cfg).items() if k in fields})


@dataclass
class Proposal:
    """One point to evaluate, and where it came from."""

    unit_x: list[float]
    tag: str  # "sobol" | "qlognehvi"


class MOBOptimizer:
    def __init__(
        self,
        space: SearchSpace,
        objectives: Sequence[ObjectiveSpec],
        cfg: OptimizerConfig | None = None,
    ) -> None:
        if not objectives:
            raise ValueError("need at least one objective")
        self.space = space
        self.objectives = list(objectives)
        self.cfg = cfg or OptimizerConfig()

        self._x: list[list[float]] = []
        self._y: list[list[float]] = []
        self._yvar: list[list[float]] = []
        self._x_infeasible: list[list[float]] = []
        self._n_proposed = 0
        self._n_sobol = 0
        self._n_failed = 0
        self._consecutive_infeasible = 0
        self._model: Any = None
        self._model_stale = True

    # ── shape ────────────────────────────────────────────────────────────────

    @property
    def n_init(self) -> int:
        return self.cfg.n_init or 2 * (self.space.dim + 1)

    @property
    def n_observations(self) -> int:
        return len(self._x)

    @property
    def n_proposed(self) -> int:
        return self._n_proposed

    @property
    def model(self) -> Any:
        """The fitted GP, or None if there has not been one yet.

        Exposed so tests can assert on what the model was actually trained on
        (in particular: that failed trials are not in there).
        """
        return self._model

    # ── telling ──────────────────────────────────────────────────────────────

    def observe(self, unit_x: Sequence[float], metrics: dict[str, float]) -> bool:
        """Record one completed evaluation. False if the metrics are unusable."""
        try:
            y = [obj.signed(float(metrics[obj.name])) for obj in self.objectives]
        except (KeyError, TypeError, ValueError) as exc:
            log.warning("dropping an observation with unusable metrics: %s", exc)
            self._n_failed += 1
            return False
        if not all(map(_finite, y)):
            log.warning("dropping an observation with non-finite objectives: %s", y)
            self._n_failed += 1
            return False

        self._x.append(list(map(float, unit_x)))
        self._y.append(y)
        if self.cfg.fixed_noise:
            self._yvar.append(
                [
                    max(float(metrics.get(f"{obj.name}_stderr", 0.0)) ** 2, 1e-9)
                    for obj in self.objectives
                ]
            )
        self._consecutive_infeasible = 0
        self._model_stale = True
        return True

    def observe_failure(self, unit_x: Sequence[float]) -> None:
        """A trial that ran and produced nothing. Counted, not modelled."""
        del unit_x
        self._n_failed += 1

    def observe_infeasible(self, unit_x: Sequence[float]) -> None:
        """A proposal the evaluator refused to build. Modelled as worthless.

        See OptimizerConfig.penalize_infeasible. The point is kept aside and
        materialized at fit time, because the reference point it is pinned to
        may not exist yet (it is anchored by the baseline trial).
        """
        self._n_failed += 1
        self._consecutive_infeasible += 1
        if self.cfg.penalize_infeasible:
            self._x_infeasible.append(list(map(float, unit_x)))
            self._model_stale = True

    @property
    def n_infeasible(self) -> int:
        return len(self._x_infeasible)

    def rehydrate(self, trials: Sequence[Trial]) -> None:
        """Rebuild the optimizer's state from the store (resume).

        Order matters: the Sobol counter has to end up at the number of Sobol
        points already handed out so the sequence continues rather than repeats.
        """
        self._x.clear()
        self._y.clear()
        self._yvar.clear()
        self._x_infeasible.clear()
        self._n_proposed = 0
        self._n_sobol = 0
        self._n_failed = 0
        self._consecutive_infeasible = 0
        for t in trials:
            self._n_proposed += 1
            if t.tag == "sobol":
                self._n_sobol += 1
            if t.status is TrialStatus.COMPLETED and t.metrics:
                self.observe(t.unit_x, t.metrics)
            elif t.status is TrialStatus.INFEASIBLE:
                self.observe_infeasible(t.unit_x)
            elif t.status.is_terminal:
                self.observe_failure(t.unit_x)
        self._model_stale = True
        log.info(
            "rehydrated: %d trials, %d observations, %d infeasible, %d sobol drawn",
            self._n_proposed,
            self.n_observations,
            self.n_infeasible,
            self._n_sobol,
        )

    def set_ref_point(self, values: dict[str, float]) -> None:
        """Fill in the reference point, e.g. from the baseline trial."""
        self.objectives = [
            ObjectiveSpec(o.name, o.direction, values.get(o.name, o.ref_point))
            for o in self.objectives
        ]

    # ── asking ───────────────────────────────────────────────────────────────

    def ask(
        self, n: int | None = None, pending: Sequence[Sequence[float]] = ()
    ) -> list[Proposal]:
        """Propose `n` points, conditioned on the ones already in flight."""
        n = int(n if n is not None else self.cfg.batch_size)
        if n <= 0:
            return []

        reason = self._why_sobol()
        if reason is None:
            try:
                points = self._ask_model(n, pending)
                return self._record(
                    [Proposal(list(map(float, p)), "qlognehvi") for p in points]
                )
            except Exception as exc:  # noqa: BLE001 - the loop must never stop here
                log.exception("acquisition failed (%s); falling back to Sobol", exc)
        else:
            log.debug("Sobol phase: %s", reason)

        if self._consecutive_infeasible >= self.cfg.max_consecutive_infeasible:
            # The breaker tripped: this Sobol point is the reset. If it is
            # infeasible too the count starts again from one.
            self._consecutive_infeasible = 0
        return self._record(self._ask_sobol(n))

    def _why_sobol(self) -> str | None:
        """Why the next batch must be Sobol, or None if the GP can drive."""
        if self._n_proposed < self.n_init:
            return f"{self._n_proposed}/{self.n_init} initial points proposed"
        if self.n_observations < max(self.cfg.min_observations, 2):
            return f"only {self.n_observations} observations"
        if ref_point_tensor(self.objectives) is None:
            missing = [o.name for o in self.objectives if o.ref_point is None]
            return "no reference point for " + ", ".join(missing)
        if self._consecutive_infeasible >= self.cfg.max_consecutive_infeasible:
            log.warning(
                "%d infeasible proposals in a row; one Sobol point to break out",
                self._consecutive_infeasible,
            )
            return f"{self._consecutive_infeasible} consecutive infeasible proposals"
        return None

    def _ask_sobol(self, n: int) -> list[Proposal]:
        points = sobol_points(self.space.dim, n, seed=self.cfg.seed, skip=self._n_sobol)
        self._n_sobol += n
        return [Proposal(list(map(float, row)), "sobol") for row in points]

    def _ask_model(self, n: int, pending: Sequence[Sequence[float]]) -> torch.Tensor:
        model = self._fit()
        ref = ref_point_tensor(self.objectives)
        assert ref is not None  # guaranteed by _why_sobol
        x_pending = (
            torch.tensor([list(map(float, p)) for p in pending], dtype=DTYPE)
            if len(pending)
            else None
        )
        return propose_qlognehvi(
            model=model,
            ref_point=ref.tolist(),
            x_baseline=self._x_tensor(),
            q=n,
            x_pending=x_pending,
            num_restarts=self.cfg.num_restarts,
            raw_samples=self.cfg.raw_samples,
            mc_samples=self.cfg.mc_samples,
            sequential=self.cfg.sequential,
            prune_baseline=self.cfg.prune_baseline,
            seed=derive_seed(self.cfg.seed, self._n_proposed),
        )

    def _record(self, proposals: list[Proposal]) -> list[Proposal]:
        self._n_proposed += len(proposals)
        return proposals

    # ── the model ────────────────────────────────────────────────────────────

    def _x_tensor(self) -> torch.Tensor:
        return torch.tensor(self._x, dtype=DTYPE)

    def _y_tensor(self) -> torch.Tensor:
        return torch.tensor(self._y, dtype=DTYPE)

    def _training_set(self) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | None]:
        """Observations plus, when it can be placed, one row per infeasible point.

        The infeasible rows sit exactly on the (signed) reference point. Their
        noise, when the model has per-point noise at all, is the median of the
        real observations': a hard zero would make the GP bend sharply to pass
        through them, and they are a statement about the region, not a
        measurement.
        """
        x = self._x_tensor()
        y = self._y_tensor()
        yvar = (
            torch.tensor(self._yvar, dtype=DTYPE)
            if self.cfg.fixed_noise and len(self._yvar) == len(self._x)
            else None
        )
        ref = ref_point_tensor(self.objectives)
        if self._x_infeasible and ref is not None:
            k = len(self._x_infeasible)
            x = torch.cat([x, torch.tensor(self._x_infeasible, dtype=DTYPE)])
            y = torch.cat([y, ref.expand(k, -1)])
            if yvar is not None:
                fill = (
                    yvar.median(dim=0).values if len(yvar) else torch.full_like(ref, 1e-9)
                )
                yvar = torch.cat([yvar, fill.expand(k, -1)])
        return x, y, yvar

    def _fit(self) -> Any:
        if self._model is None or self._model_stale:
            torch.manual_seed(derive_seed(self.cfg.seed, self._n_proposed))
            x, y, yvar = self._training_set()
            self._model = build_model(x, y, yvar)
            self._model_stale = False
        return self._model

    # ── reporting ────────────────────────────────────────────────────────────

    def hypervolume(self) -> float:
        ref = ref_point_tensor(self.objectives)
        if ref is None or not self._y:
            return 0.0
        return hypervolume(self._y_tensor(), ref)

    def pareto_indices(self) -> list[int]:
        if not self._y:
            return []
        return [
            i for i, keep in enumerate(pareto_mask(self._y_tensor()).tolist()) if keep
        ]


def _finite(value: float) -> bool:
    return value == value and abs(value) != float("inf")
