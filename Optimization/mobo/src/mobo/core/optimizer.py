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
        self._n_proposed = 0
        self._n_sobol = 0
        self._n_failed = 0
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
        self._model_stale = True
        return True

    def observe_failure(self, unit_x: Sequence[float]) -> None:
        """A trial that ran and produced nothing. Counted, not modelled."""
        del unit_x
        self._n_failed += 1

    def rehydrate(self, trials: Sequence[Trial]) -> None:
        """Rebuild the optimizer's state from the store (resume).

        Order matters: the Sobol counter has to end up at the number of Sobol
        points already handed out so the sequence continues rather than repeats.
        """
        self._x.clear()
        self._y.clear()
        self._yvar.clear()
        self._n_proposed = 0
        self._n_sobol = 0
        self._n_failed = 0
        for t in trials:
            self._n_proposed += 1
            if t.tag == "sobol":
                self._n_sobol += 1
            if t.status is TrialStatus.COMPLETED and t.metrics:
                self.observe(t.unit_x, t.metrics)
            elif t.status.is_terminal:
                self.observe_failure(t.unit_x)
        self._model_stale = True
        log.info(
            "rehydrated: %d trials, %d observations, %d sobol drawn",
            self._n_proposed,
            self.n_observations,
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

    def _fit(self) -> Any:
        if self._model is None or self._model_stale:
            torch.manual_seed(derive_seed(self.cfg.seed, self._n_proposed))
            yvar = (
                torch.tensor(self._yvar, dtype=DTYPE)
                if self.cfg.fixed_noise and len(self._yvar) == len(self._x)
                else None
            )
            self._model = build_model(self._x_tensor(), self._y_tensor(), yvar)
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
