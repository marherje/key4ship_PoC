"""The steady-state asynchronous loop: keep `max_in_flight` trials running.

Not a batch loop. A batch loop (propose q, wait for all q, repeat) wastes the
cluster twice over: it idles while the slowest job of a batch finishes, and it
throws away the information from the fast ones. Here every completion frees a
slot and immediately triggers an `ask` conditioned on whatever is still running.

The store is the state. Nothing here is remembered only in memory, so killing
the driver at any point and running `mobo-resume` continues the same experiment
— including re-adopting jobs that are still running on the cluster.
"""

from __future__ import annotations

import contextlib
import logging
import signal
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from ..exec.base import Executor, Handle
from .evaluator import Evaluator
from .optimizer import MOBOptimizer
from .store import TrialStore
from .types import Result, Trial, TrialStatus, derive_seed

log = logging.getLogger(__name__)

# Store keys. The reference point is written once and re-read forever after:
# recomputing it on resume would silently change what "hypervolume" means
# halfway through an experiment.
META_REF_POINT = "ref_point"


@dataclass
class LoopConfig:
    max_trials: int = 20
    max_in_flight: int = 4
    batch_ask: int = 1
    poll_interval: float = 10.0
    time_budget_hours: float | None = None
    max_retries: int = 1
    baseline_first: bool = True
    # objective name -> multiple of the baseline value that is still acceptable.
    # e.g. {nhits_sipad: 0.5, cost_proxy: 2.0}: half the baseline's hits, twice
    # its cost. Ignored for an objective whose ref_point is given outright.
    ref_point_from_baseline: dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_config(cls, cfg: Any) -> LoopConfig:
        fields = set(cls.__dataclass_fields__)
        kwargs = {k: v for k, v in dict(cfg).items() if k in fields}
        if "ref_point_from_baseline" in kwargs and kwargs["ref_point_from_baseline"]:
            kwargs["ref_point_from_baseline"] = {
                str(k): float(v)
                for k, v in dict(kwargs["ref_point_from_baseline"]).items()
            }
        return cls(**kwargs)


class AsyncLoop:
    def __init__(
        self,
        optimizer: MOBOptimizer,
        evaluator: Evaluator,
        executor: Executor,
        store: TrialStore,
        cfg: LoopConfig | None = None,
        global_seed: int = 0,
        on_progress: Any = None,
    ) -> None:
        self.optimizer = optimizer
        self.evaluator = evaluator
        self.executor = executor
        self.store = store
        self.cfg = cfg or LoopConfig()
        self.global_seed = int(global_seed)
        self.on_progress = on_progress

        self._handles: dict[int, Handle] = {}
        self._stop = False
        self._deadline: float | None = None

    # ── entry point ──────────────────────────────────────────────────────────

    def run(self) -> None:
        """Drive the experiment until a stopping criterion fires."""
        self._install_signal_handlers()
        if self.cfg.time_budget_hours:
            self._deadline = time.time() + self.cfg.time_budget_hours * 3600.0

        self._resume()

        while True:
            worked = self._fill_slots()
            worked |= self._harvest()

            if self._should_stop():
                break
            if not worked:
                time.sleep(self.cfg.poll_interval)

        self._finish()

    # ── resume ───────────────────────────────────────────────────────────────

    def _resume(self) -> None:
        """Rebuild from the store: observations, Sobol position, live jobs."""
        trials = self.store.all()
        if not trials:
            return
        self.optimizer.rehydrate(trials)

        stored_ref = self.store.get_meta(META_REF_POINT)
        if stored_ref:
            self.optimizer.set_ref_point({k: float(v) for k, v in stored_ref.items()})
            log.info("reference point (from the store): %s", stored_ref)

        for trial in self.store.in_flight():
            handle = self.executor.adopt(trial)
            if handle is not None:
                self._handles[trial.trial_id] = handle
                log.info("re-adopted %s (job %s)", trial.name, trial.job_id)
            else:
                log.info("%s cannot be re-adopted; re-queueing it", trial.name)
                self.store.mark(trial, TrialStatus.PROPOSED)

        for trial in self.store.by_status(TrialStatus.PROPOSED):
            self._launch(trial)

    # ── the two halves of a turn ─────────────────────────────────────────────

    def _fill_slots(self) -> bool:
        """Propose and submit until the in-flight budget is spent."""
        free = self._free_slots()
        remaining = self.cfg.max_trials - self._n_started()
        n = min(free, remaining, self.cfg.batch_ask)
        if n <= 0 or self._stop:
            return False

        worked = False
        for trial in self._next_trials(n):
            worked = True
            self._launch(trial)
        return worked

    def _next_trials(self, n: int) -> list[Trial]:
        """The baseline first (once), then whatever the optimizer proposes."""
        trials: list[Trial] = []
        if self._needs_baseline():
            baseline = self.evaluator.baseline_params()
            if baseline is not None:
                trials.append(
                    self._create(self.optimizer.space.to_unit(baseline), "baseline")
                )
                n -= 1
        if n <= 0:
            return trials

        # Everything not yet finished conditions the acquisition, the baseline
        # included — `_create` has already persisted it, so one query covers it.
        # Points *within* the requested batch are handled by the sequential
        # q-batch optimization inside `ask`.
        pending = [t.unit_x for t in self.store.pending()]
        for proposal in self.optimizer.ask(n, pending=pending):
            trials.append(self._create(proposal.unit_x, proposal.tag))
        return trials

    def _harvest(self) -> bool:
        """Collect everything the executor says has finished."""
        finished = self.executor.poll()
        for handle, info in finished:
            self._handles.pop(handle.trial_id, None)
            trial = self.store.get(handle.trial_id)
            result = self._interpret(trial, info)

            if result.status is not TrialStatus.COMPLETED and self._retriable(
                trial, info
            ):
                self._retry(trial, result)
                continue

            self.store.complete(trial, result)
            self._tell(trial, result)
            log.info(
                "%s %s%s",
                trial.name,
                result.status.value,
                f" ({result.error})" if result.error else "",
            )
            self._after_completion(trial)
        return bool(finished)

    def _interpret(self, trial: Trial, info: Any) -> Result:
        """Turn "how it ended" plus "what is on disk" into one verdict.

        Which of the two sources wins depends on who knows more:

        * the payload wrote a sentinel and metrics -> believe the payload,
          even if the executor called the job lost (a sentinel that showed up
          late on AFS is still a finished trial, and re-running it would pay
          for the same simulation twice);
        * the payload said FAILED -> its own diagnosis names the step that
          died, which "failed:" from the executor cannot;
        * anything else (lost, held, timed out) -> the executor's reason, which
          is the only description that exists.
        """
        collected = self.evaluator.collect(trial)
        if info.ok or collected.status is TrialStatus.COMPLETED:
            return collected
        if info.reason == "failed" and collected.error:
            return collected
        return Result.failed(f"{info.reason}: {info.detail or ''}".strip())

    # ── one trial's lifecycle ────────────────────────────────────────────────

    def _create(self, unit_x: Sequence[float], tag: str) -> Trial:
        """Allocate a trial id and its seed, and persist it as PROPOSED."""
        next_id = self.store.next_id()
        params = self.optimizer.space.to_params(unit_x)
        return self.store.create(
            params=params,
            unit_x=list(unit_x),
            seed=derive_seed(self.global_seed, next_id),
            tag=tag,
        )

    def _launch(self, trial: Trial) -> None:
        """Feasibility gate, then prepare and submit."""
        reason = self.evaluator.validate(trial)
        if reason is not None:
            # Rejected without spending a core-hour: the point still counts as
            # explored, and the optimizer is told so it does not come back here.
            result = Result.infeasible(reason)
            self.store.complete(trial, result)
            self._tell(trial, result)
            log.info("%s INFEASIBLE (%s)", trial.name, reason)
            return

        try:
            cmd, workdir = self.evaluator.prepare(trial)
            handle = self.executor.submit(trial, cmd, workdir)
        except Exception as exc:  # noqa: BLE001 - one bad trial must not stop the run
            log.exception("could not submit %s", trial.name)
            result = Result.failed(f"submit_failed: {exc}")
            self.store.complete(trial, result)
            self._tell(trial, result)
            return

        self._handles[trial.trial_id] = handle
        self.store.mark(
            trial,
            TrialStatus.SUBMITTED,
            workdir=str(handle.workdir),
            job_id=handle.job_id,
        )

    def _retriable(self, trial: Trial, info: Any) -> bool:
        """Retry transient infrastructure failures, not bad geometries.

        The retry keeps the *same seed*: if the second attempt fails the same
        way, the failure is a property of the design and not of the cluster.
        """
        return trial.attempt < self.cfg.max_retries and info.reason in (
            "lost",
            "timeout",
            "held",
        )

    def _retry(self, trial: Trial, result: Result) -> None:
        log.warning(
            "%s failed (%s); retry %d/%d with the same seed",
            trial.name,
            result.error,
            trial.attempt + 1,
            self.cfg.max_retries,
        )
        self.store.mark(
            trial, TrialStatus.PROPOSED, attempt=trial.attempt + 1, error=result.error
        )
        self._launch(trial)

    def _tell(self, trial: Trial, result: Result) -> None:
        if result.status is TrialStatus.COMPLETED and result.metrics:
            self.optimizer.observe(trial.unit_x, result.metrics)
        else:
            self.optimizer.observe_failure(trial.unit_x)

    def _after_completion(self, trial: Trial) -> None:
        if self._anchors_the_reference(trial):
            self._anchor_ref_point(trial)
        if self.on_progress is not None:
            try:
                self.on_progress(self)
            except Exception:  # noqa: BLE001 - a failed plot must not stop the run
                log.exception("progress callback failed")

    # ── the reference point ──────────────────────────────────────────────────

    def _anchors_the_reference(self, trial: Trial) -> bool:
        """Normally the baseline; the first completion when there is no baseline.

        With `baseline_first: false` there is no baseline to anchor to, and
        without a reference point the optimizer can never leave the Sobol phase
        — a silent way to turn a Bayesian run into random search. Anchoring on
        the first completed trial instead keeps the loop working; the reference
        is arbitrary either way, and it is still fixed exactly once.
        """
        if trial.tag == "baseline":
            return True
        return not self.cfg.baseline_first and not any(
            t.tag == "baseline" for t in self.store.all()
        )

    def _anchor_ref_point(self, baseline: Trial) -> None:
        """Fix the hypervolume reference from the baseline, once and for all."""
        if self.store.get_meta(META_REF_POINT) or not self.cfg.ref_point_from_baseline:
            return
        if baseline.status is not TrialStatus.COMPLETED or not baseline.metrics:
            log.warning(
                "the baseline trial did not complete; no reference point anchored"
            )
            return

        ref: dict[str, float] = {}
        for obj in self.optimizer.objectives:
            factor = self.cfg.ref_point_from_baseline.get(obj.name)
            if factor is None or obj.name not in baseline.metrics:
                continue
            ref[obj.name] = float(baseline.metrics[obj.name]) * float(factor)
        if not ref:
            return

        self.store.set_meta(META_REF_POINT, ref)
        self.optimizer.set_ref_point(ref)
        log.info("reference point anchored on the baseline: %s", ref)

    # ── stopping ─────────────────────────────────────────────────────────────

    def _free_slots(self) -> int:
        capacity = self.cfg.max_in_flight
        if self.executor.capacity is not None:
            capacity = min(capacity, self.executor.capacity)
        return max(0, capacity - len(self._handles))

    def _n_started(self) -> int:
        """Trials that exist, i.e. what counts against `max_trials`."""
        return len(self.store)

    def _needs_baseline(self) -> bool:
        return self.cfg.baseline_first and len(self.store) == 0

    def _should_stop(self) -> bool:
        if self._stop and not self._handles:
            return True
        if self._deadline and time.time() > self._deadline:
            log.info("time budget exhausted")
            self._stop = True
        if self._stop:
            return not self._handles
        return self._n_started() >= self.cfg.max_trials and not self._handles

    def _install_signal_handlers(self) -> None:
        def handler(signum: int, _frame: Any) -> None:
            if self._stop:  # a second one means business
                raise KeyboardInterrupt
            self._stop = True
            log.warning(
                "signal %d: no new trials will be proposed; waiting for %d in flight "
                "(press again to detach immediately, the store keeps them)",
                signum,
                len(self._handles),
            )

        for sig in (signal.SIGINT, signal.SIGTERM):
            # ValueError: not the main thread, which is how the tests drive the
            # loop. Nothing to install, nothing to complain about.
            with contextlib.suppress(ValueError):
                signal.signal(sig, handler)

    def _finish(self) -> None:
        counts = self.store.counts()
        log.info("loop finished: %s", counts or "no trials")
        if self._handles:
            log.info(
                "%d trial(s) still in flight; resume with:  mobo-resume %s",
                len(self._handles),
                self.store.path.parent,
            )
        self.executor.shutdown(cancel_running=False)
