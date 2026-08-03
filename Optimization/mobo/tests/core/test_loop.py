"""The asynchronous loop, driven over a synthetic problem by real subprocesses.

No physics here, but everything else is real: workdirs on disk, payloads in
their own processes, sentinel files, the store, resume. If this passes, what is
left to break in production is the physics payload itself.
"""

from __future__ import annotations

import json

import pytest

from mobo.core.loop import AsyncLoop, LoopConfig
from mobo.core.optimizer import MOBOptimizer, OptimizerConfig
from mobo.core.store import TrialStore
from mobo.core.types import TrialStatus
from mobo.exec.local import LocalExecutor

FAST = dict(num_restarts=4, raw_samples=64, mc_samples=32)


def make_loop(tmp_path, space2d, objectives, evaluator, **loop_kwargs):
    store = TrialStore(tmp_path / "trials.db")
    opt = MOBOptimizer(
        space2d, objectives, OptimizerConfig(n_init=4, batch_size=1, seed=0, **FAST)
    )
    executor = LocalExecutor(max_parallel=loop_kwargs.get("max_in_flight", 2))
    cfg = LoopConfig(poll_interval=0.05, **loop_kwargs)
    return AsyncLoop(opt, evaluator, executor, store, cfg, global_seed=0), store


def test_loop_runs_to_completion(tmp_path, space2d, objectives, synthetic_evaluator):
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=6,
        max_in_flight=2,
        ref_point_from_baseline={"f1": 2.0, "f2": 2.0},
    )
    loop.run()

    trials = store.all()
    assert len(trials) == 6
    assert all(t.status is TrialStatus.COMPLETED for t in trials)
    assert all(t.metrics and "f1" in t.metrics for t in trials)
    # Ids are dense and unique, and every trial got its own workdir.
    assert [t.trial_id for t in trials] == list(range(6))
    assert len({t.workdir for t in trials}) == 6


def test_trial_zero_is_the_baseline_and_anchors_the_reference_point(
    tmp_path, space2d, objectives, synthetic_evaluator
):
    """The rule from the experiment config: ref = factor x baseline."""
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=3,
        max_in_flight=1,
        ref_point_from_baseline={"f1": 2.0, "f2": 2.0},
    )
    loop.run()

    baseline = store.get(0)
    assert baseline.tag == "baseline"
    assert baseline.params == {"x0": 0.5, "x1": 0.5}

    ref = store.get_meta("ref_point")
    assert ref["f1"] == pytest.approx(baseline.metrics["f1"] * 2.0)
    assert ref["f2"] == pytest.approx(baseline.metrics["f2"] * 2.0)
    # ... and the optimizer is using it.
    assert loop.optimizer.objectives[0].ref_point == pytest.approx(ref["f1"])


def test_seeds_are_deterministic_and_distinct(
    tmp_path, space2d, objectives, synthetic_evaluator
):
    loop, store = make_loop(
        tmp_path, space2d, objectives, synthetic_evaluator, max_trials=5, max_in_flight=2
    )
    loop.run()
    seeds = [t.seed for t in store.all()]
    assert len(set(seeds)) == 5

    from mobo.core.types import derive_seed

    assert seeds == [derive_seed(0, i) for i in range(5)]


def test_infeasible_trials_cost_nothing_and_are_recorded(
    tmp_path, space2d, objectives, synthetic_evaluator
):
    synthetic_evaluator.infeasible_below = 0.5  # rejects roughly half the cube
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=8,
        max_in_flight=2,
        baseline_first=False,
    )
    loop.run()

    infeasible = store.by_status(TrialStatus.INFEASIBLE)
    assert infeasible, "the test needs at least one rejected point to prove anything"
    for t in infeasible:
        assert t.error and "feasible region" in t.error
        assert not (synthetic_evaluator.workdir(t) / "trial.json").exists()
    assert len(store.all()) == 8


def test_failed_payloads_are_retried_then_recorded(
    tmp_path, space2d, objectives, synthetic_evaluator
):
    """A payload that writes FAILED is a real failure: no retry, but classified.

    (Retries are for infrastructure — a lost or timed-out job — not for a
    payload that ran and said no.)
    """
    synthetic_evaluator.failure_rate = 1.0
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=3,
        max_in_flight=1,
        baseline_first=False,
    )
    loop.run()

    trials = store.all()
    assert all(t.status is TrialStatus.FAILED for t in trials)
    assert all("injected failure" in (t.error or "") for t in trials)
    assert all(t.attempt == 0 for t in trials)
    assert loop.optimizer.n_observations == 0


def test_lost_jobs_are_retried_with_the_same_seed(
    tmp_path, space2d, objectives, synthetic_evaluator, monkeypatch
):
    """A job that vanishes without a sentinel is infrastructure: retry it.

    Same seed on purpose — if the retry fails identically, the cause is the
    design and not the cluster.
    """
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=1,
        max_in_flight=1,
        max_retries=1,
        baseline_first=False,
    )

    calls = {"n": 0}
    real_prepare = synthetic_evaluator.prepare

    def flaky_prepare(trial):
        cmd, workdir = real_prepare(trial)
        calls["n"] += 1
        if calls["n"] == 1:
            # First attempt: a payload that dies without leaving a sentinel.
            return ["python3", "-c", "import sys; sys.exit(3)"], workdir
        return cmd, workdir

    monkeypatch.setattr(synthetic_evaluator, "prepare", flaky_prepare)
    loop.run()

    trial = store.get(0)
    assert calls["n"] == 2, "the lost job was not retried"
    assert trial.attempt == 1
    assert trial.status is TrialStatus.COMPLETED
    assert trial.seed == store.all()[0].seed


def test_resume_continues_without_losing_or_duplicating_trials(
    tmp_path, space2d, objectives, synthetic_evaluator
):
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=4,
        max_in_flight=1,
        ref_point_from_baseline={"f1": 2.0, "f2": 2.0},
    )
    loop.run()
    first_pass = {t.trial_id: t.params for t in store.all()}
    ref_before = store.get_meta("ref_point")
    store.close()

    # A fresh driver over the same store: same experiment, more budget.
    loop2, store2 = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=7,
        max_in_flight=1,
        ref_point_from_baseline={"f1": 2.0, "f2": 2.0},
    )
    loop2.run()

    trials = store2.all()
    assert [t.trial_id for t in trials] == list(range(7))
    for tid, params in first_pass.items():
        assert store2.get(tid).params == params, "a resumed run rewrote history"
    # The reference point is re-read, never recomputed: it defines the
    # hypervolume, and moving it mid-experiment would rewrite the metric.
    assert store2.get_meta("ref_point") == ref_before
    assert loop2.optimizer.n_observations == 7
    # No second baseline, and the Sobol sequence was continued rather than replayed.
    assert [t.tag for t in trials].count("baseline") == 1
    sobol_points = [tuple(t.unit_x) for t in trials if t.tag == "sobol"]
    assert len(set(sobol_points)) == len(sobol_points)


def test_max_in_flight_is_respected(tmp_path, space2d, objectives):
    """The whole point of the async loop: keep N running, not 1, not N+1."""
    from conftest import SyntheticEvaluator

    evaluator = SyntheticEvaluator(tmp_path, sleep=0.4)
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        evaluator,
        max_trials=6,
        max_in_flight=3,
        baseline_first=False,
    )

    peak = {"n": 0}
    original_fill = loop._fill_slots

    def watched_fill():
        worked = original_fill()
        peak["n"] = max(peak["n"], len(loop._handles))
        return worked

    loop._fill_slots = watched_fill
    loop.run()

    assert peak["n"] == 3, f"peak concurrency was {peak['n']}, expected 3"
    assert len(store.all()) == 6


def test_pending_trials_condition_the_ask(tmp_path, space2d, objectives):
    """With jobs in flight, the optimizer must be told about them."""
    from conftest import SyntheticEvaluator

    evaluator = SyntheticEvaluator(tmp_path, sleep=0.3)
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        evaluator,
        max_trials=8,
        max_in_flight=3,
        baseline_first=False,
    )

    seen_pending = []
    real_ask = loop.optimizer.ask

    def watched_ask(n, pending=()):
        seen_pending.append(len(pending))
        return real_ask(n, pending=pending)

    loop.optimizer.ask = watched_ask
    loop.run()

    assert max(seen_pending) >= 2, "X_pending was never populated; asks are blind"
    assert len(store.all()) == 8


def test_the_payload_is_a_real_subprocess(
    tmp_path, space2d, objectives, synthetic_evaluator
):
    """Guards against a future refactor that quietly evaluates in-process."""
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=1,
        max_in_flight=1,
        baseline_first=False,
    )
    loop.run()

    workdir = tmp_path / "trials" / "t00000"
    assert (workdir / "DONE").is_file()
    assert (workdir / "payload.log").is_file()
    metrics = json.loads((workdir / "metrics.json").read_text())
    assert set(metrics) == {"f1", "f2", "wall_time_s"}
    assert store.get(0).metrics["f1"] == pytest.approx(metrics["f1"])


def test_a_late_sentinel_still_counts_as_a_completed_trial(
    tmp_path, space2d, objectives, synthetic_evaluator, monkeypatch
):
    """The executor can call a job lost while its metrics are already on disk.

    AFS shows a file written on another machine with a delay, so "gone from the
    queue, nothing here yet" and "finished a moment ago" look identical for a
    while. Re-running such a trial pays twice for the same simulation.
    """
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=1,
        max_in_flight=1,
        baseline_first=False,
    )

    from mobo.exec.base import ExitInfo

    real_poll = loop.executor.poll

    def lying_poll():
        return [(h, ExitInfo(False, "lost", "vanished")) for h, _ in real_poll()]

    monkeypatch.setattr(loop.executor, "poll", lying_poll)
    loop.run()

    trial = store.get(0)
    assert trial.status is TrialStatus.COMPLETED, trial.error
    assert trial.metrics and "f1" in trial.metrics
    assert trial.attempt == 0, "a finished trial must not be retried"


def test_the_reference_point_is_anchored_even_without_a_baseline(
    tmp_path, space2d, objectives, synthetic_evaluator
):
    """Otherwise the optimizer never leaves the Sobol phase."""
    loop, store = make_loop(
        tmp_path,
        space2d,
        objectives,
        synthetic_evaluator,
        max_trials=2,
        max_in_flight=1,
        baseline_first=False,
        ref_point_from_baseline={"f1": 2.0, "f2": 2.0},
    )
    loop.run()

    ref = store.get_meta("ref_point")
    assert ref is not None
    first = store.get(0)
    assert ref["f1"] == pytest.approx(first.metrics["f1"] * 2.0)
    # Fixed once: the second completion must not move it.
    assert loop.optimizer.objectives[0].ref_point == pytest.approx(ref["f1"])
