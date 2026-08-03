"""Does the optimizer actually optimize?

The central test of the core is `test_beats_sobol`: with the same evaluation
budget, qLogNEHVI must dominate more hypervolume than pure quasi-random search.
Everything else in the package (executors, geometry, reports) is plumbing around
that claim, so it is worth the minute of CPU it costs.
"""

from __future__ import annotations

import math

import pytest
import torch

from mobo.core.optimizer import MOBOptimizer, OptimizerConfig
from mobo.core.pareto import hypervolume, hypervolume_trace, pareto_mask
from mobo.core.types import ObjectiveSpec

# Small but not degenerate: enough restarts to find the acquisition optimum on a
# 2D problem, few enough that the whole file runs in about a minute.
FAST = dict(num_restarts=4, raw_samples=128, mc_samples=64)


def run_campaign(braninc, objectives, space2d, n_init, n_total, seed, noise=0.0):
    """Ask/observe until `n_total` trials; returns the signed objective matrix."""
    evaluate = braninc(noise=noise, seed=seed)
    opt = MOBOptimizer(
        space2d,
        objectives,
        OptimizerConfig(n_init=n_init, batch_size=1, seed=seed, **FAST),
    )
    ys = []
    while opt.n_proposed < n_total:
        for proposal in opt.ask(1):
            metrics = evaluate(proposal.unit_x)
            if metrics is None:
                opt.observe_failure(proposal.unit_x)
                continue
            opt.observe(proposal.unit_x, metrics)
            ys.append([o.signed(metrics[o.name]) for o in objectives])
    return torch.tensor(ys, dtype=torch.double), opt


# ── the initial design ───────────────────────────────────────────────────────


def test_initial_design_is_sobol_then_switches(space2d, objectives, braninc):
    opt = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=4, seed=0, **FAST))
    evaluate = braninc(seed=0)

    tags = []
    for _ in range(6):
        (proposal,) = opt.ask(1)
        tags.append(proposal.tag)
        opt.observe(proposal.unit_x, evaluate(proposal.unit_x))

    assert tags[:4] == ["sobol"] * 4
    assert tags[4:] == ["qlognehvi"] * 2


def test_default_n_init_scales_with_dimension(space2d, objectives):
    opt = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=0))
    assert opt.n_init == 2 * (space2d.dim + 1)


def test_sobol_sequence_is_not_restarted(space2d, objectives):
    """Successive asks must continue the sequence, never replay its head."""
    opt = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=8, seed=3))
    first = [p.unit_x for p in opt.ask(4)]
    second = [p.unit_x for p in opt.ask(4)]
    assert len({tuple(p) for p in first + second}) == 8

    # ... and drawing them one at a time gives exactly the same eight points.
    one_by_one = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=8, seed=3))
    singles = [one_by_one.ask(1)[0].unit_x for _ in range(8)]
    assert singles == first + second


def test_without_a_reference_point_it_stays_on_sobol(space2d, braninc):
    opt = MOBOptimizer(
        space2d,
        [ObjectiveSpec("f1", "max"), ObjectiveSpec("f2", "max", -6.0)],
        OptimizerConfig(n_init=2, seed=0, **FAST),
    )
    evaluate = braninc(seed=0)
    tags = []
    for _ in range(5):
        (p,) = opt.ask(1)
        tags.append(p.tag)
        opt.observe(p.unit_x, evaluate(p.unit_x))
    assert tags == ["sobol"] * 5

    # Filling it in (as the loop does from the baseline trial) unblocks the GP.
    opt.set_ref_point({"f1": -18.0})
    assert opt.ask(1)[0].tag == "qlognehvi"


# ── asynchrony ───────────────────────────────────────────────────────────────


@pytest.mark.slow
def test_pending_points_are_not_reproposed(space2d, objectives, braninc):
    """X_pending is what makes the asynchronous loop work.

    Without it, every slot freed while the same observations are on the table
    gets the same proposal, and `max_in_flight` copies of one geometry get
    simulated.
    """
    opt = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=6, seed=11, **FAST))
    evaluate = braninc(seed=11)
    for _ in range(6):
        (p,) = opt.ask(1)
        opt.observe(p.unit_x, evaluate(p.unit_x))

    first = opt.ask(1)[0].unit_x
    blind = opt.ask(1)[0].unit_x  # same data, nothing declared pending
    aware = opt.ask(1, pending=[first])[0].unit_x

    def distance(a, b):
        return math.dist(a, b)

    assert distance(first, blind) < 1e-3, "sanity: without X_pending it repeats itself"
    assert distance(first, aware) > 1e-2, "X_pending did not move the proposal"


@pytest.mark.slow
def test_batch_ask_returns_distinct_points(space2d, objectives, braninc):
    opt = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=6, seed=5, **FAST))
    evaluate = braninc(seed=5)
    for _ in range(6):
        (p,) = opt.ask(1)
        opt.observe(p.unit_x, evaluate(p.unit_x))

    batch = [p.unit_x for p in opt.ask(3)]
    assert len(batch) == 3
    for i, a in enumerate(batch):
        for b in batch[i + 1 :]:
            assert math.dist(a, b) > 1e-3


# ── failures ─────────────────────────────────────────────────────────────────


@pytest.mark.slow
def test_failed_trials_are_excluded_from_the_gp(space2d, objectives, braninc):
    evaluate = braninc(failure_rate=0.3, seed=17)
    opt = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=6, seed=17, **FAST))

    n_ok = 0
    for _ in range(16):
        (p,) = opt.ask(1)
        metrics = evaluate(p.unit_x)
        if metrics is None:
            opt.observe_failure(p.unit_x)
        else:
            opt.observe(p.unit_x, metrics)
            n_ok += 1

    assert evaluate.n_failures > 0, "the fixture never failed; test proves nothing"
    assert opt.n_observations == n_ok

    # The model is fitted lazily, on ask; one more ask brings it up to date with
    # everything observed above.
    opt.ask(1)
    # ModelListGP: one GP per objective, each trained on the successes only.
    train_x = opt.model.models[0].train_inputs[0]
    assert train_x.shape[0] == n_ok


def test_unusable_metrics_are_treated_as_failures(space2d, objectives):
    opt = MOBOptimizer(space2d, objectives, OptimizerConfig(n_init=4))
    assert opt.observe([0.1, 0.1], {"f1": 1.0}) is False  # f2 missing
    assert opt.observe([0.2, 0.2], {"f1": float("nan"), "f2": 1.0}) is False
    assert opt.observe([0.3, 0.3], {"f1": 1.0, "f2": 2.0}) is True
    assert opt.n_observations == 1


# ── the cross-check that matters ─────────────────────────────────────────────


@pytest.mark.slow
def test_beats_sobol(space2d, objectives, braninc):
    """8 Sobol + 24 qLogNEHVI must dominate more volume than 32 Sobol points.

    Bayesian optimization is stochastic, so this is a majority vote over seeds
    rather than a per-seed assertion: 3 of 4 seeds must win.
    """
    ref = torch.tensor([-18.0, -6.0], dtype=torch.double)
    seeds = [0, 1, 2, 3]
    wins = []
    for seed in seeds:
        y_bo, _ = run_campaign(braninc, objectives, space2d, 8, 32, seed)
        y_sobol, _ = run_campaign(braninc, objectives, space2d, 32, 32, seed)
        hv_bo = hypervolume(y_bo, ref)
        hv_sobol = hypervolume(y_sobol, ref)
        wins.append(hv_bo > hv_sobol)
        print(f"seed {seed}: qLogNEHVI {hv_bo:.2f} vs Sobol {hv_sobol:.2f}")
    assert sum(wins) >= 3, f"qLogNEHVI won only {sum(wins)}/{len(seeds)} seeds"


# ── hypervolume bookkeeping ──────────────────────────────────────────────────


def test_hypervolume_matches_botorch_called_by_hand(objectives):
    from botorch.utils.multi_objective.hypervolume import Hypervolume

    y = torch.tensor([[-10.0, -1.0], [-2.0, -5.0], [-15.0, -5.5]], dtype=torch.double)
    ref = torch.tensor([-18.0, -6.0], dtype=torch.double)
    front = y[pareto_mask(y)]
    assert hypervolume(y, ref) == pytest.approx(Hypervolume(ref_point=ref).compute(front))


def test_hypervolume_ignores_points_worse_than_the_reference():
    y = torch.tensor([[-20.0, -7.0]], dtype=torch.double)
    ref = torch.tensor([-18.0, -6.0], dtype=torch.double)
    assert hypervolume(y, ref) == 0.0


def test_hypervolume_trace_is_non_decreasing():
    torch.manual_seed(0)
    y = torch.rand(20, 2, dtype=torch.double) * 10 - 10
    ref = torch.tensor([-11.0, -11.0], dtype=torch.double)
    trace = hypervolume_trace(y, ref)
    assert all(b >= a - 1e-9 for a, b in zip(trace, trace[1:], strict=False))


def test_pareto_mask_agrees_with_botorch(objectives):
    from botorch.utils.multi_objective.pareto import is_non_dominated

    torch.manual_seed(3)
    y = torch.rand(30, 2, dtype=torch.double)
    assert torch.equal(pareto_mask(y), is_non_dominated(y))
