"""The report, over a synthetic store.

Plots are hard to assert on, so what is tested is what can go quietly wrong:
the numbers behind them (hypervolume, which points are on the front) and the
promise that generating a report never raises — it runs after every completed
trial, and a failed plot must never take a night of cluster time with it.
"""

from __future__ import annotations

import random

import pytest
import torch

from mobo.core.pareto import hypervolume, objective_matrix, pareto_mask, ref_point_tensor
from mobo.core.store import TrialStore
from mobo.core.types import ObjectiveSpec, Result, TrialStatus
from mobo.viz.pareto import parallel_coordinates, pareto_figure, scatter_matrix
from mobo.viz.progress import (
    hypervolume_figure,
    status_figure,
    timeline_figure,
    wall_time_figure,
)
from mobo.viz.report import build_report, write_report

OBJECTIVES = [
    ObjectiveSpec("nhits_sipad", "max", 200.0),
    ObjectiveSpec("cost_proxy", "min", 3000.0),
]


@pytest.fixture
def populated_store(tmp_path):
    """20 trials with a spread of outcomes, including failures."""
    store = TrialStore(tmp_path / "trials.db")
    rng = random.Random(4)
    for i in range(20):
        trial = store.create(
            params={
                "SiPad_WThickness": rng.uniform(5, 15),
                "SiPad_dim_z": rng.uniform(250, 500),
                "sipad_fill": rng.uniform(0.3, 1.0),
            },
            unit_x=[rng.random() for _ in range(3)],
            seed=i,
            tag="baseline" if i == 0 else ("sobol" if i < 8 else "qlognehvi"),
        )
        store.mark(trial, TrialStatus.SUBMITTED, workdir=str(tmp_path / trial.name))
        if i == 5:
            store.complete(trial, Result.failed("job4_failed: k4run exited with code 1"))
        elif i == 6:
            store.complete(trial, Result.infeasible("SiPad_NLayers = 30 > max = 23"))
        else:
            store.complete(
                trial,
                Result.ok(
                    {
                        "nhits_sipad": rng.uniform(150, 700),
                        "cost_proxy": rng.uniform(400, 3500),
                        "wall_time_s": rng.uniform(60, 600),
                    }
                ),
            )
    store.set_meta("ref_point", {"nhits_sipad": 200.0, "cost_proxy": 3000.0})
    store.set_meta("versions", {"mobo": "0.1.0"})
    yield store
    store.close()


# ── the numbers behind the plots ─────────────────────────────────────────────


def test_the_front_agrees_with_botorch(populated_store):
    """What the plot calls Pareto-optimal must be what BoTorch calls it."""
    from botorch.utils.multi_objective.pareto import is_non_dominated

    y, kept = objective_matrix(populated_store.all(), OBJECTIVES)
    ours = pareto_mask(y)
    assert torch.equal(ours, is_non_dominated(y))

    # And the front really is non-dominated in physical terms: no other point
    # has both more hits and lower cost.
    front = [t for t, keep in zip(kept, ours.tolist(), strict=True) if keep]
    for a in front:
        for b in kept:
            better_hits = b.metrics["nhits_sipad"] > a.metrics["nhits_sipad"]
            better_cost = b.metrics["cost_proxy"] < a.metrics["cost_proxy"]
            assert not (better_hits and better_cost), f"{b.name} dominates {a.name}"


def test_hypervolume_matches_a_hand_call(populated_store):
    from botorch.utils.multi_objective.hypervolume import Hypervolume

    y, _ = objective_matrix(populated_store.all(), OBJECTIVES)
    ref = ref_point_tensor(OBJECTIVES)
    better = (y > ref).all(dim=-1)
    expected = Hypervolume(ref_point=ref).compute(y[better][pareto_mask(y[better])])
    assert hypervolume(y, ref) == pytest.approx(expected)


def test_failed_and_infeasible_trials_are_not_on_the_front(populated_store):
    _y, kept = objective_matrix(populated_store.all(), OBJECTIVES)
    assert all(t.status is TrialStatus.COMPLETED for t in kept)
    assert len(kept) == 18  # 20 minus one failure and one infeasible


# ── the figures ──────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "factory",
    [pareto_figure, parallel_coordinates, scatter_matrix],
    ids=["pareto", "parallel", "matrix"],
)
def test_objective_figures_build(populated_store, factory):
    fig = factory(populated_store.all(), OBJECTIVES)
    assert fig.to_html(full_html=False, include_plotlyjs=False)


@pytest.mark.parametrize(
    "factory",
    [status_figure, timeline_figure, wall_time_figure],
    ids=["status", "timeline", "walltime"],
)
def test_progress_figures_build(populated_store, factory):
    assert factory(populated_store.all()).to_html(full_html=False, include_plotlyjs=False)


def test_hypervolume_figure_builds(populated_store):
    fig = hypervolume_figure(populated_store.all(), OBJECTIVES)
    assert fig.to_html(full_html=False, include_plotlyjs=False)


def test_figures_survive_an_empty_store(tmp_path):
    """A report is written before the first trial finishes, every time."""
    store = TrialStore(tmp_path / "empty.db")
    for factory in (pareto_figure, parallel_coordinates, hypervolume_figure):
        assert factory(store.all(), OBJECTIVES) is not None
    assert status_figure(store.all()) is not None
    store.close()


def test_figures_survive_a_missing_reference_point(populated_store):
    """Before the baseline completes there is no reference point yet."""
    no_ref = [ObjectiveSpec(o.name, o.direction, None) for o in OBJECTIVES]
    assert hypervolume_figure(populated_store.all(), no_ref) is not None
    assert pareto_figure(populated_store.all(), no_ref) is not None


# ── the report ───────────────────────────────────────────────────────────────


def test_report_is_written_and_self_contained(populated_store, tmp_path):
    path = write_report(tmp_path, OBJECTIVES)
    assert path.is_file()
    html = path.read_text()

    # Inlined plotly: the report is written on lxplus and read on a laptop, so
    # nothing may be loaded from the network. (The plotly bundle does mention
    # cdn.plot.ly internally, as the default source of topojson map data — we
    # draw no maps, so nothing ever fetches it. What matters is that there is no
    # remote <script>/<link>.)
    assert "plotly" in html.lower()
    assert "<script src=" not in html.replace('<script src="data:', "")
    assert "<link" not in html
    assert "<title>" in html
    # The cards, the table and the provenance block.
    assert "hypervolume" in html
    assert "t00000" in html
    assert "job4_failed" in html  # failures are visible, not swallowed
    assert "baseline" in html


def test_report_marks_the_front_and_the_baseline(populated_store, tmp_path):
    html = build_report(tmp_path, OBJECTIVES)
    assert 'class="front' in html or 'class="baseline' in html
    assert "pareto" in html


def test_report_never_raises_on_a_broken_run(tmp_path):
    """Called after every trial: it may produce nothing, never an exception."""
    (tmp_path / "trials.db").write_text("this is not a database")
    path = write_report(tmp_path, OBJECTIVES)
    assert path.name == "report.html"
