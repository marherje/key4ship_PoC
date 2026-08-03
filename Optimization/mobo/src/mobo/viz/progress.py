"""How the campaign progressed: hypervolume, throughput, failures.

The hypervolume trace is the one plot that says whether the optimization is
working. It is computed against the reference point stored with the run, in
trial order, so it can only go up — a drop means the reference moved, which
means something rewrote history.
"""

from __future__ import annotations

from collections.abc import Sequence

import plotly.graph_objects as go

from ..core.pareto import hypervolume_trace, objective_matrix, ref_point_tensor
from ..core.types import ObjectiveSpec, Trial, TrialStatus
from .pareto import STATUS_COLOURS, _empty


def hypervolume_figure(
    trials: Sequence[Trial], objectives: Sequence[ObjectiveSpec]
) -> go.Figure:
    ref = ref_point_tensor(objectives)
    if ref is None:
        return _empty("no reference point yet (the baseline trial has not completed)")
    y, kept = objective_matrix(trials, objectives)
    if not kept:
        return _empty("no completed trials yet")

    trace = hypervolume_trace(y, ref)
    fig = go.Figure(
        go.Scatter(
            x=[t.trial_id for t in kept],
            y=trace,
            mode="lines+markers",
            line=dict(color="#4C78A8", width=2, shape="hv"),
            marker=dict(size=7),
            name="dominated hypervolume",
            text=[t.name for t in kept],
        )
    )
    fig.update_layout(
        title="Dominated hypervolume vs trial",
        xaxis_title="trial id",
        yaxis_title="hypervolume",
        template="plotly_white",
        height=420,
    )
    return fig


def status_figure(trials: Sequence[Trial]) -> go.Figure:
    """How the budget was spent: completed, failed, rejected before running."""
    counts: dict[str, int] = {}
    for t in trials:
        counts[t.status.value] = counts.get(t.status.value, 0) + 1
    if not counts:
        return _empty("no trials yet")

    order = [s for s in TrialStatus if s.value in counts]
    fig = go.Figure(
        go.Bar(
            x=[s.value for s in order],
            y=[counts[s.value] for s in order],
            marker_color=[STATUS_COLOURS.get(s, "#888") for s in order],
            text=[counts[s.value] for s in order],
            textposition="auto",
        )
    )
    fig.update_layout(
        title="Trials by status",
        template="plotly_white",
        height=340,
        yaxis_title="trials",
    )
    return fig


def timeline_figure(trials: Sequence[Trial]) -> go.Figure:
    """One bar per trial from submission to completion.

    This is where asynchrony is visible: with `max_in_flight = N` there should
    be N overlapping bars at any time, and a slot should be refilled the moment
    it frees up.
    """
    runs = [t for t in trials if t.submitted_at and t.finished_at]
    if not runs:
        return _empty("no timing information yet")

    t0 = min(t.submitted_at for t in runs)  # type: ignore[type-var]
    fig = go.Figure()
    for t in runs:
        fig.add_trace(
            go.Scatter(
                x=[(t.submitted_at - t0) / 60.0, (t.finished_at - t0) / 60.0],  # type: ignore[operator]
                y=[t.trial_id, t.trial_id],
                mode="lines",
                line=dict(color=STATUS_COLOURS.get(t.status, "#888"), width=6),
                name=t.name,
                showlegend=False,
                hovertext=f"{t.name} — {t.status.value}",
            )
        )
    fig.update_layout(
        title="Trial timeline",
        xaxis_title="minutes since the first submission",
        yaxis_title="trial id",
        template="plotly_white",
        height=420,
    )
    return fig


def wall_time_figure(trials: Sequence[Trial]) -> go.Figure:
    """Distribution of evaluation cost — what a trial actually buys."""
    times = [
        t.metrics["wall_time_s"] / 60.0
        for t in trials
        if t.metrics and t.metrics.get("wall_time_s")
    ]
    if not times:
        return _empty("no wall times recorded yet")
    fig = go.Figure(go.Histogram(x=times, marker_color="#4C78A8", nbinsx=20))
    fig.update_layout(
        title="Payload wall time",
        xaxis_title="minutes",
        yaxis_title="trials",
        template="plotly_white",
        height=340,
    )
    return fig
