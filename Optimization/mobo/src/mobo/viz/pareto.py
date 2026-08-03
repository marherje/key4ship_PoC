"""Pareto-front views of a finished (or running) experiment.

Plots are in *physical* units and physical directions — a cost axis goes up as
cost goes up — even though the optimizer works with everything maximized. The
sign flip lives in `ObjectiveSpec.signed` and stops there; a plot that silently
showed negated cost would be a very expensive kind of confusing.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import plotly.graph_objects as go

from ..core.pareto import objective_matrix, pareto_mask
from ..core.types import ObjectiveSpec, Trial, TrialStatus

STATUS_COLOURS = {
    TrialStatus.COMPLETED: "#4C78A8",
    TrialStatus.FAILED: "#E45756",
    TrialStatus.INFEASIBLE: "#B279A2",
    TrialStatus.RUNNING: "#F58518",
    TrialStatus.SUBMITTED: "#F58518",
    TrialStatus.PROPOSED: "#BAB0AC",
}


def _hover(trial: Trial, objectives: Sequence[ObjectiveSpec]) -> str:
    lines = [f"<b>{trial.name}</b> ({trial.tag or 'n/a'})"]
    for obj in objectives:
        value = (trial.metrics or {}).get(obj.name)
        if value is not None:
            lines.append(f"{obj.name}: {value:.6g}")
    for key, value in sorted(trial.params.items()):
        lines.append(
            f"{key}: {value:.4g}"
            if isinstance(value, (int, float))
            else f"{key}: {value}"
        )
    return "<br>".join(lines)


def pareto_figure(
    trials: Sequence[Trial], objectives: Sequence[ObjectiveSpec]
) -> go.Figure:
    """Objective 0 against objective 1, front highlighted, reference point shown.

    With more than two objectives this shows the first two; the scatter matrix
    below is the honest view in that case.
    """
    if len(objectives) < 2:
        return _empty("a Pareto plot needs at least two objectives")

    y, kept = objective_matrix(trials, objectives)
    if not kept:
        return _empty("no completed trials yet")
    on_front = pareto_mask(y).tolist()

    xo, yo = objectives[0], objectives[1]
    xs = [t.metrics[xo.name] for t in kept]
    ys = [t.metrics[yo.name] for t in kept]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[x for x, f in zip(xs, on_front, strict=False) if not f],
            y=[v for v, f in zip(ys, on_front, strict=False) if not f],
            mode="markers",
            name="dominated",
            marker=dict(size=8, color="#BAB0AC", line=dict(width=0.5, color="#666")),
            text=[
                _hover(t, objectives)
                for t, f in zip(kept, on_front, strict=False)
                if not f
            ],
            hoverinfo="text",
        )
    )
    front = [(x, v, t) for x, v, t, f in zip(xs, ys, kept, on_front, strict=False) if f]
    front.sort()
    fig.add_trace(
        go.Scatter(
            x=[p[0] for p in front],
            y=[p[1] for p in front],
            mode="markers+lines",
            name="Pareto front",
            marker=dict(size=11, color="#4C78A8", symbol="diamond"),
            line=dict(color="#4C78A8", width=1, dash="dot"),
            text=[_hover(p[2], objectives) for p in front],
            hoverinfo="text",
        )
    )

    baseline = [t for t in kept if t.tag == "baseline"]
    if baseline:
        fig.add_trace(
            go.Scatter(
                x=[baseline[0].metrics[xo.name]],
                y=[baseline[0].metrics[yo.name]],
                mode="markers",
                name="baseline",
                marker=dict(size=14, color="#54A24B", symbol="star"),
                text=[_hover(baseline[0], objectives)],
                hoverinfo="text",
            )
        )
    if xo.ref_point is not None and yo.ref_point is not None:
        fig.add_trace(
            go.Scatter(
                x=[xo.ref_point],
                y=[yo.ref_point],
                mode="markers",
                name="reference point",
                marker=dict(size=12, color="#E45756", symbol="x"),
            )
        )

    fig.update_layout(
        title=f"Pareto front: {xo.name} ({xo.direction}) vs {yo.name} ({yo.direction})",
        xaxis_title=xo.name,
        yaxis_title=yo.name,
        template="plotly_white",
        height=520,
    )
    return fig


def scatter_matrix(
    trials: Sequence[Trial], objectives: Sequence[ObjectiveSpec]
) -> go.Figure:
    """All objective pairs at once — the readable view beyond two objectives."""
    _y, kept = objective_matrix(trials, objectives)
    if not kept or len(objectives) < 3:
        return _empty("")
    dims = [
        dict(label=o.name, values=[t.metrics[o.name] for t in kept]) for o in objectives
    ]
    fig = go.Figure(go.Splom(dimensions=dims, marker=dict(size=6, color="#4C78A8")))
    fig.update_layout(title="Objective pairs", template="plotly_white", height=650)
    return fig


def parallel_coordinates(
    trials: Sequence[Trial], objectives: Sequence[ObjectiveSpec]
) -> go.Figure:
    """Parameters to objectives, one line per trial, coloured by objective 0.

    The quickest way to see which parameters the front actually cares about.
    """
    _y, kept = objective_matrix(trials, objectives)
    if not kept:
        return _empty("no completed trials yet")

    param_names = sorted({k for t in kept for k in t.params if _numeric(t.params.get(k))})
    dims: list[dict[str, Any]] = [
        dict(label=name, values=[float(t.params[name]) for t in kept])
        for name in param_names
    ]
    dims += [
        dict(label=o.name, values=[float(t.metrics[o.name]) for t in kept])
        for o in objectives
    ]
    colour = [float(t.metrics[objectives[0].name]) for t in kept]

    fig = go.Figure(
        go.Parcoords(
            line=dict(
                color=colour,
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title=objectives[0].name),
            ),
            dimensions=dims,
        )
    )
    fig.update_layout(
        title="Parameters -> objectives", template="plotly_white", height=520
    )
    return fig


def _numeric(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _empty(message: str) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(
        template="plotly_white",
        height=220,
        annotations=[
            dict(text=message, showarrow=False, xref="paper", yref="paper", x=0.5, y=0.5)
        ],
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    return fig
