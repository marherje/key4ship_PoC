"""One self-contained HTML file per experiment: `runs/<exp>/report.html`.

Self-contained matters more than it sounds: the report is written on lxplus and
read on a laptop, so it embeds plotly rather than linking to a CDN. It is
regenerated after every completed trial, which is why nothing here is allowed to
raise — a run must never die because a plot could not be drawn.
"""

from __future__ import annotations

import datetime as _dt
import html
import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from ..core.pareto import hypervolume, objective_matrix, pareto_mask, ref_point_tensor
from ..core.store import TrialStore
from ..core.types import ObjectiveSpec, Trial, TrialStatus
from .pareto import parallel_coordinates, pareto_figure, scatter_matrix
from .progress import hypervolume_figure, status_figure, timeline_figure, wall_time_figure

log = logging.getLogger(__name__)

CSS = """
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, sans-serif;
       margin: 0 auto; max-width: 1100px; padding: 24px; color: #222; }
h1 { margin-bottom: 4px; } h2 { margin-top: 36px; border-bottom: 1px solid #eee;
       padding-bottom: 6px; }
.sub { color: #666; margin-top: 0; }
table { border-collapse: collapse; width: 100%; font-size: 14px; }
th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid #eee; }
th { background: #fafafa; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
tr.front td { font-weight: 600; }
tr.baseline td { background: #f2f9f0; }
.cards { display: flex; flex-wrap: wrap; gap: 12px; margin: 16px 0; }
.card { border: 1px solid #eee; border-radius: 8px; padding: 12px 16px;
        min-width: 150px; }
.card .v { font-size: 22px; font-weight: 600; }
.card .k { color: #666; font-size: 12px; text-transform: uppercase; }
code, pre { background: #f6f6f6; border-radius: 4px; padding: 2px 5px;
       font-size: 13px; }
pre { padding: 10px; overflow-x: auto; }
.status-FAILED { color: #E45756; } .status-INFEASIBLE { color: #B279A2; }
"""


def write_report(
    run_dir: str | Path,
    objectives: Sequence[ObjectiveSpec],
    filename: str = "report.html",
) -> Path:
    """Render the report. Never raises: a failed plot must not kill a run."""
    run_dir = Path(run_dir)
    out = run_dir / filename
    try:
        html_text = build_report(run_dir, objectives)
    except Exception:  # noqa: BLE001
        log.exception("could not build the report for %s", run_dir)
        return out
    out.write_text(html_text)
    return out


def build_report(run_dir: Path, objectives: Sequence[ObjectiveSpec]) -> str:
    store = TrialStore(Path(run_dir) / "trials.db")
    trials = store.all()
    meta = store.all_meta()
    store.close()

    figures = [
        pareto_figure(trials, objectives),
        hypervolume_figure(trials, objectives),
        parallel_coordinates(trials, objectives),
        scatter_matrix(trials, objectives) if len(objectives) > 2 else None,
        timeline_figure(trials),
        status_figure(trials),
        wall_time_figure(trials),
    ]

    blocks = []
    first = True
    for fig in figures:
        if fig is None:
            continue
        blocks.append(
            fig.to_html(
                full_html=False,
                # plotly.js is inlined once (~3 MB) so the file works offline
                # and survives being emailed around.
                include_plotlyjs="inline" if first else False,
                default_width="100%",
            )
        )
        first = False

    return TEMPLATE.format(
        name=html.escape(str(run_dir.name)),
        generated=_dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        css=CSS,
        cards=_cards(trials, objectives),
        figures="\n".join(blocks),
        table=_table(trials, objectives),
        meta=_meta_block(meta, run_dir),
    )


# ── pieces ───────────────────────────────────────────────────────────────────


def _cards(trials: Sequence[Trial], objectives: Sequence[ObjectiveSpec]) -> str:
    counts: dict[str, int] = {}
    for t in trials:
        counts[t.status.value] = counts.get(t.status.value, 0) + 1

    y, kept = objective_matrix(trials, objectives)
    ref = ref_point_tensor(objectives)
    hv = hypervolume(y, ref) if (ref is not None and kept) else None
    front = int(pareto_mask(y).sum()) if kept else 0
    done = counts.get(TrialStatus.COMPLETED.value, 0)
    failed = counts.get(TrialStatus.FAILED.value, 0) + counts.get(
        TrialStatus.INFEASIBLE.value, 0
    )
    rate = (100.0 * failed / len(trials)) if trials else 0.0

    cards = [
        ("trials", str(len(trials))),
        ("completed", str(done)),
        ("failed / infeasible", f"{failed} ({rate:.0f}%)"),
        ("pareto front", str(front)),
        ("hypervolume", "-" if hv is None else f"{hv:.4g}"),
    ]
    return "\n".join(
        f'<div class="card"><div class="k">{html.escape(k)}</div>'
        f'<div class="v">{html.escape(v)}</div></div>'
        for k, v in cards
    )


def _table(trials: Sequence[Trial], objectives: Sequence[ObjectiveSpec]) -> str:
    y, kept = objective_matrix(trials, objectives)
    front_ids = (
        {
            t.trial_id
            for t, keep in zip(kept, pareto_mask(y).tolist(), strict=False)
            if keep
        }
        if kept
        else set()
    )

    param_names = sorted({k for t in trials for k in t.params})
    obj_names = [o.name for o in objectives]
    header = ["trial", "tag", "status", *obj_names, *param_names, "note"]

    rows = []
    for t in trials:
        classes = []
        if t.trial_id in front_ids:
            classes.append("front")
        if t.tag == "baseline":
            classes.append("baseline")
        cells = [
            f"<td>{html.escape(t.name)}</td>",
            f"<td>{html.escape(t.tag or '')}</td>",
            f'<td class="status-{t.status.value}">{t.status.value}</td>',
        ]
        cells += [
            f'<td class="num">{_num((t.metrics or {}).get(n))}</td>' for n in obj_names
        ]
        cells += [f'<td class="num">{_num(t.params.get(n))}</td>' for n in param_names]
        note = t.error or ("pareto" if t.trial_id in front_ids else "")
        cells.append(f"<td>{html.escape(str(note)[:120])}</td>")
        rows.append(f'<tr class="{" ".join(classes)}">' + "".join(cells) + "</tr>")

    head = "".join(f"<th>{html.escape(h)}</th>" for h in header)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def _meta_block(meta: dict[str, Any], run_dir: Path) -> str:
    import json

    interesting = {
        k: meta[k] for k in ("versions", "ref_point", "evaluator", "seed") if k in meta
    }
    interesting["run_dir"] = str(run_dir)
    body = html.escape(json.dumps(interesting, indent=2, sort_keys=True, default=str))
    config = meta.get("config")
    if config is not None:
        body += "\n\n" + html.escape(
            json.dumps(config, indent=2, sort_keys=True, default=str)
        )
    return f"<pre>{body}</pre>"


def _num(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f"{float(value):.6g}"
    return html.escape(str(value))


TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>mobo — {name}</title>
<style>{css}</style>
</head><body>
<h1>{name}</h1>
<p class="sub">generated {generated}</p>
<div class="cards">{cards}</div>
{figures}
<h2>Trials</h2>
{table}
<h2>Provenance</h2>
{meta}
</body></html>
"""
