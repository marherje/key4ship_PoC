#!/usr/bin/env python3
"""Pareto front of a MOBO campaign, as a figure for slides.

Efficiency against cost, with the bubble area carrying the third objective
(SiPad hits). The front is the subject and the dominated trials are context,
so the front is drawn in the accent hue and everything it beats recedes into
grey -- an emphasis chart, not a categorical one. The baseline gets its own
hue AND its own marker, so it is identifiable without relying on colour.

Usage: plot_pareto.py [--db trials.db] [--out pareto.png] [--dark]
"""
import argparse
import json
import sqlite3
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Tokens from the reference palette (validated: all-pairs CVD dE 24.7 protan,
# 33.6 normal vision, both >= 3:1 on their surface).
LIGHT = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
             grid="#e1e0d9", axis="#c3c2b7", front="#2a78d6", base="#eb6834",
             dim="#c9c8c2")
DARK = dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
            grid="#2c2c2a", axis="#383835", front="#3987e5", base="#d95926",
            dim="#4a4a46")

# Trials worth naming. Labelling all 27 would be unreadable, so these are the
# corners of the front plus the two re-processed for event displays. Every
# label reads "text (tNNNNN)" and every one is tied to its bubble by a leader
# line: with bubbles this dense -- and with the smallest bubble on the plot
# being the cheapest point -- a label floating near a marker is not the same
# as a label pointing at it.
#
# The two named by hand are fixed; the corners are FOUND, not hard-coded, so
# the labels cannot end up on the wrong trial when the campaign changes.
MANUAL = {270: ("best choice", (-96, 16)),
          237: ("max SiPad hits", (34, -30))}


def load(db):
    con = sqlite3.connect(db)
    rows = []
    for tid, params, metrics, tag in con.execute(
            "select trial_id, params, metrics, tag from trials "
            "where status='COMPLETED'"):
        m = json.loads(metrics)
        p = json.loads(params)
        if not all(k in m for k in ("eff_track", "nhits_sipad", "cost_proxy")):
            continue
        rows.append(dict(tid=tid, tag=tag, p=p, **m))
    return rows


def pareto(rows):
    def dominates(a, b):
        ge = (a["eff_track"] >= b["eff_track"]
              and a["nhits_sipad"] >= b["nhits_sipad"]
              and a["cost_proxy"] <= b["cost_proxy"])
        gt = (a["eff_track"] > b["eff_track"]
              or a["nhits_sipad"] > b["nhits_sipad"]
              or a["cost_proxy"] < b["cost_proxy"])
        return ge and gt
    return [r for r in rows if not any(dominates(o, r) for o in rows)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="Optimization/mobo/runs/snd_genie_gen_10k/trials.db")
    ap.add_argument("--out", default="pareto_snd_genie_gen.png")
    ap.add_argument("--dark", action="store_true")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    C = DARK if args.dark else LIGHT
    rows = load(args.db)
    if not rows:
        print("no completed trials", file=sys.stderr)
        return 1
    front = pareto(rows)
    front_ids = {r["tid"] for r in front}
    base = next((r for r in rows if r["tag"] == "baseline"), None)
    dominated = [r for r in rows
                 if r["tid"] not in front_ids and r is not base]

    # Bubble area proportional to hits, so the eye reads area not radius.
    hits = [r["nhits_sipad"] for r in rows]
    hmin, hmax = min(hits), max(hits)
    def area(h):
        t = (h - hmin) / (hmax - hmin) if hmax > hmin else 0.5
        return 30 + t * 340

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    fig.patch.set_facecolor(C["surface"])
    ax.set_facecolor(C["surface"])

    ax.scatter([r["cost_proxy"] for r in dominated],
               [r["eff_track"] for r in dominated],
               s=[area(r["nhits_sipad"]) for r in dominated],
               c=C["dim"], edgecolors="none", alpha=0.55, zorder=2,
               label=f"dominated ({len(dominated)})")

    front_sorted = sorted(front, key=lambda r: r["cost_proxy"])
    ax.plot([r["cost_proxy"] for r in front_sorted],
            [r["eff_track"] for r in front_sorted],
            color=C["front"], lw=1.4, alpha=0.45, zorder=3)
    ax.scatter([r["cost_proxy"] for r in front],
               [r["eff_track"] for r in front],
               s=[area(r["nhits_sipad"]) for r in front],
               c=C["front"], edgecolors=C["surface"], linewidths=2, zorder=4,
               label=f"Pareto front ({len(front)})")

    if base:
        # matplotlib's `s` is the area of the marker's bounding box, so a
        # diamond drawn at the same `s` as a circle reads far bigger -- it
        # keeps the corners the circle does not have. The baseline sits at
        # 90 % of the largest hit count, so its bubble is legitimately large
        # and the shape was making it look like the maximum. Scaling the
        # diamond's area keeps it on the same size scale as the circles.
        ax.scatter([base["cost_proxy"]], [base["eff_track"]],
                   s=area(base["nhits_sipad"]) * 0.45, c=C["base"], marker="D",
                   edgecolors=C["surface"], linewidths=2, zorder=5,
                   label="baseline")
        ax.annotate(f"baseline (t{base['tid']:05d})",
                    (base["cost_proxy"], base["eff_track"]),
                    textcoords="offset points", xytext=(18, -30),
                    color=C["base"], fontsize=10, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=C["base"], lw=0.9,
                                    shrinkA=2, shrinkB=5))

    # Corner labels, resolved from the data. The cheapest point is also the
    # one with the fewest hits, so its bubble is the smallest on the plot and a
    # bare label next to it reads as belonging to the big neighbour -- both
    # corners therefore get a leader line to the marker they name.
    best = max(front, key=lambda r: r["eff_track"])
    cheap = min(front, key=lambda r: r["cost_proxy"])
    labels = dict(MANUAL)
    labels[best["tid"]] = ("best efficiency", (10, -40))
    labels[cheap["tid"]] = ("cheapest", (30, 26))

    for r in front:
        if r["tid"] not in labels:
            continue
        text, off = labels[r["tid"]]
        ax.annotate(
            f"{text} (t{r['tid']:05d})", (r["cost_proxy"], r["eff_track"]),
            textcoords="offset points", xytext=off,
            color=C["ink2"], fontsize=9,
            arrowprops=dict(arrowstyle="-", color=C["muted"], lw=0.9,
                            shrinkA=2, shrinkB=5))

    ax.set_xlabel("cost proxy  [k€, silicon area + tungsten mass]",
                  color=C["ink2"], fontsize=11)
    ax.set_ylabel("track efficiency  [events with a global track]",
                  color=C["ink2"], fontsize=11)
    ax.set_title("SND geometry optimisation — Pareto front over 200 trials",
                 color=C["ink"], fontsize=13, fontweight="bold", pad=14,
                 loc="left")
    # NOTE: the subtitle that used to sit here explained that the front is
    # three-dimensional and this is its projection, which is why a grey point
    # can sit above the blue line without contradicting it -- it is beaten on
    # SiPad hits. Removed on request for the slides; the point still needs
    # making somewhere, so it belongs in the caption or spoken over the slide.
    ax.margins(y=0.09)
    ax.grid(True, color=C["grid"], lw=0.8, zorder=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(C["axis"])
    ax.tick_params(colors=C["muted"], labelsize=10)

    # Two legends: identity, and what the bubble area means.
    leg1 = ax.legend(loc="lower right", frameon=False, fontsize=10,
                     labelcolor=C["ink2"], scatterpoints=1,
                     markerscale=0.55)
    ax.add_artist(leg1)
    size_handles = [
        Line2D([], [], marker="o", linestyle="none",
               markersize=(area(h) ** 0.5) / 1.9,
               markerfacecolor=C["dim"], markeredgecolor="none",
               label=f"{h/1000:.0f}k")
        for h in (hmin, (hmin + hmax) / 2, hmax)]
    # Parked in the empty middle-bottom: the cheap corner of the front and the
    # identity legend already own the two lower corners.
    ax.legend(handles=size_handles, loc="lower left",
              bbox_to_anchor=(0.24, 0.005), frameon=False, ncol=3,
              fontsize=9, labelcolor=C["ink2"], title="SiPad hits (bubble area)",
              title_fontsize=9, handletextpad=1.2, columnspacing=2.6,
              borderpad=0.6)
    ax.get_legend().get_title().set_color(C["muted"])

    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, facecolor=C["surface"])
    print(f"wrote {args.out}  ({len(front)} front, {len(dominated)} dominated)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
