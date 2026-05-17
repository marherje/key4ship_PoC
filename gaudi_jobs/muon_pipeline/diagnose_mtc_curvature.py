#!/usr/bin/env python3
"""Diagnostic: truth muon trajectory vs ACTS reco trajectory in MTC.

Reads the truth muon path from MTC SciFi truth hits (timewindows.edm4hep.root)
and overlays it with the ACTS reconstructed track states (ShipHits.root).
Produces per-window (z vs x) and (z vs y) plots in plots/ and prints a
quantitative summary of x-divergence in MTC.

Notes on units (per project conventions):
  * MTC SciFi `position.x` from SND_SciFiAction is stored in ROOT cm (not mm).
    Multiply by 10 to recover mm. position.y / position.z are mm.
  * ACTSTrackStates x, y, z are all in mm.
  * MTC station centers (from parse_geometry) are in mm.
"""

import argparse
import math
import os
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
sys.path.insert(0, str(PROJECT / "simulation" / "geometry"))
from parse_geometry import SNDGeometry  # type: ignore  # noqa: E402

# ---------------------------------------------------------------------------
# Geometry constants
# ---------------------------------------------------------------------------
GEO = SNDGeometry()
MTC_Z_CENTERS = GEO.mtc_station_z_centers          # [mm]
MTC_TOTAL_LEN = GEO.mtc_n_layers * GEO.mtc_layer_thick
MTC_Z_MIN = min(MTC_Z_CENTERS) - MTC_TOTAL_LEN / 2.0 - 5.0
MTC_Z_MAX = max(MTC_Z_CENTERS) + MTC_TOTAL_LEN / 2.0 + 5.0
TAN5 = math.tan(5.0 * math.pi / 180.0)

# ---------------------------------------------------------------------------
# Truth hits reader (podio)
# ---------------------------------------------------------------------------
def read_truth_mtc_hits(edm4hep_file, collection="MTCSciFiHitsWindowed"):
    """Returns dict {frame_idx: [(x_mm, y_mm, z_mm), ...]} for the given collection."""
    try:
        from podio import reading
        reader = reading.get_reader(edm4hep_file)
    except Exception:
        import podio
        reader = podio.root_io.Reader(edm4hep_file)

    by_window = defaultdict(list)
    for iframe, frame in enumerate(reader.get("events")):
        try:
            hits = frame.get(collection)
        except Exception:
            continue
        for hit in hits:
            pos = hit.getPosition()
            # SND_SciFiAction stores hit.position.x = seg_p.x() * 10 (cm → mm),
            # y/z come from G4 step in mm. All three are in mm in the EDM4HEP file.
            by_window[iframe].append((float(pos.x), float(pos.y), float(pos.z)))
    return dict(by_window)


# ---------------------------------------------------------------------------
# Reco track states reader (RNTuple)
# ---------------------------------------------------------------------------
def read_reco_track_states(ship_hits_file):
    """Returns dict {window_id: {track_id: [(x_mm, y_mm, z_mm, loc0, tilt), ...]}}.

    Uses RDF.FromRNTuple → AsNumpy to avoid RNTupleReader teardown segfaults
    that occur in this ROOT build when the reader is held by Python.
    """
    out = defaultdict(lambda: defaultdict(list))
    try:
        df = ROOT.RDF.FromRNTuple("ACTSTrackStates", ship_hits_file)
        cols = ["window_id", "track_id", "x", "y", "z", "loc0", "tilt"]
        data = df.AsNumpy(cols)
    except Exception as e:
        print(f"[reco] cannot open ACTSTrackStates: {e}")
        return out
    n = len(data["window_id"])
    for i in range(n):
        wid = int(data["window_id"][i])
        tid = int(data["track_id"][i])
        out[wid][tid].append((
            float(data["x"][i]),
            float(data["y"][i]),
            float(data["z"][i]),
            float(data["loc0"][i]),
            float(data["tilt"][i]),
        ))
    return out


def stereo_pair_track_states(states):
    """Apply the same stereo pairing as event_display_eve.read_track_states.
    Returns list of (x_mm, y_mm, z_mm) sorted by z."""
    # U/V stereo planes within one MTC layer are 2.35 mm apart in z.
    # Sort by z so adjacent indices are siblings, then pair within a 5 mm window.
    states_sorted = sorted(states, key=lambda s: s[2])
    pts = []
    i = 0
    while i < len(states_sorted):
        x, y, z, loc0, tilt = states_sorted[i]
        if abs(tilt) > 0.01:
            if i + 1 < len(states_sorted):
                _, _, z2, loc02, tilt2 = states_sorted[i + 1]
                if abs(z - z2) < 5.0 and tilt * tilt2 < 0:
                    loc0_U = loc0  if tilt > 0 else loc02
                    loc0_V = loc02 if tilt > 0 else loc0
                    gx = (loc0_U + loc0_V) / 2.0
                    gy = (loc0_V - loc0_U) / (2.0 * TAN5)
                    gz = (z + z2) / 2.0
                    pts.append((gx, gy, gz))
                    i += 2
                    continue
            pts.append((x, y, z))
        else:
            pts.append((x, y, z))
        i += 1
    pts.sort(key=lambda p: p[2])
    return pts


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_window(window_id, truth_path, reco_tracks, outdir, raw_truth_pts=None):
    if not truth_path and not reco_tracks:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax_zx, ax_zy = axes

    if raw_truth_pts:
        rzs = [p[2] for p in raw_truth_pts]
        rxs = [p[0] for p in raw_truth_pts]
        rys = [p[1] for p in raw_truth_pts]
        ax_zx.scatter(rzs, rxs, s=6, c="lightgray", marker=".",
                      label=f"All truth hits ({len(raw_truth_pts)})", alpha=0.5)
        ax_zy.scatter(rzs, rys, s=6, c="lightgray", marker=".", alpha=0.5)

    if truth_path:
        tx = np.array([p[0] for p in truth_path])
        ty = np.array([p[1] for p in truth_path])
        tz = np.array([p[2] for p in truth_path])
        ax_zx.plot(tz, tx, "k-o", ms=3,
                   label=f"Truth muon path ({len(truth_path)} layers)")
        ax_zy.plot(tz, ty, "k-o", ms=3,
                   label=f"Truth muon path ({len(truth_path)} layers)")

    colors = ["red", "blue", "green", "magenta"]
    for it, (tid, raw_states) in enumerate(reco_tracks.items()):
        pts = stereo_pair_track_states(raw_states)
        if not pts:
            continue
        rx = np.array([p[0] for p in pts])
        ry = np.array([p[1] for p in pts])
        rz = np.array([p[2] for p in pts])
        c = colors[it % len(colors)]
        ax_zx.plot(rz, rx, "-o", ms=3, c=c, label=f"Reco track {tid}")
        ax_zy.plot(rz, ry, "-o", ms=3, c=c, label=f"Reco track {tid}")

    for ax, ylabel in [(ax_zx, "x [mm]"), (ax_zy, "y [mm]")]:
        for zc in MTC_Z_CENTERS:
            ax.axvspan(zc - MTC_TOTAL_LEN / 2.0, zc + MTC_TOTAL_LEN / 2.0,
                       alpha=0.07, color="gray")
        ax.axhline(0, lw=0.4, color="gray")
        ax.set_xlim(MTC_Z_MIN, MTC_Z_MAX)
        ax.set_xlabel("z (beam) [mm]")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    ax_zx.set_title(f"Window {window_id}: MTC z–x (curvature plane)")
    ax_zy.set_title(f"Window {window_id}: MTC z–y")

    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"mtc_curvature_w{window_id:03d}.pdf"
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


# ---------------------------------------------------------------------------
# Primary-muon trajectory extraction
# ---------------------------------------------------------------------------
def extract_primary_muon_path(truth_pts, gun_x=5.0, gun_y=-82.5):
    """Recover the muon's physical (dd_x, dd_y, dd_z) trajectory from MTC SciFi
    hits.  Each MTC layer has a U and a V plane (2.35 mm apart in z) with strips
    rotated by ±α=5°.  The hit `position.x` stored by SND_SciFiAction is the
    stereo strip centre `x_stereo = dd_x ∓ dd_y·tan α`, NOT the muon's true x.
    The hit `position.y` is the energy-weighted average step-Y written by the
    SND_SciFiAction::end() method — this IS the true muon y.

    Strategy:
      1. Group hits into layer-pairs (U+V) by z within a 5 mm window.
      2. For paired hits, dd_x = (x_U + x_V) / 2.  dd_y is taken from pos.y
         (already physical).
      3. For unpaired hits (one plane only), recover dd_x from
         pos.x + tan α·pos.y or pos.x − tan α·pos.y, picking the branch that
         lies closer to the previous trajectory point — this resolves the U/V
         ambiguity automatically.
      4. At each layer keep the hit closest to the previous (x, y) to discard
         delta-rays / secondaries.
    """
    if not truth_pts:
        return []
    # Group by 5 mm z-bin so adjacent U/V planes land in the same group.
    bins = defaultdict(list)
    for x, y, z in truth_pts:
        bins[int(round(z / 5.0))].append((x, y, z))

    path = []
    prev_x, prev_y = gun_x, gun_y
    for zk in sorted(bins.keys()):
        layer_hits = bins[zk]
        # Try to pick a hit + partner (opposite-plane companion at same z) closest to (prev_x, prev_y).
        candidates = []
        for hx, hy, hz in layer_hits:
            # Two possible stereo branches for this single hit
            ddx_plus  = hx + TAN5 * hy   # if plane has x_stereo = dd_x - tan α·dd_y
            ddx_minus = hx - TAN5 * hy   # if plane has x_stereo = dd_x + tan α·dd_y
            for ddx_try in (ddx_plus, ddx_minus):
                d2 = (ddx_try - prev_x) ** 2 + (hy - prev_y) ** 2
                candidates.append((d2, ddx_try, hy, hz))
        # Also try true pair-based recovery: pick two hits within same z-layer (diff < 5 mm)
        for i in range(len(layer_hits)):
            for j in range(i + 1, len(layer_hits)):
                hxi, hyi, hzi = layer_hits[i]
                hxj, hyj, hzj = layer_hits[j]
                if abs(hzi - hzj) < 5.0 and abs(hxi - hxj) > 1.0:
                    ddx_pair = (hxi + hxj) / 2.0
                    ddy_pair = (hyi + hyj) / 2.0
                    zp       = (hzi + hzj) / 2.0
                    d2 = (ddx_pair - prev_x) ** 2 + (ddy_pair - prev_y) ** 2
                    # Slight preference for pair recovery (more robust)
                    candidates.append((d2 * 0.5, ddx_pair, ddy_pair, zp))
        if not candidates:
            continue
        candidates.sort(key=lambda c: c[0])
        _, ddx, ddy, ddz = candidates[0]
        path.append((ddx, ddy, ddz))
        prev_x, prev_y = ddx, ddy
    return path


# ---------------------------------------------------------------------------
# Quantitative summary
# ---------------------------------------------------------------------------
def fit_curvature_in_iron(zs, xs):
    """Fit x(z) ~= x0 + slope*(z - z0) + 0.5 * curv * (z - z0)^2 (parabola).
    Returns (x0, slope, curv, R_eff_mm) where R_eff = 1/curv [mm].
    """
    if len(zs) < 3:
        return None
    z = np.asarray(zs, dtype=float)
    x = np.asarray(xs, dtype=float)
    z0 = z.mean()
    coeffs = np.polyfit(z - z0, x, 2)  # x = c2*Δz² + c1*Δz + c0
    c2, c1, c0 = coeffs
    if abs(c2) < 1e-20:
        return None
    R = 1.0 / (2.0 * c2)  # signed radius
    return {"x0": float(c0), "slope": float(c1), "curv": float(2.0 * c2),
            "radius_mm": float(R)}


def momentum_from_radius(radius_mm, B_T=1.7):
    """p [GeV] = 0.3 * B[T] * r[m] * |q|.   r in mm → /1000."""
    return 0.3 * B_T * (radius_mm / 1000.0)


def summarize(window_id, truth_path, reco_tracks):
    """Compute divergence and per-track curvature in MTC region."""
    if not truth_path:
        return None
    tz = np.array([p[2] for p in truth_path])
    tx = np.array([p[0] for p in truth_path])
    mtc_mask = (tz >= MTC_Z_MIN) & (tz <= MTC_Z_MAX)
    if mtc_mask.sum() < 3:
        return None
    tz_mtc = tz[mtc_mask]
    tx_mtc = tx[mtc_mask]
    truth_fit = fit_curvature_in_iron(tz_mtc, tx_mtc)

    results = []
    for tid, raw_states in reco_tracks.items():
        pts = stereo_pair_track_states(raw_states)
        if not pts:
            continue
        rz = np.array([p[2] for p in pts])
        rx = np.array([p[0] for p in pts])
        in_mtc = (rz >= MTC_Z_MIN) & (rz <= MTC_Z_MAX)
        if in_mtc.sum() < 3:
            continue
        rz_m = rz[in_mtc]
        rx_m = rx[in_mtc]
        tx_interp = np.interp(rz_m, tz_mtc, tx_mtc)
        dx = rx_m - tx_interp
        reco_fit = fit_curvature_in_iron(rz_m, rx_m)
        results.append({
            "track_id":      tid,
            "n_in_mtc":      int(in_mtc.sum()),
            "x_dev_mean":    float(np.mean(np.abs(dx))),
            "x_dev_max":     float(np.max(np.abs(dx))),
            "x_first_mtc":   float(rx_m[0]),
            "x_last_mtc":    float(rx_m[-1]),
            "x_first_truth": float(tx_interp[0]),
            "x_last_truth":  float(tx_interp[-1]),
            "reco_fit":      reco_fit,
            "truth_fit":     truth_fit,
        })
    return {"window_id": window_id,
            "n_truth_mtc": int(mtc_mask.sum()),
            "tracks": results}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default=str(HERE / "timewindows.edm4hep.root"),
                    help="EDM4HEP file with truth hits")
    ap.add_argument("--reco",  default=str(HERE / "ShipHits.root"),
                    help="ShipHits.root RNTuple")
    ap.add_argument("--collection", default="MTCSciFiHitsWindowed",
                    help="Truth hit collection")
    ap.add_argument("--max-windows", type=int, default=10,
                    help="Plot only the first N windows")
    ap.add_argument("--outdir", default=str(HERE / "plots_curvature"),
                    help="Output directory for plots")
    args = ap.parse_args()

    print(f"[diag] truth:      {args.truth}")
    print(f"[diag] reco:       {args.reco}")
    print(f"[diag] collection: {args.collection}")

    truth = read_truth_mtc_hits(args.truth, args.collection)
    reco  = read_reco_track_states(args.reco)
    print(f"[diag] truth windows={len(truth)}, reco windows={len(reco)}")

    all_windows = sorted(set(list(truth.keys()) + list(reco.keys())))[: args.max_windows]
    outdir = Path(args.outdir)
    summaries = []
    for w in all_windows:
        truth_raw  = truth.get(w, [])
        truth_path = extract_primary_muon_path(truth_raw)
        reco_trks  = reco.get(w, {})
        outfile = plot_window(w, truth_path, reco_trks, outdir,
                              raw_truth_pts=truth_raw)
        if outfile:
            print(f"[diag] window {w}: truth_raw={len(truth_raw)} pts, "
                  f"truth_path={len(truth_path)} layers, "
                  f"reco tracks={len(reco_trks)} -> {outfile.name}")
        s = summarize(w, truth_path, reco_trks)
        if s and s["tracks"]:
            summaries.append(s)

    print("\n===== MTC PARABOLIC-FIT CURVATURE COMPARISON =====")
    print("Reco x(z) and truth-muon x(z) inside MTC are each fit to a parabola")
    print("x = c0 + c1*Δz + 0.5*curv*Δz².  Effective curvature radius R = 1/curv [mm].")
    print("Equivalent average momentum p_eff [GeV] = 0.3 * B[T=1.7] * R[m] (only valid")
    print("for fully-immersed-in-field; here ~33% iron, so p_eff is a relative metric).")
    print()
    print(f"{'win':>4} {'tid':>3} {'nMTC':>4} "
          f"{'reco R':>10} {'truth R':>10} {'R_t/R_r':>9} "
          f"{'reco p_eff':>11} {'truth p_eff':>12} "
          f"{'|dx|mean':>9} {'|dx|max':>9}")
    for s in summaries:
        for t in s["tracks"]:
            rfit = t.get("reco_fit")  or {}
            tfit = t.get("truth_fit") or {}
            rR = rfit.get("radius_mm", float("nan"))
            tR = tfit.get("radius_mm", float("nan"))
            ratio = tR / rR if (rR and abs(rR) > 1e-6) else float("nan")
            r_p = momentum_from_radius(rR) if rR == rR else float("nan")
            t_p = momentum_from_radius(tR) if tR == tR else float("nan")
            print(f"{s['window_id']:>4d} {t['track_id']:>3d} "
                  f"{t['n_in_mtc']:>4d} "
                  f"{rR:>10.1f} {tR:>10.1f} {ratio:>9.2f} "
                  f"{r_p:>11.3f} {t_p:>12.3f} "
                  f"{t['x_dev_mean']:>9.1f} {t['x_dev_max']:>9.1f}")

    # Aggregate
    if summaries:
        all_r_R = [t["reco_fit"]["radius_mm"]   for s in summaries for t in s["tracks"]
                   if t.get("reco_fit")]
        all_t_R = [t["truth_fit"]["radius_mm"]  for s in summaries for t in s["tracks"]
                   if t.get("truth_fit")]
        if all_r_R and all_t_R:
            print()
            print(f"Median reco  R  = {np.median(np.abs(all_r_R)):8.1f} mm")
            print(f"Median truth R  = {np.median(np.abs(all_t_R)):8.1f} mm")
            print(f"Median |R_truth / R_reco| (over-curvature factor) = "
                  f"{np.median(np.abs(np.array(all_t_R) / np.array(all_r_R))):.2f}")
            print(f"Sign agreement (reco_curv * truth_curv > 0): "
                  f"{sum(1 for r, tt in zip(all_r_R, all_t_R) if r * tt > 0)}/{len(all_r_R)}")
    print(f"\nPlots in: {outdir}")


if __name__ == "__main__":
    main()
    # Work around a ROOT 6.38 RNTuple/RColumn teardown segfault when the
    # Python interpreter exits with cached RNTuple readers still alive.
    # The computation has already completed by this point; nothing in
    # interpreter finalisation needs to run.
    os._exit(0)
