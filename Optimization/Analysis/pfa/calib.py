#!/usr/bin/env python3
"""Phase 1.4: fit the global energy scales from the single-particle guns.

Reads the PFA output of every gun (produced by
gaudi_jobs/calib_pipeline/run_calib.sh) and measures, in this order:

  1. MIP peak. For the 10 GeV muon, the mode of the per-hit MIP value
     (SNDCaloHitBuilder stores it in energyError). It should sit at 1: a muon
     crossing a plane is the definition of the scale the weights use.
  2. Electromagnetic scale. Sum of E_reco over SiTarget+SiPad per event vs the
     electron energy, fitted through the origin. GlobalScaleECAL is 1/slope.
  3. Hadronic scale. Same over all four collections for pions. The ratio of
     the two slopes is e/h; GlobalScaleMTC balances the MTC against the
     silicon detectors by asking pion events to be linear when it is included.
  4. Leakage. Fraction of a pion's visible energy that lands in the MTC, per
     energy -- a design number for the MTC, and the reason E_had needs it.

Usage: calib.py [--dir /eos/.../snd_pfa/calib/reco]
"""
import math
import os
import re
import sys
from collections import defaultdict

from podio.root_io import Reader

ECAL_COLLS = ["SiTargetCaloHits", "SiPadCaloHits"]
MTC_COLLS = ["MTCSciFiCaloHits", "MTCScintCaloHits"]
DEFAULT_DIR = "/eos/experiment/ship/user/marquezh/snd_pfa/calib/reco"

# MIP scale [MeV] each plane is being reconstructed with, i.e. what job6's
# MIPEnergyOverride says. Phase 1.4 measured these from the 10 GeV muon.
CURRENT_MIP_MEV = {
    (0, 0): 0.08432, (0, 1): 0.08432,   # SiTarget X, Y
    (1, 2): 0.19531,                    # SiPad
    (2, 0): 0.19757, (2, 1): 0.19757,   # MTC SciFi U, V
    (2, 2): 3.53253,                    # MTC scintillator
}


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def mode_of(xs, lo=0.0, hi=3.0, nbins=60):
    """Peak of a histogram of xs -- the Landau most probable value, which is
    what 'the MIP scale' means; the mean of a Landau is dominated by its tail."""
    if not xs:
        return float("nan")
    w = (hi - lo) / nbins
    h = [0] * nbins
    for x in xs:
        if lo <= x < hi:
            h[int((x - lo) / w)] += 1
    if not any(h):
        return float("nan")
    k = max(range(nbins), key=lambda i: h[i])
    return lo + (k + 0.5) * w


def fit_through_origin(pairs):
    """Least squares y = a*x through the origin, plus the largest relative
    residual, which is the honest statement of linearity."""
    sxx = sum(x * x for x, _ in pairs)
    sxy = sum(x * y for x, y in pairs)
    if sxx == 0:
        return float("nan"), float("nan")
    a = sxy / sxx
    worst = max(abs(y - a * x) / (a * x) for x, y in pairs if x > 0)
    return a, worst


def sums_per_event(path):
    """Per-event (E_ecal, E_mtc) in GeV, plus every hit's MIP value."""
    events = Reader(path).get("events")
    out = []
    mips = defaultdict(list)
    for i in range(len(events)):
        frame = events[i]
        e_ecal = e_mtc = 0.0
        for coll in ECAL_COLLS + MTC_COLLS:
            try:
                hits = frame.get(coll)
            except Exception:
                continue
            for h in hits:
                if coll in ECAL_COLLS:
                    e_ecal += h.getEnergy()
                else:
                    e_mtc += h.getEnergy()
                t = h.getType()
                mips[(t // 10, t % 10)].append(h.getEnergyError())
            del hits
        out.append((e_ecal, e_mtc))
        del frame
    return out, mips


def main(argv):
    d = DEFAULT_DIR
    for a in argv[1:]:
        if a.startswith("--dir="):
            d = a.split("=", 1)[1]

    # Labels are pfa_calib_<part>_<E>GeV[_z<vertex>].root. The vertex matters:
    # the default set starts inside the SiTarget and calibrates the silicon,
    # while the _z200 set starts inside MTC40 and is the only way to calibrate
    # the MTC -- a pion fired from the SiTarget leaves under 3 % of its visible
    # energy there, and that is where most GENIE neutrinos actually interact.
    files = {}
    for fn in sorted(os.listdir(d)):
        m = re.match(r"pfa_calib_(.+?)_([0-9.]+)GeV(?:_z(-?[0-9]+))?\.root$", fn)
        if m:
            vtx = int(m.group(3)) if m.group(3) else None
            files[(m.group(1), float(m.group(2)), vtx)] = os.path.join(d, fn)
    if not files:
        print(f"no calibration files in {d}")
        return 2

    data = {}
    mip_by_plane = {}
    mip_by_plane_mtc = {}
    for key, path in sorted(files.items(), key=lambda kv: str(kv[0])):
        sums, mips = sums_per_event(path)
        data[key] = sums
        if key[0] == "mu-" and key[2] is None:
            mip_by_plane = mips
        if key[0] == "mu-" and key[2] is not None:
            mip_by_plane_mtc = mips
        print(f"read {os.path.basename(path):<40} {len(sums):>5} events")

    name = {(0, 0): "SiTarget X", (0, 1): "SiTarget Y", (1, 2): "SiPad",
            (2, 0): "MTC U", (2, 1): "MTC V", (2, 2): "MTC scint"}
    print("\n=== 1. MIP peak (mu- 10 GeV; expect ~1)")
    for k in sorted(mip_by_plane):
        xs = mip_by_plane[k]
        print(f"  {name.get(k, str(k)):<12} N={len(xs):>7}  mode={mode_of(xs):.2f}"
              f"  median={median(xs):.2f}")

    print("\n=== 2. Electromagnetic scale (e-, SiTarget+SiPad)")
    e_pairs = []
    for (part, energy, vtx), sums in sorted(data.items(), key=lambda kv: str(kv[0])):
        if part != "e-" or vtx is not None:
            continue
        vis = median([s[0] for s in sums])
        e_pairs.append((energy, vis))
        rel = [(s[0] / vis - 1.0) for s in sums] if vis > 0 else []
        res = 0.5 * (sorted(rel)[int(0.84 * len(rel))] -
                     sorted(rel)[int(0.16 * len(rel))]) if len(rel) > 3 else float("nan")
        print(f"  E={energy:>5.1f} GeV  median E_vis={vis:>9.2f}  sigma68/E={res:>6.3f}")
    slope_em, worst_em = fit_through_origin(e_pairs)
    print(f"  slope = {slope_em:.3f}  ->  GlobalScaleECAL = {1.0/slope_em:.5f}"
          f"   worst residual {100*worst_em:.1f} %")

    print("\n=== 3. Hadronic scale (pi+, all four collections)")
    h_pairs = []
    for (part, energy, vtx), sums in sorted(data.items(), key=lambda kv: str(kv[0])):
        if part != "pi+" or vtx is not None:
            continue
        vis = median([s[0] + s[1] for s in sums])
        h_pairs.append((energy, vis))
        fr = [s[1] / (s[0] + s[1]) for s in sums if (s[0] + s[1]) > 0]
        # The MEAN matters here, not the median: most pions leave nothing in
        # the MTC at all, so the median is 0 by construction and says nothing
        # about the tail that does leak.
        mean_leak = sum(fr) / len(fr) if fr else float("nan")
        frac_any = sum(1 for x in fr if x > 0) / len(fr) if fr else float("nan")
        print(f"  E={energy:>5.1f} GeV  median E_vis={vis:>9.2f}"
              f"  E_vis/E={vis/energy:>6.3f}"
              f"  mean MTC fraction={mean_leak:>6.3f}"
              f"  events reaching MTC={frac_any:>6.3f}")
    slope_h, worst_h = fit_through_origin(h_pairs)
    print(f"  slope = {slope_h:.3f}   worst residual {100*worst_h:.1f} %")
    if slope_h > 0:
        print(f"  e/h = slope_em / slope_h = {slope_em/slope_h:.3f}"
              f"   ->  EoverH = {slope_em/slope_h:.5f}")

    # The MIP scale the weights assume is the MEAN minimum ionising loss,
    # while what a thin layer actually records is the Landau most probable
    # value, which is smaller (delta rays escape). The muon measures the
    # ratio, so turn the measured modes into MIPEnergyOverride entries and the
    # MIP scale becomes empirical instead of theoretical.
    if mip_by_plane:
        print("\nMIPEnergyOverride from the muon (det:view:MeV):")
        # These are the MIP scales of the run being analysed. If job6 already
        # carries MIPEnergyOverride, the modes below are ~1 and the suggestion
        # is the same value again: the correction is multiplicative, so it
        # converges after one pass and must be CHAINED on the value in use,
        # never recomputed from the sampling scan.
        current = dict(CURRENT_MIP_MEV)
        entries = []
        for k in sorted(mip_by_plane):
            mode = mode_of(mip_by_plane[k])
            if k in current and mode == mode:
                entries.append(f"{k[0]}:{k[1]}:{current[k]*mode:.5f}")
                print(f"  {name.get(k, str(k)):<12} mode={mode:.2f}"
                      f"  {current[k]:.5f} -> {current[k]*mode:.5f} MeV")
        print("  " + repr(entries))
        print("  (chained on the scales in CURRENT_MIP_MEV; update that dict "
              "when job6's MIPEnergyOverride changes)")

    # --- 4. the MTC, from the guns that start inside it -------------------
    mtc_pairs = []
    print("\n=== 4. MTC scale (guns starting inside MTC40, z = +200 mm)")
    if mip_by_plane_mtc:
        for k in sorted(mip_by_plane_mtc):
            if k[0] == 2:
                print(f"  MIP {name.get(k, str(k)):<12} mode="
                      f"{mode_of(mip_by_plane_mtc[k]):.2f}")
    for (part, energy, vtx), sums in sorted(data.items(), key=lambda kv: str(kv[0])):
        if vtx is None or part != "pi+":
            continue
        e_mtc = median([s[1] for s in sums])
        e_all = median([s[0] + s[1] for s in sums])
        frac = median([s[1] / (s[0] + s[1]) for s in sums if (s[0] + s[1]) > 0])
        mtc_pairs.append((energy, e_mtc))
        print(f"  E={energy:>5.1f} GeV  E_vis(MTC)={e_mtc:>8.2f}"
              f"  E_vis(all)={e_all:>8.2f}  MTC fraction={frac:>6.3f}"
              f"  E_vis/E={e_all/energy:>6.3f}")
    if mtc_pairs:
        slope_mtc, worst_mtc = fit_through_origin(mtc_pairs)
        print(f"  slope(MTC only) = {slope_mtc:.3f}"
              f"   worst residual {100*worst_mtc:.1f} %")
        # Ask the MTC to answer like the silicon does for the same particle:
        # the electromagnetic scale is the common reference, since that is
        # what GlobalScaleECAL was fitted to.
        if slope_mtc > 0:
            print(f"  GlobalScaleMTC = slope_em / slope(MTC) = "
                  f"{slope_em/slope_mtc:.5f}")

    print("\nSuggested job6 defaults:")
    print(f"  GlobalScaleECAL = {1.0/slope_em:.5f}")
    print(f"  GlobalScaleMTC  = 1.0   (see note)")
    if slope_h > 0:
        print(f"  EoverH          = {slope_em/slope_h:.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
