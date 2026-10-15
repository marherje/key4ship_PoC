#!/usr/bin/env python3
"""Phase 4.2: how well SNDPFOBuilder measures E_had, E_nu and the vertex.

Reads the PFASummary that SNDPFOBuilder writes (one ParticleID per event; the
order of its parameters is the PFASummaryNames metadata entry) and compares it
with the MC truth:

  1. E_had,reco vs E_had,true -- the definition of phase 1.5, the kinetic
     energy of the primary final-state hadrons -- in bins of E_had,true AND in
     bins of the true vertex z. Next to it, in the same events, the ceiling of
     phase 1.5: E_vis,nomu / E_had,true, the calibrated energy of every hit
     the primary muon does not dominate, with the same e/h applied. The gap
     between the two is what the corridor and the thresholds cost; the rest
     is longitudinal leakage, which no algorithm gets back;
  2. nu_mu CC: E_nu,reco vs E_nu,true (the beam neutrino), in bins, split by
     how the muon momentum was obtained (pfa_p_mu_source: 0 curvature, 1 the
     iron lower bound, -1 no muon found);
  3. vertex: dx, dy, dz against the MC vertex (where the primaries start) for
     both modes, and the fraction more than 3 layers off in z;
  4. leakage: fraction of events with E_had,MTC / E_had > 0.3 per E bin;
  5. if the file has CLUE3D clusters (phase 4.3): clusters per event, how many
     point at the vertex, their energy-weighted purity (RecoMCParticleLinks
     weight), the energy left unlinked or unclustered, and E_had / E_nu with
     both definitions -- pfa_E_had (3D, linked to the vertex) against
     pfa_E_had_layers (the phase-4.1 sum over layer clusters);
  6. if the file subtracted the muon along its trajectory (pfa_E_mu_sub):
     the energy taken out against the muon's true calibrated deposit (hits
     the primary muon dominates, same e/h), on nu_mu CC with the muon found.

Criterion of the plan: |bias of E_nu| < 10 % on nu_mu CC.

Usage: enu_perf.py [pfa.root] [--max-events=N]
"""
import bisect
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from podio.root_io import Reader
from pfa_truth import TruthJoin

HIT_COLLS = ["SiTargetCaloHits", "SiPadCaloHits",
             "MTCSciFiCaloHits", "MTCScintCaloHits"]
CLUSTER_COLLS = ["SiTargetLayerClusters", "SiPadLayerClusters",
                 "MTCSciFiLayerClusters", "MTCScintLayerClusters"]
E_BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 1e9)]
# Phase 1.5 binned the true vertex z like this, and measured the ceiling
# E_vis,nomu / E_had,sum = 0.949 / 0.888 / 0.790 / 0.728 / 0.638.
Z_BINS = [(-2000, -500), (-500, 500), (500, 1500), (1500, 2500), (2500, 3400)]
BOUND_1_5 = [0.949, 0.888, 0.790, 0.728, 0.638]
NU_PDGS = {12, 14, 16}
LEPTONS = {11, 13, 15}
# Same factors SNDPFOBuilder applies (job6_pfa.py), so the ceiling is on the
# same scale as the reconstruction it is compared with.
EOVERH_ECAL = 1.07546
EOVERH_MTC = 1.0


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def sigma68(xs):
    if len(xs) < 3:
        return float("nan")
    s = sorted(xs)
    return 0.5 * (s[int(0.84 * (len(s) - 1))] - s[int(0.16 * (len(s) - 1))])


def bin_label(lo, hi):
    return f"{lo:>5}-{'inf' if hi > 1e8 else hi:<6}"


def summary_names(reader):
    try:
        md = reader.get("metadata")
        return [n.decode() if isinstance(n, bytes) else str(n)
                for n in md[0].get_parameter("PFASummaryNames")]
    except Exception:
        return []


def main(argv):
    path = "/eos/experiment/ship/user/marquezh/snd_pfa/genie/pfa_pfa_ref_a.root"
    nmax = 0
    for a in argv[1:]:
        if a.startswith("--max-events="):
            nmax = int(a.split("=", 1)[1])
        elif not a.startswith("--"):
            path = a

    reader = Reader(path)
    names = summary_names(reader)
    if not names:
        print("no PFASummaryNames metadata: was the file written by "
              "SNDPFOBuilder?")
        return 2
    events = reader.get("events")
    n_total = len(events) if not nmax else min(len(events), nmax)

    rows = []
    layer_z = {}        # (det, layer ordinal) -> z of its first plane
    for iev in range(n_total):
        frame = events[iev]
        try:
            summ = frame.get("PFASummary")
        except Exception:
            del frame
            continue
        if len(summ) == 0:
            del frame
            continue
        s = dict(zip(names, list(summ[0].getParameters())))
        truth = TruthJoin(frame)
        mcs = truth.mcs
        if mcs is None or len(mcs) == 0:
            del frame
            continue

        nu = truth.neutrino()
        lep = truth.primary_lepton()
        muon = truth.primary_muon()
        mu_idx = truth.index_of(muon)
        if nu is None:
            klass = "?"
        elif lep is None:
            klass = "NC"
        else:
            klass = {13: "numu CC", 11: "nue CC"}.get(abs(lep.getPDG()),
                                                       "nutau CC")

        e_had_true = 0.0
        vtx = None
        for p in mcs:
            if p.getGeneratorStatus() != 1:
                continue
            if vtx is None:
                v = p.getVertex()
                vtx = (v.x, v.y, v.z)
            pdg = abs(p.getPDG())
            if pdg in LEPTONS or pdg in NU_PDGS:
                continue
            e_had_true += max(0.0, p.getEnergy() - p.getMass())

        # Ceiling of phase 1.5 on the same events and the same scale.
        e_nomu = 0.0
        e_mu = 0.0
        for coll in HIT_COLLS:
            try:
                hits = frame.get(coll)
            except Exception:
                continue
            k = EOVERH_MTC if coll.startswith("MTC") else EOVERH_ECAL
            for h in hits:
                idx, _ = truth.dominant(truth.cell_truth(h.getCellID()))
                if idx != mu_idx:
                    e_nomu += k * h.getEnergy()
                else:
                    e_mu += k * h.getEnergy()
            del hits

        for coll in CLUSTER_COLLS:
            try:
                for cl in frame.get(coll):
                    sp = list(cl.getShapeParameters())
                    if len(sp) > 3:
                        key = (cl.getType() // 10, int(sp[0]))
                        layer_z[key] = min(layer_z.get(key, sp[3]), sp[3])
            except Exception:
                pass

        # 3D clusters: purity of each hadronic PFO (link weight) and its energy
        pur3d = []
        if "pfa_n_3d" in s and s["pfa_n_3d"] > 0:
            try:
                for l in frame.get("RecoMCParticleLinks"):
                    pf = l.getFrom()
                    if pf.getPDG() == 0:
                        pur3d.append((pf.getEnergy(), l.getWeight(),
                                      pf.getGoodnessOfPID() > 0.5))
            except Exception:
                pass

        rows.append({
            "pur3d": pur3d,
            "klass": klass,
            "e_nu_true": nu.getEnergy() if nu is not None else float("nan"),
            "e_had_true": e_had_true,
            "e_nomu": e_nomu,
            "e_mu": e_mu,
            "vtx": vtx,
            "has_mu_true": mu_idx >= 0,
            **s,
        })
        del frame

    # The layer sequence along the beam, recovered from the clusters (the
    # first plane of every (detector, layer) that ever held one): a z
    # difference becomes a number of layers by counting those in between.
    layers = sorted(layer_z.values())

    def layers_between(z1, z2):
        a, b = sorted((z1, z2))
        return bisect.bisect_left(layers, b) - bisect.bisect_left(layers, a)

    print(f"file   {path}")
    print(f"events {len(rows)}   ({len(layers)} layers seen along the beam)")
    classes = {}
    for r in rows:
        classes[r["klass"]] = classes.get(r["klass"], 0) + 1
    print("event classes:", classes)
    print()

    # --- 1. E_had ----------------------------------------------------------
    print("--- 1. E_had,reco / E_had,true (true = sum of E_kin of the primary "
          "hadrons; numu CC + NC)")
    print(f"{'E_had,true':>14} {'N':>5} {'ratio':>7} {'sigma68/E':>10} "
          f"{'ceiling':>8} {'gap':>7}")
    # Only events without a primary electron or tau: theirs would count as
    # "hadronic" both in the reconstruction and in the ceiling.
    good = [r for r in rows if r["e_had_true"] > 0.05
            and r["klass"] in ("numu CC", "NC")]
    for lo, hi in E_BINS:
        sel = [r for r in good if lo <= r["e_had_true"] < hi]
        if not sel:
            continue
        ratio = [r["pfa_E_had"] / r["e_had_true"] for r in sel]
        ceil_ = [r["e_nomu"] / r["e_had_true"] for r in sel]
        print(f"{bin_label(lo, hi):>14} {len(sel):>5} {median(ratio):>7.3f} "
              f"{sigma68(ratio):>10.3f} {median(ceil_):>8.3f} "
              f"{median(ratio) - median(ceil_):>+7.3f}")
    print()
    print(f"{'z_vtx,true':>14} {'N':>5} {'ratio':>7} {'sigma68/E':>10} "
          f"{'ceiling':>8} {'gap':>7} {'1.5':>6}")
    for (lo, hi), b15 in zip(Z_BINS, BOUND_1_5):
        sel = [r for r in good if r["vtx"] and lo <= r["vtx"][2] < hi]
        if not sel:
            continue
        ratio = [r["pfa_E_had"] / r["e_had_true"] for r in sel]
        ceil_ = [r["e_nomu"] / r["e_had_true"] for r in sel]
        print(f"{bin_label(lo, hi):>14} {len(sel):>5} {median(ratio):>7.3f} "
              f"{sigma68(ratio):>10.3f} {median(ceil_):>8.3f} "
              f"{median(ratio) - median(ceil_):>+7.3f} {b15:>6.3f}")
    # Where the excess over the ceiling comes from: an NC event has no muon
    # to lose, a numu CC whose muon was not found keeps the muon's own
    # deposit (a MIP through metres of sampled tungsten and iron, calibrated
    # as if it were shower energy) inside E_had.
    print()
    print(f"{'class':>22} {'N':>5} {'ratio':>7} {'sigma68/E':>10} "
          f"{'ceiling':>8} {'gap':>7}")
    for label, cut in (("NC", lambda r: r["klass"] == "NC"),
                       ("numu CC, mu found",
                        lambda r: r["klass"] == "numu CC"
                        and r["pfa_flavour"] == 14),
                       ("numu CC, mu not found",
                        lambda r: r["klass"] == "numu CC"
                        and r["pfa_flavour"] != 14)):
        sel = [r for r in good if cut(r)]
        if not sel:
            continue
        ratio = [r["pfa_E_had"] / r["e_had_true"] for r in sel]
        ceil_ = [r["e_nomu"] / r["e_had_true"] for r in sel]
        print(f"{label:>22} {len(sel):>5} {median(ratio):>7.3f} "
              f"{sigma68(ratio):>10.3f} {median(ceil_):>8.3f} "
              f"{median(ratio) - median(ceil_):>+7.3f}")
    all_ratio = [r["pfa_E_had"] / r["e_had_true"] for r in good]
    bias_ehad = median(all_ratio) - 1.0
    res_ehad = sigma68(all_ratio)
    all_ceil = [r["e_nomu"] / r["e_had_true"] for r in good]
    print(f"  all: N={len(good)}  bias={bias_ehad:+.3f}  sigma68/E={res_ehad:.3f}"
          f"  ceiling={median(all_ceil):.3f} (sigma68/E {sigma68(all_ceil):.3f})")
    print()

    # --- 2. E_nu on nu_mu CC ------------------------------------------------
    print("--- 2. E_nu,reco / E_nu,true on numu CC")
    cc = [r for r in rows if r["klass"] == "numu CC"
          and r["e_nu_true"] == r["e_nu_true"] and r["e_nu_true"] > 0]
    print(f"{'E_nu,true':>14} {'N':>5} {'ratio':>7} {'sigma68/E':>10} "
          f"{'mu found':>9}")
    for lo, hi in E_BINS:
        sel = [r for r in cc if lo <= r["e_nu_true"] < hi]
        if not sel:
            continue
        ratio = [r["pfa_E_nu"] / r["e_nu_true"] for r in sel]
        found = sum(1 for r in sel if r["pfa_flavour"] == 14) / len(sel)
        print(f"{bin_label(lo, hi):>14} {len(sel):>5} {median(ratio):>7.3f} "
              f"{sigma68(ratio):>10.3f} {found:>9.2f}")
    print()
    print(f"{'p_mu source':>14} {'N':>5} {'ratio':>7} {'sigma68/E':>10}")
    for src, label in ((0, "curvature"), (1, "iron bound"), (-1, "no muon")):
        sel = [r for r in cc if round(r["pfa_p_mu_source"]) == src]
        if not sel:
            continue
        ratio = [r["pfa_E_nu"] / r["e_nu_true"] for r in sel]
        print(f"{label:>14} {len(sel):>5} {median(ratio):>7.3f} "
              f"{sigma68(ratio):>10.3f}")
    cc_ratio = [r["pfa_E_nu"] / r["e_nu_true"] for r in cc]
    bias_enu = median(cc_ratio) - 1.0
    res_enu = sigma68(cc_ratio)
    found = [r for r in cc if r["pfa_flavour"] == 14]
    f_ratio = [r["pfa_E_nu"] / r["e_nu_true"] for r in found]
    print(f"  all numu CC: N={len(cc)}  bias={bias_enu:+.3f}  "
          f"sigma68/E={res_enu:.3f}")
    print(f"  with the muon found: N={len(found)}  "
          f"bias={median(f_ratio) - 1.0:+.3f}  sigma68/E={sigma68(f_ratio):.3f}")
    print()

    # --- 3. vertex ---------------------------------------------------------
    print("--- 3. vertex (reco - MC)")
    print(f"{'mode':>9} {'N':>5} {'found':>6} {'med dx':>8} {'med dy':>8} "
          f"{'med dz':>8} {'s68 dx':>7} {'s68 dy':>7} {'s68 dz':>7} "
          f"{'>3 lay':>7}")
    for mode, pre in (("clusters", "pfa_vtx_clu_"), ("track", "pfa_vtx_trk_")):
        sel = [r for r in rows if r["vtx"]]
        have = [r for r in sel if not (r[pre + "x"] == 0.0 and r[pre + "y"] == 0.0
                                       and r[pre + "z"] == 0.0)]
        if not have:
            print(f"{mode:>9} {len(sel):>5} {0:>6.2f}")
            continue
        dx = [r[pre + "x"] - r["vtx"][0] for r in have]
        dy = [r[pre + "y"] - r["vtx"][1] for r in have]
        dz = [r[pre + "z"] - r["vtx"][2] for r in have]
        far = sum(1 for r in have
                  if layers_between(r[pre + "z"], r["vtx"][2]) > 3) / len(have)
        print(f"{mode:>9} {len(sel):>5} {len(have)/len(sel):>6.2f} "
              f"{median(dx):>8.1f} {median(dy):>8.1f} {median(dz):>8.1f} "
              f"{sigma68(dx):>7.1f} {sigma68(dy):>7.1f} {sigma68(dz):>7.1f} "
              f"{far:>7.3f}")
    print("  (mm; '>3 lay' = fraction more than 3 layers away in z)")
    print()

    # --- 4. leakage --------------------------------------------------------
    print("--- 4. fraction of events with E_had,MTC / E_had > 0.3 (numu CC + NC)")
    for lo, hi in E_BINS:
        sel = [r for r in good if lo <= r["e_had_true"] < hi
               and r["pfa_E_had"] > 0]
        if not sel:
            continue
        f = sum(1 for r in sel
                if r["pfa_E_had_mtc"] / r["pfa_E_had"] > 0.3) / len(sel)
        print(f"  {bin_label(lo, hi)}  N={len(sel):>5}  {f:.3f}")
    print()

    # --- 5. CLUE3D --------------------------------------------------------
    if rows and "pfa_n_3d" in rows[0]:
        print("--- 5. CLUE3D hadronic clusters")
        n3 = [r["pfa_n_3d"] for r in rows]
        nl = [r["pfa_n_3d_linked"] for r in rows]
        print(f"  3D clusters/event: mean {sum(n3)/len(n3):.2f}  median "
              f"{median(n3):.0f};  linked to the vertex: mean "
              f"{sum(nl)/len(nl):.2f}")
        allp = [x for r in rows for x in r["pur3d"]]
        if allp:
            etot = sum(e for e, _, _ in allp)
            ew = sum(e * w for e, w, _ in allp) / etot if etot > 0 else float("nan")
            hi = sum(e for e, w, _ in allp if w > 0.9) / etot if etot > 0 else float("nan")
            lk = [x for x in allp if x[2]]
            elk = sum(e for e, _, _ in lk)
            ewl = sum(e * w for e, w, _ in lk) / elk if elk > 0 else float("nan")
            print(f"  purity (energy-weighted): all {ew:.3f}, linked {ewl:.3f}; "
                  f"energy in clusters with purity > 0.9: {hi:.3f}")
        tot = [r["pfa_E_had_layers"] for r in rows if r["pfa_E_had_layers"] > 0]
        fu = [r["pfa_E_had_unlinked"] / r["pfa_E_had_layers"] for r in rows
              if r["pfa_E_had_layers"] > 0]
        fc = [r["pfa_E_had_unclustered"] / r["pfa_E_had_layers"] for r in rows
              if r["pfa_E_had_layers"] > 0]
        print(f"  share of the 4.1 sum left out: unlinked mean {sum(fu)/len(fu):.3f}"
              f" (median {median(fu):.3f}), unclustered mean "
              f"{sum(fc)/len(fc):.3f} (median {median(fc):.3f})")
        print()
        print(f"{'':>24} {'4.1 layers':>22} {'4.3 CLUE3D':>22}")
        print(f"{'class':>24} {'ratio':>7} {'s68':>6} {'gap':>7} "
              f"{'ratio':>7} {'s68':>6} {'gap':>7}")
        for label, cut in (("E_had: NC", lambda r: r["klass"] == "NC"),
                           ("E_had: numu CC mu found",
                            lambda r: r["klass"] == "numu CC"
                            and r["pfa_flavour"] == 14),
                           ("E_had: numu CC no mu",
                            lambda r: r["klass"] == "numu CC"
                            and r["pfa_flavour"] != 14),
                           ("E_had: all", lambda r: True)):
            sel = [r for r in good if cut(r)]
            if not sel:
                continue
            cells = []
            ceil_ = median([r["e_nomu"] / r["e_had_true"] for r in sel])
            for key in ("pfa_E_had_layers", "pfa_E_had"):
                ratio = [r[key] / r["e_had_true"] for r in sel]
                cells.append(f"{median(ratio):>7.3f} {sigma68(ratio):>6.3f} "
                             f"{median(ratio) - ceil_:>+7.3f}")
            print(f"{label:>24} {cells[0]} {cells[1]}")
        cells = []
        for key in ("pfa_E_had_layers", "pfa_E_had"):
            ratio = [(r["pfa_E_nu"] - r["pfa_E_had"] + r[key]) / r["e_nu_true"]
                     for r in cc]
            cells.append(f"{median(ratio):>7.3f} {sigma68(ratio):>6.3f} "
                         f"{median(ratio) - 1.0:>+7.3f}")
        print(f"{'E_nu: numu CC':>24} {cells[0]} {cells[1]}")
        print("  (gap = against the ceiling for E_had, against 1 for E_nu)")
        print()

    # --- 6. muon subtraction ------------------------------------------------
    if rows and "pfa_E_mu_sub" in rows[0]:
        print("--- 6. muon deposit taken out along the trajectory vs the truth "
              "(numu CC, muon found)")
        sel = [r for r in rows if r["klass"] == "numu CC"
               and r["pfa_flavour"] == 14 and r["e_mu"] > 0.05]
        if sel:
            ratio = [r["pfa_E_mu_sub"] / r["e_mu"] for r in sel]
            match = [r["pfa_n_mu_matched"] / r["pfa_n_mu_planes"] for r in sel
                     if r["pfa_n_mu_planes"] > 0]
            print(f"  N={len(sel)}  E_mu,sub / E_mu,true: median "
                  f"{median(ratio):.3f}  sigma68 {sigma68(ratio):.3f};  "
                  f"true deposit median {median([r['e_mu'] for r in sel]):.2f} "
                  f"GeV, taken {median([r['pfa_E_mu_sub'] for r in sel]):.2f} "
                  f"GeV;  planes matched {median(match):.2f}")
        nc = [r for r in rows if r["klass"] == "NC" and r["pfa_flavour"] == 14]
        if nc:
            print(f"  NC with a fake muon: N={len(nc)}  taken median "
                  f"{median([r['pfa_E_mu_sub'] for r in nc]):.2f} GeV")
        print()

    ok = abs(bias_enu) < 0.10
    print(f"bias_ehad={bias_ehad:+.3f} res_ehad={res_ehad:.3f} "
          f"bias_enu={bias_enu:+.3f} res_enu={res_enu:.3f}")
    print(("OK" if ok else "FAIL") +
          "  (criterion: |bias of E_nu| < 0.10 on numu CC)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
