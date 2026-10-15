#!/usr/bin/env python3
"""Phase 2.3: how good the per-layer clusters are.

Measures, on the PFA output of the GENIE reference sample:

  1. occupancy -- clusters per event, hits and width per cluster, by detector;
  2. purity -- the fraction of a cluster's energy left by its dominant MC
     particle. The guion asks for a median >= 0.9 in the SiPad and the MTC
     scintillator, the two 2D readouts; the 1D strip views are reported
     without a criterion, since strips inside a hadronic shower are mixed by
     construction;
  3. granularity -- clusters per layer against the number of distinct MC
     particles depositing in it, whose ratio says whether the clusterer is
     splitting what it should;
  4. muon isolation -- for the primary muon of nu_mu CC events, the fraction
     of the layers it crosses where its hits form a cluster of their own:
     at most MaxIsoMIP MIPs, at most MaxIsoChannels channels and purity above
     MinIsoPurity. This is the number phase 3's corridor depends on, so it is
     split by whether another cluster sits within NearbyMM of it -- inside a
     hadronic shower or clear of it.

Usage: layer_perf.py [pfa.root] [--max-events N]
"""
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from podio.root_io import Reader
from pfa_truth import TruthJoin

CLUSTER_COLLS = {
    "SiTargetLayerClusters": "SiTarget",
    "SiPadLayerClusters": "SiPad",
    "MTCSciFiLayerClusters": "MTC SciFi",
    "MTCScintLayerClusters": "MTC scint",
}
MAX_ISO_MIP = 3.0
MAX_ISO_CHANNELS = 3
MIN_ISO_PURITY = 0.8
NEARBY_MM = 20.0


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main(argv):
    path = "/eos/experiment/ship/user/marquezh/snd_pfa/genie/pfa_pfa_ref_a.root"
    nmax = 0
    for a in argv[1:]:
        if a.startswith("--max-events="):
            nmax = int(a.split("=", 1)[1])
        elif not a.startswith("--"):
            path = a

    events = Reader(path).get("events")
    n_total = len(events)
    if nmax:
        n_total = min(n_total, nmax)

    per_event = defaultdict(list)
    nhits = defaultdict(list)
    widths = defaultdict(list)
    purity = defaultdict(list)
    ratio_layers = defaultdict(list)
    iso_in_shower = [0, 0]      # [isolated, total]
    iso_clear = [0, 0]

    for iev in range(n_total):
        frame = events[iev]
        truth = TruthJoin(frame)
        muon = truth.primary_muon()
        mu_idx = truth.index_of(muon)

        all_clusters = []
        for coll, label in CLUSTER_COLLS.items():
            try:
                clusters = frame.get(coll)
            except Exception:
                continue
            per_event[label].append(len(clusters))
            # particles per layer, to compare with clusters per layer
            by_layer = defaultdict(lambda: [0, set()])
            for cl in clusters:
                sp = list(cl.getShapeParameters())
                layer = int(sp[0]) if sp else -1
                n = int(sp[1]) if len(sp) > 1 else 0
                w = sp[2] if len(sp) > 2 else 0.0
                mips = sp[4] if len(sp) > 4 else 0.0
                nhits[label].append(n)
                widths[label].append(w)

                acc = {}
                for h in cl.getHits():
                    for idx, e in truth.cell_truth(h.getCellID()).items():
                        acc[idx] = acc.get(idx, 0.0) + e
                idx, pur = truth.dominant(acc)
                purity[label].append(pur)

                by_layer[layer][0] += 1
                by_layer[layer][1].update(acc.keys())
                all_clusters.append((cl, layer, n, mips, idx, pur))

            for layer, (ncl, parts) in by_layer.items():
                if parts:
                    ratio_layers[label].append(ncl / max(1, len(parts)))
            del clusters

        # muon isolation
        if mu_idx >= 0:
            for cl, layer, n, mips, idx, pur in all_clusters:
                if idx != mu_idx:
                    continue
                isolated = (mips <= MAX_ISO_MIP and n <= MAX_ISO_CHANNELS
                            and pur >= MIN_ISO_PURITY)
                p = cl.getPosition()
                nearby = False
                for other, olayer, _, _, _, _ in all_clusters:
                    if olayer != layer or other.getObjectID() == cl.getObjectID():
                        continue
                    q = other.getPosition()
                    if math.dist((p.x, p.y, p.z), (q.x, q.y, q.z)) < NEARBY_MM:
                        nearby = True
                        break
                bucket = iso_in_shower if nearby else iso_clear
                bucket[1] += 1
                bucket[0] += 1 if isolated else 0
        del frame

    print(f"file   {path}")
    print(f"events {n_total}\n")
    # Medians alone mislead here: with the flux correctly placed, most
    # neutrinos interact downstream of the silicon, so more than half the
    # events leave nothing in the SiTarget and its median cluster count is 0
    # while its mean is in the hundreds. Both are printed, and so is the
    # fraction of clusters wider than one channel -- a detector whose clusters
    # are all single channels is not really being clustered.
    print(f"{'detector':<11}{'cl/ev mean':>11}{'median':>8}{'hits mean':>10}"
          f"{'>1 ch':>7}{'width mean':>11}{'purity med':>11}{'pur=1':>7}"
          f"{'cl/part':>8}")
    for label in ["SiTarget", "SiPad", "MTC SciFi", "MTC scint"]:
        if label not in per_event:
            continue
        p = purity[label]
        n = nhits[label]
        frac1 = sum(1 for x in p if x > 0.999) / len(p) if p else float("nan")
        multi = sum(1 for x in n if x > 1) / len(n) if n else float("nan")
        mean = lambda xs: (sum(xs) / len(xs)) if xs else float("nan")
        print(f"{label:<11}{mean(per_event[label]):>11.1f}"
              f"{median(per_event[label]):>8.0f}{mean(n):>10.2f}"
              f"{multi:>7.2f}{mean(widths[label]):>11.2f}"
              f"{median(p):>11.2f}{frac1:>7.2f}"
              f"{median(ratio_layers[label]):>8.2f}")

    print(f"\nMuon isolation (cluster of its own: <= {MAX_ISO_MIP} MIP, "
          f"<= {MAX_ISO_CHANNELS} channels, purity >= {MIN_ISO_PURITY})")
    for name, b in (("inside a shower", iso_in_shower), ("clear", iso_clear)):
        frac = b[0] / b[1] if b[1] else float("nan")
        print(f"  {name:<16} {b[0]:>6}/{b[1]:<6} = {frac:.3f}")

    ok = (median(purity.get("SiPad", [0])) >= 0.9 and
          median(purity.get("MTC scint", [0])) >= 0.9)
    print("\n" + ("OK" if ok else "FAIL")
          + "  (criterion: median purity >= 0.90 in SiPad and MTC scint)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
