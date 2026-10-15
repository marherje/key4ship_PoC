#!/usr/bin/env python3
"""Where the nu_mu CC muons that SNDMuonFinder misses are lost.

Every nu_mu CC event goes to the FIRST stage it fails:

  1. acceptance -- the primary muon lit fewer than --min-layers MTC layers
     (truth, SciFi U view), so no corridor could ever reach MinMTCLayers;
  2. tracking   -- it did light them, but no GlobalTracks / MTCTracks track
     is dominated by it;
  3. finder     -- it has a track but no MuonCandidate is the muon. The cut
     that rejects it is read from a second PFA file made with the finder's
     cuts opened up (MinMTCLayers 0, MinFracLayers 0, MaxMeanMIP 1e9), where
     every track that reaches the MTC becomes a candidate with its corridor
     quantities in MuonPID;
  4. found.

Usage: muon_loss.py nominal_pfa.root loose_pfa.root [--max-events=N]
                    [--min-layers=8] [--min-frac=0.6] [--max-mean-mip=2.0]
"""
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from podio.root_io import Reader
from pfa_truth import TruthJoin

TRACKS = ("GlobalTracks", "MTCTracks")
P_BINS = [(0, 5), (5, 20), (20, 1e9)]


def muon_track_ids(frame, truth, mu_idx):
    """(collection, index) of the tracks the primary muon dominates."""
    ids = set()
    for name in TRACKS:
        try:
            coll = frame.get(name)
        except Exception:
            continue
        for t in coll:
            idx, _, _ = truth.track_truth(t)
            if idx == mu_idx:
                ids.add((name, t.getObjectID().index,
                         t.getObjectID().collectionID))
    return ids


def candidate_is_muon(cand, truth, mu_idx):
    acc = {}
    for cl in cand.getClusters():
        for h in cl.getHits():
            for idx, e in truth.cell_truth(h.getCellID()).items():
                acc[idx] = acc.get(idx, 0.0) + e
    return truth.dominant(acc)[0] == mu_idx


def mtc_layers_of(frame, truth, mu_idx):
    """MTC layers where the muon dominates a SciFi U hit (one per layer)."""
    zs = set()
    try:
        hits = frame.get("MTCSciFiCaloHits")
    except Exception:
        return 0
    for h in hits:
        if h.getType() != 20:
            continue
        if truth.dominant(truth.cell_truth(h.getCellID()))[0] == mu_idx:
            zs.add(round(h.getPosition().z))
    return len(zs)


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) for a in argv[1:] if a.startswith("--"))
    if len(args) != 2:
        print(__doc__)
        return 2
    nmax = int(opts.get("max-events", 0))
    min_layers = int(opts.get("min-layers", 8))
    min_frac = float(opts.get("min-frac", 0.6))
    max_mean_mip = float(opts.get("max-mean-mip", 2.0))

    nominal = Reader(args[0]).get("events")
    loose = Reader(args[1]).get("events")
    n = min(len(nominal), len(loose))
    if nmax:
        n = min(n, nmax)

    stage = Counter()
    stage_by_p = {b: Counter() for b in P_BINS}
    cut = Counter()          # finder losses: which cut failed first
    wrong_track = 0          # finder losses with no loose candidate either
    loose_vals = []          # (nLayersMTC, frac, meanMip) of lost muons

    for iev in range(n):
        frame = nominal[iev]
        truth = TruthJoin(frame)
        muon = truth.primary_muon()
        if muon is None:
            continue
        mu_idx = truth.index_of(muon)
        pm = muon.getMomentum()
        p = math.sqrt(pm.x ** 2 + pm.y ** 2 + pm.z ** 2)
        pbin = next(b for b in P_BINS if b[0] <= p < b[1])

        found = any(candidate_is_muon(c, truth, mu_idx)
                    for c in frame.get("MuonCandidates"))
        if found:
            s = "4 found"
        elif mtc_layers_of(frame, truth, mu_idx) < min_layers:
            s = "1 acceptance"
        elif not muon_track_ids(frame, truth, mu_idx):
            s = "2 tracking"
        else:
            s = "3 finder"
            lf = loose[iev]
            ltruth = TruthJoin(lf)
            pids = lf.get("MuonPID")
            hit = None
            for ic, c in enumerate(lf.get("MuonCandidates")):
                if any(ltruth.track_truth(t)[0] == mu_idx
                       for t in c.getTracks()):
                    hit = list(pids[ic].getParameters())
                    break
            if hit is None:
                # The muon's track never reached the MTC (type bit 4) or a
                # longer track took its clusters first.
                wrong_track += 1
                cut["track does not reach MTC / clusters taken"] += 1
            else:
                n_mtc, span, frac, mean_mip = hit[0], hit[1], hit[2], hit[3]
                loose_vals.append((n_mtc, frac, mean_mip))
                if n_mtc < min_layers:
                    cut[f"MinMTCLayers ({min_layers})"] += 1
                elif frac < min_frac:
                    cut[f"MinFracLayers ({min_frac})"] += 1
                elif mean_mip > max_mean_mip:
                    cut[f"MaxMeanMIP ({max_mean_mip})"] += 1
                else:
                    cut["passes all cuts loose (lost to claim order)"] += 1
        stage[s] += 1
        stage_by_p[pbin][s] += 1

    tot = sum(stage.values())
    print(f"nominal {args[0]}")
    print(f"loose   {args[1]}")
    print(f"nu_mu CC events: {tot}  (MTC layers needed: {min_layers})")
    print()
    hdr = "".join(f"{f'p {lo}-{hi if hi < 1e8 else 'inf'}':>12}"
                  for lo, hi in P_BINS)
    print(f"{'stage':<16}{'all':>12}{hdr}")
    for s in sorted(stage):
        row = f"{s:<16}{stage[s] / tot:>12.3f}"
        for b in P_BINS:
            nb = sum(stage_by_p[b].values())
            row += f"{(stage_by_p[b][s] / nb if nb else 0):>12.3f}"
        print(row)
    print(f"{'N':<16}{tot:>12}" + "".join(
        f"{sum(stage_by_p[b].values()):>12}" for b in P_BINS))
    print()
    print("finder losses, first cut that fails:")
    for k, v in cut.most_common():
        print(f"  {k:<46} {v:>5}")
    if loose_vals:
        med = [sorted(v)[len(v) // 2] for v in zip(*loose_vals)]
        print(f"  median (nLayersMTC, frac, meanMIP) of those: "
              f"({med[0]:.0f}, {med[1]:.2f}, {med[2]:.2f})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
