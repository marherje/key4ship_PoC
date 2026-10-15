#!/usr/bin/env python3
"""Phase 3.3: how well SNDMuonFinder finds and measures the muon.

Measures, on the PFA output of the GENIE reference sample:

  1. efficiency -- the fraction of nu_mu CC events with a MuonCandidate whose
     corridor is dominated by the PRIMARY muon, split by how much of the muon
     there was to find (hits it left, whether it reached the MTC at all). The
     reference to compare against is not 1 but "the muon had a GlobalTrack":
     the finder cannot identify what the tracker never reconstructed;
  2. fakes -- NC and nu_e CC events with a candidate, and what MC particle
     actually built those corridors. A pion crossing the whole MTC is physics,
     not a bug, and the table says how often it happens;
  3. momentum -- bias and sigma68 of (p_reco - p_true)/p_true in bins of
     p_true, separately for the candidates measured by range and by curvature
     (MuonPID parameter 7 is 1 for range, 0 for curvature);
  4. charge -- the fraction with the right sign.

Usage: muon_perf.py [pfa.root] [--max-events N]
"""
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from podio.root_io import Reader
from pfa_truth import TruthJoin

BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 1e9)]
NU_PDGS = {12, 14, 16}


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

    cls = Counter()          # event class -> count
    matched = Counter()      # event class -> events with a matched muon
    withcand = Counter()     # event class -> events with any candidate
    had_track = 0            # nu_mu CC where the muon had a GlobalTrack
    had_track_matched = 0
    by_len = {"all": [0, 0], ">=20 hits": [0, 0], "reaches MTC": [0, 0]}
    rows = []                # (p_true, p_reco, sign_ok, by_range)
    fake_pdgs = Counter()

    for iev in range(n_total):
        frame = events[iev]
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
        elif abs(lep.getPDG()) == 13:
            klass = "numu CC"
        elif abs(lep.getPDG()) == 11:
            klass = "nue CC"
        else:
            klass = "nutau CC"
        cls[klass] += 1

        try:
            cands = frame.get("MuonCandidates")
            pids = frame.get("MuonPID")
            tracks = frame.get("GlobalTracks")
        except Exception:
            del frame
            continue

        if len(cands):
            withcand[klass] += 1

        # Did the tracker give the muon a global track at all?
        mu_has_track = False
        if mu_idx >= 0:
            for t in tracks:
                idx, _, _ = truth.track_truth(t)
                if idx == mu_idx:
                    mu_has_track = True
                    break

        # Which candidate, if any, is the primary muon?
        best = None
        for ic, c in enumerate(cands):
            acc = {}
            for cl in c.getClusters():
                for h in cl.getHits():
                    for idx, e in truth.cell_truth(h.getCellID()).items():
                        acc[idx] = acc.get(idx, 0.0) + e
            idx, pur = truth.dominant(acc)
            if klass in ("NC", "nue CC") and idx >= 0 and idx < len(mcs):
                fake_pdgs[mcs[idx].getPDG()] += 1
            if mu_idx >= 0 and idx == mu_idx:
                best = (ic, c, pur)

        if klass == "numu CC" and mu_idx >= 0:
            nhits_mu = sum(1 for cid in truth.by_cell
                           if truth.dominant(truth.by_cell[cid])[0] == mu_idx)
            by_len["all"][1] += 1
            if nhits_mu >= 20:
                by_len[">=20 hits"][1] += 1
            if mu_has_track:
                had_track += 1
            if best is not None:
                matched[klass] += 1
                by_len["all"][0] += 1
                if nhits_mu >= 20:
                    by_len[">=20 hits"][0] += 1
                if mu_has_track:
                    had_track_matched += 1

                # momentum
                pm = muon.getMomentum()
                p_true = math.sqrt(pm.x ** 2 + pm.y ** 2 + pm.z ** 2)
                pr = best[1].getMomentum()
                p_reco = math.sqrt(pr.x ** 2 + pr.y ** 2 + pr.z ** 2)
                by_range = True
                if best[0] < len(pids):
                    pars = list(pids[best[0]].getParameters())
                    if len(pars) > 7:
                        by_range = pars[7] > 0.5
                sign_ok = (best[1].getCharge() < 0) == (muon.getCharge() < 0)
                if p_true > 0:
                    rows.append((p_true, p_reco, sign_ok, by_range))
        del cands, pids, tracks, frame

    print(f"file   {path}")
    print(f"events {n_total}")
    print()
    print("event classes:", dict(cls))
    print()
    print("--- efficiency on numu CC")
    for k, (num, den) in by_len.items():
        if den:
            print(f"  {k:<12} {num:>5}/{den:<5} = {num/den:.3f}")
    if had_track:
        print(f"  muon had a GlobalTrack: {had_track}/{by_len['all'][1]} = "
              f"{had_track/by_len['all'][1]:.3f}")
        print(f"  identified among those:  {had_track_matched}/{had_track} = "
              f"{had_track_matched/had_track:.3f}   <- what the finder costs")
    print()
    print("--- fakes (events with a candidate)")
    for k in ("NC", "nue CC", "nutau CC"):
        if cls[k]:
            print(f"  {k:<10} {withcand[k]:>5}/{cls[k]:<5} = "
                  f"{withcand[k]/cls[k]:.3f}")
    if fake_pdgs:
        top = ", ".join(f"{pdg}:{n}" for pdg, n in fake_pdgs.most_common(6))
        print(f"  dominant MC particle of those corridors: {top}")
    print()
    print("--- momentum")
    print(f"{'p_true bin':>12} {'N':>5} {'bias':>8} {'sigma68':>9} "
          f"{'sign ok':>8} {'by range':>9}")
    for lo, hi in BINS:
        sel = [r for r in rows if lo <= r[0] < hi]
        if not sel:
            continue
        rel = [(r[1] - r[0]) / r[0] for r in sel]
        lbl = f"{lo:>4}-{'inf' if hi > 1e8 else hi:<5}"
        print(f"{lbl:>12} {len(sel):>5} {median(rel):>8.3f} "
              f"{sigma68(rel):>9.3f} "
              f"{sum(1 for r in sel if r[2])/len(sel):>8.2f} "
              f"{sum(1 for r in sel if r[3])/len(sel):>9.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
