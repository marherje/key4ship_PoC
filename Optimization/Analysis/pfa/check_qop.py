#!/usr/bin/env python3
"""Step 0.3 check: the q/p written at AtLastHit measures the muon momentum.

For every GlobalTrack that reaches the MTC (type & 4), finds the dominant MC
particle of its hits; when that is the event's primary muon, compares
p_reco = 1/|q/p| (from the AtLastHit state, where omega carries q/p in 1/GeV)
against the truth, in bins of the muon's momentum at the vertex.

WHAT THE COMPARISON CAN AND CANNOT SHOW. AtLastHit sits deep inside the MTC,
after the muon has crossed the SiTarget tungsten, the SiPad and part of the
iron -- of order 4 GeV of ionisation loss. So p_reco is BELOW the vertex
momentum by construction and ABOVE the momentum the muon leaves the detector
with (MCParticle::momentumAtEndpoint). Asking it to match the vertex momentum
would be asking the fit to undo an energy loss nobody has corrected yet; that
correction belongs to SNDMuonFinder, which gets the material budget from phase
1. What this script therefore checks is that p_reco is BRACKETED by the two
truth momenta, and that the charge sign is right.

Acceptance (step 0.3, as amended on 2026-09-22 -- see
docs/snd_pfa_progress.md):
  * > 80 % correct charge sign for p_true > 5 GeV;
  * the median relative difference to p_endpoint within 15 % for p > 10 GeV,
    where the remaining energy loss past the last hit is small.

The fraction of tracks bracketed by (p_endpoint, p_vertex) is printed but is
NOT a criterion: measured sigma68 against the endpoint is 0.39-0.98, i.e. the
MTC spectrometer resolution is dominated by multiple scattering in 2.2 m of
iron, so a large share of tracks falls outside the bracket by resolution
alone. That is why the muon momentum is taken from the iron RANGE in phase 3
(decision of 2026-09-22), with curvature used for the charge and for muons
that leave the detector.

Usage: check_qop.py [tracks_seq.edm4hep.root] [--max-events N]
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from podio.root_io import Reader
from pfa_truth import TruthJoin

AT_LAST_HIT = 3
BINS = [(2, 5), (5, 10), (10, 20), (20, 1e9)]


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def sigma68(xs):
    """Half-width of the central 68 % interval."""
    if len(xs) < 3:
        return float("nan")
    s = sorted(xs)
    lo = s[int(0.16 * (len(s) - 1))]
    hi = s[int(0.84 * (len(s) - 1))]
    return 0.5 * (hi - lo)


def main(argv):
    path = argv[1] if len(argv) > 1 and not argv[1].startswith("--") \
        else "tracks_seq.edm4hep.root"
    nmax = 0
    for a in argv[1:]:
        if a.startswith("--max-events="):
            nmax = int(a.split("=", 1)[1])

    events = Reader(path).get("events")
    n_total = len(events)
    if nmax:
        n_total = min(n_total, nmax)

    rows = []            # (p_vtx, p_reco, sign_ok, rel_var, p_end)
    n_mtc_tracks = 0
    n_muon_tracks = 0
    n_no_state = 0
    n_unmeasured = 0     # q/p left at the seed value (no field seen)

    for iev in range(n_total):
        frame = events[iev]
        truth = TruthJoin(frame)
        muon = truth.primary_muon()
        mu_idx = truth.index_of(muon)
        try:
            tracks = frame.get("GlobalTracks")
        except Exception:
            continue
        for trk in tracks:
            if not (trk.getType() & 4):
                continue
            n_mtc_tracks += 1
            if mu_idx < 0:
                continue
            idx, purity, _ = truth.track_truth(trk)
            if idx != mu_idx:
                continue
            n_muon_tracks += 1
            state = None
            for st in trk.getTrackStates():
                if st.location == AT_LAST_HIT:
                    state = st
                    break
            if state is None:
                n_no_state += 1
                continue
            qop = state.omega
            if qop == 0.0:
                n_no_state += 1
                continue
            # covMatrix index 5 = (omega, omega) in the lower-triangular
            # 6x6 layout (d0, phi, omega, z0, tanLambda, time).
            var_qop = state.covMatrix[5]
            rel_var = var_qop / (qop * qop) if qop else float("inf")
            p_reco = 1.0 / abs(qop)
            pm = muon.getMomentum()
            pe = muon.getMomentumAtEndpoint()
            p_true = math.sqrt(pm.x ** 2 + pm.y ** 2 + pm.z ** 2)
            p_end = math.sqrt(pe.x ** 2 + pe.y ** 2 + pe.z ** 2)
            if p_true <= 0:
                continue
            sign_ok = (qop < 0) == (muon.getCharge() < 0)
            if abs(rel_var - 0.36) < 1e-3:   # seed value, never updated
                n_unmeasured += 1
            rows.append((p_true, p_reco, sign_ok, rel_var, p_end))
        del tracks, frame

    print(f"file                      {path}")
    print(f"events                    {n_total}")
    print(f"GlobalTracks with MTC     {n_mtc_tracks}")
    print(f"  dominated by primary mu {n_muon_tracks}")
    print(f"  usable AtLastHit q/p    {len(rows)}   (no state/zero q/p: {n_no_state})")
    print()
    print(f"{'p_vtx bin':>13} {'N':>5} {'dp/p vs vtx':>12} {'dp/p vs end':>12} "
          f"{'sigma68(end)':>13} {'sign ok':>8} {'bracketed':>10} {'relVar':>8}")
    ok_end = True
    ok_sign = True
    for lo, hi in BINS:
        sel = [r for r in rows if lo <= r[0] < hi]
        if not sel:
            continue
        d_vtx = median([(r[1] - r[0]) / r[0] for r in sel])
        with_end = [r for r in sel if r[4] > 0.05]
        d_end = median([(r[1] - r[4]) / r[4] for r in with_end]) if with_end \
            else float("nan")
        sig_end = sigma68([(r[1] - r[4]) / r[4] for r in with_end]) \
            if with_end else float("nan")
        frac_sign = sum(1 for r in sel if r[2]) / len(sel)
        # p_reco should sit between the endpoint and the vertex momentum,
        # allowing 10 % on each side for the fit resolution.
        frac_brk = sum(1 for r in sel
                       if 0.9 * r[4] <= r[1] <= 1.1 * r[0]) / len(sel)
        relv = median([r[3] for r in sel])
        label = f"{lo:>6}-{'inf' if hi > 1e8 else hi:<6}"
        print(f"{label} {len(sel):>5} {d_vtx:>12.3f} {d_end:>12.3f} "
              f"{sig_end:>13.3f} {frac_sign:>8.2f} {frac_brk:>10.2f} "
              f"{relv:>8.4f}")
        if lo >= 5:
            ok_sign &= frac_sign > 0.80
        if lo >= 10:
            ok_end &= abs(d_end) < 0.15
    print()
    print(f"tracks whose q/p stayed at the seed value: {n_unmeasured}")
    print(f"charge sign > 0.80 (p>5 GeV):        {'pass' if ok_sign else 'FAIL'}")
    print("bracketed fraction: informative only, see the module docstring")
    print(f"|dp/p vs endpoint| < 0.15 (p>10 GeV): {'pass' if ok_end else 'FAIL'}")
    ok = bool(rows) and ok_sign and ok_end
    print("OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
