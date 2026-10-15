#!/usr/bin/env python3
"""Phase 1.5: the ceiling on E_had, measured on GENIE events.

Before any algorithm tries to reconstruct the hadronic energy of a neutrino
interaction, this says how well it could possibly do. For every event it
compares

  E_had,sum   sum of the kinetic energy of the primary final-state hadrons
              (MCParticles, generatorStatus 1, not a lepton or a neutrino);
  E_had,nu    the energy transferred to the hadronic system, E_nu - E_lepton;
  E_vis,nomu  sum of the calibrated energy of every reconstructed hit whose
              dominant MC particle is NOT the primary muon.

TWO TRUTHS, because they are not the same number and only one of them is a
fair target. E_had,sum counts what the generator listed, and a good part of
that never reaches the calorimeter: on this sample 6 % of it is neutrons, 5.7 %
is K0 (PDG 311, which Geant4 does not even propagate as such), and more goes
into binding energy and into the neutrinos from pion and kaon decays.
E_had,nu is the physical quantity a neutrino experiment wants and does not
depend on how the generator chose to list the final state, so it is the one
the criterion uses.

The difference is what the calorimeter and the calibration lose before the
clustering, the corridor or the vertex finding get a chance to lose anything
more. Phase 4 cannot beat this number, and knowing it is the only way to tell
a bad E_had algorithm from a detector that never saw the energy.

Usage: ehad_bound.py [pfa.root] [--max-events N]
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from podio.root_io import Reader
from pfa_truth import TruthJoin

HIT_COLLS = ["SiTargetCaloHits", "SiPadCaloHits",
             "MTCSciFiCaloHits", "MTCScintCaloHits"]
BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 1e9)]
NU_PDGS = {12, 14, 16}
LEPTONS = {11, 13, 15}


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


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


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

    rows = []       # (E_had_sum, E_vis_nomu, E_had_nu, has_muon)
    n_no_truth = 0

    for iev in range(n_total):
        frame = events[iev]
        truth = TruthJoin(frame)
        mcs = truth.mcs
        if mcs is None or len(mcs) == 0:
            n_no_truth += 1
            del frame
            continue
        muon = truth.primary_muon()
        mu_idx = truth.index_of(muon)

        # Truth 1: kinetic energy of the primary hadrons, as listed.
        e_had_sum = 0.0
        for p in mcs:
            if p.getGeneratorStatus() != 1:
                continue
            pdg = abs(p.getPDG())
            if pdg in LEPTONS or pdg in NU_PDGS:
                continue
            e_had_sum += max(0.0, p.getEnergy() - p.getMass())

        # Truth 2: energy transferred to the hadronic system.
        nu = truth.neutrino()
        lep = truth.primary_lepton()
        e_had_nu = float("nan")
        if nu is not None:
            e_had_nu = nu.getEnergy() - (lep.getEnergy() if lep is not None else 0.0)

        e_nomu = 0.0
        e_all = 0.0
        for coll in HIT_COLLS:
            try:
                hits = frame.get(coll)
            except Exception:
                continue
            for h in hits:
                e = h.getEnergy()
                e_all += e
                idx, _ = truth.dominant(truth.cell_truth(h.getCellID()))
                if idx != mu_idx:
                    e_nomu += e
            del hits
        rows.append((e_had_sum, e_nomu, e_had_nu, mu_idx >= 0))
        del frame

    print(f"file   {path}")
    print(f"events {n_total}   (without usable truth: {n_no_truth})")
    print()
    for which, idx, label_t in ((1, 0, "E_had,sum (generator list)"),
                               (2, 2, "E_had,nu = E_nu - E_lepton")):
        print(f"--- against {label_t}")
        print(f"{'bin [GeV]':>14} {'N':>6} {'median ratio':>13} {'bias':>8} "
              f"{'sigma68':>9}")
        for lo, hi in BINS:
            sel = [r for r in rows
                   if r[idx] == r[idx] and lo <= r[idx] < hi and r[idx] > 0.05]
            if not sel:
                continue
            ratio = [r[1] / r[idx] for r in sel]
            rel = [(r[1] - r[idx]) / r[idx] for r in sel]
            lbl = f"{lo:>5}-{'inf' if hi > 1e8 else hi:<6}"
            print(f"{lbl:>14} {len(sel):>6} {median(ratio):>13.3f} "
                  f"{median(rel):>8.3f} {sigma68(rel):>9.3f}")
        good = [r for r in rows if r[idx] == r[idx] and r[idx] > 0.05]
        corr = pearson([r[idx] for r in good], [r[1] for r in good])
        bias = median([(r[1] - r[idx]) / r[idx] for r in good])
        print(f"  N={len(good)}  correlation={corr:.3f}  median bias={bias:+.3f}")
        print()
        if which == 2:
            corr_nu, bias_nu = corr, bias

    print(f"fraction of events with a primary muon = "
          f"{sum(1 for r in rows if r[3]) / len(rows):.3f}")
    ok = corr_nu > 0.9 and abs(bias_nu) < 0.10
    print()
    print(("OK" if ok else "FAIL") +
          "  (criterion: correlation > 0.90 and |bias| < 0.10)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
