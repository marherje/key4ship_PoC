#!/usr/bin/env python3
"""Step 1.2 check: the CalorimeterHits SNDCaloHitBuilder writes are sound.

Verifies, over the events of a job6 output file:
  * one link per written hit, each pointing at a sim hit with the same cellID
    (the link, not the index, is what later stages use to reach MC truth);
  * every hit is above the MIP threshold and carries its MIP value in
    energyError;
  * the reconstructed position comes from the segmentation, not from the sim
    hit -- checked on the MTC, where SND_SciFiAction overwrites the sim hit's y;
  * energy scales: prints the median weight actually applied (E_reco / E_dep)
    per detector and view, to compare against the phase-1.1 table.

Usage: check_calohits.py [pfa_seq.edm4hep.root] [--threshold 0.5]
"""
import sys
from collections import defaultdict

from podio.root_io import Reader

COLLS = [
    ("SiTargetCaloHits", "SiTargetCaloHitLinks", "SiTargetHitsWindowed"),
    ("SiPadCaloHits", "SiPadCaloHitLinks", "SiPadHitsWindowed"),
    ("MTCSciFiCaloHits", "MTCSciFiCaloHitLinks", "MTCSciFiHitsWindowed"),
    ("MTCScintCaloHits", "MTCScintCaloHitLinks", "MTCScintHitsWindowed"),
]


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main(argv):
    path = "pfa_seq.edm4hep.root"
    thr = 0.5
    for a in argv[1:]:
        if a.startswith("--threshold="):
            thr = float(a.split("=", 1)[1])
        elif not a.startswith("--"):
            path = a

    events = Reader(path).get("events")
    n_hits = 0
    n_links = 0
    bad_cell = 0
    below_thr = 0
    zero_mip = 0
    ratios = defaultdict(list)     # (det, view) -> E_reco/E_dep
    mtc_y_from_sim = 0
    mtc_checked = 0

    for iev in range(len(events)):
        frame = events[iev]
        for coll, linkColl, simColl in COLLS:
            try:
                hits = frame.get(coll)
                links = frame.get(linkColl)
            except Exception:
                continue
            n_hits += len(hits)
            n_links += len(links)
            for link in links:
                rec = link.getFrom()
                sim = link.getTo()
                if rec.getCellID() != sim.getCellID():
                    bad_cell += 1
                    continue
                mip = rec.getEnergyError()
                if mip <= 0:
                    zero_mip += 1
                elif mip < thr - 1e-6:
                    below_thr += 1
                edep = sim.getEnergy()
                if edep > 0:
                    t = rec.getType()
                    ratios[(t // 10, t % 10)].append(rec.getEnergy() / edep)
                if coll.startswith("MTC"):
                    mtc_checked += 1
                    # SND_SciFiAction writes an energy-weighted step y into the
                    # sim hit; a cell-centre y must differ from it in general.
                    if abs(rec.getPosition().y - sim.getPosition().y) < 1e-6:
                        mtc_y_from_sim += 1
            del hits, links
        del frame

    print(f"file                 {path}")
    print(f"events               {len(events)}")
    print(f"hits / links         {n_hits} / {n_links}")
    print(f"cellID mismatches    {bad_cell}")
    print(f"hits below threshold {below_thr}   (threshold {thr} MIP)")
    print(f"hits with MIP <= 0   {zero_mip}")
    print(f"MTC hits whose y equals the sim hit's: {mtc_y_from_sim}/{mtc_checked}")
    print()
    print(f"{'det':>4} {'view':>5} {'N':>8} {'median E_reco/E_dep':>20}")
    for (det, view), rs in sorted(ratios.items()):
        print(f"{det:>4} {view:>5} {len(rs):>8} {median(rs):>20.3f}")

    ok = (n_hits == n_links and bad_cell == 0 and below_thr == 0
          and zero_mip == 0 and n_hits > 0)
    print()
    print("OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
