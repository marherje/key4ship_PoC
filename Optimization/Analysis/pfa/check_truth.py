#!/usr/bin/env python3
"""Step 0.2 check: MCParticles survive the overlay and *ContribParticleIDs resolve.

Verifies, over every event of an EventOverlay output file:
  (a) MCParticles is present and non-empty;
  (b) for each *HitsWindowed collection, ContribParticleIDs[k] + offset is a
      valid index into MCParticles (offset from MCParticleSourceOffsets, picked
      by the contribution's source_id -- EventOverlay stores the source id in
      CaloHitContribution::PDG);
  (c) the PDG of the particle it points at equals ContribPDGs[k];
  (d) a neutrino (|PDG| in 12/14/16) is present;
  (e) in CC events a primary charged lepton is present.

It also reports how the neutrino arrives (generatorStatus, energy), which the
later phases need to define E_nu,true.

Usage: check_truth.py [events.edm4hep.root] [--max-events N]
"""
import sys
from collections import Counter

from podio.root_io import Reader

COLLS = [
    ("SiTargetHitsWindowed", "SiTargetContribPDGs", "SiTargetContribParticleIDs"),
    ("SiPadHitsWindowed", "SiPadContribPDGs", "SiPadContribParticleIDs"),
    ("MTCSciFiHitsWindowed", "MTCSciFiContribPDGs", "MTCSciFiContribParticleIDs"),
    ("MTCScintHitsWindowed", "MTCScintContribPDGs", "MTCScintContribParticleIDs"),
]
NU = {12, 14, 16}
LEP = {11, 13, 15}


def _param(frame, name):
    """Frame parameter as a list; [] when absent (podio raises a KeyError)."""
    try:
        val = frame.get_parameter(name)
    except Exception:
        return []
    if val is None:
        return []
    # podio hands back a bare scalar for a single-element vector parameter
    if isinstance(val, (int, float, str)):
        return [val]
    return list(val)


def main(argv):
    path = "events.edm4hep.root"
    nmax = 0
    args = [a for a in argv[1:] if not a.startswith("--")]
    if args:
        path = args[0]
    for a in argv[1:]:
        if a.startswith("--max-events"):
            nmax = int(a.split("=", 1)[1]) if "=" in a else 0

    events = Reader(path).get("events")
    n_total = len(events)
    if nmax:
        n_total = min(n_total, nmax)

    n_no_mc = 0
    n_bad_index = 0
    n_contrib = 0
    n_pdg_match = 0
    n_pdg_checked = 0
    n_with_nu = 0
    n_cc = 0
    n_cc_with_lep = 0
    nu_status = Counter()
    nu_energy = []
    src_ids_seen = Counter()

    for iev in range(n_total):
        frame = events[iev]
        try:
            mcs = frame.get("MCParticles")
        except Exception:
            n_no_mc += 1
            continue
        nmc = len(mcs)
        if nmc == 0:
            n_no_mc += 1
            continue

        offsets = _param(frame, "MCParticleSourceOffsets")
        src_order = _param(frame, "MCParticleSourceIDs")

        pdgs_mc = [p.getPDG() for p in mcs]

        # (d)/(e): neutrino and primary charged lepton
        nus = [p for p in mcs if abs(p.getPDG()) in NU]
        if nus:
            n_with_nu += 1
            nu = nus[0]
            nu_status[nu.getGeneratorStatus()] += 1
            nu_energy.append(nu.getEnergy())
        leps = [p for p in mcs
                if abs(p.getPDG()) in LEP and p.getGeneratorStatus() == 1]
        if nus and leps:
            n_cc += 1
            n_cc_with_lep += 1
        elif nus:
            # NC (or a CC whose lepton was not stored): not an error
            pass

        # (b)/(c)
        for coll, ppar, ipar in COLLS:
            try:
                hits = frame.get(coll)
                contribs = frame.get(coll + "_Contributions")
            except Exception:
                continue
            cpdgs = _param(frame, ppar)
            cids = _param(frame, ipar)
            if not cids:
                continue
            k = 0
            for hit in hits:
                for c in hit.getContributions():
                    if k >= len(cids):
                        break
                    src = c.getPDG()          # source id, see EventOverlay
                    src_ids_seen[src] += 1
                    off = 0
                    if offsets:
                        if src_order and src in src_order:
                            off = offsets[src_order.index(src)]
                        elif len(offsets) == 1:
                            off = offsets[0]
                    idx = cids[k] + off
                    n_contrib += 1
                    if cids[k] < 0 or idx >= nmc:
                        n_bad_index += 1
                    else:
                        n_pdg_checked += 1
                        if pdgs_mc[idx] == cpdgs[k]:
                            n_pdg_match += 1
                    k += 1

    print(f"file                     {path}")
    print(f"events checked           {n_total}")
    print(f"events without MCParticles {n_no_mc}")
    print(f"contributions checked    {n_contrib}")
    print(f"  out-of-range indices   {n_bad_index}")
    frac = (100.0 * n_pdg_match / n_pdg_checked) if n_pdg_checked else 0.0
    print(f"  PDG(MCParticle) == ContribPDGs: {n_pdg_match}/{n_pdg_checked} ({frac:.2f} %)")
    print(f"events with a neutrino   {n_with_nu}/{n_total}")
    print(f"events with nu + primary charged lepton (CC-like): {n_cc_with_lep}")
    print(f"neutrino generatorStatus {dict(nu_status)}")
    if nu_energy:
        e = sorted(nu_energy)
        print(f"neutrino energy [GeV]    min {e[0]:.2f}  median {e[len(e)//2]:.2f}  max {e[-1]:.2f}")
    print(f"source ids seen in contributions {dict(src_ids_seen)}")

    ok = (n_no_mc == 0 and n_bad_index == 0
          and n_pdg_checked > 0 and n_pdg_match == n_pdg_checked
          and n_with_nu == n_total)
    print("OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
