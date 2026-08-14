#!/usr/bin/env python3
"""Convert a GENIE "gst" summary tree into a HepMC file ddsim can read.

GENIE's ``gevgen_fnal`` (as used for the SND neutrino samples, see
``launch_genie.sh``) is run once against the detector geometry
(``SND_compact_key4ship.gdml``, a GDML export of the same geometry as
``simulation/geometry/SND_compact.xml``) and writes, among other formats, a
flat "gst" ROOT TTree (one row per neutrino interaction) via GENIE's
``gntpc -f gst`` summary converter. Each row already carries:

  * the interaction vertex (``vtxx``/``vtxy``/``vtxz``, cm) in the SAME
    reference frame as our own DD4hep geometry -- no FairShip<->key4ship
    offset needs to be re-applied here, that conversion only concerned the
    neutrino *flux* geometry GENIE was run against, not this output frame;
  * the outgoing primary lepton (``fspl``, ``pxl``/``pyl``/``pzl``/``El``);
  * the post-FSI hadronic/nuclear final state (``nf`` particles,
    ``pdgf``/``pxf``/``pyf``/``pzf``/``Ef``).

These are exactly the "stuff Geant4 needs to propagate through the
detector" -- the incoming neutrino itself is not tracked. This script turns
each gst row into one HepMC GenEvent: a single vertex at the GENIE-reported
position, with the incoming neutrino as an (untracked) beam particle and the
lepton + hadronic final state as status-1 outgoing particles, and writes the
whole batch to a HepMC ascii file. ddsim picks up the right reader
(HepMC3FileReader) automatically from the ``.hepmc`` extension -- see
``genie_base.py``.

Usage:
    ./genie_gst_to_hepmc.py --input gntp.0.gst.root --start 0 --nevents 500 \\
        --output data/genie_allflavours_0_500.hepmc

    # restrict to CC events only, or to one neutrino flavour:
    ./genie_gst_to_hepmc.py --cc-only --flavour 14 --nevents 200 -o numu_cc.hepmc
"""

import argparse
import sys

import pyhepmc


def convert(input_path, tree_name, start, nevents, output_path,
            cc_only=False, flavour=None, min_final_state_energy=0.0):
    import ROOT

    f = ROOT.TFile.Open(input_path)
    if not f or f.IsZombie():
        raise RuntimeError("Cannot open GENIE file: " + input_path)
    t = f.Get(tree_name)
    if not t:
        raise RuntimeError("Tree '%s' not found in %s" % (tree_name, input_path))

    total = t.GetEntries()
    if start >= total:
        raise RuntimeError("--start %d >= %d entries in tree" % (start, total))
    stop = total if nevents < 0 else min(total, start + nevents)

    written = 0
    skipped = 0
    with pyhepmc.open(output_path, "w") as writer:
        for i in range(start, stop):
            t.GetEntry(i)

            if cc_only and not t.cc:
                continue
            if flavour is not None and t.neu != flavour:
                continue

            evt = pyhepmc.GenEvent(
                momentum_unit=pyhepmc.Units.MomentumUnit.GEV,
                length_unit=pyhepmc.Units.LengthUnit.MM,
            )
            evt.event_number = written

            # cm -> mm; GENIE vertex time (vtxt, s) is not used downstream
            # (the SND pipelines are event-by-event, no absolute timing), so
            # it is dropped rather than converted.
            vtx = pyhepmc.GenVertex((t.vtxx * 10.0, t.vtxy * 10.0, t.vtxz * 10.0, 0.0))

            # Incoming neutrino: status 4 ("beam"), no outgoing edge from
            # this vertex on its own -> Geant4's HepMC primary generator
            # (which only tracks status-1 particles) correctly ignores it.
            nu_in = pyhepmc.GenParticle((t.pxv, t.pyv, t.pzv, t.Ev), int(t.neu), 4)
            vtx.add_particle_in(nu_in)

            n_final = 0
            # Outgoing primary lepton.
            if t.El > min_final_state_energy:
                lep = pyhepmc.GenParticle((t.pxl, t.pyl, t.pzl, t.El), int(t.fspl), 1)
                vtx.add_particle_out(lep)
                n_final += 1

            # Post-FSI hadronic/nuclear final state.
            pdgf = t.pdgf
            pxf, pyf, pzf, Ef = t.pxf, t.pyf, t.pzf, t.Ef
            for k in range(t.nf):
                if Ef[k] <= min_final_state_energy:
                    continue
                p = pyhepmc.GenParticle((pxf[k], pyf[k], pzf[k], Ef[k]), int(pdgf[k]), 1)
                vtx.add_particle_out(p)
                n_final += 1

            if n_final == 0:
                # Nothing left to propagate (can happen for very soft NC
                # events after the min-energy cut) -- skip, ddsim has
                # nothing to simulate for it anyway.
                skipped += 1
                continue

            evt.add_vertex(vtx)
            writer.write(evt)
            written += 1

    print("[genie_gst_to_hepmc] read entries [%d, %d) of %d, wrote %d events "
          "(%d skipped: no trackable final state) -> %s"
          % (start, stop, total, written, skipped, output_path))
    return written


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i", required=True,
                    help="GENIE gst.root file (local path or xrootd URL)")
    p.add_argument("--tree", default="gst", help="tree name (default: gst)")
    p.add_argument("--start", type=int, default=0, help="first entry (default: 0)")
    p.add_argument("--nevents", "-n", type=int, default=-1,
                    help="number of entries to convert, -1 = all remaining (default: -1)")
    p.add_argument("--output", "-o", required=True, help="output .hepmc file")
    p.add_argument("--cc-only", action="store_true", help="keep only CC interactions")
    p.add_argument("--flavour", type=int, default=None,
                    help="keep only this incoming neutrino PDG code "
                         "(e.g. 12/-12 nue, 14/-14 numu, 16/-16 nutau)")
    p.add_argument("--min-final-state-energy", type=float, default=0.0,
                    help="drop final-state particles at or below this energy [GeV] "
                         "(default: 0, keep everything GENIE reports)")
    args = p.parse_args()

    n = convert(args.input, args.tree, args.start, args.nevents, args.output,
                cc_only=args.cc_only, flavour=args.flavour,
                min_final_state_energy=args.min_final_state_energy)
    if n == 0:
        print("ERROR: no events written", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
