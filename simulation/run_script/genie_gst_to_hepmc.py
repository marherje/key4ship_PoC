#!/usr/bin/env python3
"""Convert a GENIE "gst" summary tree into a HepMC file ddsim can read.

GENIE's ``gevgen_fnal`` (as used for the SND neutrino samples, see
``launch_genie.sh``) is run once against the detector geometry
(``SND_compact_key4ship.gdml``, a GDML export of the same geometry as
``simulation/geometry/SND_compact.xml``) and writes, among other formats, a
flat "gst" ROOT TTree (one row per neutrino interaction) via GENIE's
``gntpc -f gst`` summary converter. Each row already carries:

  * the interaction vertex (``vtxx``/``vtxy``/``vtxz``, GENIE's native
    METERS -- NOT ROOT/TGeo's cm, see the conversion note below) in the
    reference frame of the GDML GENIE was run against. For the 2026_07_30
    EOS production that is the key4ship frame as-is (``--z-shift-m 0``).
    For samples made by ``launch_genie_gen.sh`` the GDML is shifted to the
    FairShip frame so the flux rays land where they should (see
    ``simulation/geometry/gdml/shift_gdml_frame.py``), and the vertices
    have to be brought back with ``--z-shift-m -25.76``;
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
import math
import sys

import pyhepmc


def on_shell(pdg, px, py, pz, E, pdg_db):
    """(px, py, pz, E) with E recomputed from the PDG mass.

    GENIE writes some final-state nucleons OFF-SHELL: a bound nucleon keeps
    its removal energy, so a neutron can come out with m = 935.6 MeV instead
    of 939.6. Geant4 derives the mass from (E, p) and, for one with p ~ 0,
    ends up with a neutron of zero kinetic energy and the wrong mass, tries
    to decay it, finds no channel and aborts the run (G4Exception DECAY003,
    fatal -- seen on trial t00000 of the first snd_genie_gen test). With p > 0
    the same mismatch is only a warning (part0005), which is why the shared
    2000-event sample never tripped it. Direction and |p| are kept; only E
    moves, by a few MeV. Nuclei and anything TDatabasePDG does not know keep
    GENIE's energy.
    """
    # PyROOT hands back a null proxy (falsy, and a ReferenceError on use)
    # rather than None for a code the database does not know, e.g. nuclei.
    part = pdg_db.GetParticle(int(pdg))
    if not part:
        return px, py, pz, E
    m = part.Mass()
    return px, py, pz, math.sqrt(px * px + py * py + pz * pz + m * m)


def convert(input_path, tree_name, start, nevents, output_path,
            cc_only=False, flavour=None, min_final_state_energy=0.0,
            z_shift_m=0.0, world_half_size_m=20.0):
    import ROOT

    pdg_db = ROOT.TDatabasePDG.Instance()

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
    outside = 0
    with pyhepmc.open(output_path, "w") as writer:
        for i in range(start, stop):
            t.GetEntry(i)

            if cc_only and not t.cc:
                continue
            if flavour is not None and t.neu != flavour:
                continue

            # GENIE also interacts in the air of the GDML world, which is
            # bigger than the DD4hep one once shifted (+-40 m vs +-20 m). A
            # primary vertex outside Geant4's world is a segfault in ddsim
            # (seen at z = -20.2 m), so those are dropped here and counted.
            vz = t.vtxz + z_shift_m
            if max(abs(t.vtxx), abs(t.vtxy), abs(vz)) >= world_half_size_m:
                outside += 1
                continue

            evt = pyhepmc.GenEvent(
                momentum_unit=pyhepmc.Units.MomentumUnit.GEV,
                length_unit=pyhepmc.Units.LengthUnit.MM,
            )
            evt.event_number = written

            # m -> mm. GENIE's GeomAnalyzerI reports the interaction vertex
            # in GENIE's own natural length unit (meter -- see
            # genie::units::meter), NOT in ROOT/TGeo's native cm, even
            # though the geometry it navigated (SND_compact_key4ship.gdml)
            # is a GDML/TGeo file. Confirmed against the raw gst tree: vtxz
            # spans [-20, +20] and vtxx/vtxy span roughly [-10, +16] /
            # [-8, +8], which matches the DD4hep world volume's
            # world_side = 20*m (SND_compact.xml) -- i.e. these are metres,
            # not centimetres. Using *10.0 here (cm -> mm) previously
            # shrank every vertex by a further 10x on top of that, squeezing
            # the whole ~20 m detector into a ~20 cm cube around the origin
            # (a ddsim-external unit bug -- not filtered/cut logic).
            # GENIE vertex time (vtxt, s) is not used downstream (the SND
            # pipelines are event-by-event, no absolute timing), so it is
            # dropped rather than converted.
            # z_shift_m undoes the frame the GDML was exported in (0 for a
            # key4ship-frame GDML, -25.76 for the FairShip-frame one).
            vtx = pyhepmc.GenVertex((t.vtxx * 1000.0, t.vtxy * 1000.0,
                                     vz * 1000.0, 0.0))

            # Incoming neutrino: status 4 ("beam"), no outgoing edge from
            # this vertex on its own -> Geant4's HepMC primary generator
            # (which only tracks status-1 particles) correctly ignores it.
            nu_in = pyhepmc.GenParticle((t.pxv, t.pyv, t.pzv, t.Ev), int(t.neu), 4)
            vtx.add_particle_in(nu_in)

            n_final = 0
            # Outgoing primary lepton.
            if t.El > min_final_state_energy:
                lep = pyhepmc.GenParticle(
                    on_shell(t.fspl, t.pxl, t.pyl, t.pzl, t.El, pdg_db), int(t.fspl), 1)
                vtx.add_particle_out(lep)
                n_final += 1

            # Post-FSI hadronic/nuclear final state.
            pdgf = t.pdgf
            pxf, pyf, pzf, Ef = t.pxf, t.pyf, t.pzf, t.Ef
            for k in range(t.nf):
                if Ef[k] <= min_final_state_energy:
                    continue
                p = pyhepmc.GenParticle(
                    on_shell(pdgf[k], pxf[k], pyf[k], pzf[k], Ef[k], pdg_db),
                    int(pdgf[k]), 1)
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
          "(%d skipped: no trackable final state; %d dropped: vertex outside "
          "the +-%g m world) -> %s"
          % (start, stop, total, written, skipped, outside, world_half_size_m,
             output_path))
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
    p.add_argument("--z-shift-m", type=float, default=0.0,
                    help="add this to every vertex z [m] before writing: 0 for a "
                         "gst made against the key4ship-frame GDML, -25.76 for one "
                         "made against SND_compact_fairship_frame.gdml (default: 0)")
    p.add_argument("--world-half-size-m", type=float, default=20.0,
                    help="drop events whose (shifted) vertex lies outside this "
                         "half-size cube: the DD4hep world, world_side in "
                         "SND_compact.xml (default: 20)")
    args = p.parse_args()

    n = convert(args.input, args.tree, args.start, args.nevents, args.output,
                cc_only=args.cc_only, flavour=args.flavour,
                min_final_state_energy=args.min_final_state_energy,
                z_shift_m=args.z_shift_m,
                world_half_size_m=args.world_half_size_m)
    if n == 0:
        print("ERROR: no events written", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
