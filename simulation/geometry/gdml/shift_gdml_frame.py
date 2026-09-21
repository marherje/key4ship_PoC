#!/usr/bin/env python3
"""
Translate a DD4hep-exported GDML along z, for GENIE.

Why this exists. ``gevgen_fnal`` reads the neutrino flux from a GSimple
ntuple whose ray positions are in the FairShip reference frame, in metres,
and it offers no way to translate them: ``-z`` only slides every ray along
its own direction (GSimpleNtpFlux::MoveToZ0), which leaves the line of
flight -- and therefore where it crosses the detector -- exactly where it
was. The 2026_07_30 production was made that way and its neutrino cone
opens from z ~ 0 of the key4ship frame, i.e. from the MTC40 front face
instead of from the SHiP target 25.76 m upstream.

So instead of moving the flux to the geometry we move the geometry to the
flux: this script rewrites the world-level placements of a GDML so that the
detector sits where it does in FairShip, and enlarges the world box so it
still contains it. The GENIE output then carries vertices in the FairShip
frame, and ``genie_gst_to_hepmc.py --z-shift`` takes them back to the
key4ship frame before ddsim sees them.

The offset between the two frames (from the EOS production's README):

    SiPad front face, key4ship frame:   z0 =   -37.0 cm
    SiPad front face, FairShip frame:   z1 =  2539.0 cm
    z_fairship = z_key4ship + 2576.0 cm

Only the physvols placed directly in the world volume are touched; the
detector internals are untouched, so the exported material budget is
exactly the one ddsim simulates.

    ./shift_gdml_frame.py SND_compact.gdml SND_compact_fairship_frame.gdml \\
        --dz-cm 2576.0
"""

import argparse
import xml.etree.ElementTree as ET


FAIRSHIP_DZ_CM = 2576.0


def shift(input_path, output_path, dz_cm, world_z_cm):
    # The GDML header declares the xsi namespace only for the schema
    # location attribute; without registering it ElementTree would rewrite
    # the prefix as ns0, which ROOT's TGDMLParse does not mind but which
    # makes diffs against the original unreadable.
    ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")
    tree = ET.parse(input_path)
    root = tree.getroot()

    define = root.find("define")
    solids = root.find("solids")
    structure = root.find("structure")
    world_name = root.find("setup/world").get("ref")
    world = next(v for v in structure.findall("volume") if v.get("name") == world_name)

    positions = {p.get("name"): p for p in define.findall("position")}
    unit_to_cm = {"mm": 0.1, "cm": 1.0, "m": 100.0}

    shifted = []
    for pv in world.findall("physvol"):
        ref = pv.find("positionref")
        if ref is not None:
            pos = positions[ref.get("ref")]
            scale = unit_to_cm[pos.get("unit", "mm")]
            z_old = float(pos.get("z", "0"))
            pos.set("z", repr(z_old + dz_cm / scale))
            shifted.append((pv.get("name"), z_old * scale, z_old * scale + dz_cm))
        else:
            # Placed at the origin implicitly: give it an explicit position.
            name = pv.get("name") + "_frame_shift_pos"
            pos = ET.SubElement(define, "position",
                                name=name, x="0", y="0", z=repr(dz_cm), unit="cm")
            ET.SubElement(pv, "positionref", ref=name)
            shifted.append((pv.get("name"), 0.0, dz_cm))

    box = next(s for s in solids if s.get("name") == world.find("solidref").get("ref"))
    scale = unit_to_cm[box.get("lunit", "mm")]
    z_old = float(box.get("z")) * scale
    if world_z_cm > z_old:
        box.set("z", repr(world_z_cm / scale))

    tree.write(output_path, xml_declaration=True, encoding="unicode")

    for name, z0, z1 in shifted:
        print(f"  {name:24s} z = {z0:9.3f} cm -> {z1:9.3f} cm")
    print(f"  world box z: {z_old:.0f} cm -> {max(z_old, world_z_cm):.0f} cm")
    print(f"Wrote {output_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--dz-cm", type=float, default=FAIRSHIP_DZ_CM,
                    help="z offset to ADD to every world-level placement "
                         f"[default: {FAIRSHIP_DZ_CM}, key4ship -> FairShip]")
    ap.add_argument("--world-z-cm", type=float, default=8000.0,
                    help="full z length of the world box after the shift; only "
                         "grown, never shrunk [default: 8000 = +-40 m, enough "
                         "for the FairShip frame's ~29 m]")
    args = ap.parse_args()
    shift(args.input, args.output, args.dz_cm, args.world_z_cm)


if __name__ == "__main__":
    main()
