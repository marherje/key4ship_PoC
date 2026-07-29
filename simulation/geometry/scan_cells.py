#!/usr/bin/env python3
"""Raster-scan a SiPad sensitive plane and count the distinct (x,y) cells.

Shoots a regular grid of points across one SiPad layer, navigates the real
TGeo tree for each of them and, whenever the point lands inside a sensitive
pad array, asks the DD4hep readout segmentation for the cellID. The (x,y)
fields of those cellIDs are then counted.

This is a geometry *measurement*, not arithmetic: nothing here assumes how
many wafers there are or how the segmentation offsets are wired, so it also
catches two wafers colliding onto the same global cell, an offset that is off
by one pad, or dead area that should not be dead.

Usage (with the key4hep environment and the repo's install/ on the paths):
    python3 simulation/geometry/scan_cells.py [--layer N] [--step MM]
                                              [--compact FILE] [--map]

    --step   sampling pitch in mm (default 1.0, ~5 samples per 5.53 mm pad).
             Use a value well below the pad pitch or cells will be missed.
    --map    also print an ASCII occupancy map of the plane (one char per
             wafer-sized bin), handy to eyeball the dead cross between ASUs.
"""
import argparse
import sys

import ROOT
from dd4hep import Detector


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--compact", default="simulation/geometry/SND_compact.xml")
    ap.add_argument("--layer", type=int, default=0, help="SiPad layer to scan")
    ap.add_argument("--step", type=float, default=1.0, help="sampling pitch [mm]")
    ap.add_argument("--map", action="store_true", help="print occupancy map")
    args = ap.parse_args()

    desc = Detector.getInstance()
    desc.fromXML(args.compact)

    readout = desc.readout("SiPadHits")
    seg     = readout.segmentation()
    idspec  = readout.idSpec()
    decoder = idspec.decoder()
    mgr     = ROOT.gGeoManager

    # ---- locate the sensitive slice of the requested layer -----------------
    det   = desc.detector("SiPad")
    layer = det.child("SiPad_layer_%d" % args.layer)
    lvol  = layer.placement().GetVolume()

    slice_node = None
    for i in range(lvol.GetNdaughters()):
        nd = lvol.GetNode(i)
        if nd.GetVolume().GetNdaughters() and "_wafer" in \
                nd.GetVolume().GetNode(0).GetVolume().GetName():
            slice_node = nd
            break
    if slice_node is None:
        sys.exit("no tiled sensitive slice found in SiPad_layer_%d" % args.layer)

    slice_vol = slice_node.GetVolume()
    half_x = slice_vol.GetShape().GetDX()      # cm
    half_y = slice_vol.GetShape().GetDY()

    # global z of the slice: layer placement * slice placement
    z_glob = (layer.placement().GetMatrix().GetTranslation()[2]
              + det.placement().GetMatrix().GetTranslation()[2]
              + slice_node.GetMatrix().GetTranslation()[2])

    print("Scanning %s" % slice_vol.GetName())
    print("  plane   : %.2f x %.2f mm at z = %.3f mm"
          % (2 * half_x * 10, 2 * half_y * 10, z_glob * 10))
    print("  step    : %.3f mm" % args.step)

    # ---- raster ------------------------------------------------------------
    step = args.step / 10.0                     # mm -> cm
    n_x  = int(2 * half_x / step)
    n_y  = int(2 * half_y / step)
    print("  samples : %d x %d = %d" % (n_x, n_y, n_x * n_y))

    cells      = {}          # (x, y) -> set of wafer ids that produced it
    n_sens     = 0
    n_dead     = 0
    NBIN       = 54          # bins per side for --map (10 mm bins on a 540 mm plane)
    live       = {}          # map bin -> sensitive samples
    seen       = {}          # map bin -> total samples
    import array
    lpos = array.array("d", [0.0, 0.0, 0.0])
    gpos = array.array("d", [0.0, 0.0, 0.0])

    for iy in range(n_y):
        y = -half_y + (iy + 0.5) * step
        for ix in range(n_x):
            x = -half_x + (ix + 0.5) * step

            if args.map:
                bin_key = (min(int((x + half_x) / (2 * half_x) * NBIN), NBIN - 1),
                           min(int((y + half_y) / (2 * half_y) * NBIN), NBIN - 1))
                seen[bin_key] = seen.get(bin_key, 0) + 1

            node = mgr.FindNode(x, y, z_glob)
            if node is None or not mgr.GetCurrentVolume().GetName().endswith("_wafer_pads"):
                n_dead += 1
                continue
            n_sens += 1
            if args.map:
                live[bin_key] = live.get(bin_key, 0) + 1

            # volumeID: accumulate the physVolIDs along the current path
            ids = []
            for lvl in range(mgr.GetLevel() + 1):
                pv = ROOT.dd4hep.PlacedVolume(mgr.GetMother(mgr.GetLevel() - lvl)
                                              if lvl < mgr.GetLevel()
                                              else mgr.GetCurrentNode())
                for vid in pv.volIDs():
                    ids.append((vid.first, vid.second))
            vol_id = idspec.encode(ids)

            gpos[0], gpos[1], gpos[2] = x, y, z_glob
            mgr.MasterToLocal(gpos, lpos)
            cid = seg.cellID(
                ROOT.dd4hep.Position(lpos[0], lpos[1], lpos[2]),
                ROOT.dd4hep.Position(x, y, z_glob), vol_id)

            cx = decoder.get(cid, "x")
            cy = decoder.get(cid, "y")
            wf = decoder.get(cid, "wafer")
            cells.setdefault((cx, cy), set()).add(wf)

    # ---- exact active area, straight from the TGeo shapes -------------------
    # The raster fraction above is quantisation-biased (an 88.48 mm active band
    # sampled on a 1 mm grid yields 88 or 89 hits, never 88.48), so the real
    # number is summed from the placed pad volumes instead.
    area_pads = 0.0
    area_si   = 0.0
    for i in range(slice_vol.GetNdaughters()):
        wv = slice_vol.GetNode(i).GetVolume()
        area_si += 4 * wv.GetShape().GetDX() * wv.GetShape().GetDY()
        for j in range(wv.GetNdaughters()):
            pv_shape = wv.GetNode(j).GetVolume().GetShape()
            area_pads += 4 * pv_shape.GetDX() * pv_shape.GetDY()
    area_plane = 4 * half_x * half_y

    # ---- report ------------------------------------------------------------
    xs = [c[0] for c in cells]
    ys = [c[1] for c in cells]
    wafers = set()
    for w in cells.values():
        wafers |= w
    shared = {c: w for c, w in cells.items() if len(w) > 1}

    total = n_sens + n_dead
    print()
    print("  sensitive samples : %d / %d (%.2f%%, quantisation-biased)"
          % (n_sens, total, 100.0 * n_sens / total))
    print("  active area (TGeo): %.1f / %.1f mm2 = %.2f%%   [silicon incl. rim:"
          " %.2f%%]"
          % (area_pads * 100, area_plane * 100, 100.0 * area_pads / area_plane,
             100.0 * area_si / area_plane))
    print("  distinct wafers   : %d" % len(wafers))
    print("  DISTINCT (x,y)    : %d" % len(cells))
    print("  x index range     : %d .. %d" % (min(xs), max(xs)))
    print("  y index range     : %d .. %d" % (min(ys), max(ys)))
    span_x = max(xs) - min(xs) + 1
    span_y = max(ys) - min(ys) + 1
    print("  index grid        : %d x %d = %d (%s)"
          % (span_x, span_y, span_x * span_y,
             "complete" if span_x * span_y == len(cells) else "WITH HOLES"))
    if shared:
        print("  !! %d cells reached from more than one wafer (offset clash):"
              % len(shared))
        for c, w in list(shared.items())[:5]:
            print("       (x=%d,y=%d) <- wafers %s" % (c[0], c[1], sorted(w)))
    else:
        print("  cell/wafer clash  : none")

    if args.map:
        bin_mm = 2 * half_x * 10 / NBIN
        print("\n  active-area map, %.1f mm per char: '#' fully active,"
              " '+' partly dead, ' ' fully dead" % bin_mm)
        for by in range(NBIN - 1, -1, -1):
            row = ""
            for bx in range(NBIN):
                tot = seen.get((bx, by), 0)
                act = live.get((bx, by), 0)
                if tot == 0:
                    row += "?"
                elif act == tot:
                    row += "#"
                elif act == 0:
                    row += " "
                else:
                    row += "+"
            print("    |" + row + "|")

    return 0


if __name__ == "__main__":
    sys.exit(main())
