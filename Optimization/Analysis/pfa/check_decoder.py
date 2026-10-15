#!/usr/bin/env python3
"""Step 0.4 check: SNDCellDecoder agrees with an independent Python decode.

Runs over the dump that `job6_pfa.py` produces with
`SNDCaloHitBuilder.DumpDecodedCells = True` (lines tagged
[SNDCaloHitBuilder][cell]) and re-decodes every cellID here, from the same
descriptors in SND_compact.xml but with a separate implementation, comparing
det / station / layer / plane / view / ordinal / channel.

Usage: check_decoder.py <job6.log>
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]
                       / "simulation" / "geometry"))
from parse_geometry import SNDGeometry

LINE = re.compile(
    r"\[SNDCaloHitBuilder\]\[cell\] coll=(?P<coll>\S+) cellID=(?P<cell>\d+) "
    r"det=(?P<det>-?\d+) station=(?P<station>-?\d+) layer=(?P<layer>-?\d+) "
    r"plane=(?P<plane>-?\d+) view=(?P<view>-?\d+) ordinal=(?P<ordinal>-?\d+) "
    r"channel=(?P<channel>-?\d+)")

READOUT_OF = {
    "SiTargetHitsWindowed": "SiTargetHits",
    "SiPadHitsWindowed": "SiPadHits",
    "MTCSciFiHitsWindowed": "MTCDetHits",
    "MTCScintHitsWindowed": "MTCDetHits",
}
DET_OF = {"SiTargetHitsWindowed": 0, "SiPadHitsWindowed": 1,
          "MTCSciFiHitsWindowed": 2, "MTCScintHitsWindowed": 2}


def parse_descriptor(desc):
    """'system:8,layer:8,...' -> {name: (offset, width)}; widths are unsigned
    here because none of the fields SND-PFA uses is signed."""
    out = {}
    off = 0
    for part in desc.split(","):
        name, width = part.split(":")
        out[name] = (off, int(width))
        off += int(width)
    return out


def get(fields, name, cell):
    if name not in fields:
        return -1
    off, width = fields[name]
    return (cell >> off) & ((1 << width) - 1)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    geo = SNDGeometry()
    n_mtc_layers = geo.mtc_n_layers
    fields = {coll: parse_descriptor(geo.bitfields[ro])
              for coll, ro in READOUT_OF.items()}

    n = 0
    bad = 0
    per_coll = {}
    for line in open(argv[1], errors="replace"):
        m = LINE.search(line)
        if not m:
            continue
        n += 1
        coll = m.group("coll")
        cell = int(m.group("cell"))
        f = fields[coll]
        det = DET_OF[coll]
        layer = get(f, "layer", cell)
        plane = get(f, "plane", cell) if "plane" in f else -1
        station = get(f, "station", cell) if "station" in f else -1
        if det == 0:
            view = 1 if plane == 1 else 0
            ordinal = layer
        elif det == 1:
            view = 2
            ordinal = layer
        else:
            view = 0 if plane == 0 else (1 if plane == 1 else 2)
            ordinal = station * n_mtc_layers + layer
        strip = get(f, "strip", cell) if "strip" in f else -1
        if det == 1 or (det == 2 and plane == 2):
            channel = get(f, "x", cell) * 4096 + get(f, "y", cell)
        else:
            channel = strip

        got = (int(m.group("det")), int(m.group("station")),
               int(m.group("layer")), int(m.group("plane")),
               int(m.group("view")), int(m.group("ordinal")),
               int(m.group("channel")))
        exp = (det, station, layer, plane, view, ordinal, channel)
        c = per_coll.setdefault(coll, [0, 0])
        c[0] += 1
        if got != exp:
            bad += 1
            c[1] += 1
            if bad <= 10:
                print(f"MISMATCH {coll} cell={cell}\n  algo   {got}\n  python {exp}")

    print(f"cells compared {n}")
    for coll, (tot, nbad) in sorted(per_coll.items()):
        print(f"  {coll:<24} {tot:>6} compared, {nbad:>4} mismatches")
    print(f"mismatches {bad}")
    ok = n > 0 and bad == 0
    print("OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
