#!/usr/bin/env python3
"""First figure of merit of a variant: total number of SiPad hits.

Post-processes the ShipHits.root of one variant. The SiPad RNTuple holds one
row per hit, so the FOM is just its entry count. Written as a bare integer to
fom.txt, which is what the controller reads back.

Usage (needs `source init_key4ship.sh` for ROOT):
    python3 compute_fom.py [--input ShipHits.root] [--output fom.txt]
"""

import argparse
import sys
from pathlib import Path


def count_entries(filename, container):
    import ROOT
    return ROOT.RNTupleReader.Open(container, filename).GetNEntries()


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", default="ShipHits.root", help="RNTuple file")
    ap.add_argument("--output", default="fom.txt", help="where to write the FOM")
    ap.add_argument("--container", default="SiPad", help="RNTuple to count")
    args = ap.parse_args()

    path = Path(args.input)
    if not path.is_file():
        print("  !! ERROR: no existe %s" % path)
        return 1

    try:
        n = count_entries(str(path), args.container)
    except Exception as exc:                    # missing ntuple, corrupt file…
        print("  !! ERROR leyendo %s de %s: %s" % (args.container, path, exc))
        return 1

    Path(args.output).write_text("%d\n" % n)
    print("FOM (%s nhits): %d" % (args.container, n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
