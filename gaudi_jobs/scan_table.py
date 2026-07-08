#!/usr/bin/env python3
"""Diagnostic table for multi-PG geometry-scan reconstructions.

Usage (with the key4hep environment):
    python3 scan_table.py multi_pg_geoA_pipeline [multi_pg_geoB_pipeline ...]

Reads every tracks_<pair>_angle_<a>.edm4hep.root in the given pipeline
directories and prints, per sample: % of events with the expected track
multiplicity (mu_pi: 2, e_pi: 1), fraction of events with any track, median
chi2/ndf and the multiplicity distribution.
"""

import glob
import os
import re
import sys
from collections import Counter

from podio import root_io

EXPECTED = {"mu_pi": 2, "e_pi": 1}


def stats(path, nexp):
    counts, chi2s = [], []
    for frame in root_io.Reader(path).get("events"):
        tracks = frame.get("ACTSTracks")
        counts.append(len(tracks))
        for t in tracks:
            chi2s.append(t.getChi2() / max(1, t.getNdf()))
    n = len(counts)
    match = 100.0 * sum(1 for x in counts if x == nexp) / n
    anytrk = sum(1 for x in counts if x > 0)
    med = sorted(chi2s)[len(chi2s) // 2] if chi2s else float("nan")
    dist = " ".join(f"{k}:{v}" for k, v in sorted(Counter(counts).items()))
    return match, anytrk, n, med, dist


def main():
    dirs = sys.argv[1:]
    if not dirs:
        print(__doc__)
        return 1
    for d in dirs:
        print(f"== {d}")
        files = sorted(glob.glob(os.path.join(d, "tracks_*_angle_*.edm4hep.root")))
        for f in files:
            m = re.match(r"tracks_(\w+?)_angle_([\d.]+)\.edm4hep\.root",
                         os.path.basename(f))
            pair, angle = m.group(1), m.group(2)
            nexp = EXPECTED.get(pair, 1)
            match, anytrk, n, med, dist = stats(f, nexp)
            print(f"  {pair:6s} angle={angle:>4}: n=={nexp}: {match:5.1f}%  "
                  f"con_track={anytrk}/{n}  chi2/ndf med={med:.2f}  dist[{dist}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
