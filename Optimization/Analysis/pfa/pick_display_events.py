#!/usr/bin/env python3
"""Rank events by how much there is to see in an event display.

An event is interesting when the reconstruction did something structurally
visible, not merely when it has many hits. The score therefore rewards, in
order:

  * spliced global tracks -- Track::type is a bitmask (1 = SiTarget,
    2 = SiPad, 4 = MTC), so a type with more than one bit set is a track the
    sequential tracker joined ACROSS detectors, which is the thing worth
    looking at;
  * local tracks surviving next to the globals, i.e. segments the splicer
    chose not to join: the two collections side by side show what the
    algorithm decided;
  * a muon leaving a long clean trail through the iron, with hadronic
    activity around the vertex to contrast it with;
  * several tracks, but not hundreds -- a shower that fragments into fifty
    segments is noise, not a picture.

With --clean the ranking is inverted: the aim is a picture that READS, not one
that is busy. It asks for one or two global tracks, at least one of them
spliced across detectors, few or no leftover local segments, a long track, and
a small hadronic shower around the vertex -- the event you put on a slide when
the audience has ten seconds to look at it.

Usage: pick_display_events.py <tracks.root> [--top N] [--max-events N] [--clean]
"""
import sys
from collections import Counter

from podio.root_io import Reader

LOCAL_COLLS = ["SiTargetTracks", "SiPadTracks", "MTCTracks"]
DET_NAME = {1: "SiT", 2: "SiP", 4: "MTC"}


def bits(mask):
    return [DET_NAME[b] for b in (1, 2, 4) if mask & b]


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    path = argv[1]
    top = 10
    nmax = 0
    clean = False
    for a in argv[2:]:
        if a.startswith("--top="):
            top = int(a.split("=", 1)[1])
        elif a.startswith("--max-events="):
            nmax = int(a.split("=", 1)[1])
        elif a == "--clean":
            clean = True

    events = Reader(path).get("events")
    n_total = len(events)
    if nmax:
        n_total = min(n_total, nmax)

    scored = []
    stats = Counter()
    for iev in range(n_total):
        frame = events[iev]
        try:
            globals_ = frame.get("GlobalTracks")
        except Exception:
            del frame
            continue

        n_spliced = 0
        spliced_kinds = []
        longest = 0
        n_global = len(globals_)
        for t in globals_:
            mask = t.getType()
            nstates = t.trackerHits_size()
            longest = max(longest, nstates)
            if bin(mask & 7).count("1") > 1:
                n_spliced += 1
                spliced_kinds.append("+".join(bits(mask)))

        n_local = 0
        local_by_det = {}
        for coll in LOCAL_COLLS:
            try:
                loc = frame.get(coll)
            except Exception:
                continue
            local_by_det[coll.replace("Tracks", "")] = len(loc)
            n_local += len(loc)
            del loc

        nhits = 0
        for coll in ("SiTargetHitsWindowed", "SiPadHitsWindowed",
                     "MTCSciFiHitsWindowed", "MTCScintHitsWindowed"):
            try:
                nhits += len(frame.get(coll))
            except Exception:
                pass

        stats["events"] += 1
        if n_spliced:
            stats["with_spliced"] += 1
        if n_global and n_local:
            stats["global_and_local"] += 1

        if clean:
            # Legibility: one or two long spliced tracks, nothing else. Every
            # loose local segment and every extra global costs, and so does a
            # big shower -- 3000 hits of hadronic spray hide the tracks the
            # picture is about.
            if n_spliced < 1 or n_global > 2:
                score = -1e9
            else:
                score = (20 * n_spliced
                         + longest / 4.0
                         - 6 * n_local
                         - nhits / 400.0)
        else:
            # A picture needs structure, not volume: splices weigh most, a
            # long track next, and both track counts are capped so a
            # fragmenting shower cannot win on sheer number.
            score = (10 * n_spliced
                     + 3 * min(n_local, 4)
                     + 2 * min(n_global, 4)
                     + longest / 10.0)
        scored.append((score, iev, n_global, n_local, n_spliced,
                       spliced_kinds, longest, nhits, dict(local_by_det)))
        del globals_, frame

    scored = [s for s in scored if s[0] > -1e8]
    scored.sort(reverse=True)
    print(f"file   {path}")
    print(f"events {stats['events']}   with a spliced global: "
          f"{stats['with_spliced']}   with globals AND locals: "
          f"{stats['global_and_local']}")
    print()
    print(f"{'event':>6} {'score':>6} {'glob':>5} {'loc':>4} {'spliced':>8} "
          f"{'longest':>8} {'hits':>7}  detail")
    for s in scored[:top]:
        (score, iev, ng, nl, nsp, kinds, longest, nhits, locd) = s
        detail = ", ".join(kinds) if kinds else "-"
        loc = " ".join(f"{k}:{v}" for k, v in locd.items() if v)
        print(f"{iev:>6} {score:>6.1f} {ng:>5} {nl:>4} {nsp:>8} "
              f"{longest:>8} {nhits:>7}  {detail} | {loc}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
