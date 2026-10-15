"""Shared MC-truth join for the SND-PFA analysis scripts.

The chain keeps MC truth in two pieces, and this module is the one place that
knows how to put them back together:

  * `*HitsWindowed` + their `*_Contributions` hold the energy deposits;
  * the frame parameters `*ContribPDGs` / `*ContribParticleIDs` hold, one entry
    per contribution in collection order, the PDG and the index of the
    MCParticle it came from -- that index being local to the source file it was
    read from, so `MCParticleSourceOffsets` (paired with `MCParticleSourceIDs`)
    shifts it into the merged `MCParticles` collection. EventOverlay stores the
    source id of a contribution in `CaloHitContribution::PDG`.

Since job4 links every written Track to the TrackerHit3Ds it was fitted to, a
track's MC particle is just the dominant particle of its hits, which is what
`TruthJoin.track_truth` returns.
"""

HIT_COLLECTIONS = {
    # measurement collection -> (hit collection, PDG param, particle-id param)
    "SiTargetMeasurements": ("SiTargetHitsWindowed", "SiTargetContribPDGs",
                             "SiTargetContribParticleIDs"),
    "SiPadMeasurements": ("SiPadHitsWindowed", "SiPadContribPDGs",
                          "SiPadContribParticleIDs"),
    "MTCSciFiMeasurements": ("MTCSciFiHitsWindowed", "MTCSciFiContribPDGs",
                             "MTCSciFiContribParticleIDs"),
}
# Hit collections in their own right (the MTC scint has no measurements).
ALL_HIT_COLLECTIONS = [
    ("SiTargetHitsWindowed", "SiTargetContribPDGs", "SiTargetContribParticleIDs"),
    ("SiPadHitsWindowed", "SiPadContribPDGs", "SiPadContribParticleIDs"),
    ("MTCSciFiHitsWindowed", "MTCSciFiContribPDGs", "MTCSciFiContribParticleIDs"),
    ("MTCScintHitsWindowed", "MTCScintContribPDGs", "MTCScintContribParticleIDs"),
]

NU_PDGS = {12, 14, 16}
LEPTON_PDGS = {11, 13, 15}


def param(frame, name):
    """Frame parameter as a list; [] when absent or unreadable."""
    try:
        val = frame.get_parameter(name)
    except Exception:
        return []
    if val is None:
        return []
    if isinstance(val, (int, float, str)):
        return [val]
    return list(val)


class TruthJoin:
    """Per-event map from cellID to the MC particles that lit it."""

    def __init__(self, frame):
        self.frame = frame
        self.offsets = param(frame, "MCParticleSourceOffsets")
        self.source_ids = param(frame, "MCParticleSourceIDs")
        try:
            self.mcs = frame.get("MCParticles")
        except Exception:
            self.mcs = None
        # cellID -> {mc index: energy}
        self.by_cell = {}
        for coll, ppar, ipar in ALL_HIT_COLLECTIONS:
            try:
                hits = frame.get(coll)
            except Exception:
                continue
            ids = param(frame, ipar)
            if not ids:
                continue
            k = 0
            for hit in hits:
                acc = self.by_cell.setdefault(hit.getCellID(), {})
                for c in hit.getContributions():
                    if k >= len(ids):
                        break
                    idx = self._global_index(ids[k], c.getPDG())
                    if idx is not None:
                        acc[idx] = acc.get(idx, 0.0) + c.getEnergy()
                    k += 1

    def _global_index(self, local_index, source_id):
        if local_index is None or local_index < 0:
            return None
        off = 0
        if self.offsets:
            if self.source_ids and source_id in self.source_ids:
                off = self.offsets[self.source_ids.index(source_id)]
            elif len(self.offsets) == 1:
                off = self.offsets[0]
        idx = local_index + off
        if self.mcs is None or idx >= len(self.mcs):
            return None
        return idx

    # ---- event-level truth -------------------------------------------------

    def neutrino(self):
        """The incoming neutrino (generatorStatus 4 beam particle), or None."""
        if self.mcs is None:
            return None
        for p in self.mcs:
            if abs(p.getPDG()) in NU_PDGS:
                return p
        return None

    def primary_lepton(self):
        """Outgoing primary charged lepton of a CC event, or None."""
        if self.mcs is None:
            return None
        best = None
        for p in self.mcs:
            if p.getGeneratorStatus() != 1:
                continue
            if abs(p.getPDG()) not in LEPTON_PDGS:
                continue
            if best is None or p.getEnergy() > best.getEnergy():
                best = p
        return best

    def primary_muon(self):
        lep = self.primary_lepton()
        return lep if (lep is not None and abs(lep.getPDG()) == 13) else None

    def index_of(self, particle):
        return particle.getObjectID().index if particle is not None else -1

    # ---- object-level truth ------------------------------------------------

    def cell_truth(self, cell_id):
        """{mc index: energy} for one cellID."""
        return self.by_cell.get(cell_id, {})

    def dominant(self, weights):
        """(mc index, purity) of the largest entry of a {index: energy} map."""
        if not weights:
            return -1, 0.0
        total = sum(weights.values())
        idx = max(weights, key=weights.get)
        return idx, (weights[idx] / total if total > 0 else 0.0)

    def track_truth(self, track):
        """(mc index, purity, n_hits) of a Track, from its linked hits."""
        acc = {}
        n = 0
        for hit in track.getTrackerHits():
            n += 1
            for idx, e in self.cell_truth(hit.getCellID()).items():
                acc[idx] = acc.get(idx, 0.0) + e
        idx, pur = self.dominant(acc)
        return idx, pur, n
