# ---------------------------------------------------------------------------
# MOBO port of gaudi_jobs/genie_pipeline/job5_rntuple_seq.py.
#
# Identical physics and identical tuning -- the file is generated from that one
# and every comment explaining a value is preserved verbatim. Only the three
# things an optimization run needs are different:
#   * the geometry comes from $SND_COMPACT (the variant being evaluated), not
#     from the committed baseline,
#   * verbosity comes from $SND_OUTPUT_LEVEL, because at DEBUG one trial writes
#     hundreds of MB of log,
#   * the repo root is parents[2], one level less deep than the original.
# Keep the two in sync: this is the GENIE analysis chain the optimizer scores.
# ---------------------------------------------------------------------------
import os as _os
from pathlib import Path as _Path
from Gaudi.Configuration import DEBUG as _DEBUG, INFO as _INFO, \
                                WARNING as _WARNING, ERROR as _ERROR
_LEVELS = {"DEBUG": _DEBUG, "INFO": _INFO, "WARNING": _WARNING, "ERROR": _ERROR}
OUTPUT_LEVEL = _LEVELS.get(_os.environ.get("SND_OUTPUT_LEVEL", "INFO").upper(), _INFO)
_REPO = _Path(__file__).resolve().parents[2]
COMPACT = _os.environ.get(
    "SND_COMPACT", str(_REPO / "simulation" / "geometry" / "SND_compact.xml"))

# RNTuple export for the sequential tracker (GENIE).
#
# Same as job5_rntuple.py, but reads tracks_seq.edm4hep.root and exports all
# FOUR track collections instead of one: GlobalTracks (spliced) plus the three
# local ones. Each becomes a pair of RNTuples, <name> and <name>States, with
# exactly the schema the event display already reads for ACTSTracks/
# ACTSTrackStates.
#
#   k4run job5_rntuple_seq.py      # tracks_seq.edm4hep.root -> ShipHits_seq.root

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EDM4HEP2RNTuple
import sys
from pathlib import Path

sys.path.insert(0, str(_REPO / "simulation" / "geometry"))
from parse_geometry import SNDGeometry
geo = SNDGeometry(COMPACT)

iosvc = IOSvc()
iosvc.Input = ["tracks_seq.edm4hep.root"]

converter = EDM4HEP2RNTuple("EDM4HEP2RNTuple")
converter.InputFile      = "tracks_seq.edm4hep.root"
converter.OutputFile     = "ShipHits_seq.root"
converter.NTupleNames    = ["SiTarget", "SiPad", "MTCSciFi", "MTCScint"]
converter.Collections    = [
    "SiTargetHitsWindowed",
    "SiPadHitsWindowed",
    "MTCSciFiHitsWindowed",
    "MTCScintHitsWindowed",
]
converter.BitFields      = [
    geo.bitfields["SiTargetHits"],
    geo.bitfields["SiPadHits"],
    geo.bitfields["MTCDetHits"],   # MTCSciFi hits (plane=0,1)
    geo.bitfields["MTCDetHits"],   # MTCScint hits (plane=2)
]
converter.SourceIDParams = ["SiTargetSourceIDs", "SiPadSourceIDs",
                            "MTCSciFiSourceIDs", "MTCScintSourceIDs"]
converter.ContribPDGParams = ["SiTargetContribPDGs", "SiPadContribPDGs",
                              "MTCSciFiContribPDGs", "MTCScintContribPDGs"]
converter.DetectorIDs    = [0, 1, 3, 3]
converter.TrackFile      = "tracks_seq.edm4hep.root"
# GlobalTracks first so it is the one a reader defaulting to "the first track
# collection" picks up.
converter.TrackCollections = [
    "GlobalTracks",
    "SiTargetTracks",
    "SiPadTracks",
    "MTCTracks",
]
converter.MeasCollections = [
    "SiTargetMeasurements",
    "SiPadMeasurements",
    "MTCSciFiMeasurements",
]
converter.MeasNtupleNames = [
    "SiTargetMeas",
    "SiPadMeas",
    "MTCSciFiMeas",
]
converter.MeasBitFields = [
    geo.bitfields["SiTargetHits"],
    geo.bitfields["SiPadHits"],
    geo.bitfields["MTCDetHits"],
]

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [converter],
    ExtSvc  = [iosvc]
)
