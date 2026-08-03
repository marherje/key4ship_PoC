from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EDM4HEP2RNTuple
import os
import sys
from pathlib import Path

# Optimization/Analysis/jobX.py -> parents[2] is the repo root (one level less
# deep than the gaudi_jobs/<pipeline>/ originals, which use parent.parent.parent).
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "simulation" / "geometry"))
from parse_geometry import SNDGeometry

# The variant being analysed; falls back to the committed baseline.
COMPACT = os.environ.get(
    "SND_COMPACT", str(REPO / "simulation" / "geometry" / "SND_compact.xml"))
geo = SNDGeometry(COMPACT)

iosvc = IOSvc()
iosvc.Input = ["tracks.edm4hep.root"]

converter = EDM4HEP2RNTuple("EDM4HEP2RNTuple")
converter.InputFile      = "tracks.edm4hep.root"
converter.OutputFile     = "ShipHits.root"
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
    geo.bitfields["MTCDetHits"],
    geo.bitfields["MTCDetHits"],
]
converter.SourceIDParams = ["SiTargetSourceIDs", "SiPadSourceIDs",
                            "MTCSciFiSourceIDs", "MTCScintSourceIDs"]
converter.ContribPDGParams = ["SiTargetContribPDGs", "SiPadContribPDGs",
                              "MTCSciFiContribPDGs", "MTCScintContribPDGs"]
converter.DetectorIDs         = [0, 1, 3, 3]
converter.TrackFile           = "tracks.edm4hep.root"
converter.TrackCollectionName = "ACTSTracks"
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
