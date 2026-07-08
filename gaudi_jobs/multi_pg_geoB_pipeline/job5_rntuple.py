from k4FWCore import ApplicationMgr, IOSvc
import os
from Configurables import EDM4HEP2RNTuple
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "simulation" / "geometry"))
from parse_geometry import SNDGeometry
geo = SNDGeometry(str(Path(__file__).resolve().parent.parent.parent / "simulation" / "geometry" / "SND_compact_geoB.xml"))

iosvc = IOSvc()
tracks_file   = os.environ.get("TRACKS_FILE", "tracks.edm4hep.root")
shiphits_file = os.environ.get("SHIPHITS_FILE", "ShipHits.root")

iosvc.Input = [tracks_file]

converter = EDM4HEP2RNTuple("EDM4HEP2RNTuple")
converter.InputFile      = tracks_file
converter.OutputFile     = shiphits_file
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
converter.DetectorIDs         = [0, 1, 3, 3]
converter.TrackFile           = tracks_file
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
