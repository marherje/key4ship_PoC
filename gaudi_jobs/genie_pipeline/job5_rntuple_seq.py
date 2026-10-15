# RNTuple export for the sequential tracker (GENIE).
#
# Same as job5_rntuple.py, but reads tracks_seq.edm4hep.root and exports all
# FOUR track collections instead of one: GlobalTracks (spliced) plus the three
# local ones. Each becomes a pair of RNTuples, <name> and <name>States, with
# exactly the schema the event display already reads for ACTSTracks/
# ACTSTrackStates.
#
#   k4run job5_rntuple_seq.py      # tracks_seq.edm4hep.root -> ShipHits_seq.root

import os

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EDM4HEP2RNTuple
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "simulation" / "geometry"))
from parse_geometry import SNDGeometry
geo = SNDGeometry()

iosvc = IOSvc()
iosvc.Input = [os.environ.get("INPUT_FILE",
    "/eos/experiment/ship/user/marquezh/snd_pfa/genie/tracks_pfa_ref_a.root")]

converter = EDM4HEP2RNTuple("EDM4HEP2RNTuple")
# EDM4HEP2RNTuple opens the files itself (it bypasses IOSvc), so these have to
# follow INPUT_FILE too or the RNTuple silently describes a different sample
# than the one IOSvc was given.
converter.InputFile      = iosvc.Input[0]
converter.OutputFile     = os.environ.get(
    "HITS_FILE",
    "/eos/experiment/ship/user/marquezh/snd_pfa/genie/ShipHits_pfa_ref_a.root")
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
converter.TrackFile      = iosvc.Input[0]
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
