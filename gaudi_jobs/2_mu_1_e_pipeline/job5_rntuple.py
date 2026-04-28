from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EDM4HEP2RNTuple

iosvc = IOSvc()
iosvc.Input = ["tracks.edm4hep.root"]

converter = EDM4HEP2RNTuple("EDM4HEP2RNTuple")
converter.InputFile      = "tracks.edm4hep.root"
converter.OutputFile     = "ShipHits.root"
converter.NTupleNames    = ["SiTarget", "SiPad"]
converter.Collections    = ["SiTargetHitsWindowed", "SiPadHitsWindowed"]
converter.BitFields      = [
    "system:8,layer:8,slice:4,plane:1,column:2,row:2,strip:14",
    "system:8,layer:8,slice:4,x:9,y:9"
]

converter.SourceIDParams = ["SiTargetSourceIDs", "SiPadSourceIDs"]
converter.DetectorIDs         = [0, 1]
converter.TrackFile           = "tracks.edm4hep.root"   # path to tracking output
converter.TrackCollectionName = "ACTSTracks"    # default, can omit
converter.MeasCollections = [
    "SiTargetMeasurements",
    "SiPadMeasurements",
]
converter.MeasNtupleNames = [
    "SiTargetMeas",
    "SiPadMeas",
]
converter.MeasBitFields = [
    "system:8,layer:8,slice:4,plane:1,column:2,row:2,strip:14",
    "system:8,layer:8,slice:4,x:9,y:9",
]

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [converter],
    ExtSvc  = [iosvc]
)
