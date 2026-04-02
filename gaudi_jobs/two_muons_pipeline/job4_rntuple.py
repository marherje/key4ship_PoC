from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EDM4HEP2RNTuple

iosvc = IOSvc()
iosvc.Input = "digitized.edm4hep.root"

converter = EDM4HEP2RNTuple("EDM4HEP2RNTuple")
converter.InputFile      = "digitized.edm4hep.root"
converter.OutputFile     = "ShipHits.root"
converter.NTupleNames    = ["SiTarget", "SiPad"]
converter.Collections    = ["SiTargetHitsDigi", "SiPadHitsDigi"]
converter.BitFields      = [
    "system:8,layer:8,slice:4,plane:1,strip:14",
    "system:8,layer:8,slice:4,x:9,y:9"
]
converter.SourceIDParams = ["SiTargetSourceIDs", "SiPadSourceIDs"]
converter.DetectorIDs    = [0, 1]

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [converter],
    ExtSvc  = [iosvc]
)
