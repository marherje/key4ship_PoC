from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventWindowSplitter

iosvc = IOSvc()
iosvc.Input = "superevent.edm4hep.root"
# No Output in IOSvc: EventWindowSplitter writes directly via podio::ROOTWriter

splitter = EventWindowSplitter("EventWindowSplitter")
splitter.InputCollectionSiTarget  = "SiTargetHitsMerged"
splitter.InputCollectionSiPixel   = "SiPixelHitsMerged"
splitter.OutputFile               = "timewindows.edm4hep.root"
splitter.OutputCollectionSiTarget = "SiTargetHitsWindowed"
splitter.OutputCollectionSiPixel  = "SiPixelHitsWindowed"
splitter.WindowSize = 25.0   # ns, adjust as needed

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [splitter],
    ExtSvc  = [iosvc]
)
