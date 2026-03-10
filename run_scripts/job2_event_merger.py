from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventMerger

iosvc = IOSvc()
iosvc.Input  = "delayed.edm4hep.root"
# No Output in IOSvc: EventMerger writes directly via podio::ROOTWriter

merger = EventMerger("EventMerger")
merger.InputCollectionSiTarget  = "SiTargetHitsDelayed"
merger.InputCollectionSiPixel   = "SiPixelHitsDelayed"
merger.OutputFile               = "superevent.edm4hep.root"
merger.OutputCollectionSiTarget = "SiTargetHitsMerged"
merger.OutputCollectionSiPixel  = "SiPixelHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [merger],
    ExtSvc  = [iosvc]
)
