from k4FWCore import ApplicationMgr, IOSvc
from Configurables import DelayTagger

iosvc = IOSvc()
iosvc.Input  = "output.edm4hep.root"
iosvc.Output = "delayed.edm4hep.root"

tagger = DelayTagger("DelayTagger")
tagger.InputCollectionSiTarget  = "SiTargetHits"
tagger.InputCollectionSiPixel   = "SiPixelHits"
tagger.OutputCollectionSiTarget = "SiTargetHitsDelayed"
tagger.OutputCollectionSiPixel  = "SiPixelHitsDelayed"
tagger.EventDelay = 25.0   # ns, adjust as needed

ApplicationMgr(
    EvtSel = "NONE",
    EvtMax = -1,
    TopAlg = [tagger],
    ExtSvc = [iosvc]
)
