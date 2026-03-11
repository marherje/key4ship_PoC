from k4FWCore import ApplicationMgr
from Configurables import EventWindowSplitter

splitter = EventWindowSplitter("EventWindowSplitter")
splitter.InputFile                = "shuffled.root"
splitter.InputCollectionSiTarget  = "SiTargetHitsMerged"
splitter.InputCollectionSiPixel   = "SiPixelHitsMerged"
splitter.OutputFile               = "timewindows.root"
splitter.OutputCollectionSiTarget = "SiTargetHitsWindowed"
splitter.OutputCollectionSiPixel  = "SiPixelHitsWindowed"
splitter.WindowSize = 25.0   # ns

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [splitter],
    ExtSvc  = []
)
