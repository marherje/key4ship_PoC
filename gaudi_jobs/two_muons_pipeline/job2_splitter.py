from k4FWCore import ApplicationMgr
from Configurables import EventWindowSplitter

splitter = EventWindowSplitter("EventWindowSplitter")
splitter.InputFile                = "shuffled.edm4hep.root"
splitter.InputCollectionSiTarget  = "SiTargetHitsMerged"
splitter.InputCollectionSiPad   = "SiPadHitsMerged"
splitter.OutputFile               = "timewindows.edm4hep.root"
splitter.OutputCollectionSiTarget = "SiTargetHitsWindowed"
splitter.OutputCollectionSiPad  = "SiPadHitsWindowed"
splitter.WindowSize = 25.0   # ns

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [splitter],
    ExtSvc  = []
)
