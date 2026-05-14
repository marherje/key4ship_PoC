from k4FWCore import ApplicationMgr
from Configurables import EventWindowSplitter

splitter = EventWindowSplitter("EventWindowSplitter")
splitter.InputFile                = "shuffled.edm4hep.root"
splitter.InputCollectionSiTarget  = "SiTargetHitsMerged"
splitter.InputCollectionSiPad     = "SiPadHitsMerged"
splitter.InputCollectionMTCSciFi  = "MTCSciFiHitsMerged"
splitter.InputCollectionMTCScint  = "MTCScintHitsMerged"
splitter.OutputFile               = "timewindows.edm4hep.root"
splitter.OutputCollectionSiTarget = "SiTargetHitsWindowed"
splitter.OutputCollectionSiPad    = "SiPadHitsWindowed"
splitter.OutputCollectionMTCSciFi = "MTCSciFiHitsWindowed"
splitter.OutputCollectionMTCScint = "MTCScintHitsWindowed"
splitter.WindowSize = 25.0   # ns

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [splitter],
    ExtSvc  = []
)
