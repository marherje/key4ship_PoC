from k4FWCore import ApplicationMgr
from Configurables import EventShuffler

# EventShuffler reads input files directly via podio::ROOTReader.
# IOSvc is NOT used here: do not set iosvc.Input.
# All work happens in finalize() after Gaudi processes exactly 1 dummy event.

shuffler = EventShuffler("EventShuffler")
shuffler.InputFiles = [
    "output.edm4hep.root",
    "output2.edm4hep.root"
]
shuffler.SourceIDs = [1, 2]
shuffler.Delays    = [25.0001, 15.0]   # ns, adjust per source
shuffler.CollectionsSiTarget = [
    "SiTargetHits",
    "SiTargetHits"
]
shuffler.CollectionsSiPixel = [
    "SiPixelHits",
    "SiPixelHits"
]
shuffler.OutputFile               = "shuffled.root"
shuffler.OutputCollectionSiTarget = "SiTargetHitsMerged"
shuffler.OutputCollectionSiPixel  = "SiPixelHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,       # Only 1 Gaudi event: execute() is a no-op, all work is in finalize()
    TopAlg  = [shuffler],
    ExtSvc  = []
)
