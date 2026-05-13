from k4FWCore import ApplicationMgr
from Configurables import EventShuffler

# EventShuffler reads input files directly via podio::ROOTReader.
# IOSvc is NOT used here: do not set iosvc.Input.
# All work happens in finalize() after Gaudi processes exactly 1 dummy event.

shuffler = EventShuffler("EventShuffler")
shuffler.InputFiles = [
    "../../output_mu-_xyz_1_-42.5_-1000_dir_0_0_1_E75_10events.edm4hep.root",
    # "../../output_e-_xyz_1_-82.5_-1000_dir_0_0_1_E75.edm4hep.root"
]
shuffler.SourceIDs = [
    1,
                    #   2
                      ]
shuffler.Delays    = [
    25.0001,
    # 15.0
    ]   # ns, adjust per source
shuffler.CollectionsSiTarget = [
    "SiTargetHits",
    # "SiTargetHits"
]
shuffler.CollectionsSiPad = [
    "SiPadHits",
    # "SiPadHits"
]
shuffler.CollectionsMTC = [
    "MTCDetHits",
    # "MTCDetHits"
]
shuffler.MaxEventsPerSource        = 10000
shuffler.OutputFile                = "shuffled.edm4hep.root"
shuffler.OutputCollectionSiTarget  = "SiTargetHitsMerged"
shuffler.OutputCollectionSiPad     = "SiPadHitsMerged"
shuffler.OutputCollectionMTCSciFi  = "MTCSciFiHitsMerged"
shuffler.OutputCollectionMTCScint  = "MTCScintHitsMerged"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,       # Only 1 Gaudi event: execute() is a no-op, all work is in finalize()
    TopAlg  = [shuffler],
    ExtSvc  = []
)
