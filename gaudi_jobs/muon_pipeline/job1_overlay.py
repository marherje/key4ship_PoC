from k4FWCore import ApplicationMgr
from Configurables import EventOverlay

# EventOverlay reads the input files directly via podio::ROOTReader and
# overlays event i of every source into output event i: event-by-event
# processing, no time shuffling and no window splitting. Output collections
# keep the *Windowed names so job3/job4/job5 run unchanged.

overlay = EventOverlay("EventOverlay")
overlay.InputFiles = [
    "../../simulation/run_script/data/output_mu-_xyz_1_-82.5_-1000_dir_0_0_1_E5.edm4hep.root"
]
overlay.SourceIDs           = [1]
overlay.CollectionsSiTarget = ['SiTargetHits']
overlay.CollectionsSiPad    = ['SiPadHits']
overlay.CollectionsMTC      = ['MTCDetHits']
overlay.MaxEventsPerSource  = 50
overlay.OutputFile          = "events.edm4hep.root"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,       # Only 1 Gaudi event: all work happens in execute()
    TopAlg  = [overlay],
    ExtSvc  = []
)
