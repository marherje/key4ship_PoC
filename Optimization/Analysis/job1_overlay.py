from k4FWCore import ApplicationMgr
from Configurables import EventOverlay

# Single-source "overlay": converts the multi-PG ddsim output (both primaries
# already in the same event) to the *Windowed collections consumed by
# job3/job4/job5. Event-by-event, no time manipulation.

overlay = EventOverlay("EventOverlay")
import os

infile = os.environ.get(
    "INPUT_FILE",
    "../../simulation/run_script/data/output_QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0.edm4hep.root")

overlay.InputFiles = [infile]
overlay.SourceIDs           = [1]
overlay.CollectionsSiTarget = ['SiTargetHits']
overlay.CollectionsSiPad    = ['SiPadHits']
overlay.CollectionsMTC      = ['MTCDetHits']
overlay.MaxEventsPerSource  = 0     # all events
overlay.OutputFile          = "events.edm4hep.root"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [overlay],
    ExtSvc  = []
)
