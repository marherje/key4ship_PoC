from k4FWCore import ApplicationMgr
from Configurables import EventOverlay

# Single-source "overlay": converts GENIE-driven ddsim output (each event
# already carries the full GENIE final state -- outgoing lepton + post-FSI
# hadronic system, see simulation/run_script/genie_gst_to_hepmc.py and
# genie_base.py) to the *Windowed collections consumed by job3/job4/job5.
# Event-by-event, no time manipulation -- same pattern as e_pi_pipeline.

overlay = EventOverlay("EventOverlay")
import os

infile = os.environ.get(
    "INPUT_FILE",
    "../../simulation/run_script/data/output_genie_allflavours.edm4hep.root")

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
