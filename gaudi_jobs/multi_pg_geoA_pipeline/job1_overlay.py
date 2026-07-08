from k4FWCore import ApplicationMgr
from Configurables import EventOverlay
import os
import sys

# Single-source overlay for the geoA multi-PG samples: converts the ddsim
# output (both primaries already in the same event) to the *Windowed
# collections consumed by job3/job4. Select the sample with INPUT_FILE.

infile = os.environ.get("INPUT_FILE", "")
if not infile:
    print("Use: INPUT_FILE=<data/output_geoA_*.edm4hep.root> k4run job1_overlay.py")
    sys.exit(1)

overlay = EventOverlay("EventOverlay")
overlay.InputFiles          = [infile]
overlay.SourceIDs           = [1]
overlay.CollectionsSiTarget = ['SiTargetHits']
overlay.CollectionsSiPad    = ['SiPadHits']
overlay.CollectionsMTC      = ['MTCDetHits']
overlay.MaxEventsPerSource  = 0
overlay.OutputFile          = "events.edm4hep.root"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [overlay],
    ExtSvc  = []
)
