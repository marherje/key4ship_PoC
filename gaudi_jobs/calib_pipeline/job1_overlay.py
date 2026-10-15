# Overlay for the single-particle calibration guns (phase 1.3/1.4 of
# docs/snd_pfa_plan.md). Same single-source conversion the GENIE pipeline
# does, but reading and writing on EOS: these samples are ~370 MB of input
# and the AFS quota is 10 GB.
#
#   INPUT_FILE=<gun.edm4hep.root> OUTPUT_FILE=<events.root> k4run job1_overlay.py

import os

from k4FWCore import ApplicationMgr
from Configurables import EventOverlay

overlay = EventOverlay("EventOverlay")
overlay.InputFiles = [os.environ["INPUT_FILE"]]
overlay.SourceIDs           = [1]
overlay.CollectionsSiTarget = ['SiTargetHits']
overlay.CollectionsSiPad    = ['SiPadHits']
overlay.CollectionsMTC      = ['MTCDetHits']
overlay.MaxEventsPerSource  = int(os.environ.get("MAX_EVENTS", "0"))
overlay.OutputFile          = os.environ["OUTPUT_FILE"]

ApplicationMgr(EvtSel="NONE", EvtMax=1, TopAlg=[overlay], ExtSvc=[])
