from k4FWCore import ApplicationMgr
from Configurables import EventOverlay
import os

# mu+ 10 GeV + pi+ 5 GeV from a common vertex at the FRONT of the SiTarget.
#
# Same pair and opening angle as mu_pi_pipeline, but fired from z = -2150 mm
# instead of -1000 mm. The SiTarget spans [-2200, -370] mm, so the pion now
# crosses its full ~1.8 m — some six nuclear interaction lengths of tungsten —
# rather than the last 0.6 m. It therefore usually interacts before reaching
# the SiPad, and that is the point of this sample: it produces tracks that
# begin and end inside the SiTarget, which the sequential tracker must emit as
# LOCAL SiTargetTracks rather than force into a global chain or discard.
#
# Single-source "overlay": the multi-PG ddsim output already carries both
# primaries in the same event, so this only converts them to the *Windowed
# collections consumed by job3/job4/job5.
#
# Regenerate the input with the committed launcher (pair "mu+:10:pi+:5:-2150"):
#
#     cd simulation/run_script && RUN_LOCAL=1 ./launch_multiplePG.sh
#     # fewer events:  NEVENTS=30 RUN_LOCAL=1 ./launch_multiplePG.sh

infile = os.environ.get(
    "INPUT_FILE",
    "../../simulation/run_script/data/"
    "output_QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0_z-2150.edm4hep.root")

overlay = EventOverlay("EventOverlay")
overlay.InputFiles          = [infile]
overlay.SourceIDs           = [1]
overlay.CollectionsSiTarget = ['SiTargetHits']
overlay.CollectionsSiPad    = ['SiPadHits']
overlay.CollectionsMTC      = ['MTCDetHits']
overlay.MaxEventsPerSource  = 0     # all events
overlay.OutputFile          = "events.edm4hep.root"

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,       # Only 1 Gaudi event: all work happens in execute()
    TopAlg  = [overlay],
    ExtSvc  = []
)
