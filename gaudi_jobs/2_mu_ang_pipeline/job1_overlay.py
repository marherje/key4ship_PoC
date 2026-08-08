from k4FWCore import ApplicationMgr
from Configurables import EventOverlay
import os

# Two muons from a COMMON vertex, separated by an opening angle.
#
# Single-source "overlay": the multi-PG ddsim output already carries both
# primaries in the same event, so this only converts them to the *Windowed
# collections consumed by job3/job4/job5. Event-by-event, no time manipulation.
#
# Regenerate the input with the committed launcher (pair "mu-:10:mu-:10",
# angles 1.0 2.0 5.0 10.0):
#
#     cd simulation/run_script && ./launch_multiplePG.sh
#     # or, without HTCondor:
#     RUN_LOCAL=1 ./launch_multiplePG.sh
#
# WHY AN ANGLE. This pipeline replaces the old 2_mu_pipeline, which overlaid two
# single-gun files both fired with dir=(0,0,1) — exactly parallel to the beam
# and to each other. That is degenerate for the SiTarget: its planes measure X
# and Y separately, so with two parallel tracks the swapped pairing
# (track A's StripX, track B's StripY) is ALSO a perfect straight line and no
# fit can tell it from the real one. Two truly parallel tracks essentially do
# not occur in physics, and the reconstruction was being judged on an ambiguity
# that only that sample has. With an opening angle the ghost is not straight and
# the fit rejects it on its own.

infile = os.environ.get(
    "INPUT_FILE",
    "../../simulation/run_script/data/"
    "output_QGSP_BERT_SND_mu-_10GeV_mu-_10GeV_angle_1.0.edm4hep.root")

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
