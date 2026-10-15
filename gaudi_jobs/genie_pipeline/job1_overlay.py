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
    # SND-PFA reference run, generated 2026-09-22 with the current geometry
    # and the corrected flux frame. It replaces
    # output_genie_allflavours_0_2000.edm4hep.root, which must NOT be used any
    # more: simulated in August against the 166-layer SiTarget compact, its
    # hits reach layer 165 where today's geometry stops at 119, so 13 % of its
    # SiTarget hits address volumes that no longer exist and the VolumeManager
    # cannot place them. Its neutrino cone is also centred at z ~ 0 instead of
    # 25.76 m upstream. See docs/snd_pfa_progress.md.
    "../../simulation/run_script/data/output_genie_pfa_ref_a.edm4hep.root")

overlay.InputFiles = [infile]
overlay.SourceIDs           = [1]
overlay.CollectionsSiTarget = ['SiTargetHits']
overlay.CollectionsSiPad    = ['SiPadHits']
overlay.CollectionsMTC      = ['MTCDetHits']
overlay.MaxEventsPerSource  = 0     # all events
# Outputs live on EOS (docs/snd_pfa_plan.md section B); OUTPUT_FILE overrides.
overlay.OutputFile          = os.environ.get(
    "OUTPUT_FILE",
    "/eos/experiment/ship/user/marquezh/snd_pfa/genie/events_pfa_ref_a.root")

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = 1,
    TopAlg  = [overlay],
    ExtSvc  = []
)
