# SND-PFA output as flat RNTuples (phase 5.1).
#
#   k4run job7_pfa_rntuple.py      # pfa (job6_pfa.py) -> ShipPFA.root
#
# EventPFA (one row per event: truth + PFASummary), PFOs and, with
# PFA_WRITE_CLUSTERS=1, LayerClusters. Read by mobo.ship.metrics.pfa_metrics.

import os

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import PFASummaryWriter
from Gaudi.Configuration import INFO

iosvc = IOSvc()
iosvc.Input = [os.environ.get("INPUT_FILE",
    "/eos/experiment/ship/user/marquezh/snd_pfa/genie/pfa_pfa_ref_a.root")]
# Only the event count matters here: the writer reads the file itself.
iosvc.CollectionNames = ["PFASummary"]

writer = PFASummaryWriter("PFASummaryWriter")
# PFASummaryWriter opens the file itself (the truth is in frame parameters and
# the PFASummary column names in the run metadata, which IOSvc does not hand
# to algorithms), so it must follow INPUT_FILE.
writer.InputFile = iosvc.Input[0]
writer.OutputFile = os.environ.get(
    "PFA_RNTUPLE_FILE",
    "/eos/experiment/ship/user/marquezh/snd_pfa/genie/ShipPFA_pfa_ref_a.root")
writer.WriteClusters = os.environ.get("PFA_WRITE_CLUSTERS", "0") == "1"
# Same e/h as SNDPFOBuilder in job6_pfa.py, for the truth ceiling E_vis_nomu.
writer.EoverHECAL = 1.07546
writer.EoverHMTC = 1.0
writer.OutputLevel = INFO

ApplicationMgr(
    EvtSel="NONE",
    EvtMax=-1,
    TopAlg=[writer],
    ExtSvc=[iosvc],
)
