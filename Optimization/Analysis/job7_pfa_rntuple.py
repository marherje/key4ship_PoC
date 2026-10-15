# ---------------------------------------------------------------------------
# MOBO port of gaudi_jobs/genie_pipeline/job7_pfa_rntuple.py.
#
# Identical physics and identical tuning -- the file is generated from that one
# and every comment explaining a value is preserved verbatim. Only the things
# an optimization run needs are different:
#   * the geometry comes from $SND_COMPACT (the variant being evaluated), not
#     from the committed baseline -- the sampling weights, the iron profile and
#     the layer sequence of the PFA all follow it,
#   * verbosity comes from $SND_OUTPUT_LEVEL, because at DEBUG one trial writes
#     hundreds of MB of log,
#   * the files are the bare names of the trial directory
#     (pfa_seq.edm4hep.root -> ShipPFA.root).
# Keep the two in sync: this is the PFA chain the optimizer scores.
# ---------------------------------------------------------------------------
import os as _os
from pathlib import Path as _Path
from Gaudi.Configuration import DEBUG as _DEBUG, INFO as _INFO, \
                                WARNING as _WARNING, ERROR as _ERROR
_LEVELS = {"DEBUG": _DEBUG, "INFO": _INFO, "WARNING": _WARNING, "ERROR": _ERROR}
OUTPUT_LEVEL = _LEVELS.get(_os.environ.get("SND_OUTPUT_LEVEL", "INFO").upper(), _INFO)
_REPO = _Path(__file__).resolve().parents[2]
COMPACT = _os.environ.get(
    "SND_COMPACT", str(_REPO / "simulation" / "geometry" / "SND_compact.xml"))

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
iosvc.Input = ["pfa_seq.edm4hep.root"]
# Only the event count matters here: the writer reads the file itself.
iosvc.CollectionNames = ["PFASummary"]

writer = PFASummaryWriter("PFASummaryWriter")
# PFASummaryWriter opens the file itself (the truth is in frame parameters and
# the PFASummary column names in the run metadata, which IOSvc does not hand
# to algorithms), so it must follow INPUT_FILE.
writer.InputFile = iosvc.Input[0]
writer.OutputFile = "ShipPFA.root"
writer.WriteClusters = os.environ.get("PFA_WRITE_CLUSTERS", "0") == "1"
# Same e/h as SNDPFOBuilder in job6_pfa.py, for the truth ceiling E_vis_nomu.
writer.EoverHECAL = 1.07546
writer.EoverHMTC = 1.0
writer.OutputLevel = OUTPUT_LEVEL

ApplicationMgr(
    EvtSel="NONE",
    EvtMax=-1,
    TopAlg=[writer],
    ExtSvc=[iosvc],
)
