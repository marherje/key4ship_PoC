# SND-PFA hit building for the calibration guns (phase 1.4).
#
#   INPUT_FILE=<events.root> OUTPUT_FILE=<pfa.root> k4run job6_pfa.py
#
# The weights come from the committed baseline compact, which is the geometry
# the guns were simulated in.

import os
import sys
from pathlib import Path

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import SNDCaloHitBuilder
from Gaudi.Configuration import INFO

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "simulation" / "geometry"))
from parse_geometry import SNDGeometry
geo = SNDGeometry()

iosvc = IOSvc()
iosvc.Input = [os.environ["INPUT_FILE"]]
iosvc.Output = os.environ["OUTPUT_FILE"]

builder = SNDCaloHitBuilder("SNDCaloHitBuilder")
builder.InputCollections = ["SiTargetHitsWindowed", "SiPadHitsWindowed",
                            "MTCSciFiHitsWindowed", "MTCScintHitsWindowed"]
builder.ReadoutDescriptors = [geo.bitfields["SiTargetHits"],
                              geo.bitfields["SiPadHits"],
                              geo.bitfields["MTCDetHits"],
                              geo.bitfields["MTCDetHits"]]
builder.DetectorIDs = [0, 1, 2, 2]
builder.OutputCollections = ["SiTargetCaloHits", "SiPadCaloHits",
                             "MTCSciFiCaloHits", "MTCScintCaloHits"]
builder.OutputLinkCollections = ["SiTargetCaloHitLinks", "SiPadCaloHitLinks",
                                 "MTCSciFiCaloHitLinks", "MTCScintCaloHitLinks"]
builder.MTCLayersPerStation = geo.mtc_n_layers
builder.CompactFile = str(REPO / "simulation" / "geometry" / "SND_compact.xml")
builder.WeightMode = "sampling"
# Fitted in phase 1.4 on the single-particle guns.
#
# ECAL: E_vis / E_true came out 0.660 for electrons, i.e. the sampling weights
# recover two thirds of an EM shower and this closes the rest. With it applied
# the response is exact across 1-20 GeV (0.99, 1.97, 5.01, 10.00, 20.00 GeV).
#
# MTC: fitted on the guns fired from INSIDE MTC40 (z = +200 mm), which is
# where most GENIE neutrinos interact -- the iron outweighs the SiTarget
# tungsten about five to one. Pions fired from the SiTarget leave under 3 % of
# their energy in the MTC, so they cannot constrain it. Its response is also
# more linear than the silicon's (5.1 % against 15.1 % worst residual over
# 2-40 GeV), which is worth knowing when E_had is built in phase 4.
builder.GlobalScaleECAL = 1.51601
builder.GlobalScaleMTC = 1.31372
builder.ThresholdMIP = float(os.environ.get("THRESHOLD_MIP", "0.5"))
builder.MIPEnergyOverride = [
    # 1 MIP = what the 10 GeV muon actually deposits (Landau most probable
    # value), not the mean minimum-ionising loss the sampling scan computes:
    # delta rays escape a thin layer, so the MPV is ~0.73 of the mean in
    # 300 um of silicon and ~0.88 in the 20 mm scintillator, exactly the
    # trend measured in phase 1.4. The MTC scintillator lands on 3.53 MeV,
    # which the tracking report had measured independently as ~3.5 MeV/MIP.
    "0:0:0.08432", "0:1:0.08432",   # SiTarget X, Y
    "1:2:0.19531",                  # SiPad
    "2:0:0.19757", "2:1:0.19757",   # MTC SciFi U, V
    "2:2:3.53253",                  # MTC scintillator
]
builder.OutputLevel = INFO

ApplicationMgr(EvtSel="NONE", EvtMax=-1, TopAlg=[builder], ExtSvc=[],
               OutputLevel=INFO)
