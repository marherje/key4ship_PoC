# SND-PFA over the sequential tracking output (docs/snd_pfa_plan.md).
#
#   k4run job6_pfa.py       # tracks_seq.edm4hep.root -> pfa_seq.edm4hep.root
#
# Reads the windowed hits, the measurements and the tracks that job4 wrote, and
# runs the particle-flow stages on them: hits -> layer clusters -> muon ->
# PFOs, with the truth clusters alongside as the ceiling.

import os

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import SNDCaloHitBuilder, SNDLayerClusterer, \
                          SNDTruthClusterer, SNDMuonFinder, SNDPFOBuilder, \
                          SNDHadronicClusterer
from Gaudi.Configuration import INFO
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent
                       / "simulation" / "geometry"))
from parse_geometry import SNDGeometry
geo = SNDGeometry()

iosvc = IOSvc()
iosvc.Input = [os.environ.get("INPUT_FILE",
    "/eos/experiment/ship/user/marquezh/snd_pfa/genie/tracks_pfa_ref_a.root")]
iosvc.Output = os.environ.get("OUTPUT_FILE",
    "/eos/experiment/ship/user/marquezh/snd_pfa/genie/pfa_pfa_ref_a.root")

builder = SNDCaloHitBuilder("SNDCaloHitBuilder")
builder.InputCollections = [
    "SiTargetHitsWindowed",
    "SiPadHitsWindowed",
    "MTCSciFiHitsWindowed",
    "MTCScintHitsWindowed",
]
# Same descriptors the measurement converters use, straight from the compact.
builder.ReadoutDescriptors = [
    geo.bitfields["SiTargetHits"],
    geo.bitfields["SiPadHits"],
    geo.bitfields["MTCDetHits"],
    geo.bitfields["MTCDetHits"],
]
builder.DetectorIDs = [0, 1, 2, 2]
builder.OutputCollections = [
    "SiTargetCaloHits",
    "SiPadCaloHits",
    "MTCSciFiCaloHits",
    "MTCScintCaloHits",
]
builder.OutputLinkCollections = [
    "SiTargetCaloHitLinks",
    "SiPadCaloHitLinks",
    "MTCSciFiCaloHitLinks",
    "MTCScintCaloHitLinks",
]
builder.ThresholdMIP = 0.5
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
builder.MTCLayersPerStation = geo.mtc_n_layers
# The sampling weights come from the material the compact really has in front
# of each sensitive plane, so they follow a MOBO variant's tungsten thickness
# instead of being frozen at the baseline's.
builder.CompactFile = "../../simulation/geometry/SND_compact.xml"
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
builder.OutputLevel = INFO

# One clusterer per collection: the cuts differ between 75 um strips, 5.5 mm
# pads and 50 mm tiles, and keeping the outputs separate lets a later stage
# ask for one detector's clusters without filtering.
CLUSTERERS = [
    # name,          input,               output,                 det, gap, seed
    ("SiTarget",  "SiTargetCaloHits",  "SiTargetLayerClusters",  0, 1, 2.0),
    ("SiPad",     "SiPadCaloHits",     "SiPadLayerClusters",     1, 0, 2.0),
    ("MTCSciFi",  "MTCSciFiCaloHits",  "MTCSciFiLayerClusters",  2, 1, 2.0),
    ("MTCScint",  "MTCScintCaloHits",  "MTCScintLayerClusters",  2, 0, 2.0),
]
READOUT = {0: geo.bitfields["SiTargetHits"], 1: geo.bitfields["SiPadHits"],
           2: geo.bitfields["MTCDetHits"]}

clusterers = []
for name, inColl, outColl, det, gap, seed in CLUSTERERS:
    alg = SNDLayerClusterer(f"LayerClusterer_{name}")
    alg.InputCollection = inColl
    alg.OutputCollection = outColl
    alg.ReadoutDescriptor = READOUT[det]
    alg.DetectorID = det
    alg.MTCLayersPerStation = geo.mtc_n_layers
    alg.MaxGapChannels = gap
    alg.SeedMIP = seed
    alg.OutputLevel = INFO
    clusterers.append(alg)

# One cluster per MC particle: the ceiling every later stage is measured
# against (CyberPFA's TruthClusteringAlg).
truth = SNDTruthClusterer("SNDTruthClusterer")
truth.HitCollections = ["SiTargetCaloHits", "SiPadCaloHits",
                        "MTCSciFiCaloHits", "MTCScintCaloHits"]
truth.SimCollections = ["SiTargetHitsWindowed", "SiPadHitsWindowed",
                        "MTCSciFiHitsWindowed", "MTCScintHitsWindowed"]
truth.ParticleIDParameters = ["SiTargetContribParticleIDs",
                              "SiPadContribParticleIDs",
                              "MTCSciFiContribParticleIDs",
                              "MTCScintContribParticleIDs"]
truth.OutputCollection = "TruthClusters"
truth.InputFile = iosvc.Input[0]
truth.OutputLevel = INFO

# Muon finding: a corridor along every track that reaches the MTC, with the
# momentum from the curvature and the iron range as a lower bound (see the
# algorithm's header for why).
muon = SNDMuonFinder("SNDMuonFinder")
muon.TrackCollections = ["GlobalTracks", "MTCTracks"]
muon.ClusterCollections = ["SiTargetLayerClusters", "SiPadLayerClusters",
                           "MTCSciFiLayerClusters", "MTCScintLayerClusters"]
muon.CompactFile = "../../simulation/geometry/SND_compact.xml"
# Chosen on the 3x3 scan of phase 3.3 (2026-09-27). Efficiency and fakes move
# together and no point reaches both targets at once: 6 layers / 25 mm gives
# 0.905 efficiency with 13.7 % NC fakes, 10 layers / 10 mm gives 0.744 with
# 5.1 %. This is the middle of that curve. What limits it is not the cuts but
# how many MTC layers the tracking reconstructs -- a median of 8 at 2 GeV and
# 34 above 20 -- so buying both at once means improving the MTC tracking.
#
# Retuned on 2026-10-14 after the corridor learned the stereo term (U/V
# clusters compared at s = x -+ y tan(a)): the window no longer has to be
# wide to absorb the missing |y| tan(a), and a narrow one keeps pions out. A
# corridor that runs to the last MTC layer -- a muon leaving through the
# back -- needs only 6 layers. Efficiency on nu_mu CC 0.636 -> 0.647, NC
# fakes 0.076 -> 0.049, nu_e CC fakes 0.093 -> 0.069 (docs/snd_pfa_progress.md,
# "Muon efficiency study").
muon.MinMTCLayers = 8
muon.MinMTCLayersAtExit = 6
muon.CorridorHalfWidthMTCSciFi = 2.0
muon.OutputLevel = INFO

# How the muon's own deposit leaves E_had. "trajectory" (default since
# 2026-10-14): the muon guns' mean deposit at the reconstructed p_mu is taken
# out plane by plane along the muon track, from whichever cluster holds it,
# and what a cluster cannot give comes off the event total -- unbiased, E_had
# in nu_mu CC on its ceiling. "trajectory-capped" stops at what each cluster
# holds (85 % of the muon taken), "corridor" is phase 4.1's rule (drop the
# corridor clusters, keep what they hold above one MIP), which left E_had
# 16 % high in nu_mu CC.
MUON_SUBTRACTION = os.environ.get("PFA_MUON_SUBTRACTION", "trajectory")

# 3D clusters of the hadronic system with CLUE3D (phase 4.3): the views of
# each layer are paired into (x, y) points and clustered by density across
# neighbouring layers, with one critical distance and density per detector.
hadclu = SNDHadronicClusterer("SNDHadronicClusterer")
# Along the trajectory the muon's clusters go into CLUE3D like any other and
# the PFO builder takes the muon's share back out of whichever 3D cluster
# holds them; with the corridor rule they stay out.
hadclu.ExcludeMuonCorridor = MUON_SUBTRACTION == "corridor"
hadclu.ClusterCollections = muon.ClusterCollections
hadclu.CompactFile = "../../simulation/geometry/SND_compact.xml"
hadclu.MTCStereoAngleDeg = geo.mtc_fiber_angle_deg
hadclu.CriticalDistance = [10.0, 15.0, 40.0]   # SiTarget, SiPad, MTC [mm]
hadclu.CriticalDensity = [0.02, 0.04, 0.2]     # [GeV]
hadclu.OutputLevel = INFO

# The event: the muon, a hadronic system made of every cluster less the
# muon's expected deposit, the vertex and E_nu. The summary numbers go to the PFASummary
# collection, in the order the PFASummaryNames metadata gives.
pfo = SNDPFOBuilder("SNDPFOBuilder")
pfo.ClusterCollections = muon.ClusterCollections
pfo.VertexTrackCollections = ["SiTargetTracks", "GlobalTracks"]
pfo.CompactFile = "../../simulation/geometry/SND_compact.xml"
pfo.MTCStereoAngleDeg = geo.mtc_fiber_angle_deg
# e/h of the silicon, measured on the phase-1.4 guns with the scales above in
# place (slope e 1.000 / slope pi 0.930). The MTC needs none: GlobalScaleMTC
# was itself fitted on pions.
pfo.EoverHECAL = 1.07546
pfo.EoverHMTC = 1.0
pfo.VertexMode = "clusters"
pfo.MuonSubtraction = MUON_SUBTRACTION
# E_had from the CLUE3D clusters is opt-in (PFA_HADRONIC_CLUSTERS=
# HadronicClusters); by default it stays the phase-4.1 sum over layer
# clusters, which measured better on 2000 GENIE events (phase 4.3,
# 2026-10-13): CLUE3D's outliers only remove energy, and in NC events the
# layer sum already sits on the ceiling.
pfo.HadronicClusters = os.environ.get("PFA_HADRONIC_CLUSTERS", "")
# The plan's 30 mm axis-to-vertex link keeps only 20 % of the 3D clusters'
# energy: a straight fit through a piece of a hadron shower in iron does not
# point back to the vertex that well (47 % within 100 mm, 82 % within 300 mm),
# and in a GENIE-only sample there is nothing else for it to reject. Off
# until a sample with backgrounds gives it something to do.
pfo.VertexLinkMM = 1e9
pfo.OutputLevel = INFO

ApplicationMgr(
    EvtSel="NONE",
    EvtMax=-1,
    TopAlg=[builder] + clusterers + [truth, muon, hadclu, pfo],
    ExtSvc=[],
    OutputLevel=INFO,
)
