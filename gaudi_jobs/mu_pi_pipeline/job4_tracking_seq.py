# Sequential per-subdetector tracking (ACTSSequentialTracker).
#
# The alternative to job4_tracking.py, which runs the global CKF
# (ACTSProtoTracker) with the SiPad up-weighted. Both are built; this steering
# picks the sequential one and writes to its own output file, so the two can be
# run back to back on the same events and compared.
#
#   k4run job4_tracking_seq.py     # events.edm4hep.root -> tracks_seq.edm4hep.root
#
# Everything up to the measurement converters is identical to job4_tracking.py
# — the difference starts at the tracker.

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, SiTargetMeasConverter, \
                          SiPadMeasConverter, MTCSciFiMeasConverter, \
                          ACTSSequentialTracker
from Gaudi.Configuration import DEBUG, INFO
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "simulation" / "geometry"))
from parse_geometry import SNDGeometry
geo = SNDGeometry()

# ── MTC outer iron magnetic field map ─────────────────────────────────────────
# All geometry constants come from SND_compact.xml via SNDGeometry.
MTC_BFIELD_Y    = geo.mtc_bfield_y
MTC_FE_THICK    = geo.mtc_outer_fe_thick
MTC_LAYER_THICK = geo.mtc_layer_thick
MTC_N_LAYERS    = geo.mtc_n_layers
MTC_TOTAL_LEN   = MTC_N_LAYERS * MTC_LAYER_THICK

# Station (z_center, ACTS-Y half-size, ACTS-Z half-size) in mm + 1 mm margin.
# ACTS coords: x = beam axis (= DD4hep Z), y = DD4hep Y, z = DD4hep X.
_station_params = [
    (z, h + 1.0, h + 1.0)
    for z, h in zip(geo.mtc_station_z_centers, geo.mtc_station_env_half_heights)
]
# Flat list: [xlo, xhi, ylo, yhi, zlo, zhi, by] per outer iron slab
_iron_ranges = []
for _z_ctr, _yhalf, _zhalf in _station_params:
    _front = _z_ctr - MTC_TOTAL_LEN / 2.0
    for _layer in range(MTC_N_LAYERS):
        _lo = _front + _layer * MTC_LAYER_THICK
        _hi = _lo + MTC_FE_THICK
        _iron_ranges.extend([_lo, _hi, -_yhalf, _yhalf, -_zhalf, _zhalf, MTC_BFIELD_Y])
# ──────────────────────────────────────────────────────────────────────────────

iosvc = IOSvc()
iosvc.Input          = "events.edm4hep.root"
iosvc.Output         = "tracks_seq.edm4hep.root"
# keep * forwards all input collections to the output so job5 can write them.
iosvc.outputCommands = ["keep *"]

geosvc = ACTSGeoSvc("ACTSGeoSvc")
geosvc.CompactFile = "../../simulation/geometry/SND_compact.xml"
geosvc.OutputLevel = INFO

sitarget_conv = SiTargetMeasConverter("SiTargetMeasConverter")
sitarget_conv.GeoSvc           = "ACTSGeoSvc"
sitarget_conv.InputCollection  = "SiTargetHitsWindowed"
sitarget_conv.OutputCollection = "SiTargetMeasurements"
sitarget_conv.BitField         = geo.bitfields["SiTargetHits"]
sitarget_conv.StripPitch       = geo.sitarget_strip_pitch
sitarget_conv.NSensorCols  = geo.sitarget_sensor_ncols
sitarget_conv.NSensorRows  = geo.sitarget_sensor_nrows
sitarget_conv.SensorWidth  = geo.sitarget_sensor_width
sitarget_conv.SensorHeight = geo.sitarget_sensor_height
sitarget_conv.SensorGap    = geo.sitarget_sensor_gap
sitarget_conv.OutputLevel      = INFO

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = "SiPadHitsWindowed"
sipad_conv.OutputCollection = "SiPadMeasurements"
sipad_conv.BitField         = geo.bitfields["SiPadHits"]
sipad_conv.PixelSizeX       = geo.ecal_cell_size_x
sipad_conv.PixelSizeY       = geo.ecal_cell_size_y
sipad_conv.OutputLevel      = INFO

mtcscifi_conv = MTCSciFiMeasConverter("MTCSciFiMeasConverter")
mtcscifi_conv.GeoSvc           = "ACTSGeoSvc"
mtcscifi_conv.InputCollection  = "MTCSciFiHitsWindowed"
mtcscifi_conv.OutputCollection = "MTCSciFiMeasurements"
mtcscifi_conv.BitField         = geo.bitfields["MTCDetHits"]
mtcscifi_conv.StripPitch       = geo.mtc_scifi_channel_size
mtcscifi_conv.StereoAngleDeg   = geo.mtc_fiber_angle_deg
mtcscifi_conv.PairMethod       = "AllPairs"
mtcscifi_conv.PairMaxDz        = 25.0   # mm — accept U/V partners within this Δz
mtcscifi_conv.OutputLevel      = INFO

seq = ACTSSequentialTracker("ACTSSequentialTracker")
seq.GeoSvc         = "ACTSGeoSvc"
seq.InputSiTarget  = "SiTargetMeasurements"
seq.InputSiPad     = "SiPadMeasurements"
seq.InputMTC       = "MTCSciFiMeasurements"
seq.MTCStereoAngle = geo.mtc_fiber_angle_deg
seq.SiTargetBitField = geo.bitfields["SiTargetHits"]
seq.MTCBitField      = geo.bitfields["MTCDetHits"]
seq.IronFieldRanges  = _iron_ranges  # field inside the MTC outer iron slabs

# Outputs: three local collections (segments that spliced with nothing) and the
# global one (spliced chains). A segment absorbed into a global track does NOT
# appear in its local collection.
seq.OutputSiTargetTracks = "SiTargetTracks"
seq.OutputSiPadTracks    = "SiPadTracks"
seq.OutputMTCTracks      = "MTCTracks"
seq.OutputGlobalTracks   = "GlobalTracks"

# ── Seeding (Hough, run once per detector per iteration) ──────────────────────
seq.MaxSeeds         = 5
seq.HoughBinSize     = 5.0    # mm
seq.HoughHalfSize    = 200.0  # mm
seq.HoughMinVotes    = 3
seq.SeedCompatRadius = 8.0    # mm — radius for centroid refinement
seq.SeedStripPitch   = geo.sitarget_strip_pitch
seq.SeedMomentum     = 3.0    # GeV

# ── Segment finding ──────────────────────────────────────────────────────────
# Two iterations per detector. The first takes the clean tracks; its hits are
# then removed and the second sweeps what is left with looser cuts, which is
# what recovers a track masked by a better one.
seq.Chi2CutOff       = 70.0   # CKF per-surface chi2, first iteration
seq.NumMeasCutOff    = 1
seq.LooseChi2Scale   = 2.0    # second iteration: 140
seq.SegMaxChi2Tight  = 10.0   # segment chi2/ndf, first iteration
seq.SegMaxChi2Loose  = 20.0   # ... second
seq.MinMeasPerSegment = [4, 3, 4]   # [SiTarget, SiPad, MTC]; SiPad hits are 2D
seq.DuplicateOverlapFraction = 0.7

# ── Splicing ─────────────────────────────────────────────────────────────────
# Only neighbours: SiTarget<->SiPad and SiPad<->MTC. A full chain forms through
# both links; there is deliberately no SiTarget->MTC shortcut. The geometric
# window only proposes candidates — the combined refit decides.
seq.SpliceMaxDist       = 15.0   # mm at the inter-detector boundary
seq.SpliceMaxAngle      = 0.1    # rad
seq.SpliceAngleWeight   = 100.0  # mm^2/rad^2
seq.GlobalMaxChi2PerNdf = 10.0   # combined refit acceptance

# ── Shower cleaning (same values as the global tracker) ──────────────────────
seq.HitPurgeWindow        = [1.0, 8.0, 0.0]   # mm, per [SiTarget, SiPad, MTC]
seq.HitPurgeMaxNeighbors  = [8, 4, 0]
seq.IsolationWindow       = 5.0   # mm — 2D radius for neighbour counting
seq.IsolationMaxNeighbors = 2
seq.HoughMaxMultiplicity  = 10.0  # safety net after the isolation filter
# SiTarget ghost resolution. With N tracks the 1D strips give N^2 crossings, of
# which N are real; a ghost seed can otherwise build a hybrid segment out of the
# strips of two different tracks, which then extrapolates onto neither and costs
# BOTH tracks their splice. A SiPad 2D point within this window confirms a
# crossing as real; only the strip-swapped ghosts of confirmed crossings are
# dropped, so a track that stops inside the SiTarget keeps its seed.
seq.GhostResolveWindow    = 10.0  # mm

# NOTE: no SiPadWeight and no AnnealingIterations. Those properties do not
# exist on this algorithm: the per-detector fits have no cross-detector
# ambiguity to resolve, so every measurement enters at its nominal variance.

seq.OutputLevel = DEBUG

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [sitarget_conv, sipad_conv, mtcscifi_conv, seq],
    ExtSvc  = [iosvc, geosvc]
)
