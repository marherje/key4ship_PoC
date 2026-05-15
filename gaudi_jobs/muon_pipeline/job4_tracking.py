from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, SiTargetMeasConverter, \
                          SiPadMeasConverter, MTCSciFiMeasConverter, ACTSProtoTracker
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
MTC_INTER_GAP   = geo.mtc_inter_gap

MTC40_Z, MTC50_Z, MTC60_Z = geo.mtc_station_z_centers

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
iosvc.Input          = "timewindows.edm4hep.root"
iosvc.Output         = "tracks.edm4hep.root"
# keep * forwards all input collections (including MTCSciFiHitsWindowed and
# MTCScintHitsWindowed) to the output so job5 can write them to the RNTuple.
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
sitarget_conv.OutputLevel      = DEBUG

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = "SiPadHitsWindowed"
sipad_conv.OutputCollection = "SiPadMeasurements"
sipad_conv.BitField         = geo.bitfields["SiPadHits"]
sipad_conv.PixelSizeX       = geo.ecal_cell_size_x
sipad_conv.PixelSizeY       = geo.ecal_cell_size_y
sipad_conv.OutputLevel      = DEBUG

mtcscifi_conv = MTCSciFiMeasConverter("MTCSciFiMeasConverter")
mtcscifi_conv.GeoSvc           = "ACTSGeoSvc"
mtcscifi_conv.InputCollection  = "MTCSciFiHitsWindowed"
mtcscifi_conv.OutputCollection = "MTCSciFiMeasurements"
mtcscifi_conv.BitField         = geo.bitfields["MTCDetHits"]
mtcscifi_conv.StripPitch       = geo.mtc_scifi_channel_size
mtcscifi_conv.StereoAngleDeg   = geo.mtc_fiber_angle_deg
mtcscifi_conv.OutputLevel      = DEBUG

proto = ACTSProtoTracker("ACTSProtoTracker")
proto.GeoSvc           = "ACTSGeoSvc"
proto.InputSiTarget    = "SiTargetMeasurements"
proto.InputSiPad       = "SiPadMeasurements"
proto.InputMTC         = "MTCSciFiMeasurements"
proto.MTCStereoAngle   = geo.mtc_fiber_angle_deg
proto.OutputCollection = "ACTSTracks"
proto.IronFieldRanges  = _iron_ranges  # MultiRangeBField in MTC outer iron slabs
# Hough Transform automatic seeding
proto.AutoSeed         = True
proto.MaxSeeds         = 3
proto.HoughBinSize     = 5.0    # mm — coarser ok since we use 2D crossings
proto.HoughHalfSize    = 200.0  # mm
proto.HoughMinVotes    = 3      # crossings: each station contributes 1 crossing
proto.SeedCompatRadius = 8.0    # mm — radius for centroid refinement
proto.SeedStripPitch   = geo.sitarget_strip_pitch
proto.SeedMomentum     = 10.0   # GeV
proto.MaxChi2PerMeas        = 500.0
proto.HoughMaxMultiplicity  = 10.0  # safety net after isolation filter
# Crossing isolation filter for shower rejection
# Filters crossings (SiTarget) and positions (SiPad) by 2D density
# within each station/layer. Isolated = few neighbors within IsolationWindow.
proto.IsolationWindow       = 5.0   # mm — 2D radius for neighbor counting
proto.IsolationMaxNeighbors = 2     # max neighbors to be considered isolated
                                    # muon: 0 neighbors → always passes
                                    # shower: hundreds → always rejected

# Disable manual seeding (AutoSeed=True overrides these)
# proto.SeedPositions  = [...]
# proto.SeedDirections = [...]
proto.OutputLevel      = DEBUG

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [sitarget_conv, sipad_conv, mtcscifi_conv, proto],
    ExtSvc  = [iosvc, geosvc]
)
