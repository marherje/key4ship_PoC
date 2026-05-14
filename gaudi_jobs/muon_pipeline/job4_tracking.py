from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, SiTargetMeasConverter, \
                          SiPadMeasConverter, MTCSciFiMeasConverter, ACTSProtoTracker
from Gaudi.Configuration import DEBUG, INFO

# ── MTC outer iron magnetic field map ─────────────────────────────────────────
# MTC_BFIELD_Y must match MTC_BFieldY in simulation/geometry/SND_compact.xml.
# All geometry constants below must mirror the SND_compact.xml <define> block.
MTC_BFIELD_Y    = 0.0    # Tesla
MTC_FE_THICK    = 50.0   # mm — outer iron absorber (slice 0 per layer)
MTC_LAYER_THICK = 74.7   # mm — 50+3+1.35+1+1.35+3+15
MTC_N_LAYERS    = 15
ECAL_DIM_Z      = 20 * (3.5 + 11.0)   # 290 mm  (Ecal_NLayers × Ecal_LayerThickness)
MTC_TOTAL_LEN   = MTC_N_LAYERS * MTC_LAYER_THICK   # 1120.5 mm
MTC_INTER_GAP   = 1.0    # mm

MTC40_Z = ECAL_DIM_Z / 2.0 + 600.0 + MTC_TOTAL_LEN / 2.0
MTC50_Z = MTC40_Z + MTC_TOTAL_LEN + MTC_INTER_GAP
MTC60_Z = MTC50_Z + MTC_TOTAL_LEN + MTC_INTER_GAP

# Station (z_center, ACTS-Y half-size, ACTS-Z half-size) in mm + 1 mm margin.
# ACTS coords: x = beam axis (= DD4hep Z), y = DD4hep Y, z = DD4hep X.
_station_params = [
    (MTC40_Z, 201.0, 201.0),
    (MTC50_Z, 251.0, 251.0),
    (MTC60_Z, 301.0, 301.0),
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
sitarget_conv.BitField         = "system:8,layer:8,slice:4,plane:1,column:2,row:2,strip:14"
sitarget_conv.StripPitch       = 0.0755
sitarget_conv.NSensorCols  = 4
sitarget_conv.NSensorRows  = 2
sitarget_conv.SensorWidth  = 99.25
sitarget_conv.SensorHeight = 199.5
sitarget_conv.SensorGap    = 1.0
sitarget_conv.OutputLevel      = DEBUG

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = "SiPadHitsWindowed"
sipad_conv.OutputCollection = "SiPadMeasurements"
sipad_conv.BitField         = "system:8,layer:8,slice:4,x:9,y:9"
sipad_conv.PixelSizeX       = 5.5
sipad_conv.PixelSizeY       = 5.5
sipad_conv.OutputLevel      = DEBUG

mtcscifi_conv = MTCSciFiMeasConverter("MTCSciFiMeasConverter")
mtcscifi_conv.GeoSvc           = "ACTSGeoSvc"
mtcscifi_conv.InputCollection  = "MTCSciFiHitsWindowed"
mtcscifi_conv.OutputCollection = "MTCSciFiMeasurements"
mtcscifi_conv.BitField         = "system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9"
mtcscifi_conv.StripPitch       = 1.0   # mm
mtcscifi_conv.StereoAngleDeg   = 5.0
mtcscifi_conv.OutputLevel      = DEBUG

proto = ACTSProtoTracker("ACTSProtoTracker")
proto.GeoSvc           = "ACTSGeoSvc"
proto.InputSiTarget    = "SiTargetMeasurements"
proto.InputSiPad       = "SiPadMeasurements"
proto.InputMTC         = "MTCSciFiMeasurements"
proto.MTCStereoAngle   = 5.0
proto.OutputCollection = "ACTSTracks"
proto.IronFieldRanges  = _iron_ranges  # MultiRangeBField in MTC outer iron slabs
# Hough Transform automatic seeding
proto.AutoSeed         = True
proto.MaxSeeds         = 3
proto.HoughBinSize     = 5.0    # mm — coarser ok since we use 2D crossings
proto.HoughHalfSize    = 200.0  # mm
proto.HoughMinVotes    = 3      # crossings: each station contributes 1 crossing
proto.SeedCompatRadius = 8.0    # mm — radius for centroid refinement
proto.SeedStripPitch   = 0.0755 # mm — SiTarget strip pitch
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
