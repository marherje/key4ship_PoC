# Sequential per-subdetector tracking (ACTSSequentialTracker) for GENIE.
#
# The alternative to job4_tracking.py, which runs the global CKF
# (ACTSProtoTracker) with the SiPad up-weighted. Both are built; this steering
# picks the sequential one and writes to its own output file, so the two can
# be run back to back on the same events and compared.
#
#   k4run job4_tracking_seq.py     # events.edm4hep.root -> tracks_seq.edm4hep.root
#
# Everything up to the measurement converters is identical to job4_tracking.py
# — the difference starts at the tracker. See gaudi_jobs/e_pi_pipeline's
# job4_tracking_seq.py for the pattern this mirrors.

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
# (z_center, ACTS-Y half-size, ACTS-Z half-size). ACTS z is DD4hep X, and THAT
# is the axis along which the magnetised iron core is narrower than the
# envelope: MTCDetector.cpp registers each slab as
#     {zlo, zhi, xhalf = iron_width/2, yhalf = env_height/2 + 1}
# and its MTCIronField returns By ONLY inside those slabs -- the steel filler
# flanking the core in X carries no field. Using the envelope half-height for
# both axes, as this used to, magnetised a 100 mm band per side in MTC40 and
# 50 mm in MTC50 that Geant4 leaves field-free. A track in that band is
# simulated straight and propagated curved: the prediction misses, the CKF
# rejects the hit (chi2 > Chi2CutOff) and the segment dies. The event-1 muon
# of the GENIE sample crosses that band in MTC50 (it reaches x = -300). This
# fix alone did not recover it -- the loc0 sign bug in the calibrator did --
# but the band is real and the map must match what Geant4 simulated.
_iron_x_half = [
    geo.constants[f"MTC{w}_iron_width"] / 2.0 for w in (40, 50, 60)
]
_station_params = [
    (z, h + 1.0, xh + 1.0)
    for z, h, xh in zip(geo.mtc_station_z_centers,
                        geo.mtc_station_env_half_heights, _iron_x_half)
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
# PairMaxDz is a sanity guard against pairing ACROSS layers, not a tuning knob:
# measured this session, U and V of one layer sit 2.65 mm apart (z=54.675 vs
# 57.325 in station 0) while consecutive layers are 80 mm apart, so every value
# between ~3 and ~79 selects exactly the same pairs. Scans at 25/15/8 mm gave
# bit-identical output; going below 2.65 would drop ALL pairs at once.
mtcscifi_conv.PairMaxDz        = 25.0   # mm — accept U/V partners within this Δz
# U/V ghost veto (MTCSciFiMeasConverter.ScintCollection/ScintConfirmWindow):
# implemented and measured this session, left OFF because it trades the metric
# that matters for one that does not.
#
# The idea: with N particles crossing a layer, AllPairs builds N^2 crossings of
# which N are real; a ghost takes U from one particle and V from another, so it
# lands where no particle was and drags the CKF along a fictitious trajectory.
# The MTC scintillator plane (50 mm cells, otherwise unused by the
# reconstruction) can confirm the real ones -- and it does separate particles:
# a particle lights a median of 2 cells, and 53% of multi-particle layers have
# completely disjoint cell sets.
#
# But on the full 2000-event GENIE chunk, at its best window (30 mm) it gives
# GlobalTracks 1204 -> 1224 (+20) while events with >=1 GlobalTracks fall
# 618 -> 606 (-12): it sometimes vetoes the only 2D pair holding up an event's
# only track. Since the goal is that a track gets reconstructed at all, events
# is the metric that rules, and it loses on that one. (12 mm is far worse:
# segments -47%, tracks -40%. 50 mm is untested at full statistics.)
# To re-enable, set both properties -- setting only one warns and does nothing.
#
# RE-MEASURED 2026-09-11 after the loc0 sign fix, judged by MUON efficiency
# (events whose MC muon comes out as a muon-dominated segment; 1117 events):
#   scint off              : muon eff 65.0%  meas/ev 1086  >5 tracks/ev 6.2%
#   scint 30 mm            : muon eff 64.9%  meas/ev  836  >5 tracks/ev 5.3%
#   scint 30 + MTC purge   : muon eff 65.4%  meas/ev  309  >5 tracks/ev 0.6%
# The -12 events above were shower fragments, not muons. ON since then, together
# with the MTC density purge (see HitPurgeWindow).
mtcscifi_conv.ScintCollection    = "MTCScintHitsWindowed"
mtcscifi_conv.ScintConfirmWindow = 30.0   # mm
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

# ── Seeding (local straight-line tracklets, per detector per iteration) ───────
# The seeder finds short straight tracklets in sliding windows of layers and
# emits position AND direction. It replaced a transverse-position accumulator
# that could only find beam-parallel tracks — see the header comment above
# findSeeds() in ACTSSequentialTracker.cpp for the measurements behind these
# values.
seq.MaxSeeds         = 20     # windows produce more candidates than the old
                              # single-pass accumulator; 5 was throwing tracks away
seq.HoughBinSize     = 5.0    # mm — intercept bin
seq.HoughHalfSize    = 350.0  # mm — MUST exceed the detector half-size (MTC
                              # station envelope is 300 mm). At the old 200 mm
                              # every track outside that window was lost.
seq.HoughMinVotes    = 3
seq.SeedWindowLayers = 4      # ~240 mm in the MTC. 6 was tried first from the
                              # MC-truth residual table alone and is WORSE
                              # (33 vs 39 events with a track in 60): that table
                              # averages over muons, and a 3 GeV muon bending in
                              # the MTC iron changes slope by dtx~0.24 over 10
                              # layers, so a 6-layer straight fit returns a
                              # chord, not the tangent at the seed point, and
                              # the CKF loses the track after a few surfaces.
seq.SeedWindowStride = 2      # overlapping windows, so a track crossing a
                              # boundary is still found whole by its neighbour
seq.SeedMinLayers    = 3
seq.SeedSlopeMax     = 1.5    # dc/dz — ~56 deg, covers the GENIE angular range
seq.SeedSlopeBin     = 0.01   # 0.01 * 200 mm half-window = 2 mm < HoughBinSize
seq.SeedResidWindow  = 8.0    # mm — measured p90 of the 4-layer straight-line
                              # residual is 8.0 mm, and using the p90 rather
                              # than the median is deliberate: a seed that
                              # misses is free, a track that is never seeded is
                              # lost for good.

# Per-detector overrides. The scalars above are all in units that only mean
# something relative to a detector's LAYER PITCH, and the three differ by
# almost an order of magnitude:
#
#     SiTarget  11.0 mm/layer   166 layers   1826 mm deep
#     SiPad     16.8 mm/layer    22 layers    370 mm deep
#     MTC       80.0 mm/layer    42 layers   3360 mm deep
#
# so the scalar window of 4 layers is 240 mm in the MTC but only 33 mm in the
# SiTarget, where it cannot determine a slope at all. These pick a comparable
# window LENGTH (~230-240 mm) in each detector instead.
#
# CAVEAT, do not read these as equals: the MTC column is MEASURED (MC-truth
# straight-line residuals of 199 muons, see findSeeds()). The SiTarget and
# SiPad columns are derived from geometry and sensor pitch only — the same
# measurement could not be repeated for them because ddsim leaves
# CaloHitContribution::stepPosition at zero and the strip cell position
# degenerates to the strip centre in the unmeasured coordinate. They need the
# same scan the MTC window got, on a pipeline that actually populates them
# (muon/e_pi), before anyone quotes an efficiency from them.
#                              [SiTarget, SiPad,  MTC]
seq.SeedWindowLayersPerDet = [      22,     15,     4]   # ~231 / 235 / 240 mm
seq.SeedWindowStridePerDet = [      11,      7,     2]   # half a window
seq.SeedResidWindowPerDet  = [     2.0,    5.0,   8.0]   # mm; SiPad cell is 5.53
seq.HoughBinSizePerDet     = [     1.0,    5.0,   5.0]   # mm intercept bin
seq.SeedSlopeBinPerDet     = [   0.005,   0.01,  0.01]   # keep bin*halfWindow < intercept bin
seq.SeedCompatRadius = 8.0    # unused by the tracklet seeder, kept for load
seq.SeedStripPitch   = geo.sitarget_strip_pitch   # idem
seq.SeedMomentum     = 3.0    # GeV — still a prior: the lever arm short enough
                              # for a reliable association is too short to
                              # measure curvature against multiple scattering

# ── Segment finding ──────────────────────────────────────────────────────────
# Two iterations per detector. The first takes the clean tracks; its hits are
# then removed and the second sweeps what is left with looser cuts, which is
# what recovers a track masked by a better one.
seq.Chi2CutOff       = 70.0   # CKF per-surface chi2, first iteration
seq.NumMeasCutOff    = 1
seq.LooseChi2Scale   = 2.0    # second iteration: 140
seq.SegMaxChi2Tight  = 10.0   # segment chi2/ndf, first iteration
seq.SegMaxChi2Loose  = 20.0   # ... second
# Per-detector override of the two above: SiTarget and SiPad keep 10/20, the
# MTC gets 50/100. THIS IS THE LARGEST SINGLE EFFECT MEASURED ON GENIE.
#
# Why: the scalar SegMaxChi2Tight=10 was by far the most restrictive gate in
# the whole reconstruction. Of the segment candidates the cuts reject on a
# 500-event chunk, 4048 are MTC killed by chi2/ndf, versus 740 SiTarget and
# 384 SiPad -- and the rejected MTC candidates are not wild (chi2/ndf p10=14.5,
# median=32.6, p90=87.3). One scalar tuned for the precision detectors was
# throttling the MTC, which supplies ~77% of all segments in GENIE.
#
# Measured on 500 events (events with >=1 GlobalTracks / GlobalTracks):
#   10/20 everywhere (before)  : 139 (28.4%) / 279
#   50/100 everywhere          : 205 (41.9%) / 806
#   50/100 on MTC only (this)  : 190 (38.9%) / 787
# Loosening MTC alone keeps ~93% of the gain, and the fits are independent per
# detector so SiTarget/SiPad cannot be degraded by it. Loosening them too buys
# 15 more events but inflates their unspliced local segments from 4/25 to
# 27/81 -- cost without benefit on the detectors that fix the vertex position.
#
# THE TRADE-OFF IS REAL AND DELIBERATE: MTC segment purity (fraction of a
# segment's measurements belonging to its dominant MC particle) drops 0.78 ->
# 0.60 and hybrid segments rise 71.5% -> 81.4%. Accepted because the objective
# is that a track is reconstructed at all. If downstream analysis ever uses the
# TRACK MOMENTUM rather than just its existence and position, revisit this:
# at 0.60 purity, 40% of a segment's hits come from other particles, so the
# fitted momentum averages more than one trajectory.
# MTC TIGHTENED 50 -> 5. The 50 was set when the seeder could only produce
# beam-parallel seeds: a fit started pointing the wrong way has a large chi2 by
# construction, so the gate had to be opened to keep anything at all. With
# directed seeds that reason is gone, and the loose gate is now the main source
# of fake tracks. Measured on 197 GENIE events, MTC gate 50/25/15/10/5:
#   efficiency  64.0 / 62.4 / 61.9 / 60.9 / 60.4 %
#   tracks/event 7.19 / 6.30 / 5.01 / 4.00 / 2.13
# Cutting on fit quality is strictly better than cutting on length: the length
# route reaches 4.05 tracks/event at 52.8% efficiency, this one reaches 4.00 at
# 60.9% and 2.13 at 60.4%. SiTarget and SiPad keep 10/20 -- 100% of their
# accepted segments already sit below 20, so their gate is not the binding one.
#
# SIPAD TIGHTENED 10 -> 3 after the same MC-truth study, on 500 GENIE events.
# The SiPad fragments differently from the MTC -- its pieces are contiguous and
# complete (median completeness 1.00 over 6 of 22 layers) rather than
# overlapping -- but the chi2 gate separates them just as well, and far more
# cheaply than it did in the MTC:
#   gate 10 / 5 / 3 : segments per MC particle 2.45 / 2.25 / 1.72
#                     median purity            0.80 / 0.86 / 1.00
#                     global efficiency        55.1 / 54.4 / 54.0 %
# At 3 the median SiPad segment is built from ONE MC particle, for 1.1 points.
#
# SITARGET LEFT AT 10 deliberately. Tightening it to 5 costs 2 points of global
# efficiency and 31% of the events holding a SiTarget segment (52 -> 36), and
# buys purity 0.88 -> 0.92 on a sample of 68 MC particles -- inside the noise.
# SiTarget is already the healthy detector: 1.51 segments per particle against
# the SiPad's 2.45, with 26-measurement segments. Nothing to fix here, and the
# GENIE statistics could not settle it anyway; confirm on muon/e_pi if it ever
# matters.
# RETUNED 2026-09-11 after the loc0 sign fix (SeqCalibrator): every number
# above was measured with seed and measurements in mirror-image frames, so the
# MTC gates were really selecting near-axis tracks. Scan on 2000 GENIE events
# (MTC MinMeas x tight/loose gate; simulation/run_script/steer/scan_mtc.sub),
# with the metric that matters -- events whose MC muon is reconstructed as a
# muon-dominated segment (1117 events with a muon leaving hits):
#
#   MTC min / gates   ev>=1 track  tracks/ev  >5 tracks  muon eff
#      3   /  5,10       63.9%       2.30      16.2%      65.8%
#      3   /  8,16       65.2%       2.72      19.9%        --
#      5   /  5,10       53.8%       1.07       3.1%      61.2%
#      5   /  8,16       56.4%       1.37       6.2%      65.0%   <- chosen
#      6   /  8,16       52.0%       1.00       2.8%      63.4%
#      8   /  8,16       45.7%       0.76       1.4%      60.2%
#
# MinMeas 3 -> 5 in the MTC removes 8 points of "events with a track" but only
# 0.8 points of muon efficiency: what goes is 3-4 hit shower fragments (67% of
# all MTC segments at min 3, 40-46% MC-pure), not muons. The looser gate 8/16
# is what keeps the muons: with the frame right a pure 11-hit muon fits at
# chi2/ndf ~6 in 50 mm Fe per layer, and 5 was cutting it.
seq.SegMaxChi2TightPerDet = [10.0, 3.0, 8.0]
seq.SegMaxChi2LoosePerDet = [20.0, 6.0, 16.0]
# Scanned per detector on 400 GENIE events, one detector at a time with the
# other two held here, measuring efficiency AND MC-truth purity at every point.
# The three answers came out different, which is why this is a vector:
#
#   MTC 3 (and 3 is the floor the code allows: runCKFPass discards candidates
#   below 3 measurements). Raising it is STRICTLY DOMINATED -- 3/4/5/6 give
#   60.4/57.3/50.6/44.2% of events with a track while the median purity stays
#   0.67/0.67/0.60/0.67 and the fraction of segments built from >=2 particles
#   gets steadily WORSE, 72.3/74.7/75.8/75.4%. Length does not discriminate
#   quality in the MTC: the dirty population sits at 5-9 measurements, not in
#   the short tail, so cutting from below removes the clean stubs too.
#
#   SiTarget 5, raised from 4. Here the scan is monotone in the good direction:
#   4/5/6/8 give purity 0.88/0.92/0.92/0.94 and completeness 0.73/0.76/0.80/
#   0.82, paid for in events holding a SiTarget segment, 40/37/34/33. 5 is the
#   knee. CAVEAT: only 87 segments in 400 GENIE events carry this arm of the
#   scan, so it is indicative, not settled -- confirm it on muon/e_pi, which
#   actually populate the silicon.
#
#   SiPad 3, unchanged. Its two purity measures disagree in direction: the
#   dominant-particle fraction improves (0.76 -> 0.83 from 3 to 6) while the
#   fraction of segments containing >=2 particles gets worse (62.9 -> 67.4%).
#   That is the contamination changing shape, not going away, and it costs 13
#   of 69 events. Not a trade worth making.
# Left LOW on purpose. Raising it to 5 everywhere does cut the track count
# (8.79 -> 4.05 per event) but costs 11 points of efficiency, because a real
# track can start anywhere in the detector and a flat length floor kills the
# short and the late-starting ones indiscriminately. Cutting on FIT QUALITY
# instead reaches the same track count for 8 points more efficiency -- see
# SegMaxChi2TightPerDet below. Measured on 197 GENIE events.
# MTC 3 -> 5 on 2026-09-11, see the scan table above SegMaxChi2TightPerDet.
# SITARGET RETUNED 2026-09-11 (night), first scan with inclined tracks reachable
# (sign + navigator fixes). 2000 GENIE events, muon eff / muon events with
# exactly one global / muonless events with a global / SiTarget locals per MC
# particle:
#   min 5,  compl 0.01 (ref)      74.2% / 64.2% / 30.4% / 1.65
#   min 8                         73.9% / 65.5% / 27.9% / 1.52
#   min 12                        73.9% / 65.4% / 27.1% / 1.48
#   compl 0.5                     74.1% / 65.4% / 27.9% / 1.53
#   chi2 5/10                     73.1% / 65.1% / 24.6% / 1.41   (costs 1 pt of muon)
#   min 8 + compl 0.5             73.9% / 65.9% / 26.2% / 1.48   <- chosen
#   min 8 + compl 0.3 + chi2 5/10 73.0% / 65.6% / 22.8% / 1.30
seq.MinMeasPerSegment = [8, 4, 5]   # [SiTarget, SiPad, MTC]
# Duplicates counted in PHYSICAL HITS, not measurement indices. The per-seed
# masking already makes segments disjoint in measurement index, so a threshold
# on those can never fire; with AllPairs one U strip yields one measurement per
# V partner, so disjoint measurement sets routinely rest on shared hits.
# Measured on events 0-5: 148 pairs share 1 cellID, 18 share 2, 4 share 3.
seq.DuplicateSharedCellIDs = 1
# Fragments of one particle are disjoint by construction, so no duplicate
# filter can catch them. This rebuilds them: same operation as the
# inter-detector splice, applied within a detector.
seq.MergeSegments = True
# Fragments of one particle OVERLAP in z, they are not consecutive pieces: the
# per-seed masking makes the next seed rebuild the particle from the hits the
# first one left, interleaved along the same z range. With the strict
# end-to-end rule, 15012 of 15210 candidate pairs were rejected on the z gap
# and only 12 merges happened in 1000 events. Allowing overlap (3500 mm covers
# the deepest detector, so any two segments sharing z are tested) is the only
# change today that costs NO efficiency: 53.5% before and after, tracks/event
# 1.54 -> 1.42, merges 1 -> 79 per 389 events. It joins instead of discarding.
seq.MergeMaxOverlap = 3500.0
# The floor above is flat, and a flat floor is not neutral in where a track
# starts: a segment entering at MTC layer 39 of 42 can hold at most 4 layers of
# measurements, one entering at layer 0 has 42 ahead of it. GENIE vertices are
# spread over the full instrumented depth, so this asymmetry is systematic
# here, not a corner case. With this on, the floor is capped by the layers
# still available downstream of where the segment actually starts. It can only
# LOWER the floor, so a track crossing the whole detector is judged unchanged.
seq.ScaleMinMeasByDepth = True
# Completeness floor per detector (layers hit / layers spanned by the segment).
# Left permissive on purpose: MEASURED this session, cutting on it costs
# reconstructed events and buys nothing.
#
# MC truth does say the low-completeness SiTarget/SiPad segments are dirty
# (median 3 distinct MC particles at ~0.5 dominant purity), and [0.5, 0.5, 0.01]
# duly removes them -- but on the full 2000-event GENIE chunk that dropped
# events with >=1 GlobalTracks from 618 to 589 while the track count stayed
# flat (1204 -> 1206). Segment purity and track yield pull in OPPOSITE
# directions here: a dirty SiPad segment still BRIDGES a chain, and removing it
# forbids BOTH of its neighbour links, so the event is lost entirely.
# Impure hits inside a track are acceptable as long as the track is
# reconstructed at all -- that is the objective, not per-hit purity.
#
# For MTC the relation inverts on top of that: its lowest-completeness bin is
# its PUREST (median 1 particle, 53% single-particle, dominant purity 1.00) --
# real single-particle tracks with U/V pairing gaps, not junk.
seq.MinMeasFracPerSegment = [0.5, 0.01, 0.01]   # SiTarget 0.01 -> 0.5, 2026-09-11 night scan
seq.DuplicateOverlapFraction = 0.7

# ── Splicing ─────────────────────────────────────────────────────────────────
# Only neighbours: SiTarget<->SiPad and SiPad<->MTC. A full chain forms through
# both links; there is deliberately no SiTarget->MTC shortcut. The geometric
# window only proposes candidates — the combined refit decides.
seq.SpliceMaxDist       = 15.0   # mm at the inter-detector boundary
seq.SpliceMaxAngle      = 0.1    # rad
seq.SpliceAngleWeight   = 100.0  # mm^2/rad^2
seq.GlobalMaxChi2PerNdf = 10.0   # combined refit acceptance

# GENIE vertices are spread across the full ~5.6 m instrumented depth (unlike
# the single-vertex PG samples this tracker was first tuned on), so a
# genuine track that starts inside MTC or SiPad -- with nothing upstream to
# splice against, not a failed splice -- is expected to be common. Promote
# those to GlobalTracks (own single-detector Track::type bit) instead of
# losing them in the local-only collections. See ACTSSequentialTracker.cpp's
# PromoteUnspliced doc: only promotes when the neighbouring detector(s)
# found ZERO segments this event -- a segment whose neighbour had a real
# candidate that simply failed the splice window/refit is never promoted.
seq.PromoteUnspliced    = True
# Evening of 2026-09-11, validated on the 2000 GENIE events (muon efficiency /
# muon events with exactly one global / muonless events with a global):
#   production 17:35                     72.9% / 56.6% / 27.3%
#   + extension, proxdup, bounds, X&Y    71.9% / 61.0% / 28.4%
#   + backward fill (conditional)        74.2% / 64.2% / 30.4%   <- this
seq.ExtendSegments      = True

# ── Shower cleaning ────────────────────────────────────────────────────────
# HitPurgeWindow/MaxNeighbors: validated this session against the GENIE
# baseline (2000-event chunk) — [1.0, 8.0, 0.0] mm purges 4.0% of the
# measurement pool (up to 20-31% in individual dense-shower events), unlike
# the pitch-sized windows tried first ([0.075,5.5,0] / [0.15,5.5,0]) which
# purged ~0%/0.09%. Same values as the other 5 ACTSSequentialTracker
# pipelines (muon/e_pi/2_mu_ang/mu_pi/mu_pi_early). Note: even at this scale
# the purge only touches SiTarget+SiPad (15% of the GENIE pool) and measured
# ZERO net improvement in accepted-track yield on the global tracker, because
# MTC (85% of the pool, deliberately excluded — window[2]=0) is the real
# bottleneck for GENIE's dense hadronic showers. Kept anyway: it's the
# correct, proven-safe setting, just not expected to be the fix by itself.
# MTC PURGE ON since 2026-09-11 (was window 0 = off). Scan on 2000 GENIE
# events with the sign fix in place, MTC min-meas 5 / gates 8,16:
#   MTC purge      ev>=1 track  tracks/ev  =1 track  >5 tracks  muon eff
#   off               56.4%       1.37      32.5%      6.2%      65.0%
#   20 mm / 4         50.8%       0.76      36.0%      0.7%      64.0%
#   30 mm / 8         50.6%       0.76      36.1%      0.6%      64.5%
#   30 mm / 8 + scint 49.7%       0.72      36.4%      0.6%      65.4%
# Per hit the purge hardly discriminates (it removes ~9% of the muon's own MTC
# hits against ~17% of the rest, the muon leaves the vertex INSIDE the shower),
# but per track it does: the CKF rides through the holes while the 3-5 hit
# shower fragments lose the cluster they were built from. Muon efficiency is
# flat within +-1%; what goes is the multiplicity tail (>5 tracks: 6.2 -> 0.6%).
# SiPad 8 -> 25 mm and MTC 30/8 -> 40/6 on 2026-09-11 (evening), after the
# navigator fix. EM showers in the SiPad (nu_e CC: one MC particle, ~12 hits per
# layer over 11 layers) were yielding 2-11 "tracks" of 3-4 hits each (events
# 36, 159, 1338) because 8 mm only reaches the adjacent pad. Scan on 2000
# events, muon efficiency / globals in muonless events:
#   SiPad 8/4,  min 3 (ref)      73.9% / 33.2%
#   SiPad 25/4, min 3            73.9% / 32.3%
#   SiPad 25/4, min 4            73.7% / 29.6%
#   + MTC 40/6                   72.9% / 27.3%   <- chosen (=1 global: 56.6%)
#   + MTC 40/4                   72.6% / 26.2%
#   + MTC 40/6, MTC min 6        71.5% / 24.9%
seq.HitPurgeWindow        = [1.0, 25.0, 40.0]  # mm, per [SiTarget, SiPad, MTC]
seq.HitPurgeMaxNeighbors  = [8, 4, 6]
# Border mask (2026-09-11): also drop the survivors within the window of a
# purged core hit -- the DBSCAN border -- so a blob goes whole. Same 2000
# events, on top of the purge above: muon eff 65.4 -> 66.2%, MTC pool 309 -> 272
# measurements/event, >5 tracks/event 0.6 -> 0.2%. NOTE: k4run cannot switch
# this off from the command line -- `--ACTSSequentialTracker.HitPurgeMaskBorder
# False` still yields True (non-empty string); edit the steering instead.
seq.HitPurgeMaskBorder    = True
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
