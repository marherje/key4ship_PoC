// ---------------------------------------------------------------------------
// ACTSSequentialTracker — sequential per-subdetector tracking + splicing.
//
// The alternative to ACTSProtoTracker's single global CKF. Instead of running
// one combinatorial fit across SiTarget + SiPad + MTC at once (and inflating
// the SiPad's information weight so its unambiguous 2D points anchor the
// trajectory in y), this algorithm:
//
//   1. finds track SEGMENTS independently inside each subdetector, with a CKF
//      that navigates only that subdetector's surfaces — two search iterations
//      per detector, each seed refined by a re-seed pass, with the hits of the
//      first iteration's winners removed before the second;
//   2. SPLICES segments of ADJACENT detectors (SiTarget<->SiPad, SiPad<->MTC)
//      whose ends extrapolate onto each other, and accepts the splice only if
//      a combined Kalman refit over the union of their hits fits;
//   3. writes the spliced chains to GlobalTracks and every segment that was
//      not absorbed into a chain to its own detector's local collection.
//
// There are no per-detector information weights and no annealing here, by
// design: those existed to resolve the SiTarget strip / MTC U-V ambiguity
// against the SiPad inside ONE fit. Per-detector fits do not have that
// ambiguity to resolve, so every measurement enters at its nominal variance.
//
// A track may be born anywhere — at the front, middle or back of the SiTarget,
// inside the SiPad, or inside the MTC. Nothing here requires a SiTarget segment
// to exist: each detector seeds itself, and a single-detector segment that
// splices with nothing is a legitimate local track, not a failure.
//
// ACTSProtoTracker is untouched and still selectable; the two live side by side
// and the job4 steering picks one. Everything below is in an anonymous
// namespace because both translation units define these helper types and both
// are linked into the same module (libSND_reco).
// ---------------------------------------------------------------------------

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/DataHandle.h"

// DD4hep segmentation
#include "DDSegmentation/BitFieldCoder.h"

// edm4hep
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MutableTrack.h"

// ACTS
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"

#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "Acts/EventData/Types.hpp"

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"

// SND geometry service
#include "ISNDGeoSvc.h"

// Standard
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Measurement and source-link plumbing.
//
// Same shapes as ACTSProtoTracker's, deliberately re-declared here rather than
// shared through a header: the two trackers are independent by construction, so
// tuning one can never silently move the other. The anonymous namespace keeps
// the duplicate names from colliding at link time.
// ---------------------------------------------------------------------------

struct SeqMeasurement {
  const Acts::Surface* surface     = nullptr;
  double               localCoord  = 0.0;
  double               localCoord2 = 0.0;
  double               variance    = 0.0;
  double               variance2   = 0.0;
  bool                 is2D        = false;
  int                  detectorID  = -1;   // 0=SiTarget, 1=SiPad, 2=MTC
  int                  plane       = -1;
  float                time        = 0.0f;
  float                eDep        = 0.0f;

  SeqMeasurement(const Acts::Surface* sf, double lc, double lc2,
                 double var, double var2, bool twod, int detID, int pl,
                 float t, float e)
      : surface(sf), localCoord(lc), localCoord2(lc2),
        variance(var), variance2(var2), is2D(twod),
        detectorID(detID), plane(pl), time(t), eDep(e) {}
};

struct SeqSourceLink {
  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
  std::size_t index = 0;
  void setGeometryId(Acts::GeometryIdentifier gid) { m_geometryId = gid; }
private:
  Acts::GeometryIdentifier m_geometryId;
};

// Range of source links on one surface. The vector must be sorted by geometryId.
struct SeqSourceLinkAccessor {
  const std::vector<Acts::SourceLink>* slinks = nullptr;
  std::pair<std::vector<Acts::SourceLink>::const_iterator,
            std::vector<Acts::SourceLink>::const_iterator>
  operator()(const Acts::Surface& surface) const {
    const auto geoId = surface.geometryId();
    auto lo = std::lower_bound(
        slinks->begin(), slinks->end(), geoId,
        [](const Acts::SourceLink& sl, const Acts::GeometryIdentifier& id) {
          return sl.get<SeqSourceLink>().geometryId() < id;
        });
    auto hi = std::upper_bound(
        lo, slinks->end(), geoId,
        [](const Acts::GeometryIdentifier& id, const Acts::SourceLink& sl) {
          return id < sl.get<SeqSourceLink>().geometryId();
        });
    return {lo, hi};
  }
};

// Resolves a source link back to its surface; required by KalmanFitterExtensions.
struct SeqSurfaceAccessor {
  const std::vector<SeqMeasurement>* meas = nullptr;
  const Acts::Surface* operator()(const Acts::SourceLink& sl) const {
    return (*meas)[sl.get<SeqSourceLink>().index].surface;
  }
};

// Wraps DirectNavigator and injects the surface list at makeState() time, so
// the CKF's setPlainOptions() (which only copies the NavigatorPlainOptions base
// fields) cannot erase it. Here the list is ONE detector's surfaces — that is
// what makes the per-detector search independent.
struct SeqFixedNavigator {
  std::vector<const Acts::Surface*> surfaces;

  struct Options : public Acts::NavigatorPlainOptions {
    explicit Options(const Acts::GeometryContext& gctx)
        : Acts::NavigatorPlainOptions(gctx) {}
    void setPlainOptions(const Acts::NavigatorPlainOptions& opts) {
      static_cast<Acts::NavigatorPlainOptions&>(*this) = opts;
    }
  };

  using State = Acts::DirectNavigator::State;

  State makeState(const Options& opts) const {
    Acts::DirectNavigator::Options dirOpts(opts.geoContext);
    static_cast<Acts::NavigatorPlainOptions&>(dirOpts) =
        static_cast<const Acts::NavigatorPlainOptions&>(opts);
    dirOpts.surfaces = surfaces;
    Acts::DirectNavigator inner;
    return inner.makeState(dirOpts);
  }

  Acts::Result<void> initialize(State& state, const Acts::Vector3& pos,
                                const Acts::Vector3& dir,
                                Acts::Direction propDir) const {
    Acts::DirectNavigator inner;
    return inner.initialize(state, pos, dir, propDir);
  }
  Acts::NavigationTarget nextTarget(State& state, const Acts::Vector3& pos,
                                    const Acts::Vector3& dir) const {
    Acts::DirectNavigator inner;
    return inner.nextTarget(state, pos, dir);
  }
  bool checkTargetValid(const State& state, const Acts::Vector3& pos,
                        const Acts::Vector3& dir) const {
    Acts::DirectNavigator inner;
    return inner.checkTargetValid(state, pos, dir);
  }
  void handleSurfaceReached(State& state, const Acts::Vector3& pos,
                            const Acts::Vector3& dir,
                            const Acts::Surface& sf) const {
    Acts::DirectNavigator inner;
    inner.handleSurfaceReached(state, pos, dir, sf);
  }
  const Acts::Surface* currentSurface(const State& s) const { return s.currentSurface; }
  const Acts::TrackingVolume* currentVolume(const State&) const { return nullptr; }
  const Acts::IVolumeMaterial* currentVolumeMaterial(const State&) const { return nullptr; }
  const Acts::Surface* startSurface(const State& s) const { return s.options.startSurface; }
  const Acts::Surface* targetSurface(const State& s) const { return s.options.targetSurface; }
  bool endOfWorldReached(State&) const { return false; }
  bool navigationBreak(const State& s) const { return s.navigationBreak; }
};

using SeqStepper       = Acts::EigenStepper<>;
using SeqCKFPropagator = Acts::Propagator<SeqStepper, SeqFixedNavigator>;

using SeqTrackContainer = Acts::TrackContainer<
    Acts::VectorTrackContainer,
    Acts::VectorMultiTrajectory,
    std::shared_ptr>;

using SeqCKF = Acts::CombinatorialKalmanFilter<SeqCKFPropagator, SeqTrackContainer>;

// Refit fitter: the surface-sequence fit() overload needs a plain
// DirectNavigator, since the sequence is handed over per fit() call. Used both
// for the per-segment refit (that detector's surfaces) and for the combined
// refit of a spliced chain (all surfaces).
using SeqKFPropagator = Acts::Propagator<SeqStepper, Acts::DirectNavigator>;
using SeqKF           = Acts::KalmanFitter<SeqKFPropagator, Acts::VectorMultiTrajectory>;

// By inside the registered iron slabs, exactly zero everywhere else — so the
// same provider is correct for all three detectors and no special-casing is
// needed for the field-free SiTarget/SiPad region.
class SeqIronSlabBField : public Acts::MagneticFieldProvider {
public:
  struct Cache { explicit Cache(const Acts::MagneticFieldContext&) {} };
  struct Slab { double xlo, xhi, ylo, yhi, zlo, zhi, by; };

  explicit SeqIronSlabBField(std::vector<Slab> slabs) : m_slabs(std::move(slabs)) {}

  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override {
    return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& pos,
      Acts::MagneticFieldProvider::Cache&) const override {
    for (const auto& s : m_slabs) {
      if (pos.x() >= s.xlo && pos.x() <= s.xhi &&
          pos.y() >= s.ylo && pos.y() <= s.yhi &&
          pos.z() >= s.zlo && pos.z() <= s.zhi) {
        return Acts::Result<Acts::Vector3>::success(Acts::Vector3(0.0, s.by, 0.0));
      }
    }
    return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
  }

private:
  std::vector<Slab> m_slabs;
};

// Measurement calibrator. Note what is NOT here: the per-detector weight
// division of ACTSProtoTracker's calibrator. Every measurement enters at its
// nominal variance, always — see the file header for why the weights are gone.
struct SeqCalibrator {
  const std::vector<SeqMeasurement>* meas = nullptr;

  void operator()(const Acts::GeometryContext& /*gctx*/,
                  const Acts::CalibrationContext& /*cctx*/,
                  const Acts::SourceLink& sl,
                  Acts::VectorMultiTrajectory::TrackStateProxy ts) const {
    const auto& ssl = sl.get<SeqSourceLink>();
    const auto& m   = (*meas)[ssl.index];
    // TrackStateCreator leaves this to the calibrator (ACTS convention).
    // Without it the hit fingerprints — and therefore the frozen-hit refit,
    // the duplicate filter and the whole splice stage — come out empty.
    ts.setUncalibratedSourceLink(Acts::SourceLink{sl});

    if (m.is2D) {
      // SiPad (pixel) or MTC combined U+V surface: (loc0, loc1).
      constexpr std::array<Acts::BoundIndices, 2> indices = {
          Acts::eBoundLoc0, Acts::eBoundLoc1};
      ts.allocateCalibrated(2);
      ts.template calibrated<2>() = Acts::ActsVector<2>(m.localCoord, m.localCoord2);
      ts.template calibratedCovariance<2>() =
          Acts::ActsSquareMatrix<2>{{m.variance, 0.0}, {0.0, m.variance2}};
      ts.setProjectorSubspaceIndices(indices);
      return;
    }
    if (m.detectorID == 2) {
      // MTC SciFi unpaired: 1D on the stereo-rotated U/V surface -> eBoundLoc0.
      ts.allocateCalibrated(1);
      ts.template calibrated<1>() = Acts::ActsVector<1>(m.localCoord);
      ts.template calibratedCovariance<1>() =
          Acts::ActsSquareMatrix<1>{{m.variance}};
      constexpr std::array<Acts::BoundIndices, 1> mtcIdx = {Acts::eBoundLoc0};
      ts.setProjectorSubspaceIndices(mtcIdx);
      return;
    }
    // SiTarget strip: plane=0 (StripX) -> eBoundLoc0, plane=1 (StripY) -> eBoundLoc1.
    ts.allocateCalibrated(1);
    ts.template calibrated<1>() = Acts::ActsVector<1>(m.localCoord);
    ts.template calibratedCovariance<1>() = Acts::ActsSquareMatrix<1>{{m.variance}};
    if (m.plane == 1) {
      constexpr std::array<Acts::BoundIndices, 1> indices = {Acts::eBoundLoc1};
      ts.setProjectorSubspaceIndices(indices);
    } else {
      constexpr std::array<Acts::BoundIndices, 1> indices = {Acts::eBoundLoc0};
      ts.setProjectorSubspaceIndices(indices);
    }
  }
};

// One end of a segment: the fitted state at its most upstream or most
// downstream measurement, in ACTS global coordinates (x = beam).
struct SeqEndState {
  Acts::Vector3 pos{Acts::Vector3::Zero()};
  Acts::Vector3 dir{Acts::Vector3::UnitX()};
  double        qOverP = 0.0;
};

// A track found inside ONE detector.
struct SeqSegment {
  int                   det       = -1;
  int                   iteration = 0;
  Acts::TrackIndexType  idx       = 0;    // into the event's track container
  std::set<std::size_t> fp;               // global measurement indices
  double                chi2      = 0.0;
  int                   ndf       = 1;
  SeqEndState           front;            // most upstream fitted state
  SeqEndState           back;             // most downstream fitted state
  double                seedX     = 0.0;  // DD4hep transverse X of its seed
  double                seedY     = 0.0;
  bool                  spliced   = false;  // absorbed into a global track
};

// A candidate splice between an upstream and a downstream segment.
struct SeqLink {
  std::size_t up = 0, dn = 0;   // indices into the segment vector
  double dist  = 0.0;           // transverse mismatch at the boundary [mm]
  double angle = 0.0;           // 3D angle between the two directions [rad]
  double score = 0.0;           // dist^2 + w * angle^2 — lower is better
};

constexpr int kHelixParams = 5;  // loc0, loc1, phi, theta, q/p

}  // namespace

// ---------------------------------------------------------------------------
// ACTSSequentialTracker
// ---------------------------------------------------------------------------

class ACTSSequentialTracker : public Gaudi::Algorithm {
public:
  ACTSSequentialTracker(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  // ---- Inputs -------------------------------------------------------------
  Gaudi::Property<std::string> m_inputSiTarget{
      this, "InputSiTarget", "SiTargetMeasurements",
      "SiTarget TrackerHit3DCollection from SiTargetMeasConverter"};
  Gaudi::Property<std::string> m_inputSiPad{
      this, "InputSiPad", "SiPadMeasurements",
      "SiPad TrackerHit3DCollection from SiPadMeasConverter"};
  Gaudi::Property<std::string> m_inputMTC{
      this, "InputMTC", "MTCSciFiMeasurements",
      "MTC SciFi TrackerHit3DCollection from MTCSciFiMeasConverter"};

  // ---- Outputs ------------------------------------------------------------
  Gaudi::Property<std::string> m_outSiTarget{
      this, "OutputSiTargetTracks", "SiTargetTracks",
      "Local SiTarget segments — those NOT absorbed into a global track"};
  Gaudi::Property<std::string> m_outSiPad{
      this, "OutputSiPadTracks", "SiPadTracks",
      "Local SiPad segments — those NOT absorbed into a global track"};
  Gaudi::Property<std::string> m_outMTC{
      this, "OutputMTCTracks", "MTCTracks",
      "Local MTC segments — those NOT absorbed into a global track"};
  Gaudi::Property<std::string> m_outGlobal{
      this, "OutputGlobalTracks", "GlobalTracks",
      "Spliced tracks (2 or 3 segments) after the combined Kalman refit. "
      "Track::type is a composition bitmask: 1=SiTarget, 2=SiPad, 4=MTC, "
      "so 3=SiT+SiP, 6=SiP+MTC, 7=full chain."};

  // ---- cellID decoding ----------------------------------------------------
  Gaudi::Property<std::string> m_siTargetBitFieldStr{
      this, "SiTargetBitField",
      "system:8,layer:8,slice:4,plane:1,column:2,row:2,strip:14",
      "BitField string for SiTarget cellID decoding. The layer and plane are "
      "decoded from the cellID and fed to ISNDGeoSvc::surfaceByAddress, which "
      "is exact — unlike a nearest-z scan, which needs a correct surface "
      "partition to even be well posed."};
  Gaudi::Property<std::string> m_mtcBitFieldStr{
      this, "MTCBitField",
      "system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9",
      "BitField string for MTC cellID decoding"};
  Gaudi::Property<double> m_mtcStereoAngle{
      this, "MTCStereoAngle", 5.0, "MTC SciFi stereo angle [degrees]"};

  // ---- Magnetic field -----------------------------------------------------
  Gaudi::Property<double> m_bFieldX{this, "BFieldX", 0.0, "BField X [T]"};
  Gaudi::Property<double> m_bFieldY{this, "BFieldY", 0.0, "BField Y [T]"};
  Gaudi::Property<double> m_bFieldZ{this, "BFieldZ", 0.0, "BField Z [T]"};
  Gaudi::Property<std::vector<double>> m_ironFieldRanges{
      this, "IronFieldRanges", {},
      "Per-slab field: [xlo,xhi, ylo,yhi, zlo,zhi, by] x N slabs in ACTS coords [mm, T]"};

  // ---- Propagation --------------------------------------------------------
  Gaudi::Property<double> m_seedMomentum{
      this, "SeedMomentum", 10.0, "Seed momentum magnitude [GeV]"};

  // The seed covariance is much more important here than in a global fit, and
  // for a reason worth stating. The Hough seeder votes in a purely transverse
  // (x, y) histogram, i.e. it already assumes a beam-parallel track, and it
  // seeds from a crossing found in THIS detector — so the seed position is
  // good. Left as loose as a global fit's seed (sigma_loc = 10 mm, sigma_angle
  // = 1 rad), the prior is so uninformative that over the SiTarget's ~2 m the
  // CKF can wander from one track onto another: with two tracks 114 mm apart it
  // picks track A's StripX and track B's StripY and fits the hybrid. The
  // SiTarget cannot resolve that pairing on its own — 1D strips, one
  // measurement per surface — so the prior has to do it.
  Gaudi::Property<double> m_seedVarLoc{
      this, "SeedVarLoc", 100.0,
      "Seed variance on loc0/loc1 [mm^2]. 100 = 10 mm, about two Hough bins; "
      "measured to work better than a tighter 25 on two-track events."};
  Gaudi::Property<double> m_seedVarAngle{
      this, "SeedVarAngle", 0.01,
      "Seed variance on phi/theta [rad^2]. 0.01 = 0.1 rad, which is wide "
      "compared to the tracks the transverse Hough can seed at all, and narrow "
      "enough that a fit cannot cross to a neighbouring track."};
  Gaudi::Property<double> m_seedVarQOverP{
      this, "SeedVarQOverP", 0.04, "Seed variance on q/p [(GeV^-1)^2]"};
  Gaudi::Property<int> m_maxPropSteps{
      this, "MaxPropSteps", 100000, "Maximum propagation steps per fit"};
  Gaudi::Property<double> m_maxStepSize{
      this, "MaxStepSize", 100.0, "Maximum single step size [mm]"};

  // ---- Seeding (Hough, per detector) --------------------------------------
  Gaudi::Property<int> m_maxSeeds{
      this, "MaxSeeds", 5, "Maximum number of Hough seeds per detector per iteration"};
  Gaudi::Property<double> m_houghBinSize{
      this, "HoughBinSize", 5.0, "Hough histogram bin size [mm]"};
  Gaudi::Property<double> m_houghHalfSize{
      this, "HoughHalfSize", 200.0, "Hough histogram transverse half-size [mm]"};
  Gaudi::Property<int> m_houghMinVotes{
      this, "HoughMinVotes", 3, "Minimum votes for a Hough peak to seed"};
  Gaudi::Property<double> m_seedCompatRadius{
      this, "SeedCompatRadius", 10.0, "Radius [mm] for seed centroid refinement"};
  Gaudi::Property<double> m_stripPitch{
      this, "SeedStripPitch", 0.0755, "Strip pitch [mm] for seed refinement"};

  // ---- Shower cleaning (carried over unchanged) ---------------------------
  Gaudi::Property<double> m_houghMaxMultiplicity{
      this, "HoughMaxMultiplicity", 1e9,
      "Maximum crossing multiplicity per station for a Hough peak to be a "
      "track candidate. Peaks above this are classified as showers and "
      "skipped. 1e9 disables."};
  Gaudi::Property<double> m_isolationWindow{
      this, "IsolationWindow", 0.0,
      "2D distance [mm] for the seed-level crossing isolation filter. A "
      "crossing is kept only if fewer than IsolationMaxNeighbors other "
      "crossings of the same station lie within it. 0.0 disables."};
  Gaudi::Property<int> m_isolationMaxNeighbors{
      this, "IsolationMaxNeighbors", 2,
      "Neighbour budget for the isolation filter. Muon: 0. Shower: hundreds."};
  Gaudi::Property<std::vector<double>> m_hitPurgeWindow{
      this, "HitPurgeWindow", {0.0, 0.0, 0.0},
      "Distance [mm] per detector [SiTarget, SiPad, MTC] for the density-based "
      "shower-hit purge of the measurement pool, applied per surface before "
      "seeding and fitting. 0.0 disables that detector."};
  Gaudi::Property<std::vector<int>> m_hitPurgeMaxNeighbors{
      this, "HitPurgeMaxNeighbors", {8, 4, 0},
      "Same-surface neighbour budget within HitPurgeWindow, per detector. "
      "MIP + delta rays: a few. Shower core: tens."};

  // ---- Per-detector segment finding ---------------------------------------
  Gaudi::Property<double> m_chi2CutOff{
      this, "Chi2CutOff", 70.0,
      "CKF MeasurementSelector: maximum local chi2 to attach a measurement to "
      "a surface, in the FIRST (tight) iteration."};
  Gaudi::Property<int> m_numMeasCutOff{
      this, "NumMeasCutOff", 1,
      "CKF MeasurementSelector: maximum measurements accepted per surface."};
  Gaudi::Property<double> m_looseChi2Scale{
      this, "LooseChi2Scale", 2.0,
      "Multiplier applied to Chi2CutOff in the SECOND iteration. The second "
      "pass runs on the hits the first one did not take, so it may be more "
      "permissive without competing with the tracks already found."};
  Gaudi::Property<double> m_segMaxChi2Tight{
      this, "SegMaxChi2Tight", 10.0,
      "Maximum chi2/ndf for a segment found in the first iteration."};
  Gaudi::Property<double> m_segMaxChi2Loose{
      this, "SegMaxChi2Loose", 20.0,
      "Maximum chi2/ndf for a segment found in the second iteration."};
  Gaudi::Property<std::vector<int>> m_minMeasPerSegment{
      this, "MinMeasPerSegment", {4, 3, 4},
      "Minimum measurements for a segment, per detector [SiTarget, SiPad, MTC]. "
      "SiPad hits are 2D, so 3 of them already carry 6 degrees of freedom."};
  Gaudi::Property<double> m_dupOverlapFraction{
      this, "DuplicateOverlapFraction", 0.7,
      "Two segments of one detector are duplicates when they share more than "
      "this fraction of the smaller one's hits; the worse chi2/ndf is dropped."};

  // ---- Splicing -----------------------------------------------------------
  Gaudi::Property<double> m_spliceMaxDist{
      this, "SpliceMaxDist", 15.0,
      "Maximum transverse distance [mm] at the inter-detector boundary between "
      "the straight-line extrapolations of an upstream segment's downstream end "
      "and a downstream segment's upstream end."};
  Gaudi::Property<double> m_spliceMaxAngle{
      this, "SpliceMaxAngle", 0.1,
      "Maximum 3D angle [rad] between those two segment directions."};
  Gaudi::Property<double> m_spliceAngleWeight{
      this, "SpliceAngleWeight", 100.0,
      "Weight [mm^2/rad^2] of the angular term in the splice score "
      "(score = dist^2 + w * angle^2); lower score = better match."};
  Gaudi::Property<double> m_globalMaxChi2PerNdf{
      this, "GlobalMaxChi2PerNdf", 10.0,
      "Maximum chi2/ndf of the combined Kalman refit for a chain to be "
      "accepted as a global track. This, not the geometric window, is what "
      "actually decides whether a splice is real."};

  // ---- Hough seeding ------------------------------------------------------
  struct SeedCandidate {
    double x;
    double y;
    double z_start;
    int    nVotes;
    double multiplicity;
  };

  Gaudi::Property<double> m_ghostResolveWindow{
      this, "GhostResolveWindow", 0.0,
      "Transverse distance [mm] within which a SiTarget StripX x StripY "
      "crossing is 'confirmed' by an unambiguous SiPad 2D point. The "
      "strip-swapped ghosts of two confirmed crossings are then dropped from "
      "SiTarget seeding. With N tracks the SiTarget alone produces N^2 "
      "crossings of which N are real, and a ghost seed can build a hybrid "
      "segment out of strips belonging to two different tracks — which then "
      "extrapolates onto neither. Crossings that no SiPad point confirms are "
      "never dropped, so a track stopping inside the SiTarget keeps its "
      "crossing. This is a SEEDING aid only: it changes no variance and no "
      "fit. 0.0 disables."};

  // `pool` is the detector's own (possibly hit-masked) measurements; `allMeas`
  // is the full event, needed only for the SiPad points that resolve SiTarget
  // ghosts.
  std::vector<SeedCandidate> findSeeds(
      const std::vector<SeqMeasurement>& pool,
      const std::vector<SeqMeasurement>& allMeas,
      const Acts::GeometryContext& gctx,
      int detMask) const;

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siTargetHandle, m_siPadHandle, m_mtcHandle;
  mutable std::array<
      std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackCollection>>, 3>
      m_localHandles;                         // index = detID
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackCollection>>
      m_globalHandle;

  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_siTargetBitField;
  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_mtcBitField;

  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  mutable std::atomic<long long> m_eventCount{0};
  mutable Acts::MagneticFieldContext m_mctx;
  mutable Acts::CalibrationContext   m_cctx;
};

// ---------------------------------------------------------------------------
// findSeeds() — Hough transform over unambiguous 2D points of ONE detector.
//
// detMask selects the source: 0 = SiTarget StripX x StripY crossings,
// 1 = SiPad 2D hits, 2 = MTC U x V stereo crossings. Unlike the global tracker
// this is never called with -1: every seed here belongs to exactly one
// detector, which is the whole point of the sequential design.
// ---------------------------------------------------------------------------

std::vector<ACTSSequentialTracker::SeedCandidate>
ACTSSequentialTracker::findSeeds(
    const std::vector<SeqMeasurement>& measurements,
    const std::vector<SeqMeasurement>& allMeas,
    const Acts::GeometryContext& gctx,
    int detMask) const
{
  std::vector<SeedCandidate> seeds;
  if (measurements.empty()) return seeds;

  // SiPad reference points for SiTarget ghost resolution. Taken from the whole
  // event, not from the SiTarget pool, which by construction has none.
  const double ghostWin  = m_ghostResolveWindow.value();
  const double ghostWin2 = ghostWin * ghostWin;
  std::vector<std::pair<double, double>> sipadXY;
  if (detMask == 0 && ghostWin > 0.0) {
    for (const auto& m : allMeas)
      if (m.is2D && m.detectorID == 1) sipadXY.emplace_back(m.localCoord, m.localCoord2);
  }

  struct Point2D { double x, y, z; int weight; };
  std::vector<Point2D> points2D;

  const double isolWin      = m_isolationWindow.value();
  const int    isolMaxNeigh = m_isolationMaxNeighbors.value();
  const bool   doIsolation  = (isolWin > 0.0);
  const double isolWin2     = isolWin * isolWin;

  const double stationTolerance = 10.0;  // mm — max z gap grouping one station

  // --- Source A: SiPad 2D hits, with per-layer position isolation ----------
  if (detMask == 1) {
    std::map<const Acts::Surface*, std::vector<const SeqMeasurement*>> bySurface;
    for (const auto& m : measurements) {
      if (!m.is2D || m.detectorID != 1) continue;
      bySurface[m.surface].push_back(&m);
    }
    for (const auto& [surf, layerHits] : bySurface) {
      const double beamZ = surf->center(gctx).x();
      for (std::size_t i = 0; i < layerHits.size(); ++i) {
        const auto* mp = layerHits[i];
        if (doIsolation) {
          int nNeigh = 0;
          for (std::size_t j = 0; j < layerHits.size(); ++j) {
            if (j == i) continue;
            const double dx = layerHits[j]->localCoord  - mp->localCoord;
            const double dy = layerHits[j]->localCoord2 - mp->localCoord2;
            if (dx * dx + dy * dy < isolWin2) ++nNeigh;
            if (nNeigh > isolMaxNeigh) break;
          }
          if (nNeigh > isolMaxNeigh) continue;
        }
        points2D.push_back({mp->localCoord, mp->localCoord2, beamZ, 1});
      }
    }
  }

  // --- Source B: SiTarget StripX x StripY crossings ------------------------
  if (detMask == 0) {
    std::map<int, std::vector<const SeqMeasurement*>> stationGroups;
    for (const auto& m : measurements) {
      if (m.is2D || m.detectorID != 0) continue;
      const int key = static_cast<int>(
          std::round(m.surface->center(gctx).x() / stationTolerance));
      stationGroups[key].push_back(&m);
    }
    for (const auto& [key, stationMeas] : stationGroups) {
      std::vector<const SeqMeasurement*> stripX, stripY;
      for (const auto* m : stationMeas) {
        if (m->plane == 0)      stripX.push_back(m);
        else if (m->plane == 1) stripY.push_back(m);
      }
      if (stripX.empty() || stripY.empty()) continue;

      // The originating strips are kept so ghost resolution can tell a real
      // crossing from the strip-swapped combination of two real ones.
      struct Crossing2D {
        double x, y, z;
        const SeqMeasurement* mx;
        const SeqMeasurement* my;
      };
      std::vector<Crossing2D> crossings;
      for (const auto* mx : stripX) {
        for (const auto* my : stripY) {
          const double zX = mx->surface->center(gctx).x();
          const double zY = my->surface->center(gctx).x();
          crossings.push_back({mx->localCoord, my->localCoord,
                               0.5 * (zX + zY), mx, my});
        }
      }
      // Hard cap: a station this busy is a shower core, not a set of tracks.
      const int maxCrossingsPerStation = 500;
      if (static_cast<int>(crossings.size()) > maxCrossingsPerStation) continue;

      // Confirm pass: a crossing lying on a SiPad point is real, and its two
      // strips are marked. A crossing that is NOT confirmed but whose BOTH
      // strips are confirmed is the swapped-strip ghost of two real tracks and
      // is dropped below. Nothing else is ever dropped.
      std::set<const SeqMeasurement*> confirmedX, confirmedY;
      std::vector<char> isConfirmed(crossings.size(), 0);
      if (!sipadXY.empty()) {
        for (std::size_t ci = 0; ci < crossings.size(); ++ci) {
          const auto& c = crossings[ci];
          for (const auto& sp : sipadXY) {
            const double dx = sp.first - c.x, dy = sp.second - c.y;
            if (dx * dx + dy * dy <= ghostWin2) {
              isConfirmed[ci] = 1;
              confirmedX.insert(c.mx);
              confirmedY.insert(c.my);
              break;
            }
          }
        }
      }

      for (std::size_t ci = 0; ci < crossings.size(); ++ci) {
        const auto& c = crossings[ci];
        if (!isConfirmed[ci] &&
            confirmedX.count(c.mx) && confirmedY.count(c.my)) {
          continue;  // strip-swapped ghost of two confirmed crossings
        }
        if (doIsolation) {
          int nNeigh = 0;
          for (std::size_t cj = 0; cj < crossings.size(); ++cj) {
            if (cj == ci) continue;
            const double dx = crossings[cj].x - c.x;
            const double dy = crossings[cj].y - c.y;
            if (dx * dx + dy * dy < isolWin2) ++nNeigh;
            if (nNeigh > isolMaxNeigh) break;
          }
          if (nNeigh > isolMaxNeigh) continue;
        }
        points2D.push_back({c.x, c.y, c.z, 2});
      }
    }
  }

  // --- Source C: MTC ------------------------------------------------------
  // Two shapes reach us, depending on MTCSciFiMeasConverter.PairMethod:
  //   * "AllPairs" (the configured default) already pairs U with V and emits
  //     ONE 2D measurement per combination on the flat combined surface
  //     (plane=3). Those are unambiguous (x, y) points — used directly, like
  //     SiPad hits. Missing this is what left the MTC unseeded: with AllPairs
  //     there are no 1D planes left to cross, so a crossing-only source finds
  //     nothing at all and the detector silently produces no segments.
  //   * "None" (legacy) emits 1D hits on the U/V stereo surfaces, which still
  //     have to be crossed here.
  if (detMask == 2) {
    const double tanStereo = std::tan(m_mtcStereoAngle.value() * M_PI / 180.0);
    const int maxCrossingsPerLayer = 200;

    std::map<int, std::vector<const SeqMeasurement*>> mtcGroups;
    for (const auto& m : measurements) {
      if (m.detectorID != 2) continue;
      if (m.is2D) {
        // Already a resolved U x V point; weight 2, since two strips made it.
        points2D.push_back({m.localCoord, m.localCoord2,
                            m.surface->center(gctx).x(), 2});
        continue;
      }
      const int key = static_cast<int>(std::round(m.surface->center(gctx).x() / 5.0));
      mtcGroups[key].push_back(&m);
    }

    for (const auto& [key, meas] : mtcGroups) {
      std::vector<const SeqMeasurement*> uPl, vPl;
      for (const auto* m : meas) {
        if      (m->plane == 0) uPl.push_back(m);
        else if (m->plane == 1) vPl.push_back(m);
      }
      if (uPl.empty() || vPl.empty()) continue;
      if (static_cast<int>(uPl.size() * vPl.size()) > maxCrossingsPerLayer) continue;
      for (const auto* mu : uPl) {
        for (const auto* mv : vPl) {
          const double xc = 0.5 * (mu->localCoord + mv->localCoord);
          const double yc = (mv->localCoord - mu->localCoord) / (2.0 * tanStereo);
          const double zc = 0.5 * (mu->surface->center(gctx).x()
                                 + mv->surface->center(gctx).x());
          points2D.push_back({xc, yc, zc, 2});
        }
      }
    }
  }

  if (points2D.empty()) return seeds;

  // --- Hough histogram -----------------------------------------------------
  const double halfSize = m_houghHalfSize.value();
  const double binSize  = m_houghBinSize.value();
  const int    nBins    = static_cast<int>(2.0 * halfSize / binSize) + 1;

  std::vector<std::vector<int>> histo(nBins, std::vector<int>(nBins, 0));
  auto toBin = [&](double c) {
    return std::max(0, std::min(nBins - 1,
                                static_cast<int>((c + halfSize) / binSize)));
  };
  auto fromBin = [&](int b) { return -halfSize + (b + 0.5) * binSize; };

  for (const auto& p : points2D) histo[toBin(p.x)][toBin(p.y)] += p.weight;

  const int minVotes    = m_houghMinVotes.value();
  const int suppressRad = static_cast<int>(std::ceil(15.0 / binSize));

  struct Peak { int ix, iy, votes; };
  std::vector<Peak> peaks;
  for (int ix = 0; ix < nBins; ++ix) {
    for (int iy = 0; iy < nBins; ++iy) {
      const int v = histo[ix][iy];
      if (v < minVotes) continue;
      bool isMax = true;
      for (int dx = -1; dx <= 1 && isMax; ++dx) {
        for (int dy = -1; dy <= 1 && isMax; ++dy) {
          if (dx == 0 && dy == 0) continue;
          const int nx = ix + dx, ny = iy + dy;
          if (nx < 0 || nx >= nBins || ny < 0 || ny >= nBins) continue;
          if (histo[nx][ny] > v) isMax = false;
        }
      }
      if (isMax) peaks.push_back({ix, iy, v});
    }
  }
  std::sort(peaks.begin(), peaks.end(),
            [](const Peak& a, const Peak& b) { return a.votes > b.votes; });

  // --- Track/shower classification by crossing multiplicity ----------------
  const double compatR  = m_seedCompatRadius.value();
  const double compatR2 = compatR * compatR;
  const double maxMult  = m_houghMaxMultiplicity.value();

  struct PeakWithMult { int ix, iy, votes; double multiplicity; };
  std::vector<PeakWithMult> peaksWithMult;
  for (const auto& pk : peaks) {
    const double peakX = fromBin(pk.ix), peakY = fromBin(pk.iy);
    int nCompat = 0;
    std::set<int> compatStations;
    for (const auto& p : points2D) {
      const double dx = p.x - peakX, dy = p.y - peakY;
      if (dx * dx + dy * dy < compatR2) {
        nCompat += p.weight;
        compatStations.insert(static_cast<int>(std::round(p.z / stationTolerance)));
      }
    }
    const int nStations = static_cast<int>(compatStations.size());
    const double mult = (nStations > 0)
        ? static_cast<double>(nCompat) / nStations : 0.0;
    peaksWithMult.push_back({pk.ix, pk.iy, pk.votes, mult});
  }

  std::vector<bool> suppressed(peaksWithMult.size(), false);
  for (std::size_t i = 0; i < peaksWithMult.size(); ++i) {
    if (suppressed[i]) continue;
    for (std::size_t j = i + 1; j < peaksWithMult.size(); ++j) {
      if (suppressed[j]) continue;
      if (std::abs(peaksWithMult[i].ix - peaksWithMult[j].ix) <= suppressRad &&
          std::abs(peaksWithMult[i].iy - peaksWithMult[j].iy) <= suppressRad)
        suppressed[j] = true;
    }
  }

  const int    maxS  = m_maxSeeds.value();
  const double pitch = m_stripPitch.value();

  for (std::size_t pi = 0;
       pi < peaksWithMult.size() && static_cast<int>(seeds.size()) < maxS; ++pi) {
    if (suppressed[pi]) continue;
    const auto& pk = peaksWithMult[pi];
    const double peakX = fromBin(pk.ix), peakY = fromBin(pk.iy);

    if (pk.multiplicity > maxMult) {
      debug() << "[ACTSSequentialTracker] det=" << detMask
              << " Hough peak x=" << peakX << " y=" << peakY
              << " -> SHOWER (multiplicity=" << pk.multiplicity
              << " > " << maxMult << "), skipped." << endmsg;
      continue;
    }

    std::map<int, int> xFreq, yFreq;
    double firstZ = std::numeric_limits<double>::max();
    int nPts = 0;
    for (const auto& p : points2D) {
      const double dx = p.x - peakX, dy = p.y - peakY;
      if (dx * dx + dy * dy < compatR2) {
        xFreq[static_cast<int>(std::round(p.x / pitch))] += p.weight;
        yFreq[static_cast<int>(std::round(p.y / pitch))] += p.weight;
        nPts += p.weight;
        if (p.z < firstZ) firstZ = p.z;
      }
    }
    if (nPts == 0) continue;

    double refinedX = peakX; int maxXF = 0;
    for (const auto& [strip, freq] : xFreq)
      if (freq > maxXF) { maxXF = freq; refinedX = strip * pitch; }
    double refinedY = peakY; int maxYF = 0;
    for (const auto& [strip, freq] : yFreq)
      if (freq > maxYF) { maxYF = freq; refinedY = strip * pitch; }

    seeds.push_back({refinedX, refinedY, firstZ, pk.votes, pk.multiplicity});
  }

  return seeds;
}

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode ACTSSequentialTracker::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    m_siTargetHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiTarget.value(), Gaudi::DataHandle::Reader, this);
    m_siPadHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiPad.value(), Gaudi::DataHandle::Reader, this);
    m_mtcHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputMTC.value(), Gaudi::DataHandle::Reader, this);

    m_localHandles[0] = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outSiTarget.value(), Gaudi::DataHandle::Writer, this);
    m_localHandles[1] = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outSiPad.value(), Gaudi::DataHandle::Writer, this);
    m_localHandles[2] = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outMTC.value(), Gaudi::DataHandle::Writer, this);
    m_globalHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outGlobal.value(), Gaudi::DataHandle::Writer, this);

    m_siTargetBitField = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
        m_siTargetBitFieldStr.value());
    m_mtcBitField = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
        m_mtcBitFieldStr.value());

    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[ACTSSequentialTracker] Failed to retrieve ACTSGeoSvc." << endmsg;
      return StatusCode::FAILURE;
    }

    // The per-detector surface lists come from the geometry service, which
    // recorded them during the TGeo walk. There is deliberately no local
    // largest-z-gap heuristic here: with 40/20/90 surfaces the MTC
    // inter-station gaps exceed the inter-detector ones and such a heuristic
    // silently mis-partitions the geometry.
    std::size_t nTot = 0;
    for (int det = 0; det < 3; ++det) {
      const auto& sf = m_geoSvc->surfacesOf(det);
      if (sf.empty()) {
        error() << "[ACTSSequentialTracker] Geometry service reports no "
                << "surfaces for detector " << det
                << ". Sequential tracking needs all three." << endmsg;
        return StatusCode::FAILURE;
      }
      nTot += sf.size();
    }
    info() << "[ACTSSequentialTracker] Initialized. Surfaces: SiTarget="
           << m_geoSvc->surfacesOf(0).size()
           << " SiPad=" << m_geoSvc->surfacesOf(1).size()
           << " MTC="   << m_geoSvc->surfacesOf(2).size()
           << " (total " << nTot << " of " << m_geoSvc->allSurfaces().size()
           << ")" << endmsg;

    if (m_minMeasPerSegment.value().size() != 3) {
      error() << "[ACTSSequentialTracker] MinMeasPerSegment must have exactly "
              << "3 entries [SiTarget, SiPad, MTC]." << endmsg;
      return StatusCode::FAILURE;
    }

    return sc;
  } catch (const std::exception& e) {
    error() << "[ACTSSequentialTracker] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSSequentialTracker] Unknown exception in initialize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode ACTSSequentialTracker::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    // Create all four collections unconditionally, so downstream always finds
    // them even in an event where nothing was reconstructed.
    std::array<edm4hep::TrackCollection*, 3> localOut{};
    for (int det = 0; det < 3; ++det) localOut[det] = m_localHandles[det]->createAndPut();
    auto* globalOut = m_globalHandle->createAndPut();

    const auto& allSurfaces = m_geoSvc->allSurfaces();
    const auto& gctx        = m_geoSvc->geometryContext();

    // =====================================================================
    // PHASE 0 — assemble the measurement pool
    // =====================================================================
    // ONE vector for all three detectors, so a measurement index is globally
    // unique: fingerprints from different detectors are then directly
    // comparable, which is what makes the splice stage's set union trivial.
    std::vector<SeqMeasurement> measurements;

    const auto* stHits = m_siTargetHandle->get();
    if (stHits) {
      for (std::size_t i = 0; i < stHits->size(); ++i) {
        const auto& hit = (*stHits)[i];
        const auto& pos = hit.getPosition();
        const auto& cov = hit.getCovMatrix();
        const uint64_t cellID = hit.getCellID();
        const int layer = static_cast<int>((*m_siTargetBitField)["layer"].value(cellID));
        const int plane = static_cast<int>((*m_siTargetBitField)["plane"].value(cellID));

        const Acts::Surface* surf = m_geoSvc->surfaceByAddress(0, -1, layer, plane);
        if (!surf) {
          warning() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " SiTarget hit " << i << " layer=" << layer
                    << " plane=" << plane << " has no matching surface." << endmsg;
          continue;
        }
        // plane=0 (StripX) measures DD4hep X, plane=1 (StripY) measures Y.
        const double localCoord = (plane == 0) ? pos.x : pos.y;
        const double var        = (plane == 0) ? cov[0] : cov[3];
        measurements.emplace_back(surf, localCoord, 0.0, var, 0.0, false, 0, plane,
                                  hit.getTime(), hit.getEDep());
      }
    }

    const auto* spHits = m_siPadHandle->get();
    if (spHits) {
      for (std::size_t i = 0; i < spHits->size(); ++i) {
        const auto& hit   = (*spHits)[i];
        const auto& pos   = hit.getPosition();
        const auto& cov   = hit.getCovMatrix();
        const int   layer = hit.getQuality();
        const Acts::Surface* surf = m_geoSvc->surfaceByAddress(1, -1, layer, -1);
        if (!surf) {
          warning() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " SiPad hit " << i << " layer=" << layer
                    << " has no matching surface." << endmsg;
          continue;
        }
        measurements.emplace_back(surf, pos.x, pos.y, cov[0], cov[3], true, 1, -1,
                                  hit.getTime(), hit.getEDep());
      }
    }

    const auto* mtcHits = m_mtcHandle->get();
    if (mtcHits) {
      for (std::size_t i = 0; i < mtcHits->size(); ++i) {
        const auto& hit   = (*mtcHits)[i];
        const auto& pos   = hit.getPosition();
        const auto& cov   = hit.getCovMatrix();
        const int   plane = hit.getQuality();  // 0=U, 1=V, 3=combined 2D
        const uint64_t cellID = hit.getCellID();
        const int station = static_cast<int>((*m_mtcBitField)["station"].value(cellID));
        const int layer   = static_cast<int>((*m_mtcBitField)["layer"].value(cellID));

        const Acts::Surface* surf =
            m_geoSvc->surfaceByAddress(2, station, layer, plane);
        if (!surf) {
          warning() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " MTC hit " << i << " station=" << station
                    << " layer=" << layer << " plane=" << plane
                    << " has no matching surface." << endmsg;
          continue;
        }
        const bool is2D = (plane == 3);
        measurements.emplace_back(surf, pos.x, is2D ? pos.y : 0.0,
                                  cov[0], is2D ? cov[3] : 0.0, is2D, 2, plane,
                                  hit.getTime(), hit.getEDep());
      }
    }

    if (measurements.empty()) {
      debug() << "[ACTSSequentialTracker] evt=" << evtNum << " no measurements." << endmsg;
      return StatusCode::SUCCESS;
    }

    std::sort(measurements.begin(), measurements.end(),
              [&](const SeqMeasurement& a, const SeqMeasurement& b) {
                return a.surface->center(gctx).x() < b.surface->center(gctx).x();
              });

    // =====================================================================
    // PHASE 1 — density purge of the measurement pool (shower cleaning)
    // =====================================================================
    // Per surface, a measurement with more than HitPurgeMaxNeighbors others
    // within HitPurgeWindow is shower-like and is removed BEFORE seeding and
    // fitting — so unlike the seed-level isolation filter, this keeps shower
    // cores out of the fit itself. Per detector, because the physical density
    // scales differ (75 um strips vs 5.5 mm pads vs U/V combinatorics).
    const std::size_t nBeforePurge = measurements.size();
    {
      const auto& purgeWin  = m_hitPurgeWindow.value();
      const auto& purgeMaxN = m_hitPurgeMaxNeighbors.value();
      const bool purgeActive =
          purgeWin.size() == 3 && purgeMaxN.size() == 3 &&
          (purgeWin[0] > 0.0 || purgeWin[1] > 0.0 || purgeWin[2] > 0.0);
      if (purgeActive) {
        std::unordered_map<const Acts::Surface*, std::vector<std::size_t>> bySurf;
        for (std::size_t i = 0; i < measurements.size(); ++i) {
          const int det = measurements[i].detectorID;
          if (det < 0 || det > 2 || purgeWin[det] <= 0.0) continue;
          bySurf[measurements[i].surface].push_back(i);
        }
        std::vector<bool> keep(measurements.size(), true);
        for (const auto& [surf, idxs] : bySurf) {
          const int    det  = measurements[idxs.front()].detectorID;
          const double win  = purgeWin[det];
          const int    maxN = purgeMaxN[det];
          if (static_cast<int>(idxs.size()) <= maxN) continue;
          for (std::size_t i : idxs) {
            const auto& a = measurements[i];
            int n = 0;
            for (std::size_t j : idxs) {
              if (j == i) continue;
              const auto& b = measurements[j];
              const double d = a.is2D
                  ? std::hypot(a.localCoord - b.localCoord,
                               a.localCoord2 - b.localCoord2)
                  : std::abs(a.localCoord - b.localCoord);
              if (d < win && ++n > maxN) break;
            }
            if (n > maxN) keep[i] = false;
          }
        }
        std::vector<SeqMeasurement> purged;
        purged.reserve(measurements.size());
        for (std::size_t i = 0; i < measurements.size(); ++i)
          if (keep[i]) purged.push_back(measurements[i]);
        measurements.swap(purged);
        if (measurements.empty()) {
          debug() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " all measurements purged as shower-like." << endmsg;
          return StatusCode::SUCCESS;
        }
      }
    }

    debug() << "[ACTSSequentialTracker] evt=" << evtNum
            << " SiTarget=" << (stHits  ? stHits->size()  : 0)
            << " SiPad="    << (spHits  ? spHits->size()  : 0)
            << " MTC="      << (mtcHits ? mtcHits->size() : 0)
            << " measurements=" << measurements.size()
            << " (purged " << (nBeforePurge - measurements.size()) << ")" << endmsg;

    // =====================================================================
    // Shared fitting infrastructure (once per event)
    // =====================================================================
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
    {
      const auto& ironRanges = m_ironFieldRanges.value();
      if (!ironRanges.empty() && ironRanges.size() % 7 == 0) {
        std::vector<SeqIronSlabBField::Slab> slabs;
        slabs.reserve(ironRanges.size() / 7);
        for (std::size_t i = 0; i < ironRanges.size(); i += 7) {
          slabs.push_back({ironRanges[i+0] * Acts::UnitConstants::mm,
                           ironRanges[i+1] * Acts::UnitConstants::mm,
                           ironRanges[i+2] * Acts::UnitConstants::mm,
                           ironRanges[i+3] * Acts::UnitConstants::mm,
                           ironRanges[i+4] * Acts::UnitConstants::mm,
                           ironRanges[i+5] * Acts::UnitConstants::mm,
                           ironRanges[i+6] * Acts::UnitConstants::T});
        }
        bField = std::make_shared<SeqIronSlabBField>(std::move(slabs));
      } else {
        bField = std::make_shared<Acts::ConstantBField>(
            Acts::Vector3(m_bFieldX.value() * Acts::UnitConstants::T,
                          m_bFieldY.value() * Acts::UnitConstants::T,
                          m_bFieldZ.value() * Acts::UnitConstants::T));
      }
    }

    SeqCalibrator calibrator;
    calibrator.meas = &measurements;

    Acts::GainMatrixUpdater  gainMatrixUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    SeqSurfaceAccessor kfSurfaceAccessor;
    kfSurfaceAccessor.meas = &measurements;

    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> kfExtensions;
    kfExtensions.calibrator.template connect<&SeqCalibrator::operator()>(&calibrator);
    kfExtensions.updater
        .connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
            &gainMatrixUpdater);
    kfExtensions.smoother
        .connect<&Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
            &kfSmoother);
    kfExtensions.surfaceAccessor
        .connect<&SeqSurfaceAccessor::operator()>(&kfSurfaceAccessor);

    // One KF for the whole event: the surface sequence is a fit() argument, so
    // the same fitter serves the per-detector segment refits (that detector's
    // surfaces) and the combined refit of a spliced chain (all surfaces).
    SeqStepper kfStepper(bField);
    Acts::DirectNavigator kfNavigator;
    SeqKFPropagator kfPropagator(std::move(kfStepper), std::move(kfNavigator));
    SeqKF kf(std::move(kfPropagator),
             Acts::getDefaultLogger("SeqKF", Acts::Logging::WARNING));

    Acts::PropagatorPlainOptions pOptions(gctx, m_mctx);
    pOptions.direction = Acts::Direction::Forward();
    pOptions.stepping.maxStepSize = m_maxStepSize.value();
    pOptions.maxSteps = static_cast<std::size_t>(m_maxPropSteps.value());

    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(
        gctx, m_mctx, std::cref(m_cctx), kfExtensions, pOptions);

    // One track container for the whole event; only selected indices are written.
    auto trackBackend = std::make_shared<Acts::VectorTrackContainer>();
    auto trajBackend  = std::make_shared<Acts::VectorMultiTrajectory>();
    SeqTrackContainer tracks(trackBackend, trajBackend);

    // Seed covariance — see the SeedVar* properties for why this is not simply
    // "as loose as possible".
    Acts::BoundSquareMatrix seedCov = Acts::BoundSquareMatrix::Zero();
    seedCov(Acts::eBoundLoc0,   Acts::eBoundLoc0)   = m_seedVarLoc.value();
    seedCov(Acts::eBoundLoc1,   Acts::eBoundLoc1)   = m_seedVarLoc.value();
    seedCov(Acts::eBoundPhi,    Acts::eBoundPhi)    = m_seedVarAngle.value();
    seedCov(Acts::eBoundTheta,  Acts::eBoundTheta)  = m_seedVarAngle.value();
    seedCov(Acts::eBoundQOverP, Acts::eBoundQOverP) = m_seedVarQOverP.value();
    seedCov(Acts::eBoundTime,   Acts::eBoundTime)   = 1e9;

    // The refits, in contrast, must stay measurement-dominated: their starting
    // parameters come from a fit that already used these hits, so a tight prior
    // there would double-count the information.
    Acts::BoundSquareMatrix refitCov = Acts::BoundSquareMatrix::Zero();
    refitCov(Acts::eBoundLoc0,   Acts::eBoundLoc0)   = 1e2;
    refitCov(Acts::eBoundLoc1,   Acts::eBoundLoc1)   = 1e2;
    refitCov(Acts::eBoundPhi,    Acts::eBoundPhi)    = 1.0;
    refitCov(Acts::eBoundTheta,  Acts::eBoundTheta)  = 1.0;
    refitCov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.04;
    refitCov(Acts::eBoundTime,   Acts::eBoundTime)   = 1e9;

    const double seedQoverP =
        -1.0 / (m_seedMomentum.value() * Acts::UnitConstants::GeV);

    // ---- Shared little helpers ---------------------------------------------
    auto chi2NdfOf = [&](const auto& track) -> double {
      const int ndf = std::max(1, static_cast<int>(track.nDoF()) - kHelixParams);
      return track.chi2() / static_cast<double>(ndf);
    };

    // TrackStateCreator sets MeasurementFlag on outlier states too, so outliers
    // must be excluded explicitly: they are not part of the fit and must never
    // reach a frozen-hit refit.
    auto fingerprintOf = [](const auto& track) -> std::set<std::size_t> {
      std::set<std::size_t> fp;
      for (const auto& ts : track.trackStatesReversed()) {
        if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag) &&
            !ts.typeFlags().test(Acts::TrackStateFlag::OutlierFlag) &&
            ts.hasUncalibratedSourceLink()) {
          fp.insert(ts.getUncalibratedSourceLink().template get<SeqSourceLink>().index);
        }
      }
      return fp;
    };

    // Fitted states at the two ends of a track, in ACTS global coordinates.
    // Ordered by beam (ACTS x), so `first` is upstream and `second` downstream
    // regardless of the order ACTS happens to store the states in.
    auto endStatesOf = [&](const auto& track)
        -> std::optional<std::pair<SeqEndState, SeqEndState>> {
      std::vector<std::pair<double, SeqEndState>> states;
      for (const auto& ts : track.trackStatesReversed()) {
        if (!ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) continue;
        if (ts.typeFlags().test(Acts::TrackStateFlag::OutlierFlag)) continue;
        Acts::BoundVector bv;
        if      (ts.hasSmoothed()) bv = ts.smoothed();
        else if (ts.hasFiltered()) bv = ts.filtered();
        else continue;
        const double phi = bv[Acts::eBoundPhi], theta = bv[Acts::eBoundTheta];
        Acts::Vector3 dir(std::sin(theta) * std::cos(phi),
                          std::sin(theta) * std::sin(phi),
                          std::cos(theta));
        const auto& surf = ts.referenceSurface();
        Acts::Vector3 pos = surf.localToGlobal(
            gctx, Acts::Vector2(bv[Acts::eBoundLoc0], bv[Acts::eBoundLoc1]), dir);
        states.push_back({pos.x(), SeqEndState{pos, dir, bv[Acts::eBoundQOverP]}});
      }
      if (states.empty()) return std::nullopt;
      std::sort(states.begin(), states.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
      return std::make_pair(states.front().second, states.back().second);
    };

    // Straight-line back-extrapolation of an end state onto a surface, as
    // BoundTrackParameters at the loose seed covariance. Used both to re-seed a
    // CKF pass and to start a refit.
    auto paramsAtSurface = [&](const SeqEndState& st, const Acts::Surface* target,
                               const Acts::BoundSquareMatrix& cov)
        -> std::optional<Acts::BoundTrackParameters> {
      if (target == nullptr || std::abs(st.dir.x()) < 1e-6) return std::nullopt;
      const double xT = target->center(gctx).x();
      Acts::Vector3 pos = st.pos + st.dir * ((xT - st.pos.x()) / st.dir.x());
      Acts::Vector4 pos4(pos.x(), pos.y(), pos.z(), 0.0);
      auto res = Acts::BoundTrackParameters::create(
          gctx, target->getSharedPtr(), pos4, st.dir, st.qOverP, cov,
          Acts::ParticleHypothesis::muon());
      if (!res.ok()) return std::nullopt;
      return *res;
    };

    // Frozen-hit Kalman refit over a given surface sequence. `fp` is the hit
    // set; nothing is re-selected. Returns the new track index on success.
    auto refit = [&](const std::set<std::size_t>& fp,
                     const Acts::BoundTrackParameters& start,
                     const std::vector<const Acts::Surface*>& surfaces)
        -> std::optional<Acts::TrackIndexType> {
      if (fp.size() < 3) return std::nullopt;
      std::vector<Acts::SourceLink> fitSLinks;
      fitSLinks.reserve(fp.size());
      for (std::size_t mIdx : fp) {
        SeqSourceLink ssl;
        ssl.index = mIdx;
        ssl.setGeometryId(measurements[mIdx].surface->geometryId());
        fitSLinks.emplace_back(ssl);
      }
      auto res = kf.fit(fitSLinks.begin(), fitSLinks.end(), start, kfOptions,
                        surfaces, tracks);
      if (!res.ok()) return std::nullopt;
      if ((*res).nMeasurements() < 3) return std::nullopt;
      return (*res).index();
    };

    // =====================================================================
    // PHASE 2 — per-detector segment finding, two iterations each
    // =====================================================================
    std::vector<SeqSegment> segments;

    for (int det = 0; det < 3; ++det) {
      const auto& detSurfaces = m_geoSvc->surfacesOf(det);
      if (detSurfaces.empty()) continue;

      // Global indices of this detector's measurements, minus the ones already
      // taken by an accepted segment of the same detector.
      std::set<std::size_t> used;
      const int minMeas = m_minMeasPerSegment.value()[det];

      for (int iteration = 0; iteration < 2; ++iteration) {
        const double chi2Cut = (iteration == 0)
            ? m_chi2CutOff.value()
            : m_chi2CutOff.value() * m_looseChi2Scale.value();
        const double segMaxChi2 = (iteration == 0)
            ? m_segMaxChi2Tight.value()
            : m_segMaxChi2Loose.value();

        // Remaining pool of this detector, as a standalone vector for the
        // seeder, plus the source links the CKF may use.
        std::vector<SeqMeasurement> pool;
        std::vector<Acts::SourceLink> slinks;
        for (std::size_t i = 0; i < measurements.size(); ++i) {
          if (measurements[i].detectorID != det) continue;
          if (used.count(i)) continue;
          pool.push_back(measurements[i]);
          SeqSourceLink ssl;
          ssl.index = i;
          ssl.setGeometryId(measurements[i].surface->geometryId());
          slinks.push_back(Acts::SourceLink(ssl));
        }
        if (static_cast<int>(pool.size()) < minMeas) {
          debug() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " det=" << det << " iter=" << iteration
                  << " pool=" << pool.size() << " below MinMeasPerSegment="
                  << minMeas << ", stopping this detector." << endmsg;
          break;
        }

        std::sort(slinks.begin(), slinks.end(),
                  [](const Acts::SourceLink& a, const Acts::SourceLink& b) {
                    return a.get<SeqSourceLink>().geometryId() <
                           b.get<SeqSourceLink>().geometryId();
                  });

        auto seeds = findSeeds(pool, measurements, gctx, det);
        if (seeds.empty()) continue;

        // ---- CKF for this (detector, iteration) --------------------------
        SeqSourceLinkAccessor slAccessor;
        slAccessor.slinks = &slinks;

        Acts::MeasurementSelectorCuts measCuts;
        measCuts.chi2CutOff            = {chi2Cut};
        measCuts.numMeasurementsCutOff =
            {static_cast<std::size_t>(m_numMeasCutOff.value())};
        Acts::MeasurementSelector measSelector(measCuts);

        using SLinkIter = std::vector<Acts::SourceLink>::const_iterator;
        Acts::TrackStateCreator<SLinkIter, SeqTrackContainer> tsc;
        tsc.sourceLinkAccessor.connect<&SeqSourceLinkAccessor::operator()>(&slAccessor);
        tsc.calibrator.template connect<&SeqCalibrator::operator()>(&calibrator);
        tsc.measurementSelector
            .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
                &measSelector);

        Acts::CombinatorialKalmanFilterExtensions<SeqTrackContainer> ckfExtensions;
        ckfExtensions.updater
            .connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
                &gainMatrixUpdater);
        ckfExtensions.createTrackStates
            .connect<&Acts::TrackStateCreator<SLinkIter, SeqTrackContainer>::createTrackStates>(
                &tsc);

        // THE point of this algorithm: the navigator sees only this detector.
        SeqFixedNavigator ckfNavigator;
        ckfNavigator.surfaces.assign(detSurfaces.begin(), detSurfaces.end());

        SeqStepper       ckfStepper(bField);
        SeqCKFPropagator ckfPropagator(std::move(ckfStepper), std::move(ckfNavigator));
        SeqCKF ckf(std::move(ckfPropagator),
                   Acts::getDefaultLogger("SeqCKF", Acts::Logging::WARNING));

        Acts::CombinatorialKalmanFilterOptions<SeqTrackContainer> ckfOptions(
            gctx, m_mctx, std::cref(m_cctx), ckfExtensions, pOptions, true, true);

        auto runCKFPass = [&](const Acts::BoundTrackParameters& sp)
            -> std::optional<Acts::TrackIndexType> {
          auto res = ckf.findTracks(sp, ckfOptions, tracks);
          if (!res.ok()) return std::nullopt;
          std::optional<Acts::TrackIndexType> bestIdx;
          double best = std::numeric_limits<double>::max();
          for (const auto& t : *res) {
            if (t.nMeasurements() < 3) continue;
            const double c2n = chi2NdfOf(t);
            if (c2n < best) { best = c2n; bestIdx = t.index(); }
          }
          return bestIdx;
        };

        // The CKF always starts at this detector's most upstream surface. A
        // track born in the middle of the detector simply leaves holes in front
        // of its first hit, which the loose seed covariance absorbs.
        const Acts::Surface* sfSeed = detSurfaces.front();

        std::vector<SeqSegment> found;

        for (std::size_t iSeed = 0; iSeed < seeds.size(); ++iSeed) {
          // DD4hep (x=transverse X, y=transverse Y, z=beam) -> ACTS, where the
          // geometry's 90-degree Y rotation makes ACTS x the beam axis:
          // ePos1 = dd_y, ePos2 = dd_x, direction = (dd_dz, dd_dx, dd_dy).
          const double ddx = seeds[iSeed].x;
          const double ddy = seeds[iSeed].y;

          Acts::Vector4 seedPos4;
          seedPos4[Acts::ePos0] = sfSeed->center(gctx).x();
          seedPos4[Acts::ePos1] = ddy;
          seedPos4[Acts::ePos2] = ddx;
          seedPos4[Acts::eTime] = 0.0;
          Acts::Vector3 seedDir(1.0, 0.0, 0.0);  // along the beam

          auto seedParamsRes = Acts::BoundTrackParameters::create(
              gctx, sfSeed->getSharedPtr(), seedPos4, seedDir, seedQoverP,
              seedCov, Acts::ParticleHypothesis::muon());
          if (!seedParamsRes.ok()) continue;

          // ---- Search, then refine by re-seeding ------------------------
          auto candIdx = runCKFPass(*seedParamsRes);
          if (!candIdx) {
            debug() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " det=" << det << " iter=" << iteration
                    << " seed=" << iSeed << " at (" << ddx << "," << ddy
                    << ") mm: CKF returned no candidate with >=3 measurements"
                    << endmsg;
            continue;
          }

          double bestC2n = chi2NdfOf(tracks.getTrack(*candIdx));
          auto   bestIdx = *candIdx;

          if (auto ends = endStatesOf(tracks.getTrack(*candIdx))) {
            if (auto reseed = paramsAtSurface(ends->first, sfSeed, seedCov)) {
              if (auto refIdx = runCKFPass(*reseed)) {
                const double c2n = chi2NdfOf(tracks.getTrack(*refIdx));
                if (c2n < bestC2n) { bestC2n = c2n; bestIdx = *refIdx; }
              }
            }
          }

          // ---- Frozen-hit refit at nominal variances --------------------
          auto fp = fingerprintOf(tracks.getTrack(bestIdx));
          if (auto ends = endStatesOf(tracks.getTrack(bestIdx))) {
            if (auto start = paramsAtSurface(ends->first, sfSeed, refitCov)) {
              if (auto rIdx = refit(fp, *start, detSurfaces)) {
                bestIdx = *rIdx;
                fp      = fingerprintOf(tracks.getTrack(bestIdx));
              }
            }
          }

          auto finalTrack = tracks.getTrack(bestIdx);
          const int    ndf  = std::max(
              1, static_cast<int>(finalTrack.nDoF()) - kHelixParams);
          const double chi2 = finalTrack.chi2();
          const double c2n  = chi2 / static_cast<double>(ndf);
          const int    nMeas = static_cast<int>(finalTrack.nMeasurements());

          // Logged rather than silently dropped: "detector X found no segment"
          // is the single most common reason a global track fails to form
          // (a missing SiPad segment forbids BOTH of its neighbour links), so
          // the cut that rejected the candidate has to be visible.
          if (nMeas < minMeas || c2n > segMaxChi2) {
            debug() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " det=" << det << " iter=" << iteration
                    << " seed=" << iSeed << " candidate rejected: nMeas="
                    << nMeas << " (min " << minMeas << ") chi2/ndf=" << c2n
                    << " (max " << segMaxChi2 << ")" << endmsg;
            continue;
          }

          auto ends = endStatesOf(finalTrack);
          if (!ends) continue;

          // Mask this segment's hits from the pool the CKF may still use, so
          // the seeds after it cannot rebuild the same particle.
          //
          // This is not an optimisation, it is what makes multi-track events
          // work. The seed covariance is deliberately loose (phi/theta free to
          // ~1 rad), so a seed's transverse position barely constrains the fit:
          // without masking, a seed sitting on track B happily converges onto
          // track A and the two come out with an identical hit set — one of
          // them is then discarded as a duplicate and a real track is lost.
          // Masking between ITERATIONS is not enough; it has to happen between
          // seeds. `slinks` stays sorted under erase and the accessor holds a
          // pointer to it, so nothing needs re-wiring.
          slinks.erase(
              std::remove_if(slinks.begin(), slinks.end(),
                             [&seg_fp = std::as_const(fp)](const Acts::SourceLink& sl) {
                               return seg_fp.count(
                                          sl.get<SeqSourceLink>().index) > 0;
                             }),
              slinks.end());

          SeqSegment seg;
          seg.det       = det;
          seg.iteration = iteration;
          seg.idx       = bestIdx;
          seg.fp        = std::move(fp);
          seg.chi2      = chi2;
          seg.ndf       = ndf;
          seg.front     = ends->first;
          seg.back      = ends->second;
          seg.seedX     = ddx;
          seg.seedY     = ddy;
          found.push_back(std::move(seg));
        }

        // ---- Best-first duplicate filter within this iteration ----------
        std::sort(found.begin(), found.end(),
                  [](const SeqSegment& a, const SeqSegment& b) {
                    return a.chi2 / std::max(1, a.ndf) < b.chi2 / std::max(1, b.ndf);
                  });
        const double dupFrac = m_dupOverlapFraction.value();
        std::vector<SeqSegment> survivors;
        for (auto& cand : found) {
          bool dup = false;
          for (const auto& keep : survivors) {
            std::size_t nShared = 0;
            for (std::size_t h : cand.fp) if (keep.fp.count(h)) ++nShared;
            const double smaller =
                static_cast<double>(std::min(cand.fp.size(), keep.fp.size()));
            if (smaller > 0 && nShared / smaller > dupFrac) {
              // Two seeds converging on one particle is normal; two seeds on
              // two DIFFERENT particles landing here means one real track is
              // being thrown away, so say which seeds were involved.
              debug() << "[ACTSSequentialTracker] evt=" << evtNum
                      << " det=" << det << " iter=" << iteration
                      << " segment from seed (" << cand.seedX << ","
                      << cand.seedY << ") dropped as duplicate of the one from ("
                      << keep.seedX << "," << keep.seedY << "): "
                      << nShared << "/" << smaller << " hits shared" << endmsg;
              dup = true;
              break;
            }
          }
          if (!dup) survivors.push_back(std::move(cand));
        }

        // Mask the survivors' hits so the next iteration works on what is left.
        for (const auto& s : survivors) {
          used.insert(s.fp.begin(), s.fp.end());
          segments.push_back(s);
        }

        debug() << "[ACTSSequentialTracker] evt=" << evtNum
                << " det=" << det << " iter=" << iteration
                << " pool=" << pool.size() << " seeds=" << seeds.size()
                << " segments=" << survivors.size() << endmsg;
      }
    }

    // =====================================================================
    // PHASE 3 — splice adjacent segments and refit the chains
    // =====================================================================
    // Boundaries in the beam coordinate, halfway between the last surface of
    // one detector and the first of the next. Derived from the geometry, never
    // hard-coded: a variant geometry moves them on its own.
    auto boundaryBetween = [&](int detUp, int detDn) -> double {
      const auto& up = m_geoSvc->surfacesOf(detUp);
      const auto& dn = m_geoSvc->surfacesOf(detDn);
      return 0.5 * (up.back()->center(gctx).x() + dn.front()->center(gctx).x());
    };

    // Straight-line extrapolation to a beam plane. Everything upstream of the
    // MTC iron is field-free, so this is exact there; from the MTC's first
    // SciFi plane back to the boundary it crosses one iron slab, which biases
    // it slightly at low momentum. That bias lives inside the tolerance
    // window — the combined refit below is what actually judges the splice.
    auto extrapTo = [](const SeqEndState& st, double xB) -> std::optional<Acts::Vector3> {
      if (std::abs(st.dir.x()) < 1e-6) return std::nullopt;
      return st.pos + st.dir * ((xB - st.pos.x()) / st.dir.x());
    };

    const double maxDist  = m_spliceMaxDist.value();
    const double maxAngle = m_spliceMaxAngle.value();
    const double angleW   = m_spliceAngleWeight.value();

    std::vector<SeqLink> links;
    for (int detUp = 0; detUp < 2; ++detUp) {         // 0->1 and 1->2 only
      const int detDn = detUp + 1;
      const double xB = boundaryBetween(detUp, detDn);
      for (std::size_t i = 0; i < segments.size(); ++i) {
        if (segments[i].det != detUp) continue;
        for (std::size_t j = 0; j < segments.size(); ++j) {
          if (segments[j].det != detDn) continue;
          auto pUp = extrapTo(segments[i].back,  xB);
          auto pDn = extrapTo(segments[j].front, xB);
          if (!pUp || !pDn) continue;
          // Transverse plane at the boundary: ACTS y and z (= DD4hep y and x).
          const double dist = std::hypot(pUp->y() - pDn->y(), pUp->z() - pDn->z());
          const double cosA = std::clamp(
              segments[i].back.dir.dot(segments[j].front.dir), -1.0, 1.0);
          const double angle = std::acos(cosA);
          const bool ok = (dist <= maxDist && angle <= maxAngle);
          // Every considered pair is logged, passing or not. A missing global
          // track is almost always a pair that just failed one of these two
          // windows, and without the numbers there is no way to tell which.
          debug() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " splice det" << detUp << "seg" << i
                  << " -> det" << detDn << "seg" << j
                  << " dist=" << dist << " (max " << maxDist << ")"
                  << " angle=" << angle << " (max " << maxAngle << ")"
                  << (ok ? " CANDIDATE" : " rejected") << endmsg;
          if (!ok) continue;
          links.push_back({i, j, dist, angle,
                           dist * dist + angleW * angle * angle});
        }
      }
    }

    // Greedy by score. Each segment has one upstream and one downstream port;
    // a link is taken only if both of its ports are still free, so a SiPad
    // segment can pick up a SiTarget upstream AND an MTC downstream and form a
    // 3-chain, while two SiTarget segments can never claim the same SiPad one.
    std::sort(links.begin(), links.end(),
              [](const SeqLink& a, const SeqLink& b) { return a.score < b.score; });

    const std::size_t kNone = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> nextOf(segments.size(), kNone);
    std::vector<std::size_t> prevOf(segments.size(), kNone);
    std::vector<SeqLink>     takenLinks;
    for (const auto& lk : links) {
      if (nextOf[lk.up] != kNone || prevOf[lk.dn] != kNone) continue;
      nextOf[lk.up] = lk.dn;
      prevOf[lk.dn] = lk.up;
      takenLinks.push_back(lk);
    }

    // Chains: start from every segment with a successor and no predecessor.
    std::vector<std::vector<std::size_t>> chains;
    for (std::size_t i = 0; i < segments.size(); ++i) {
      if (prevOf[i] != kNone || nextOf[i] == kNone) continue;
      std::vector<std::size_t> chain{i};
      std::size_t cur = i;
      while (nextOf[cur] != kNone) { cur = nextOf[cur]; chain.push_back(cur); }
      chains.push_back(std::move(chain));
    }

    // Combined refit of one chain over ALL surfaces. Accepts or rejects.
    struct AcceptedGlobal {
      Acts::TrackIndexType idx;
      double chi2; int ndf;
      double seedX, seedY;
      int typeMask;
      std::vector<std::size_t> members;
    };
    std::vector<AcceptedGlobal> globals;

    auto tryChain = [&](const std::vector<std::size_t>& chain)
        -> std::optional<AcceptedGlobal> {
      if (chain.size() < 2) return std::nullopt;
      std::set<std::size_t> unionFp;
      int mask = 0;
      for (std::size_t s : chain) {
        unionFp.insert(segments[s].fp.begin(), segments[s].fp.end());
        mask |= (1 << segments[s].det);
      }
      // Start from the chain's most upstream state, taken back to the very
      // first surface of the geometry — the sequence handed to the fitter.
      auto start = paramsAtSurface(segments[chain.front()].front,
                                   allSurfaces.front(), refitCov);
      if (!start) return std::nullopt;
      auto idx = refit(unionFp, *start, allSurfaces);
      if (!idx) return std::nullopt;

      auto t = tracks.getTrack(*idx);
      const int    ndf  = std::max(1, static_cast<int>(t.nDoF()) - kHelixParams);
      const double chi2 = t.chi2();
      const double c2n  = chi2 / static_cast<double>(ndf);
      debug() << "[ACTSSequentialTracker] evt=" << evtNum
              << " chain of " << chain.size() << " (mask=" << mask << ")"
              << " hits=" << unionFp.size()
              << " combined refit chi2/ndf=" << c2n
              << " (max " << m_globalMaxChi2PerNdf.value() << ")"
              << (c2n > m_globalMaxChi2PerNdf.value() ? " REJECTED" : " accepted")
              << endmsg;
      if (c2n > m_globalMaxChi2PerNdf.value()) return std::nullopt;

      return AcceptedGlobal{*idx, chi2, ndf,
                            segments[chain.front()].seedX,
                            segments[chain.front()].seedY,
                            mask, chain};
    };

    for (const auto& chain : chains) {
      if (auto g = tryChain(chain)) { globals.push_back(std::move(*g)); continue; }

      // A 3-chain that does not fit as a whole may still hold a real 2-chain.
      // Try both sub-chains, better-scored link first; the odd segment out
      // stays local.
      if (chain.size() == 3) {
        auto scoreOf = [&](std::size_t up, std::size_t dn) {
          for (const auto& lk : takenLinks)
            if (lk.up == up && lk.dn == dn) return lk.score;
          return std::numeric_limits<double>::max();
        };
        std::vector<std::vector<std::size_t>> subs{
            {chain[0], chain[1]}, {chain[1], chain[2]}};
        std::sort(subs.begin(), subs.end(),
                  [&](const auto& a, const auto& b) {
                    return scoreOf(a[0], a[1]) < scoreOf(b[0], b[1]);
                  });
        for (const auto& sub : subs) {
          if (auto g = tryChain(sub)) { globals.push_back(std::move(*g)); break; }
        }
      }
    }

    for (const auto& g : globals)
      for (std::size_t s : g.members) segments[s].spliced = true;

    // =====================================================================
    // PHASE 4 — write out
    // =====================================================================
    // The TrackState encoding is deliberately identical to ACTSProtoTracker's:
    // job5 (EDM4HEP2RNTuple), the event display and diagnose_mtc_curvature.py
    // all decode it, and they are reused as-is on the new collections.
    //   state[0]      location=AtIP,    D0 = seed x,  Z0 = seed y
    //   state[1..N]   location=AtOther, referencePoint = fitted global position,
    //                 D0 = raw loc0, Z0 = per-state chi2, omega = stereo tilt.
    //   MTC U/V partners are averaged into one state with omega = 0.
    auto writeTrack = [&](edm4hep::TrackCollection* out,
                          Acts::TrackIndexType idx, double chi2, int ndf,
                          double ddx, double ddy, int typeMask) {
      auto finalTrack = tracks.getTrack(idx);
      auto track = out->create();
      track.setType(typeMask);
      track.setChi2(static_cast<float>(chi2));
      track.setNdf(ndf);

      {
        edm4hep::TrackState seedState{};
        seedState.location = edm4hep::TrackState::AtIP;
        seedState.D0       = static_cast<float>(ddx);
        seedState.Z0       = static_cast<float>(ddy);
        track.addToTrackStates(seedState);
      }

      try {
        auto tipIdx = finalTrack.tipIndex();
        auto& mutableTraj = tracks.trackStateContainer();
        std::vector<edm4hep::TrackState> collected;
        while (true) {
          auto ts = mutableTraj.getTrackState(tipIdx);
          if (ts.hasCalibrated()) {
            edm4hep::TrackState edm4ts;
            edm4ts.location = edm4hep::TrackState::AtOther;
            edm4ts.D0       = 0.0f;
            edm4ts.Z0       = static_cast<float>(ts.chi2());
            const auto& surf  = ts.referenceSurface();
            const float beamZ = static_cast<float>(surf.center(gctx).x());

            // Stereo detection: R * localY projected on ACTS z is sin(+-alpha);
            // zero on the flat SiTarget/SiPad/combined-MTC surfaces.
            const Acts::Vector3 localYinGlobal =
                surf.transform(gctx).rotation() * Acts::Vector3::UnitY();
            const double stereoTiltZ = localYinGlobal.z();
            const double cosAlpha    = localYinGlobal.y();
            edm4ts.omega = static_cast<float>(stereoTiltZ);

            auto toGlobalX = [&](double l0, double l1) {
              return static_cast<float>(l0 * cosAlpha + l1 * stereoTiltZ);
            };
            auto toGlobalY = [&](double l0, double l1) {
              return static_cast<float>(-l0 * stereoTiltZ + l1 * cosAlpha);
            };

            if (ts.hasSmoothed() || ts.hasFiltered()) {
              const Acts::BoundVector bv =
                  ts.hasSmoothed() ? ts.smoothed() : ts.filtered();
              edm4ts.phi = static_cast<float>(bv[Acts::eBoundPhi]);
              edm4ts.tanLambda = ts.hasSmoothed()
                  ? static_cast<float>(std::tan(M_PI / 2.0 - bv[Acts::eBoundTheta]))
                  : 0.0f;
              const double l0 = bv[Acts::eBoundLoc0];
              const double l1 = bv[Acts::eBoundLoc1];
              edm4ts.D0 = static_cast<float>(l0);
              edm4ts.referencePoint =
                  edm4hep::Vector3f{toGlobalX(l0, l1), toGlobalY(l0, l1), beamZ};
            } else {
              edm4ts.phi = 0.0f;
              edm4ts.tanLambda = 0.0f;
              edm4ts.referencePoint = edm4hep::Vector3f{0.f, 0.f, beamZ};
            }
            collected.push_back(edm4ts);
          }
          if (!ts.hasPrevious()) break;
          tipIdx = ts.previous();
        }

        // Pair-average MTC U/V partners. On a SciFi strip surface only loc0 is
        // updated; the bound-to-global rotation then maps the untouched loc1
        // into +-sin(alpha) in global x with the sign flipping between U and V
        // — the visible zigzag. Averaging the partners cancels it to first
        // order. omega=0 marks the result as "already paired".
        std::sort(collected.begin(), collected.end(),
                  [](const edm4hep::TrackState& a, const edm4hep::TrackState& b) {
                    return a.referencePoint.z < b.referencePoint.z;
                  });
        for (std::size_t i = 0; i < collected.size(); ) {
          bool didPair = false;
          if (i + 1 < collected.size()) {
            const auto& a = collected[i];
            const auto& b = collected[i + 1];
            const float dz = std::abs(a.referencePoint.z - b.referencePoint.z);
            if (std::abs(a.omega) > 0.01f && std::abs(b.omega) > 0.01f &&
                a.omega * b.omega < 0.0f && dz < 5.0f) {
              edm4hep::TrackState avg;
              avg.location  = edm4hep::TrackState::AtOther;
              avg.D0        = 0.0f;
              avg.Z0        = a.Z0 + b.Z0;
              avg.omega     = 0.0f;
              avg.phi       = 0.5f * (a.phi + b.phi);
              avg.tanLambda = 0.5f * (a.tanLambda + b.tanLambda);
              avg.referencePoint = edm4hep::Vector3f{
                  0.5f * (a.referencePoint.x + b.referencePoint.x),
                  0.5f * (a.referencePoint.y + b.referencePoint.y),
                  0.5f * (a.referencePoint.z + b.referencePoint.z)};
              track.addToTrackStates(avg);
              i += 2;
              didPair = true;
            }
          }
          if (!didPair) { track.addToTrackStates(collected[i]); ++i; }
        }
      } catch (const std::exception& e) {
        warning() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " trackStates iteration failed: " << e.what() << endmsg;
      }
    };

    for (const auto& g : globals)
      writeTrack(globalOut, g.idx, g.chi2, g.ndf, g.seedX, g.seedY, g.typeMask);

    std::array<int, 3> nLocal{0, 0, 0};
    for (const auto& s : segments) {
      if (s.spliced) continue;
      writeTrack(localOut[s.det], s.idx, s.chi2, s.ndf, s.seedX, s.seedY,
                 1 << s.det);
      ++nLocal[s.det];
    }

    info() << "[ACTSSequentialTracker] evt=" << evtNum
           << " measurements=" << measurements.size()
           << " segments=" << segments.size()
           << " links=" << takenLinks.size() << "/" << links.size()
           << " global=" << globals.size()
           << " local(SiT/SiP/MTC)=" << nLocal[0] << "/" << nLocal[1]
           << "/" << nLocal[2] << endmsg;

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[ACTSSequentialTracker] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSSequentialTracker] Unknown exception in execute()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode ACTSSequentialTracker::finalize() {
  try {
    m_siTargetHandle.reset();
    m_siPadHandle.reset();
    m_mtcHandle.reset();
    for (auto& h : m_localHandles) h.reset();
    m_globalHandle.reset();
    info() << "[ACTSSequentialTracker] Done. Total events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[ACTSSequentialTracker] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSSequentialTracker] Unknown exception in finalize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(ACTSSequentialTracker)
