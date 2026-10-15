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
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
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
#include <numeric>
#include <set>
#include <sstream>
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
  // Carried through purely so a segment can be joined back to MC truth
  // offline: all three measurement converters copy the originating hit's
  // cellID onto their output, so it is a stable key into the *HitsWindowed
  // collections and their parallel *ContribParticleIDs frame parameter.
  uint64_t             cellID      = 0;
  // Index of the TrackerHit3D this came from, inside its own input
  // collection (SiTarget/SiPad/MTC, i.e. `detectorID` picks the collection).
  // It is what lets a written Track point back at its hits, which is how any
  // consumer -- MC-truth studies, SND-PFA's muon corridor -- finds out which
  // particle built a track. -1 when unknown.
  int                  srcIndex    = -1;

  SeqMeasurement(const Acts::Surface* sf, double lc, double lc2,
                 double var, double var2, bool twod, int detID, int pl,
                 float t, float e, uint64_t cid, int srcIdx = -1)
      : surface(sf), localCoord(lc), localCoord2(lc2),
        variance(var), variance2(var2), is2D(twod),
        detectorID(detID), plane(pl), time(t), eDep(e), cellID(cid),
        srcIndex(srcIdx) {}
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
  // The surface sequence IS this navigator's world: once it is exhausted the
  // CKF must finish, and "end of world" is the only aborter it has for that.
  // Acts::DirectNavigator answers false here and relies on the caller's path
  // limit instead. That was fine for beam-parallel tracks, which crossed the
  // last surface and then coasted 10 m to the loop-protection limit in ~100
  // steps of 100 mm. A track leaving the detector SIDEWAYS did not: it fell
  // off the field map into sub-mm stepping, hit the 100000-step ceiling, and
  // the whole CKF pass came back as PropagatorError::StepCountLimitReached --
  // with 130 clean measurements on it, thrown away as "no candidate". GENIE
  // event 9's SiTarget track (tx = 0.35, x from -55 to +66 mm) was lost to
  // exactly this, from every one of its twelve correct seeds, and the fit
  // then settled for a beam-parallel ghost through the Y strips.
  bool endOfWorldReached(State& s) const { return s.navigationBreak; }
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
// THE SIGN OF loc0. Every flat surface ACTSGeoSvc builds carries the rot90Y
// transform, whose local axes expressed in ACTS global are
//     e0 = (0, 0, -1),  e1 = (0, 1, 0),  normal = (1, 0, 0)
// and ACTS global is (beam, DD4hep Y, DD4hep X). So on those surfaces
//     loc0 = -ddX,   loc1 = +ddY.
// A seed built with BoundTrackParameters::create() from a global position
// honours that (a point at ddX = -138 lands at loc0 = +138), and so does every
// localToGlobal() on a fitted state. The measurements did not: the converters
// hand over ddX in localCoord and this calibrator used to pass it straight
// through as loc0, i.e. with the wrong sign. Seed and measurements then sat in
// mirror-image frames, 2|ddX| apart.
//
// It went unnoticed for as long as the seeder could only find beam-parallel
// tracks near the axis, where 2|ddX| is a few mm and hides inside the 10 mm
// seed sigma. The rewritten seeder finds inclined tracks anywhere in +-300 mm,
// and there the mismatch is hundreds of mm: measured on GENIE event 1, a seed
// placed 1 mm from the muon's own hit was 274 mm from it in loc0, the selector
// took a pion hit 29 mm away instead, and no seed pointed along the muon ever
// produced a segment. That one sign is behind a large part of the MTC's
// fragmentation and impurity.
//
// Only loc0 flips; loc1 (= +ddY) was already right, which is why the SiTarget
// Y strips and the y of every 2D point never showed the problem.
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
      ts.template calibrated<2>() = Acts::ActsVector<2>(-m.localCoord, m.localCoord2);
      ts.template calibratedCovariance<2>() =
          Acts::ActsSquareMatrix<2>{{m.variance, 0.0}, {0.0, m.variance2}};
      ts.setProjectorSubspaceIndices(indices);
      return;
    }
    if (m.detectorID == 2) {
      // MTC SciFi unpaired: 1D on the stereo-rotated U/V surface -> eBoundLoc0.
      // The converter delivers localCoord = cosα·s_raw = s_U (or s_V); the
      // tilted surface's loc0 is −s_U (resp. −s_V), same sign rule as the
      // flat surfaces (see the header comment).
      ts.allocateCalibrated(1);
      ts.template calibrated<1>() = Acts::ActsVector<1>(-m.localCoord);
      ts.template calibratedCovariance<1>() =
          Acts::ActsSquareMatrix<1>{{m.variance}};
      constexpr std::array<Acts::BoundIndices, 1> mtcIdx = {Acts::eBoundLoc0};
      ts.setProjectorSubspaceIndices(mtcIdx);
      return;
    }
    // SiTarget strip: plane=0 (StripX) -> eBoundLoc0, plane=1 (StripY) -> eBoundLoc1.
    // StripX measures ddX, which on this frame is -loc0; StripY measures ddY = +loc1.
    ts.allocateCalibrated(1);
    ts.template calibrated<1>() =
        Acts::ActsVector<1>(m.plane == 1 ? m.localCoord : -m.localCoord);
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
  const Acts::Surface* surface = nullptr;   // the measurement surface it sits on
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
  bool                  extension = false;  // built by PHASE 2.7 from a
                                            // neighbour's track, never a
                                            // standalone find; dropped if it
                                            // ends up unspliced
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
      "Spliced tracks (2 or 3 segments) after the combined Kalman refit, "
      "plus -- when PromoteUnspliced=true -- single-detector segments that "
      "had no candidate to splice against. Track::type is a composition "
      "bitmask: 1=SiTarget, 2=SiPad, 4=MTC, so 3=SiT+SiP, 6=SiP+MTC, "
      "7=full chain (a promoted segment keeps its own single bit)."};

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
  // ACTS logger level for the CKF/KF, as Acts::Logging::Level (0=VERBOSE ..
  // 4=WARNING). Diagnostic only: VERBOSE prints every surface the filter
  // visits, which is what settles "why did the CKF drop this seed".
  Gaudi::Property<int> m_actsLogLevel{
      this, "ActsLogLevel", static_cast<int>(Acts::Logging::WARNING),
      "Acts::Logging::Level for the CKF and KF loggers (0=VERBOSE, 4=WARNING)."};
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
      this, "SeedCompatRadius", 10.0,
      "UNUSED since the tracklet seeder replaced the transverse-position "
      "accumulator; superseded by SeedResidWindow. Kept declared so the "
      "existing steering files that set it still load."};
  Gaudi::Property<double> m_stripPitch{
      this, "SeedStripPitch", 0.0755,
      "UNUSED since the tracklet seeder replaced the transverse-position "
      "accumulator: the seed position now comes from a least-squares line fit "
      "to the supporting hits, not from a strip-frequency centroid. Kept "
      "declared so the existing steering files that set it still load."};

  // ---- Tracklet seeding ---------------------------------------------------
  // The lever arm is the one real tuning knob here, and it is bounded from
  // BOTH sides: too short and there is no direction resolution, too long and
  // multiple scattering in the MTC iron breaks the straight-line model. See
  // the measured residual table above findSeeds().
  Gaudi::Property<int> m_seedWindowLayers{
      this, "SeedWindowLayers", 6,
      "Layers spanned by one seeding window. 6 is ~400 mm in the MTC, where a "
      "1 GeV+ muon's straight-line residual is 2.8 mm (median) — a few SciFi "
      "channels. Raise it for a detector with less material."};
  Gaudi::Property<int> m_seedWindowStride{
      this, "SeedWindowStride", 3,
      "Layers between consecutive window starts. Below SeedWindowLayers the "
      "windows overlap, so a track crossing a window boundary is still found "
      "whole by the neighbouring window; duplicates are merged."};
  Gaudi::Property<int> m_seedMinLayers{
      this, "SeedMinLayers", 3,
      "Distinct layers a tracklet must span. 2 points always define a line, so "
      "3 is the smallest value that confirms rather than just fits."};
  Gaudi::Property<double> m_seedSlopeMax{
      this, "SeedSlopeMax", 1.5,
      "Half-range of the slope axis [dc/dz]. 1.5 is ~56 deg, comfortably "
      "beyond the GENIE final-state angular distribution (median 243 mrad for "
      "the most energetic charged particle of an event)."};
  Gaudi::Property<double> m_seedSlopeBin{
      this, "SeedSlopeBin", 0.01,
      "Slope bin [dc/dz]. Coupled to HoughBinSize: over a half-window of L the "
      "binning smears the intercept by SeedSlopeBin*L, which must stay below "
      "HoughBinSize or the peak spreads across bins and never reaches "
      "HoughMinVotes. Default: 0.01*200 mm = 2 mm < 5 mm."};
  // The four knobs above are all in units that only mean something RELATIVE TO
  // A DETECTOR'S LAYER PITCH, and the three detectors differ by almost an order
  // of magnitude: SiTarget 11 mm, SiPad 16.8 mm, MTC 80 mm. A single scalar
  // window of 4 layers is 240 mm in the MTC but 33 mm in the SiTarget, where it
  // cannot determine a slope at all. Same pattern as Chi2CutOffPerDet: empty
  // means "use the scalar", so a steering that sets neither is unchanged; a
  // non-empty vector must have 3 entries [SiTarget, SiPad, MTC].
  Gaudi::Property<std::vector<int>> m_seedWindowLayersPerDet{
      this, "SeedWindowLayersPerDet", {},
      "Per-detector override of SeedWindowLayers [SiTarget, SiPad, MTC]. Pick "
      "them so the window spans a comparable LENGTH in each detector, not a "
      "comparable number of layers."};
  Gaudi::Property<std::vector<int>> m_seedWindowStridePerDet{
      this, "SeedWindowStridePerDet", {},
      "Per-detector override of SeedWindowStride [SiTarget, SiPad, MTC]."};
  Gaudi::Property<std::vector<double>> m_seedResidWindowPerDet{
      this, "SeedResidWindowPerDet", {},
      "Per-detector override of SeedResidWindow [SiTarget, SiPad, MTC]. The "
      "MTC value is measured; a detector with 0.0755 mm strips has no business "
      "using the MTC's 8 mm."};
  Gaudi::Property<std::vector<double>> m_houghBinSizePerDet{
      this, "HoughBinSizePerDet", {},
      "Per-detector override of HoughBinSize [SiTarget, SiPad, MTC]."};
  Gaudi::Property<std::vector<double>> m_seedSlopeBinPerDet{
      this, "SeedSlopeBinPerDet", {},
      "Per-detector override of SeedSlopeBin [SiTarget, SiPad, MTC]. Keep the "
      "coupling: SeedSlopeBin * (window length / 2) < HoughBinSize."};

  Gaudi::Property<double> m_seedResidWindow{
      this, "SeedResidWindow", 5.0,
      "Residual [mm] within which a hit supports a tracklet. Sized from the "
      "measured straight-line residuals, not from the sensor pitch: at 6 MTC "
      "layers those are 2.8 mm median / 9.5 mm p90."};

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
  Gaudi::Property<bool> m_hitPurgeMaskBorder{
      this, "HitPurgeMaskBorder", false,
      "Also drop every measurement within HitPurgeWindow of a purged one (the "
      "DBSCAN 'border' of the dense core), so a shower blob is removed whole "
      "and only measurements isolated in their own layer survive. Off: only "
      "the over-populated cores go, and the rim of a shower stays in the pool."};

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
  // Per-detector overrides for the four fit gates above. Empty (the default)
  // means "use the scalar", so every existing steering keeps its behaviour
  // byte-for-byte; a non-empty vector must have 3 entries [SiTarget, SiPad, MTC].
  //
  // These exist because the gates are needed at very different strengths per
  // detector. Measured on 500 GENIE events: of the segment candidates the cuts
  // reject, 4048 are MTC killed by chi2/ndf against SegMaxChi2Tight=10, versus
  // 740 for SiTarget and 384 for SiPad. Those MTC candidates are not wild —
  // their chi2/ndf runs p10=14.5, median=32.6, p90=87.3 — so the single scalar
  // is set by what SiTarget/SiPad need and silently throttles the MTC, which
  // supplies ~77% of all segments in GENIE.
  Gaudi::Property<std::vector<double>> m_chi2CutOffPerDet{
      this, "Chi2CutOffPerDet", {},
      "Per-detector override of Chi2CutOff [SiTarget, SiPad, MTC]. Empty uses "
      "the scalar."};
  Gaudi::Property<std::vector<int>> m_numMeasCutOffPerDet{
      this, "NumMeasCutOffPerDet", {},
      "Per-detector override of NumMeasCutOff [SiTarget, SiPad, MTC]. Empty "
      "uses the scalar."};
  Gaudi::Property<std::vector<double>> m_segMaxChi2TightPerDet{
      this, "SegMaxChi2TightPerDet", {},
      "Per-detector override of SegMaxChi2Tight [SiTarget, SiPad, MTC]. Empty "
      "uses the scalar."};
  Gaudi::Property<std::vector<double>> m_segMaxChi2LoosePerDet{
      this, "SegMaxChi2LoosePerDet", {},
      "Per-detector override of SegMaxChi2Loose [SiTarget, SiPad, MTC]. Empty "
      "uses the scalar."};

  Gaudi::Property<std::vector<int>> m_minMeasPerSegment{
      this, "MinMeasPerSegment", {4, 3, 4},
      "Minimum measurements for a segment, per detector [SiTarget, SiPad, MTC]. "
      "SiPad hits are 2D, so 3 of them already carry 6 degrees of freedom."};
  Gaudi::Property<int> m_minMeasPerProjection{
      this, "MinMeasPerProjection", 2,
      "For the 1D-strip detector (SiTarget): a candidate must hold at least "
      "this many measurements on StripX planes AND on StripY planes. A "
      "segment measured in one projection only is a line with one free "
      "coordinate -- the fit puts it anywhere (GENIE event 1678: 34 Y strips, "
      "no X strip, fitted x wandering from -66 to +748 mm, drawn as a zigzag "
      "outside the target). 0 disables."};
  Gaudi::Property<double> m_backfillMaxChi2Growth{
      this, "BackfillMaxChi2Growth", 1.5,
      "The backward fill (hits upstream of the seed, see the PHASE 2 code) is "
      "kept only if the refit chi2/ndf of the union stays below this factor "
      "times the segment's own, or below the tight gate."};
  Gaudi::Property<double> m_stateBoundsTolerance{
      this, "StateBoundsTolerance", 5.0,
      "A candidate is rejected when any of its fitted measurement states lies "
      "outside its surface's bounds by more than this [mm]. On a 1D stereo "
      "surface the fit is free along the fibre, and a shower ghost of a few "
      "crossings plus one unpaired fibre can be 'fitted' with y = 411 mm "
      "in a 300 mm plane (GENIE event 228) -- drawn outside the detector. "
      "Negative disables."};
  Gaudi::Property<std::vector<double>> m_minMeasFracPerSegment{
      this, "MinMeasFracPerSegment", {0.01, 0.01, 0.01},
      "Minimum fraction of layers a segment must hold a measurement on, per "
      "detector [SiTarget, SiPad, MTC], counted between its OWN first and last "
      "hit layer: completeness = layersHit / (lastLayer - firstLayer + 1). "
      "Unlike MinMeasPerSegment this is a completeness measure, so it does NOT "
      "penalise a track that stops inside the detector (common in MTC for "
      "hadronic showers) — it penalises gaps, and hybrid segments stitched "
      "from two different particles' hits. A segment must pass BOTH cuts. "
      "0 does not mean 'no cut': it means DO NOT TRACK this detector at all "
      "(a 0% requirement is meaningless, which frees the value as a bypass "
      "switch); the detectors it disables are listed at initialize(). The 0.01 "
      "default is below 1/nLayers for every detector, i.e. deliberately unable "
      "to reject anything — the real per-pipeline values are meant to be set "
      "from the measured completeness distribution, which is logged for every "
      "candidate at DEBUG."};
  // MinMeasPerSegment is a FLAT floor, and a flat floor is not neutral with
  // respect to where a track starts. A segment entering at MTC layer 39 of 42
  // can hold at most 4 layers of measurements, so a floor of 3 eats most of
  // what is achievable; one entering at layer 0 has 42 layers ahead and the
  // same floor is irrelevant. GENIE vertices are spread over the full
  // instrumented depth, so this asymmetry is not a corner case here.
  //
  // With this on, the floor is capped by the layers still AVAILABLE downstream
  // of where the thing being judged actually starts -- the event's entry layer
  // for the detector bypass, the segment's own first layer for the segment
  // acceptance. It can only ever LOWER the floor, never raise it, so a track
  // that crosses the whole detector is judged exactly as before.
  //
  // Default false: the five other pipelines wiring this tracker keep their
  // behaviour byte-for-byte until they opt in.
  Gaudi::Property<bool> m_scaleMinMeasByDepth{
      this, "ScaleMinMeasByDepth", false,
      "Cap MinMeasPerSegment by the layers available downstream of the entry "
      "point, so a track that starts near the end of a detector is not judged "
      "against a floor it could never reach. Never raises the floor."};

  Gaudi::Property<bool> m_dumpSegmentCellIDs{
      this, "DumpSegmentCellIDs", false,
      "Log the cellID of every measurement of each accepted segment, so a "
      "segment can be joined offline against MC truth (the *HitsWindowed "
      "collections and their parallel *ContribParticleIDs frame parameter) to "
      "count how many distinct MC particles it is built from. That is the "
      "measurement that says whether a low-completeness segment is a genuine "
      "track with gaps or a hybrid of several particles — which the "
      "completeness value alone cannot distinguish. Diagnostic only: it "
      "changes no reconstruction, and it is verbose, so default false."};
  Gaudi::Property<int> m_entryPointMinHits{
      this, "EntryPointMinHits", 1,
      "Measurements a layer needs before it can be taken as this event's entry "
      "point into a detector. Above 1 it stops a single stray upstream hit "
      "from dragging the entry point back to the front of the detector, which "
      "would undo the navigator trim."};
  Gaudi::Property<int> m_dupSharedCellIDs{
      this, "DuplicateSharedCellIDs", 0,
      "Two segments of one detector are duplicates when they share at least "
      "this many CELL IDs -- physical hits, not measurement indices. 0 "
      "disables. The distinction is the whole point: the per-seed hit masking "
      "already makes the segments of one iteration DISJOINT in measurement "
      "index, so a threshold on those can never fire. But with "
      "MTCSciFiMeasConverter's AllPairs one physical U strip produces one "
      "measurement per V partner, so two disjoint measurement sets routinely "
      "rest on the same physical hits. Measured on GENIE events 0-5: 148 "
      "segment pairs share 1 cellID, 18 share 2, 4 share 3, and zero share a "
      "measurement index."};
  Gaudi::Property<std::vector<double>> m_dupProxDistPerDet{
      this, "DuplicateProxDistPerDet", {0.3, 6.0, 10.0},
      "PHASE 2.4, per detector [mm]: two segments of one detector are "
      "duplicates when, on at least DuplicateProxSurfaces surfaces, each holds "
      "a measurement within this distance of the other's. Catches what the "
      "cellID test cannot: the SAME particle rebuilt from its OTHER fibres. A "
      "muon crossing a SciFi plane at an angle lights two adjacent fibres, "
      "and AllPairs turns each into its own family of crossings; the first "
      "segment masks one family, the loose iteration rebuilds the track from "
      "the second (GENIE event 66: one muon, three MTC segments, two of them "
      "sharing zero cellIDs with the first). Unlike the per-iteration filter "
      "this runs ACROSS iterations, before the merge. For the MTC a pair also "
      "counts when it sits on the same or an adjacent fibre in EITHER view "
      "(|ds_U| or |ds_V| < 1.5 mm), whatever the transverse distance; the "
      "10 mm transverse window is below the plane's own y resolution "
      "(1 mm / 2 tan 5 deg = 5.7 mm per fibre). Empty disables."};
  Gaudi::Property<int> m_dupProxSurfaces{
      this, "DuplicateProxSurfaces", 3,
      "Minimum number of coincident surfaces for DuplicateProxDistPerDet."};
  Gaudi::Property<double> m_dupOverlapFraction{
      this, "DuplicateOverlapFraction", 0.7,
      "Two segments of one detector are duplicates when they share more than "
      "this fraction of the smaller one's hits; the worse chi2/ndf is dropped."};

  // ---- Intra-detector merging ---------------------------------------------
  // The per-seed hit masking makes the segments of one iteration disjoint in
  // hits, which is right, but it also means a particle whose track the CKF
  // drops halfway is rebuilt by the NEXT seed from the hits that are left --
  // as a second, non-overlapping segment of the same particle. No duplicate
  // filter can see that: there is nothing shared to compare. Measured on GENIE
  // event 3, MTC: 34 segments belonging to just 6 MC particles, one of them
  // holding 14 segments of 3-7 measurements each.
  //
  // This is the same operation PHASE 3 performs between detectors, applied
  // WITHIN one: propose by a geometric window at the junction, let a combined
  // refit decide. Default off so the five other pipelines are unchanged.
  Gaudi::Property<bool> m_mergeSegments{
      this, "MergeSegments", false,
      "Merge segments of the same detector that are contiguous in beam z and "
      "compatible at their junction. Off by default."};
  Gaudi::Property<double> m_mergeMaxGap{
      this, "MergeMaxGap", 500.0,
      "Maximum beam-direction gap [mm] between the end of the upstream segment "
      "and the start of the downstream one for a merge to be proposed."};
  // Fragments of one particle OVERLAP in z; they are not consecutive pieces.
  // Measured on 1000 GENIE events with the contiguity-only rule: 15012 pairs
  // rejected on the gap against 198 reaching the geometric test -- the rule
  // could only ever see 1% of the problem. The per-seed hit masking is what
  // produces this: the second seed rebuilds the particle from the hits the
  // first one left, which are interleaved along the same z range, not
  // downstream of it.
  Gaudi::Property<double> m_mergeMaxOverlap{
      this, "MergeMaxOverlap", 0.0,
      "How far [mm] the downstream segment may START UPSTREAM of the other's "
      "end and still be proposed for merging. 0 keeps the strict "
      "end-to-end rule; a value comparable to a detector's depth lets two "
      "segments covering the same z range be tested against each other. The "
      "combined refit still decides."};
  Gaudi::Property<double> m_mergeMaxDist{
      this, "MergeMaxDist", 30.0,
      "Transverse miss distance [mm] at the junction, straight-line "
      "extrapolated. Looser than the inter-detector SpliceMaxDist because the "
      "extrapolation crosses MTC iron."};
  Gaudi::Property<double> m_mergeMaxAngle{
      this, "MergeMaxAngle", 0.3,
      "Maximum angle [rad] between the two segments' directions at the "
      "junction. Generous: a 3 GeV muon changes slope by dtx~0.24 over ten MTC "
      "layers, so a tight cut would forbid merging exactly the tracks that "
      "fragment."};
  Gaudi::Property<double> m_mergeMaxChi2PerNdf{
      this, "MergeMaxChi2PerNdf", 50.0,
      "Combined-refit acceptance for a merge. This, not the geometric window, "
      "is what decides whether two fragments are one particle."};
  Gaudi::Property<int> m_mergeMaxPasses{
      this, "MergeMaxPasses", 4,
      "How many times the merging sweep is repeated, so a track broken into "
      "three or more pieces can be rebuilt one junction at a time."};

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

  // An unspliced segment normally stays in its own detector's local
  // collection. That is right when a splice was attempted and failed (a real
  // candidate existed in the neighbouring detector but missed the
  // dist/angle window or the combined refit) -- silently promoting THAT case
  // would hide genuine mismatches. It is wrong when the neighbouring
  // detector found no segments at all this event: there was nothing to even
  // attempt a splice against, so the segment is not "unmatched", it is the
  // whole track (e.g. a shower vertex that starts inside MTC, downstream of
  // SiTarget/SiPad's z-range -- expected to be common for GENIE interactions
  // spread across the full detector depth, unlike the single-vertex PG
  // samples this tracker was first tuned on).
  Gaudi::Property<bool> m_promoteUnspliced{
      this, "PromoteUnspliced", false,
      "When true, an unspliced segment is written to OutputGlobalTracks "
      "(Track::type = its own single-detector bit) instead of its local "
      "collection, IF EVERY detector it could have spliced through (its "
      "neighbour(s) in the SiTarget<->SiPad<->MTC chain) found zero "
      "segments this event -- i.e. there was no candidate to even attempt "
      "a splice against. A segment whose neighbour had candidates that "
      "simply failed the splice window or the combined refit is never "
      "promoted by this flag; it stays local, unchanged from today's "
      "behaviour. Default false: preserves the existing 5 pipelines "
      "wiring this tracker until they opt in."};

  // ---- Track extension into a neighbour (PHASE 2.7) -----------------------
  // A good segment whose neighbour detector holds nothing to splice with is
  // followed INTO that neighbour: its fitted end state, with the fit's own
  // (tight) covariance instead of a Hough seed's loose one, drives the
  // neighbour's CKF, and whatever it picks up becomes a small segment there,
  // which PHASE 3 then splices and refits like any other. This is what
  // recovers a muon that clips only two or three layers of the next detector
  // -- below that detector's MinMeasPerSegment, so never a standalone find,
  // but unambiguous when the track is already known.
  Gaudi::Property<bool> m_extendSegments{
      this, "ExtendSegments", false,
      "Enable PHASE 2.7: extend unmatched segments into the neighbouring "
      "detector(s) and let PHASE 3 splice the result."};
  Gaudi::Property<int> m_extendMinSourceMeas{
      this, "ExtendMinSourceMeas", 6,
      "A segment must hold at least this many measurements to be extended. "
      "Short segments carry too little direction to point into a neighbour."};
  Gaudi::Property<double> m_extendMaxSourceChi2PerNdf{
      this, "ExtendMaxSourceChi2PerNdf", 5.0,
      "... and fit at most this chi2/ndf."};
  Gaudi::Property<int> m_extendMinMeas{
      this, "ExtendMinMeas", 2,
      "Minimum measurements the neighbour's CKF must collect for the "
      "extension to be kept. Below the detector's MinMeasPerSegment on "
      "purpose: the source track is the evidence, not the hit count."};
  Gaudi::Property<double> m_extendChi2CutOff{
      this, "ExtendChi2CutOff", 30.0,
      "CKF per-measurement chi2 cut for the extension pass. Tighter than "
      "Chi2CutOff: the start parameters are a fitted track, not a seed."};
  Gaudi::Property<double> m_extendVarLoc{
      this, "ExtendVarLoc", 25.0,
      "Position variance [mm^2] given to the extension's start parameters "
      "on the neighbour's first surface. The straight-line hop across the "
      "boundary and the material there justify a few mm, not the fit's own "
      "sub-mm."};
  Gaudi::Property<double> m_extendVarAngle{
      this, "ExtendVarAngle", 4e-4,
      "Direction variance [rad^2] of the start parameters (20 mrad): the "
      "segment's own direction plus scattering at the boundary."};
  Gaudi::Property<bool> m_extendBackward{
      this, "ExtendBackward", true,
      "Also extend UPSTREAM (MTC -> SiPad, SiPad -> SiTarget), propagating "
      "backward from the segment's front state."};

  // ---- Hough seeding ------------------------------------------------------
  struct SeedCandidate {
    double x;        // DD4hep transverse x where the tracklet starts [mm]
    double y;        // DD4hep transverse y where the tracklet starts [mm]
    double z_start;  // beam coordinate the seed is defined at [mm]
    double tx;       // dx/dz — the whole point of the rewrite: seeds now
    double ty;       // dy/dz   carry a direction instead of assuming the beam
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
// findSeeds() — local straight-line tracklets, per detector.
//
// WHAT THIS REPLACED, AND WHY. The previous seeder accumulated a 2D histogram
// of the transverse positions (x, y) of the detector's points, discarding z
// entirely. A peak therefore meant "many hits at the same transverse position
// at different z", i.e. it could only find tracks PARALLEL TO THE BEAM, and
// the seed it produced carried no direction at all (the CKF was started with
// a hardcoded along-beam unit vector). That was tenable for the single-vertex
// particle-gun samples this tracker was first tuned on. It is not tenable for
// GENIE: measured over 2000 events, the most energetic charged final-state
// particle of an event has a median angle of 243 mrad to the beam and only
// 4.6% of events have one below 30 mrad — which is roughly where the old
// accumulator stopped working (a 5 mm bin over the MTC's 80 mm layer pitch
// tolerates ~60 mrad before consecutive hits fall in different bins).
//
// THE MODEL HERE IS A STRAIGHT LINE OVER A SHORT WINDOW OF LAYERS, in both
// projections, and it is deliberately NOT a curvature fit. Measured on 199
// MC-truth muons above 1 GeV crossing the MTC (median RMS residual of a
// straight-line fit, MC truth positions):
//
//     layers   lever    RMS y (median/p90)   RMS x (median/p90)
//        3     160 mm      0.51 / 5.74 mm       1.81 / 6.06 mm
//        4     240 mm      1.11 / 6.79 mm       3.27 / 8.01 mm
//        6     400 mm      2.83 / 9.52 mm       4.51 / 11.00 mm
//       10     720 mm      4.28 / 11.19 mm      6.41 / 14.87 mm
//       20    1520 mm      6.82 / 12.39 mm      9.95 / 16.16 mm
//
// y is the field-free projection (B is along DD4hep Y, so the Lorentz force
// q v x B has no y component and y-z is exactly straight): its residual is
// pure multiple scattering, and the MTC scatters hard — 50 mm of iron per
// 80 mm layer. x adds the bending on top. The two are COMPARABLE, not
// separable: at 6 layers the field's sagitta contribution (~2.4 mm at 3 GeV)
// is already below the scattering (2.8 mm). That is why there is no
// curvature term and no q/p estimate here — over any lever short enough for
// the association to be reliable, the curvature signal is buried in
// scattering, and over a lever long enough to measure it the scattering has
// destroyed the association. Momentum stays the CKF's job: it propagates in
// the real per-slab field (IronFieldRanges) starting from SeedMomentum.
//
// STRUCTURE. Per detector, per sliding window of SeedWindowLayers layers:
//   * 2D detectors (SiPad; MTC with PairMethod=AllPairs) give full (x, y, z)
//     points. A linear Hough in the field-free y-z projection finds the
//     tracklets, then a second linear Hough in x-z runs ONLY on the points
//     supporting a given y peak. Because both projections come from the same
//     point, the x<->y association is exact and there is no ghost problem.
//   * 1D-strip detectors (SiTarget; MTC with PairMethod=None) measure the two
//     projections on separate planes, so each is Hough'd independently and the
//     peaks are combined. This is where the old seeder built explicit
//     StripX x StripY crossings — N^2 of them for N tracks, needing the whole
//     GhostResolveWindow machinery and a hard 500-crossings-per-station cap
//     that silently dropped busy stations. Hough'ing the projections instead
//     moves the ambiguity from pairs of HITS to pairs of TRACK CANDIDATES,
//     which is a far smaller and better-posed problem.
//
// Cost is O(N * nSlopeBins) per window — linear in the number of points, which
// is what keeps dense hadronic showers affordable (single GENIE events reach
// 9000+ measurements).
// ---------------------------------------------------------------------------

namespace {

// One 1D projection measurement: coordinate c on a plane at beam position z.
struct ViewPoint {
  double      c;
  double      z;
  int         weight;
  int         layer;
  std::size_t src;    // index into the caller's 3D point list, when there is one
};

// A straight line c = c0 + t * (z - zRef) found in one projection.
struct LinePeak {
  double                   c0 = 0.0;
  double                   t  = 0.0;
  double                   zRef = 0.0;
  int                      votes = 0;
  int                      nLayers = 0;
  std::vector<std::size_t> support;   // indices into the ViewPoint vector
};

// Weighted least-squares c = a + b*(z - zRef) over the supporting points.
bool refineLine(const std::vector<ViewPoint>& v,
                const std::vector<std::size_t>& sup,
                double zRef, double& a, double& b) {
  double sw = 0, sz = 0, sc = 0, szz = 0, szc = 0;
  for (std::size_t i : sup) {
    const double w = static_cast<double>(v[i].weight);
    const double dz = v[i].z - zRef;
    sw += w; sz += w * dz; sc += w * v[i].c;
    szz += w * dz * dz; szc += w * dz * v[i].c;
  }
  const double det = sw * szz - sz * sz;
  if (sw <= 0.0 || std::abs(det) < 1e-9) return false;
  b = (sw * szc - sz * sc) / det;
  a = (sc - b * sz) / sw;
  return true;
}

// Linear Hough over (intercept at zRef, slope) for one projection.
//
// Out-of-range intercepts are DROPPED, not clamped into the edge bin. The old
// seeder clamped, which piled every track outside its +-HoughHalfSize window
// into the border bins and manufactured peaks there — with HoughHalfSize left
// at 200 mm against an MTC station half-height of 300 mm, that was a live
// defect, not a corner case.
std::vector<LinePeak> houghLine(const std::vector<ViewPoint>& v, double zRef,
                                double halfSize, double posBin, double slopeMax,
                                double slopeBin, int minVotes, double residWin,
                                int minLayers, int maxPeaks) {
  std::vector<LinePeak> peaks;
  if (v.size() < 2) return peaks;

  const int nT = std::max(1, static_cast<int>(2.0 * slopeMax / slopeBin) + 1);
  const int nC = std::max(1, static_cast<int>(2.0 * halfSize / posBin) + 1);
  std::vector<int> histo(static_cast<std::size_t>(nT) * nC, 0);

  auto cBin = [&](double c) -> int {
    const double f = (c + halfSize) / posBin;
    if (f < 0.0 || f >= static_cast<double>(nC)) return -1;
    return static_cast<int>(f);
  };

  for (int it = 0; it < nT; ++it) {
    const double t = -slopeMax + (it + 0.5) * slopeBin;
    int* row = histo.data() + static_cast<std::size_t>(it) * nC;
    for (const auto& p : v) {
      const int b = cBin(p.c - t * (p.z - zRef));
      if (b >= 0) row[b] += p.weight;
    }
  }

  struct Cell { int it, ic, votes; };
  std::vector<Cell> cells;
  for (int it = 0; it < nT; ++it) {
    for (int ic = 0; ic < nC; ++ic) {
      const int val = histo[static_cast<std::size_t>(it) * nC + ic];
      if (val < minVotes) continue;
      bool isMax = true;
      for (int dt = -1; dt <= 1 && isMax; ++dt) {
        for (int dc = -1; dc <= 1 && isMax; ++dc) {
          if (dt == 0 && dc == 0) continue;
          const int jt = it + dt, jc = ic + dc;
          if (jt < 0 || jt >= nT || jc < 0 || jc >= nC) continue;
          if (histo[static_cast<std::size_t>(jt) * nC + jc] > val) isMax = false;
        }
      }
      if (isMax) cells.push_back({it, ic, val});
    }
  }
  std::sort(cells.begin(), cells.end(),
            [](const Cell& a, const Cell& b) { return a.votes > b.votes; });

  std::vector<std::vector<std::size_t>> taken;
  for (const auto& cell : cells) {
    if (static_cast<int>(peaks.size()) >= maxPeaks) break;
    double t  = -slopeMax + (cell.it + 0.5) * slopeBin;
    double c0 = -halfSize + (cell.ic + 0.5) * posBin;

    // Collect support around the bin centre, refit, collect again: the bin is
    // coarse on purpose (it only has to group the hits) and the fit is what
    // actually sets the seed's position and direction.
    std::vector<std::size_t> sup;
    for (int pass = 0; pass < 2; ++pass) {
      sup.clear();
      for (std::size_t i = 0; i < v.size(); ++i)
        if (std::abs(v[i].c - (c0 + t * (v[i].z - zRef))) < residWin)
          sup.push_back(i);
      if (sup.size() < 2) break;
      double a = c0, b = t;
      if (!refineLine(v, sup, zRef, a, b)) break;
      c0 = a; t = b;
    }
    if (sup.size() < 2) continue;

    std::set<int> layers;
    int votes = 0;
    for (std::size_t i : sup) { layers.insert(v[i].layer); votes += v[i].weight; }
    if (votes < minVotes || static_cast<int>(layers.size()) < minLayers) continue;

    // Two Hough cells often describe the same tracklet; keep the first.
    bool dup = false;
    for (const auto& prev : taken) {
      std::size_t shared = 0;
      for (std::size_t i : sup)
        if (std::find(prev.begin(), prev.end(), i) != prev.end()) ++shared;
      if (2 * shared > std::min(prev.size(), sup.size())) { dup = true; break; }
    }
    if (dup) continue;

    LinePeak pk;
    pk.c0 = c0; pk.t = t; pk.zRef = zRef;
    pk.votes = votes; pk.nLayers = static_cast<int>(layers.size());
    pk.support = sup;
    taken.push_back(sup);
    peaks.push_back(std::move(pk));
  }
  return peaks;
}

}  // namespace

std::vector<ACTSSequentialTracker::SeedCandidate>
ACTSSequentialTracker::findSeeds(
    const std::vector<SeqMeasurement>& measurements,
    const std::vector<SeqMeasurement>& /*allMeas*/,
    const Acts::GeometryContext& gctx,
    int detMask) const
{
  std::vector<SeedCandidate> seeds;
  if (measurements.empty()) return seeds;

  // Empty override vector => the scalar, so existing steerings are unchanged.
  auto perDet = [detMask](const auto& vec, auto scalar) {
    return vec.empty() ? scalar : vec[detMask];
  };

  const double halfSize  = m_houghHalfSize.value();
  const double posBin    = perDet(m_houghBinSizePerDet.value(), m_houghBinSize.value());
  const double slopeMax  = m_seedSlopeMax.value();
  const double slopeBin  = perDet(m_seedSlopeBinPerDet.value(), m_seedSlopeBin.value());
  const int    minVotes  = m_houghMinVotes.value();
  const double residWin  = perDet(m_seedResidWindowPerDet.value(), m_seedResidWindow.value());
  const int    winLayers = std::max(2, perDet(m_seedWindowLayersPerDet.value(),
                                              m_seedWindowLayers.value()));
  const int    stride    = std::max(1, perDet(m_seedWindowStridePerDet.value(),
                                              m_seedWindowStride.value()));
  const int    minLayers = std::max(2, m_seedMinLayers.value());
  const int    maxS      = std::max(1, m_maxSeeds.value());
  const double maxMult   = m_houghMaxMultiplicity.value();

  // --- Points, split by how the detector measures them ---------------------
  struct P3 { double x, y, z; int weight; int layer; };
  std::vector<P3>        pts3;             // detectors giving a full (x,y,z)
  std::vector<ViewPoint> viewA, viewB;     // 1D strips: plane 0 and plane 1
  const bool stereo = (detMask == 2);      // MTC 1D planes are U/V, not X/Y

  const double isolWin  = m_isolationWindow.value();
  const int    isolMaxN = m_isolationMaxNeighbors.value();

  for (const auto& m : measurements) {
    if (m.detectorID != detMask) continue;
    const int    layer = m_geoSvc->layerOrdinalOf(m.surface);
    const double z     = m.surface->center(gctx).x();
    if (m.is2D) {
      pts3.push_back({m.localCoord, m.localCoord2, z, 2, layer});
    } else if (m.plane == 0) {
      viewA.push_back({m.localCoord, z, 1, layer, 0});
    } else if (m.plane == 1) {
      viewB.push_back({m.localCoord, z, 1, layer, 0});
    }
  }

  // Per-layer isolation filter on 2D points, carried over from the old seeder:
  // a point with too many same-layer neighbours is shower core, not a track.
  if (isolWin > 0.0 && !pts3.empty()) {
    const double w2 = isolWin * isolWin;
    std::vector<P3> kept;
    kept.reserve(pts3.size());
    for (std::size_t i = 0; i < pts3.size(); ++i) {
      int n = 0;
      for (std::size_t j = 0; j < pts3.size() && n <= isolMaxN; ++j) {
        if (i == j || pts3[j].layer != pts3[i].layer) continue;
        const double dx = pts3[j].x - pts3[i].x, dy = pts3[j].y - pts3[i].y;
        if (dx * dx + dy * dy < w2) ++n;
      }
      if (n <= isolMaxN) kept.push_back(pts3[i]);
    }
    pts3.swap(kept);
  }

  // --- Layer range the windows slide over -----------------------------------
  int loOrd = std::numeric_limits<int>::max(), hiOrd = std::numeric_limits<int>::min();
  auto noteLayer = [&](int l) {
    if (l < 0) return;
    loOrd = std::min(loOrd, l); hiOrd = std::max(hiOrd, l);
  };
  for (const auto& p : pts3)  noteLayer(p.layer);
  for (const auto& p : viewA) noteLayer(p.layer);
  for (const auto& p : viewB) noteLayer(p.layer);
  if (loOrd > hiOrd) return seeds;

  // A seed already found is not sought again in the next (overlapping) window.
  auto isDuplicate = [&](const SeedCandidate& s) {
    for (const auto& e : seeds) {
      const double dz = s.z_start - e.z_start;
      const double dx = s.x - (e.x + e.tx * dz);
      const double dy = s.y - (e.y + e.ty * dz);
      if (std::abs(s.tx - e.tx) < 2.0 * slopeBin &&
          std::abs(s.ty - e.ty) < 2.0 * slopeBin &&
          dx * dx + dy * dy < residWin * residWin) return true;
    }
    return false;
  };

  auto emit = [&](double x0, double y0, double z0, double tx, double ty,
                  int votes, int nLayers) {
    const double mult = (nLayers > 0)
        ? static_cast<double>(votes) / nLayers : 0.0;
    if (mult > maxMult) {
      debug() << "[ACTSSequentialTracker] det=" << detMask
              << " tracklet at (" << x0 << "," << y0 << "," << z0
              << ") -> SHOWER (multiplicity=" << mult << " > " << maxMult
              << "), skipped." << endmsg;
      return;
    }
    SeedCandidate s{x0, y0, z0, tx, ty, votes, mult};
    if (!isDuplicate(s)) seeds.push_back(s);
  };

  for (int wLo = loOrd; wLo <= hiOrd; wLo += stride) {
    const int wHi = wLo + winLayers - 1;

    // ---- 2D-point detectors: y first (field-free), then x on its support ---
    if (!pts3.empty()) {
      std::vector<ViewPoint> vy;
      double zSum = 0.0;
      for (std::size_t i = 0; i < pts3.size(); ++i) {
        const auto& p = pts3[i];
        if (p.layer < wLo || p.layer > wHi) continue;
        vy.push_back({p.y, p.z, p.weight, p.layer, i});
        zSum += p.z;
      }
      if (static_cast<int>(vy.size()) >= minLayers) {
        const double zRef = zSum / vy.size();
        for (const auto& yPk : houghLine(vy, zRef, halfSize, posBin, slopeMax,
                                         slopeBin, minVotes, residWin,
                                         minLayers, maxS)) {
          std::vector<ViewPoint> vx;
          for (std::size_t i : yPk.support) {
            const auto& p = pts3[vy[i].src];
            vx.push_back({p.x, p.z, p.weight, p.layer, vy[i].src});
          }
          for (const auto& xPk : houghLine(vx, zRef, halfSize, posBin, slopeMax,
                                           slopeBin, minVotes, residWin,
                                           minLayers, maxS)) {
            double z0 = std::numeric_limits<double>::max();
            for (std::size_t i : xPk.support) z0 = std::min(z0, vx[i].z);
            emit(xPk.c0 + xPk.t * (z0 - zRef), yPk.c0 + yPk.t * (z0 - zRef),
                 z0, xPk.t, yPk.t, std::min(xPk.votes, yPk.votes),
                 std::min(xPk.nLayers, yPk.nLayers));
          }
        }
      }
    }

    // ---- 1D-strip detectors: the two projections, then combined ------------
    if (!viewA.empty() && !viewB.empty()) {
      std::vector<ViewPoint> va, vb;
      double zSum = 0.0; std::size_t nz = 0;
      for (const auto& p : viewA)
        if (p.layer >= wLo && p.layer <= wHi) { va.push_back(p); zSum += p.z; ++nz; }
      for (const auto& p : viewB)
        if (p.layer >= wLo && p.layer <= wHi) { vb.push_back(p); zSum += p.z; ++nz; }
      if (nz > 0 && static_cast<int>(va.size()) >= minLayers &&
          static_cast<int>(vb.size()) >= minLayers) {
        const double zRef = zSum / nz;
        const auto aPk = houghLine(va, zRef, halfSize, posBin, slopeMax,
                                   slopeBin, minVotes, residWin, minLayers, maxS);
        const auto bPk = houghLine(vb, zRef, halfSize, posBin, slopeMax,
                                   slopeBin, minVotes, residWin, minLayers, maxS);
        // The x<->y (or U<->V) pairing a 1D-strip detector cannot resolve on
        // its own. Ranked by the weaker of the two projections and capped, so
        // a busy event cannot blow the combination up.
        const double tanStereo =
            std::tan(m_mtcStereoAngle.value() * M_PI / 180.0);
        for (const auto& a : aPk) {
          for (const auto& b : bPk) {
            double z0 = std::numeric_limits<double>::max();
            for (std::size_t i : a.support) z0 = std::min(z0, va[i].z);
            for (std::size_t i : b.support) z0 = std::min(z0, vb[i].z);
            const double ca = a.c0 + a.t * (z0 - zRef);
            const double cb = b.c0 + b.t * (z0 - zRef);
            double x0, y0, tx, ty;
            if (stereo) {
              // u = x - y*tan(theta), v = x + y*tan(theta)
              if (std::abs(tanStereo) < 1e-9) continue;
              x0 = 0.5 * (ca + cb);          tx = 0.5 * (a.t + b.t);
              y0 = (cb - ca) / (2.0 * tanStereo);
              ty = (b.t - a.t) / (2.0 * tanStereo);
            } else {
              x0 = ca; y0 = cb; tx = a.t; ty = b.t;
            }
            emit(x0, y0, z0, tx, ty, std::min(a.votes, b.votes),
                 std::min(a.nLayers, b.nLayers));
          }
        }
      }
    }
  }

  std::sort(seeds.begin(), seeds.end(),
            [](const SeedCandidate& a, const SeedCandidate& b) {
              return a.nVotes > b.nVotes;
            });
  if (static_cast<int>(seeds.size()) > maxS) seeds.resize(maxS);
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
    if (m_minMeasFracPerSegment.value().size() != 3) {
      error() << "[ACTSSequentialTracker] MinMeasFracPerSegment must have "
              << "exactly 3 entries [SiTarget, SiPad, MTC]." << endmsg;
      return StatusCode::FAILURE;
    }
    // A per-detector override of the wrong size would otherwise be indexed out
    // of bounds; refuse at startup rather than read past the end.
    {
      const std::array<std::pair<const char*, std::size_t>, 9> overrides{{
          {"Chi2CutOffPerDet",      m_chi2CutOffPerDet.value().size()},
          {"NumMeasCutOffPerDet",   m_numMeasCutOffPerDet.value().size()},
          {"SegMaxChi2TightPerDet", m_segMaxChi2TightPerDet.value().size()},
          {"SegMaxChi2LoosePerDet", m_segMaxChi2LoosePerDet.value().size()},
          {"SeedWindowLayersPerDet", m_seedWindowLayersPerDet.value().size()},
          {"SeedWindowStridePerDet", m_seedWindowStridePerDet.value().size()},
          {"SeedResidWindowPerDet",  m_seedResidWindowPerDet.value().size()},
          {"HoughBinSizePerDet",     m_houghBinSizePerDet.value().size()},
          {"SeedSlopeBinPerDet",     m_seedSlopeBinPerDet.value().size()}}};
      for (const auto& [name, sz] : overrides) {
        if (sz != 0 && sz != 3) {
          error() << "[ACTSSequentialTracker] " << name << " must be empty "
                  << "(use the scalar) or have exactly 3 entries "
                  << "[SiTarget, SiPad, MTC]; got " << sz << "." << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }

    info() << "[ACTSSequentialTracker] Layers per detector: SiTarget="
           << m_geoSvc->nLayersOf(0) << " SiPad=" << m_geoSvc->nLayersOf(1)
           << " MTC=" << m_geoSvc->nLayersOf(2) << endmsg;
    // A 0 in MinMeasFracPerSegment silently removes a whole detector from the
    // reconstruction, and "0 = no minimum" is the natural misreading. Say it
    // out loud at INFO so it can never be a surprise in the output.
    {
      static const std::array<const char*, 3> kDetNames{"SiTarget", "SiPad", "MTC"};
      for (int det = 0; det < 3; ++det) {
        if (m_minMeasFracPerSegment.value()[det] <= 0.0) {
          info() << "[ACTSSequentialTracker] " << kDetNames[det]
                 << " tracking is DISABLED (MinMeasFracPerSegment[" << det
                 << "]=0): no segments will be found in it." << endmsg;
        }
      }
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
                                  hit.getTime(), hit.getEDep(), cellID,
                                  static_cast<int>(i));
      }
    }

    const auto* spHits = m_siPadHandle->get();
    if (spHits) {
      for (std::size_t i = 0; i < spHits->size(); ++i) {
        const auto& hit   = (*spHits)[i];
        const auto& pos   = hit.getPosition();
        const auto& cov   = hit.getCovMatrix();
        const int   layer = hit.getQuality();
        const uint64_t cellID = hit.getCellID();
        const Acts::Surface* surf = m_geoSvc->surfaceByAddress(1, -1, layer, -1);
        if (!surf) {
          warning() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " SiPad hit " << i << " layer=" << layer
                    << " has no matching surface." << endmsg;
          continue;
        }
        measurements.emplace_back(surf, pos.x, pos.y, cov[0], cov[3], true, 1, -1,
                                  hit.getTime(), hit.getEDep(), cellID,
                                  static_cast<int>(i));
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
                                  hit.getTime(), hit.getEDep(), cellID,
                                  static_cast<int>(i));
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
          if (m_hitPurgeMaskBorder.value()) {
            // One pass: a survivor within `win` of any core hit goes too.
            // Borders do not propagate (no chaining), so an isolated track
            // hit two windows away from the blob is untouched.
            std::vector<std::size_t> core;
            for (std::size_t i : idxs) if (!keep[i]) core.push_back(i);
            if (core.empty()) continue;
            for (std::size_t i : idxs) {
              if (!keep[i]) continue;
              const auto& a = measurements[i];
              for (std::size_t j : core) {
                const auto& b = measurements[j];
                const double d = a.is2D
                    ? std::hypot(a.localCoord - b.localCoord,
                                 a.localCoord2 - b.localCoord2)
                    : std::abs(a.localCoord - b.localCoord);
                if (d < win) { keep[i] = false; break; }
              }
            }
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
    // PHASE 1.5 — event topology scan
    // =====================================================================
    // Where does the activity actually BEGIN in each detector this event? The
    // fit criteria are derived from that instead of assuming a vertex at the
    // front of the detector. That assumption held for the PG samples this
    // tracker was tuned on and does not hold for GENIE, whose vertices are
    // spread over the full instrumented depth: a real track often starts in
    // the middle of a detector, and the detectors upstream of it see nothing
    // at all.
    struct DetTopology {
      int  nHits        = 0;
      int  entryLayer   = -1;  // layer ordinal; -1 = no entry point this event
      int  exitLayer    = -1;
      int  entrySurfIdx = 0;   // index into surfacesOf(det) where entryLayer starts
      bool active       = false;
    };
    std::array<DetTopology, 3> topo{};
    {
      // std::map keeps the ordinals ascending, which is beam order, so the
      // first qualifying layer found is the entry point and the last is the exit.
      std::array<std::map<int, int>, 3> hitsPerLayer;
      for (const auto& m : measurements) {
        const int det = m.detectorID;
        if (det < 0 || det > 2) continue;
        const int ord = m_geoSvc->layerOrdinalOf(m.surface);
        if (ord < 0) continue;
        ++topo[det].nHits;
        ++hitsPerLayer[det][ord];
      }

      const int minLayerHits = std::max(1, m_entryPointMinHits.value());
      for (int det = 0; det < 3; ++det) {
        for (const auto& [ord, n] : hitsPerLayer[det]) {
          if (n < minLayerHits) continue;
          if (topo[det].entryLayer < 0) topo[det].entryLayer = ord;
          topo[det].exitLayer = ord;
        }

        const auto& detSurfaces = m_geoSvc->surfacesOf(det);
        if (topo[det].entryLayer >= 0) {
          // Ordinals were assigned in beam order and surfacesOf() is in beam
          // order, so the first surface reaching the entry ordinal is where the
          // navigator should start.
          for (std::size_t i = 0; i < detSurfaces.size(); ++i) {
            if (m_geoSvc->layerOrdinalOf(detSurfaces[i]) >= topo[det].entryLayer) {
              topo[det].entrySurfIdx = static_cast<int>(i);
              break;
            }
          }
        }

        // Two independent bypasses. Static: the pipeline turned this detector
        // off. Dynamic: nothing to fit here this event — no hits at all, or
        // fewer than a segment could ever need. The dynamic one is what makes
        // a shower starting inside the MTC skip SiTarget and SiPad by itself,
        // with nothing configured. Compared in MEASUREMENTS, not layers: one
        // SiTarget layer yields up to two (StripX + StripY), so a layer count
        // would wrongly disable a track entering near the detector's end.
        const double frac    = m_minMeasFracPerSegment.value()[det];
        int          minMeas = m_minMeasPerSegment.value()[det];
        // Same cap as the segment acceptance, applied at the event's entry
        // layer: a detector the track only clips at its downstream end must
        // not be bypassed for holding fewer hits than a full crossing would.
        if (m_scaleMinMeasByDepth.value() && topo[det].entryLayer >= 0) {
          const int avail = m_geoSvc->nLayersOf(det) - topo[det].entryLayer;
          minMeas = std::max(1, std::min(minMeas, avail));
        }
        topo[det].active = (frac > 0.0) && (topo[det].entryLayer >= 0) &&
                           (topo[det].nHits >= minMeas) && !detSurfaces.empty();
      }

      debug() << "[ACTSSequentialTracker] evt=" << evtNum << " topology:"
              << " SiTarget[hits=" << topo[0].nHits << " layers="
              << topo[0].entryLayer << ".." << topo[0].exitLayer
              << " of " << m_geoSvc->nLayersOf(0)
              << " active=" << topo[0].active << "]"
              << " SiPad[hits=" << topo[1].nHits << " layers="
              << topo[1].entryLayer << ".." << topo[1].exitLayer
              << " of " << m_geoSvc->nLayersOf(1)
              << " active=" << topo[1].active << "]"
              << " MTC[hits=" << topo[2].nHits << " layers="
              << topo[2].entryLayer << ".." << topo[2].exitLayer
              << " of " << m_geoSvc->nLayersOf(2)
              << " active=" << topo[2].active << "]" << endmsg;
    }

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
        states.push_back({pos.x(), SeqEndState{pos, dir, bv[Acts::eBoundQOverP], &surf}});
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
      const auto& allDetSurfaces = m_geoSvc->surfacesOf(det);
      if (allDetSurfaces.empty()) continue;

      if (!topo[det].active) {
        debug() << "[ACTSSequentialTracker] evt=" << evtNum << " det=" << det
                << " skipped: hits=" << topo[det].nHits
                << " entryLayer=" << topo[det].entryLayer
                << " frac=" << m_minMeasFracPerSegment.value()[det] << endmsg;
        continue;
      }

      // The detector as this event sees it: everything from the entry layer
      // downstream. Trimming the head off is what stops the CKF from starting
      // on a surface the track never crossed and dragging a tail of holes
      // behind its first real hit.
      const std::vector<const Acts::Surface*> detSurfaces(
          allDetSurfaces.begin() + topo[det].entrySurfIdx, allDetSurfaces.end());

      // Global indices of this detector's measurements, minus the ones already
      // taken by an accepted segment of the same detector.
      std::set<std::size_t> used;
      const int    minMeas = m_minMeasPerSegment.value()[det];
      const double minFrac = m_minMeasFracPerSegment.value()[det];

      // Per-detector override, falling back to the scalar when the vector is
      // empty (validated to be size 3 or empty at initialize()).
      auto perDet = [det](const auto& vec, auto scalar) {
        return vec.empty() ? scalar : vec[det];
      };
      const double chi2CutBase   = perDet(m_chi2CutOffPerDet.value(),
                                          m_chi2CutOff.value());
      const double segChi2Tight  = perDet(m_segMaxChi2TightPerDet.value(),
                                          m_segMaxChi2Tight.value());
      const double segChi2Loose  = perDet(m_segMaxChi2LoosePerDet.value(),
                                          m_segMaxChi2Loose.value());
      const int    numMeasCut    = perDet(m_numMeasCutOffPerDet.value(),
                                          m_numMeasCutOff.value());

      // See ScaleMinMeasByDepth. `firstLayer` is where the thing being judged
      // begins; the cap is the layer count still ahead of it. Never raises the
      // floor, and never drops below 1.
      auto effMinMeas = [&](int firstLayer) {
        if (!m_scaleMinMeasByDepth.value() || firstLayer < 0) return minMeas;
        const int avail = m_geoSvc->nLayersOf(det) - firstLayer;
        return std::max(1, std::min(minMeas, avail));
      };

      for (int iteration = 0; iteration < 2; ++iteration) {
        const double chi2Cut = (iteration == 0)
            ? chi2CutBase
            : chi2CutBase * m_looseChi2Scale.value();
        const double segMaxChi2 = (iteration == 0)
            ? segChi2Tight
            : segChi2Loose;

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
        if (seeds.empty()) {
          // Logged because this, not a rejected fit, is the most common way a
          // GENIE event ends up with no segment: it used to be a bare
          // `continue` and the detector then produced no output line at all,
          // which made the dominant failure mode invisible in the logs.
          debug() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " det=" << det << " iter=" << iteration
                  << " pool=" << pool.size()
                  << ": seeder found no tracklet, detector yields nothing."
                  << endmsg;
          continue;
        }

        // ---- CKF for this (detector, iteration) --------------------------
        SeqSourceLinkAccessor slAccessor;
        slAccessor.slinks = &slinks;

        Acts::MeasurementSelectorCuts measCuts;
        measCuts.chi2CutOff            = {chi2Cut};
        measCuts.numMeasurementsCutOff =
            {static_cast<std::size_t>(numMeasCut)};
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
                   Acts::getDefaultLogger("SeqCKF", static_cast<Acts::Logging::Level>(m_actsLogLevel.value())));

        Acts::CombinatorialKalmanFilterOptions<SeqTrackContainer> ckfOptions(
            gctx, m_mctx, std::cref(m_cctx), ckfExtensions, pOptions, true, true);
        // Same search, upstream: used to fill in the hits between the
        // detector's entry and the tracklet the seed came from.
        Acts::PropagatorPlainOptions pOptionsBwd = pOptions;
        pOptionsBwd.direction = Acts::Direction::Backward();
        Acts::CombinatorialKalmanFilterOptions<SeqTrackContainer> ckfOptionsBwd(
            gctx, m_mctx, std::cref(m_cctx), ckfExtensions, pOptionsBwd, true, true);

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

        std::vector<SeqSegment> found;

        for (std::size_t iSeed = 0; iSeed < seeds.size(); ++iSeed) {
          const auto& sd = seeds[iSeed];

          // The CKF starts where the TRACKLET starts, not at the detector's
          // entry layer. Seeds are now found in sliding windows, so a track
          // born in the middle of the detector gets a start surface in the
          // middle of it; starting every seed at the front left a tail of
          // holes in front of the first real hit and paid for it in
          // propagation and in every hole-counting quantity.
          const Acts::Surface* sfSeed = detSurfaces.front();
          for (const auto* sf : detSurfaces) {
            if (sf->center(gctx).x() >= sd.z_start - 1.0) { sfSeed = sf; break; }
          }
          const double zSeed = sfSeed->center(gctx).x();

          // DD4hep (x=transverse X, y=transverse Y, z=beam) -> ACTS, where the
          // geometry's 90-degree Y rotation makes ACTS x the beam axis:
          // ePos1 = dd_y, ePos2 = dd_x. A unit step along the beam therefore
          // moves the ACTS vector by (1, ty, tx) — same axis mapping as the
          // position, which the previous comment here got wrong while the
          // hardcoded along-beam direction made it moot.
          const double ddx = sd.x + sd.tx * (zSeed - sd.z_start);
          const double ddy = sd.y + sd.ty * (zSeed - sd.z_start);

          Acts::Vector4 seedPos4;
          seedPos4[Acts::ePos0] = zSeed;
          seedPos4[Acts::ePos1] = ddy;
          seedPos4[Acts::ePos2] = ddx;
          seedPos4[Acts::eTime] = 0.0;
          Acts::Vector3 seedDir(1.0, sd.ty, sd.tx);
          seedDir.normalize();

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
                    << ") mm tx=" << sd.tx << " ty=" << sd.ty
                    << ": CKF returned no candidate with >=3 measurements"
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
          // Refit `hits`, started from the fitted front state of `fromIdx`
          // carried to the most upstream surface among `hits`. Returns the new
          // index and its chi2/ndf.
          auto refitFrom = [&](const std::set<std::size_t>& hits,
                               Acts::TrackIndexType fromIdx)
              -> std::optional<std::pair<Acts::TrackIndexType, double>> {
            const Acts::Surface* sfFront = nullptr;
            for (std::size_t mIdx : hits) {
              const auto* sf = measurements[mIdx].surface;
              if (sfFront == nullptr ||
                  sf->center(gctx).x() < sfFront->center(gctx).x())
                sfFront = sf;
            }
            auto ends = endStatesOf(tracks.getTrack(fromIdx));
            if (!ends || sfFront == nullptr) return std::nullopt;
            auto start = paramsAtSurface(ends->first, sfFront, refitCov);
            if (!start) return std::nullopt;
            auto rIdx = refit(hits, *start, detSurfaces);
            if (!rIdx) return std::nullopt;
            return std::make_pair(*rIdx, chi2NdfOf(tracks.getTrack(*rIdx)));
          };
          double refitC2n = bestC2n;
          if (auto r = refitFrom(fp, bestIdx)) {
            bestIdx  = r->first;
            refitC2n = r->second;
            fp       = fingerprintOf(tracks.getTrack(bestIdx));
          }

          // ---- Backward fill: the hits upstream of the seed -------------
          // The CKF only ever runs downstream from the tracklet's start, so
          // a seed found deep in the detector left everything in front of it
          // on the table; the loose iteration then built a second, short
          // segment there, or nothing. GENIE event 895: the muon's MTC
          // segment began at layer 14, 1.1 m of magnetised iron after the
          // MTC's front, and its back-propagated front state missed the SiPad
          // by 58 mm -- the splice failed and a ghost took its place. Run the
          // same CKF backward from the fitted front state, refit the union,
          // and keep it only if the fit does not degrade: the upstream end
          // of a muon's segment is often where it leaves its own vertex
          // shower, and an unconditional union there cost 3.5 points of muon
          // efficiency while buying the cleanest multiplicity yet.
          if (auto ends = endStatesOf(tracks.getTrack(bestIdx))) {
            if (ends->first.surface != nullptr &&
                ends->first.surface != detSurfaces.front()) {
              if (auto bstart = paramsAtSurface(ends->first, ends->first.surface,
                                                seedCov)) {
                if (auto bres = ckf.findTracks(*bstart, ckfOptionsBwd, tracks);
                    bres.ok()) {
                  std::set<std::size_t> bestUnion = fp;
                  for (const auto& t : *bres) {
                    auto bfp = fingerprintOf(t);
                    bfp.insert(fp.begin(), fp.end());
                    if (bfp.size() > bestUnion.size()) bestUnion = std::move(bfp);
                  }
                  if (bestUnion.size() > fp.size()) {
                    auto r = refitFrom(bestUnion, bestIdx);
                    const bool take = r && r->second <=
                        std::max(refitC2n * m_backfillMaxChi2Growth.value(),
                                 segChi2Tight);
                    debug() << "[ACTSSequentialTracker] evt=" << evtNum
                            << " det=" << det << " iter=" << iteration
                            << " seed=" << iSeed << " backward fill: +"
                            << (bestUnion.size() - fp.size()) << " measurements, chi2/ndf "
                            << refitC2n << " -> " << (r ? r->second : -1.0)
                            << (take ? " taken" : " rejected") << endmsg;
                    if (take) {
                      bestIdx  = r->first;
                      refitC2n = r->second;
                      fp       = fingerprintOf(tracks.getTrack(bestIdx));
                    }
                  }
                }
              }
            }
          }

          auto finalTrack = tracks.getTrack(bestIdx);
          const int    ndf  = std::max(
              1, static_cast<int>(finalTrack.nDoF()) - kHelixParams);
          const double chi2 = finalTrack.chi2();
          const double c2n  = chi2 / static_cast<double>(ndf);
          const int    nMeas = static_cast<int>(finalTrack.nMeasurements());

          // Completeness: of the layers between this segment's OWN first and
          // last hit, how many does it actually hold a measurement on? The span
          // being the segment's own is the point — a track that stops inside
          // the detector is not penalised (common in MTC for hadronic showers),
          // while gaps are, and a hybrid segment stitched from two different
          // particles' hits is mostly gaps.
          std::set<int> segLayers;
          for (std::size_t mIdx : fp) {
            const int ord = m_geoSvc->layerOrdinalOf(measurements[mIdx].surface);
            if (ord >= 0) segLayers.insert(ord);
          }
          const int segSpan = segLayers.empty()
              ? 0 : (*segLayers.rbegin() - *segLayers.begin() + 1);
          // Left at 1.0 when the span is unavailable, so a geometry lookup that
          // failed can never reject a candidate on its own.
          const double completeness = (segSpan > 0)
              ? static_cast<double>(segLayers.size()) / segSpan
              : 1.0;

          // Logged rather than silently dropped: "detector X found no segment"
          // is the single most common reason a global track fails to form
          // (a missing SiPad segment forbids BOTH of its neighbour links), so
          // the cut that rejected the candidate has to be visible.
          // Capped by what this segment could have reached from where it
          // actually starts, not by the detector's full depth.
          const int minMeasHere = effMinMeas(
              segLayers.empty() ? -1 : *segLayers.begin());

          // See MinMeasPerProjection.
          if (det == 0 && m_minMeasPerProjection.value() > 0) {
            int nX = 0, nY = 0;
            for (std::size_t mIdx : fp)
              (measurements[mIdx].plane == 0 ? nX : nY)++;
            if (nX < m_minMeasPerProjection.value() ||
                nY < m_minMeasPerProjection.value()) {
              debug() << "[ACTSSequentialTracker] evt=" << evtNum
                      << " det=" << det << " iter=" << iteration
                      << " seed=" << iSeed << " candidate rejected: single "
                      << "projection (X " << nX << ", Y " << nY << ")" << endmsg;
              continue;
            }
          }

          // See StateBoundsTolerance.
          int nOutside = 0;
          if (m_stateBoundsTolerance.value() >= 0.0) {
            const auto tol = Acts::BoundaryTolerance::AbsoluteEuclidean(
                m_stateBoundsTolerance.value());
            for (const auto& ts : finalTrack.trackStatesReversed()) {
              if (!ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) continue;
              if (ts.typeFlags().test(Acts::TrackStateFlag::OutlierFlag)) continue;
              Acts::BoundVector bv;
              if      (ts.hasSmoothed()) bv = ts.smoothed();
              else if (ts.hasFiltered()) bv = ts.filtered();
              else continue;
              if (!ts.referenceSurface().bounds().inside(
                      Acts::Vector2(bv[Acts::eBoundLoc0], bv[Acts::eBoundLoc1]), tol))
                ++nOutside;
            }
          }
          if (nOutside > 0) {
            debug() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " det=" << det << " iter=" << iteration
                    << " seed=" << iSeed << " candidate rejected: " << nOutside
                    << " fitted state(s) outside the surface bounds" << endmsg;
            continue;
          }
          if (nMeas < minMeasHere || c2n > segMaxChi2 || completeness < minFrac) {
            debug() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " det=" << det << " iter=" << iteration
                    << " seed=" << iSeed << " candidate rejected: nMeas="
                    << nMeas << " (min " << minMeasHere << ") chi2/ndf=" << c2n
                    << " (max " << segMaxChi2 << ") completeness="
                    << completeness << " (" << segLayers.size() << "/" << segSpan
                    << ", min " << minFrac << ")" << endmsg;
            continue;
          }
          debug() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " det=" << det << " iter=" << iteration
                  << " seed=" << iSeed << " candidate accepted: nMeas=" << nMeas
                  << " chi2/ndf=" << c2n << " completeness=" << completeness
                  << " (" << segLayers.size() << "/" << segSpan << ")" << endmsg;

          if (m_dumpSegmentCellIDs.value()) {
            // One line per segment, parseable: the truth join is done offline
            // against the *ContribParticleIDs frame parameter, keyed by cellID.
            std::ostringstream cids;
            for (std::size_t mIdx : fp) cids << ' ' << measurements[mIdx].cellID;
            info() << "[SegDump] evt=" << evtNum << " det=" << det
                   << " completeness=" << completeness
                   << " nMeas=" << nMeas << " cellIDs:" << cids.str() << endmsg;
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
        // Best-first, and "best" is stated here rather than inherited from the
        // seed ordering. The duplicate filter below is first-wins, so whatever
        // order `found` is in decides which of two overlapping segments is
        // kept. That used to be the Hough peak's vote count, which was
        // harmless when one pass produced one seed per track. It is not
        // harmless now: the seeder sweeps OVERLAPPING WINDOWS, so one particle
        // is normally seeded several times, and a seed from a window deep in
        // the detector can carry more votes than the one at the entry while
        // producing a strictly shorter segment. Order by measurement count,
        // then by fit quality, so the most complete segment is the one that
        // survives and masks its hits.
        std::sort(found.begin(), found.end(),
                  [](const SeqSegment& a, const SeqSegment& b) {
                    if (a.fp.size() != b.fp.size()) return a.fp.size() > b.fp.size();
                    const double ca = a.chi2 / std::max(1, a.ndf);
                    const double cb = b.chi2 / std::max(1, b.ndf);
                    return ca < cb;
                  });

        const double dupFrac  = m_dupOverlapFraction.value();
        const int    dupCells  = m_dupSharedCellIDs.value();
        // Physical hits behind a segment. See DuplicateSharedCellIDs for why
        // the measurement index is the wrong unit to compare.
        auto cellsOf = [&](const SeqSegment& sg) {
          std::set<uint64_t> c;
          for (std::size_t h : sg.fp) c.insert(measurements[h].cellID);
          return c;
        };
        std::vector<std::set<uint64_t>> keptCells;
        std::vector<SeqSegment> survivors;
        for (auto& cand : found) {
          const auto candCells = cellsOf(cand);
          bool dup = false;
          for (std::size_t ki = 0; ki < survivors.size(); ++ki) {
            const auto& keep = survivors[ki];
            std::size_t nShared = 0;
            for (std::size_t h : cand.fp) if (keep.fp.count(h)) ++nShared;
            std::size_t nCells = 0;
            if (dupCells > 0)
              for (uint64_t c : candCells)
                if (keptCells[ki].count(c)) ++nCells;
            const double smaller =
                static_cast<double>(std::min(cand.fp.size(), keep.fp.size()));
            if ((dupCells > 0 && static_cast<int>(nCells) >= dupCells) ||
                (smaller > 0 && nShared / smaller > dupFrac)) {
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
          if (!dup) {
            keptCells.push_back(cellsOf(cand));
            survivors.push_back(std::move(cand));
          }
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

    // ---- Boundary matching through the field -----------------------------
    // How well does the upstream segment's downstream end meet the
    // downstream segment's upstream end? Judged by PROPAGATING the downstream
    // segment's front state backward, in the real field, onto the surface the
    // upstream segment ends on, and comparing there. The old straight-line
    // meeting at the boundary plane was exact for SiTarget|SiPad (no field)
    // and wrong for SiPad|MTC: an MTC segment starting a few layers in has
    // curved through that much iron, its front direction extrapolated
    // straight back misses by tens of mm (GENIE event 895: 96 mm, so the
    // muon's own 24-hit MTC segment failed the 15 mm window and a ghost
    // family of the same fibres got spliced instead). The downstream segment
    // is the one propagated because it is the one that MEASURED the momentum
    // (the MTC has the field; the silicon's q/p is a prior).
    using SeqFreePropagator = Acts::Propagator<SeqStepper, Acts::VoidNavigator>;
    SeqFreePropagator freeProp{SeqStepper(bField), Acts::VoidNavigator{}};
    auto propagateTo = [&](const SeqEndState& from, const Acts::Surface* target,
                           Acts::Direction direction) -> std::optional<SeqEndState> {
      if (from.surface == nullptr || target == nullptr) return std::nullopt;
      Acts::Vector4 pos4(from.pos.x(), from.pos.y(), from.pos.z(), 0.0);
      auto start = Acts::BoundTrackParameters::create(
          gctx, from.surface->getSharedPtr(), pos4, from.dir, from.qOverP,
          std::nullopt, Acts::ParticleHypothesis::muon());
      if (!start.ok()) return std::nullopt;
      SeqFreePropagator::Options<> opts(gctx, m_mctx);
      opts.direction = direction;
      opts.maxSteps  = static_cast<std::size_t>(m_maxPropSteps.value());
      opts.stepping.maxStepSize = m_maxStepSize.value();
      auto res = freeProp.propagate(*start, *target, opts);
      if (!res.ok() || !res->endParameters) return std::nullopt;
      const auto& ep = *res->endParameters;
      return SeqEndState{ep.position(gctx), ep.direction(), ep.qOverP(), target};
    };
    // (dist [mm] transverse, angle [rad]) at the upstream segment's last
    // surface; straight-line fallback if the propagation fails.
    auto matchAtBoundary = [&](const SeqEndState& upBack, const SeqEndState& dnFront)
        -> std::optional<std::pair<double, double>> {
      if (auto p = propagateTo(dnFront, upBack.surface, Acts::Direction::Backward())) {
        const double dist  = std::hypot(p->pos.y() - upBack.pos.y(),
                                        p->pos.z() - upBack.pos.z());
        const double angle = std::acos(std::clamp(p->dir.dot(upBack.dir), -1.0, 1.0));
        return std::make_pair(dist, angle);
      }
      if (std::abs(upBack.dir.x()) < 1e-6 || std::abs(dnFront.dir.x()) < 1e-6)
        return std::nullopt;
      const double xB = 0.5 * (upBack.pos.x() + dnFront.pos.x());
      const Acts::Vector3 a = upBack.pos  + upBack.dir  * ((xB - upBack.pos.x())  / upBack.dir.x());
      const Acts::Vector3 b = dnFront.pos + dnFront.dir * ((xB - dnFront.pos.x()) / dnFront.dir.x());
      return std::make_pair(std::hypot(a.y() - b.y(), a.z() - b.z()),
                            std::acos(std::clamp(upBack.dir.dot(dnFront.dir), -1.0, 1.0)));
    };

    // =====================================================================
    // PHASE 2.4 — duplicates across iterations, BEFORE the merge
    // =====================================================================
    // See DuplicateProxDistPerDet. Best-first within each detector: longest,
    // then best chi2/ndf, survives; anything that runs along it on
    // DuplicateProxSurfaces surfaces, or shares DuplicateSharedCellIDs
    // physical hits with it, goes. Runs BEFORE the merge: the per-iteration
    // filter cannot see a loose-pass rebuild of a tight-pass segment, and the
    // merge would happily glue the two copies into one "26-measurement"
    // segment with every layer counted twice (GENIE event 1678).
    int nProxDup = 0;
    if (!m_dupProxDistPerDet.value().empty() && segments.size() > 1) {
      const int minSurf = m_dupProxSurfaces.value();
      std::vector<std::size_t> order(segments.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        const auto& A = segments[a]; const auto& B = segments[b];
        if (A.fp.size() != B.fp.size()) return A.fp.size() > B.fp.size();
        return A.chi2 / std::max(1, A.ndf) < B.chi2 / std::max(1, B.ndf);
      });
      // Per surface, the measurements of a segment.
      auto bySurface = [&](const SeqSegment& sg) {
        std::map<const Acts::Surface*, std::vector<std::size_t>> m;
        for (std::size_t h : sg.fp) m[measurements[h].surface].push_back(h);
        return m;
      };
      // "Close" means "on the same physical channel(s)". For the MTC that is
      // NOT a transverse distance: two crossings built on the same V fibre
      // and adjacent U fibres sit 0.5 mm apart in x but 5.7 mm apart in y
      // (1 mm / 2 tan 5 deg), so the test is done in fibre space -- the
      // stereo coordinates s_U = x cos(a) - y sin(a), s_V = x cos(a) + y
      // sin(a) -- and sharing EITHER fibre within `d` counts.
      const double alpha = m_mtcStereoAngle.value() * M_PI / 180.0;
      const double ca = std::cos(alpha), sa = std::sin(alpha);
      auto close = [&](const SeqMeasurement& a, const SeqMeasurement& b, double d) {
        if (a.detectorID == 2 && a.is2D) {
          const double dU = std::abs((a.localCoord - b.localCoord) * ca -
                                     (a.localCoord2 - b.localCoord2) * sa);
          const double dV = std::abs((a.localCoord - b.localCoord) * ca +
                                     (a.localCoord2 - b.localCoord2) * sa);
          if (std::min(dU, dV) < 1.5) return true;   // same/adjacent fibre
          return std::hypot(a.localCoord - b.localCoord,
                            a.localCoord2 - b.localCoord2) < d;
        }
        if (a.is2D)
          return std::hypot(a.localCoord - b.localCoord,
                            a.localCoord2 - b.localCoord2) < d;
        return std::abs(a.localCoord - b.localCoord) < d;
      };
      std::vector<char> dropped(segments.size(), 0);
      std::vector<std::size_t> kept;
      for (std::size_t ci : order) {
        const auto& cand = segments[ci];
        const double dist = m_dupProxDistPerDet.value()[cand.det];
        const auto candBy = bySurface(cand);
        for (std::size_t ki : kept) {
          const auto& keep = segments[ki];
          if (keep.det != cand.det) continue;
          const auto keepBy = bySurface(keep);
          int nCoinc = 0;
          for (const auto& [surf, hs] : candBy) {
            auto it = keepBy.find(surf);
            if (it == keepBy.end()) continue;
            bool hit = false;
            for (std::size_t a : hs) {
              for (std::size_t b : it->second)
                if (close(measurements[a], measurements[b], dist)) { hit = true; break; }
              if (hit) break;
            }
            if (hit && ++nCoinc >= minSurf) break;
          }
          std::size_t nCells = 0;
          if (m_dupSharedCellIDs.value() > 0) {
            std::set<uint64_t> kc;
            for (std::size_t h : keep.fp) kc.insert(measurements[h].cellID);
            std::set<uint64_t> seen;
            for (std::size_t h : cand.fp) {
              const uint64_t cid = measurements[h].cellID;
              if (kc.count(cid) && seen.insert(cid).second) ++nCells;
            }
          }
          if (nCoinc >= minSurf ||
              (m_dupSharedCellIDs.value() > 0 &&
               static_cast<int>(nCells) >= m_dupSharedCellIDs.value())) {
            debug() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " det=" << cand.det << " segment (" << cand.fp.size()
                    << "m, iter " << cand.iteration << ") dropped: runs within "
                    << dist << " mm of a kept " << keep.fp.size() << "m segment on "
                    << nCoinc << " surfaces, " << nCells << " shared cellIDs" << endmsg;
            dropped[ci] = 1;
            ++nProxDup;
            break;
          }
        }
        if (!dropped[ci]) kept.push_back(ci);
      }
      if (nProxDup > 0) {
        std::vector<SeqSegment> left;
        for (std::size_t i = 0; i < segments.size(); ++i)
          if (!dropped[i]) left.push_back(std::move(segments[i]));
        segments.swap(left);
      }
    }

    // =====================================================================
    // PHASE 2.5 — merge contiguous segments of the SAME detector
    // =====================================================================
    int nMerged = 0;
    if (m_mergeSegments.value() && segments.size() > 1) {
      const double maxGap     = m_mergeMaxGap.value();
      const double maxOverlap = m_mergeMaxOverlap.value();
      const double maxDist  = m_mergeMaxDist.value();
      const double maxAngle = m_mergeMaxAngle.value();
      const double maxC2n   = m_mergeMaxChi2PerNdf.value();
      const double angleW   = m_spliceAngleWeight.value();

      for (int pass = 0; pass < m_mergeMaxPasses.value(); ++pass) {
        struct MergeCand { std::size_t up, dn; double score; };
        std::vector<MergeCand> cands;
        int nGapRejected = 0;
        for (std::size_t i = 0; i < segments.size(); ++i) {
          for (std::size_t j = 0; j < segments.size(); ++j) {
            if (i == j || segments[i].det != segments[j].det) continue;
            const auto& A = segments[i];
            const auto& B = segments[j];
            const double gap = B.front.pos.x() - A.back.pos.x();
            if (gap < -maxOverlap || gap > maxGap) { ++nGapRejected; continue; }
            // Overlapping pairs are symmetric: (i,j) and (j,i) describe the
            // same merge, so only keep the ordering that starts furthest
            // upstream and halve the work.
            if (gap < 0.0 && A.front.pos.x() > B.front.pos.x()) continue;
            if (std::abs(A.back.dir.x()) < 1e-6) continue;
            const Acts::Vector3 at =
                A.back.pos + A.back.dir * (gap / A.back.dir.x());
            const double dist =
                std::hypot(at.y() - B.front.pos.y(), at.z() - B.front.pos.z());
            const double ang = std::acos(std::clamp(
                A.back.dir.dot(B.front.dir), -1.0, 1.0));
            // Every proposed pair is logged, passing or not: "the merge did
            // not happen" is otherwise indistinguishable from "no pair was
            // ever adjacent", and those need different fixes.
            debug() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " det=" << segments[i].det << " merge? seg" << i
                    << "(" << A.fp.size() << "m) -> seg" << j
                    << "(" << B.fp.size() << "m)"
                    << " gap=" << gap << " dist=" << dist
                    << " (max " << maxDist << ") angle=" << ang
                    << " (max " << maxAngle << ")"
                    << ((dist > maxDist || ang > maxAngle) ? " rejected"
                                                           : " CANDIDATE")
                    << endmsg;
            if (dist > maxDist || ang > maxAngle) continue;
            cands.push_back({i, j, dist * dist + angleW * ang * ang});
          }
        }
        debug() << "[ACTSSequentialTracker] evt=" << evtNum
                << " merge pass " << pass << ": " << segments.size()
                << " segments, " << nGapRejected
                << " pairs rejected on the z gap, " << cands.size()
                << " candidates" << endmsg;
        if (cands.empty()) break;
        std::sort(cands.begin(), cands.end(),
                  [](const MergeCand& a, const MergeCand& b) {
                    return a.score < b.score;
                  });

        std::vector<char> consumed(segments.size(), 0);
        std::vector<SeqSegment> mergedOut;
        for (const auto& mc : cands) {
          if (consumed[mc.up] || consumed[mc.dn]) continue;
          const auto& A = segments[mc.up];
          const auto& B = segments[mc.dn];
          std::set<std::size_t> fp = A.fp;
          fp.insert(B.fp.begin(), B.fp.end());

          const auto& allSurf = m_geoSvc->surfacesOf(A.det);
          if (allSurf.empty()) continue;
          const std::vector<const Acts::Surface*> surf(
              allSurf.begin() + topo[A.det].entrySurfIdx, allSurf.end());
          if (surf.empty()) continue;

          auto start = paramsAtSurface(A.front, surf.front(), refitCov);
          if (!start) continue;
          auto idx = refit(fp, *start, surf);
          if (!idx) continue;
          auto t = tracks.getTrack(*idx);
          const double c2n = chi2NdfOf(t);
          if (c2n > maxC2n) continue;
          auto ends = endStatesOf(t);
          if (!ends) continue;

          SeqSegment m;
          m.det       = A.det;
          m.iteration = std::min(A.iteration, B.iteration);
          m.idx       = *idx;
          m.fp        = std::move(fp);
          m.chi2      = t.chi2();
          m.ndf       = std::max(1, static_cast<int>(t.nDoF()) - kHelixParams);
          m.front     = ends->first;
          m.back      = ends->second;
          m.seedX     = A.seedX;
          m.seedY     = A.seedY;
          debug() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " det=" << A.det << " merged " << A.fp.size() << "+"
                  << B.fp.size() << " -> " << m.fp.size()
                  << " measurements, chi2/ndf=" << c2n << endmsg;
          consumed[mc.up] = consumed[mc.dn] = 1;
          mergedOut.push_back(std::move(m));
          ++nMerged;
        }
        if (mergedOut.empty()) break;
        for (std::size_t i = 0; i < segments.size(); ++i)
          if (!consumed[i]) mergedOut.push_back(std::move(segments[i]));
        segments.swap(mergedOut);
      }
    }

    // =====================================================================
    // PHASE 2.7 — extend unmatched segments into their neighbours
    // =====================================================================
    // See the ExtendSegments property. Runs after the merge so a fragmented
    // segment has been made whole first, and before the splice so its output
    // is judged by the same combined refit as everything else.
    int nExtended = 0;
    if (m_extendSegments.value() && !segments.empty()) {
      const int    minSrc   = m_extendMinSourceMeas.value();
      const double maxSrcC2 = m_extendMaxSourceChi2PerNdf.value();
      const int    minExt   = m_extendMinMeas.value();
      const double maxDistE = m_spliceMaxDist.value();
      const double maxAngE  = m_spliceMaxAngle.value();

      Acts::BoundSquareMatrix extCov = Acts::BoundSquareMatrix::Zero();
      extCov(Acts::eBoundLoc0,   Acts::eBoundLoc0)   = m_extendVarLoc.value();
      extCov(Acts::eBoundLoc1,   Acts::eBoundLoc1)   = m_extendVarLoc.value();
      extCov(Acts::eBoundPhi,    Acts::eBoundPhi)    = m_extendVarAngle.value();
      extCov(Acts::eBoundTheta,  Acts::eBoundTheta)  = m_extendVarAngle.value();
      extCov(Acts::eBoundQOverP, Acts::eBoundQOverP) = m_seedVarQOverP.value();
      extCov(Acts::eBoundTime,   Acts::eBoundTime)   = 1e9;

      // Measurements already owned by any segment or by an earlier extension.
      std::set<std::size_t> owned;
      for (const auto& sg : segments) owned.insert(sg.fp.begin(), sg.fp.end());

      // Does `nb` already hold a segment that PHASE 3 could splice to `src`
      // across this boundary? Then the extension is not needed and would
      // only compete with it.
      auto hasSpliceMate = [&](const SeqSegment& src, int nb, bool forward) {
        for (const auto& other : segments) {
          if (other.det != nb) continue;
          auto m = forward ? matchAtBoundary(src.back, other.front)
                           : matchAtBoundary(other.back, src.front);
          if (m && m->first <= maxDistE && m->second <= maxAngE) return true;
        }
        return false;
      };

      // Sources in quality order, so the best track claims its hits first.
      std::vector<std::size_t> order(segments.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return segments[a].fp.size() > segments[b].fp.size();
      });

      std::vector<SeqSegment> extensions;
      for (std::size_t si : order) {
        const SeqSegment src = segments[si];   // copy: `segments` may grow
        if (src.extension) continue;
        if (static_cast<int>(src.fp.size()) < minSrc) continue;
        if (src.chi2 / std::max(1, src.ndf) > maxSrcC2) continue;

        for (int step : {+1, -1}) {
          const bool forward = (step > 0);
          if (!forward && !m_extendBackward.value()) continue;
          const int nb = src.det + step;
          if (nb < 0 || nb > 2) continue;
          const auto& nbSurfacesAll = m_geoSvc->surfacesOf(nb);
          if (nbSurfacesAll.empty()) continue;
          if (hasSpliceMate(src, nb, forward)) continue;

          // The neighbour's free measurements.
          std::vector<Acts::SourceLink> slinks;
          for (std::size_t i = 0; i < measurements.size(); ++i) {
            if (measurements[i].detectorID != nb || owned.count(i)) continue;
            SeqSourceLink ssl;
            ssl.index = i;
            ssl.setGeometryId(measurements[i].surface->geometryId());
            slinks.push_back(Acts::SourceLink(ssl));
          }
          if (static_cast<int>(slinks.size()) < minExt) continue;
          std::sort(slinks.begin(), slinks.end(),
                    [](const Acts::SourceLink& a, const Acts::SourceLink& b) {
                      return a.get<SeqSourceLink>().geometryId() <
                             b.get<SeqSourceLink>().geometryId();
                    });

          // The surface list stays in beam order: DirectNavigator walks it
          // backwards itself when the propagation direction is Backward. The
          // start surface is the neighbour's edge nearest to the source.
          std::vector<const Acts::Surface*> nbSurfaces(nbSurfacesAll.begin(),
                                                       nbSurfacesAll.end());
          const Acts::Surface* sfStart = forward ? nbSurfaces.front()
                                                 : nbSurfaces.back();

          const SeqEndState& from = forward ? src.back : src.front;
          auto start = paramsAtSurface(from, sfStart, extCov);
          if (!start) continue;

          SeqSourceLinkAccessor slAccessor;
          slAccessor.slinks = &slinks;
          Acts::MeasurementSelectorCuts measCuts;
          measCuts.chi2CutOff            = {m_extendChi2CutOff.value()};
          measCuts.numMeasurementsCutOff = {std::size_t{1}};
          Acts::MeasurementSelector measSelector(measCuts);
          using SLinkIter = std::vector<Acts::SourceLink>::const_iterator;
          Acts::TrackStateCreator<SLinkIter, SeqTrackContainer> tsc;
          tsc.sourceLinkAccessor.connect<&SeqSourceLinkAccessor::operator()>(&slAccessor);
          tsc.calibrator.template connect<&SeqCalibrator::operator()>(&calibrator);
          tsc.measurementSelector
              .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
                  &measSelector);
          Acts::CombinatorialKalmanFilterExtensions<SeqTrackContainer> ext;
          ext.updater
              .connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
                  &gainMatrixUpdater);
          ext.createTrackStates
              .connect<&Acts::TrackStateCreator<SLinkIter, SeqTrackContainer>::createTrackStates>(
                  &tsc);
          SeqFixedNavigator nav;
          nav.surfaces = nbSurfaces;
          SeqStepper stepper(bField);
          SeqCKFPropagator prop(std::move(stepper), std::move(nav));
          SeqCKF extCkf(std::move(prop),
                        Acts::getDefaultLogger("SeqCKFExt",
                            static_cast<Acts::Logging::Level>(m_actsLogLevel.value())));
          Acts::PropagatorPlainOptions extPOpts = pOptions;
          extPOpts.direction = forward ? Acts::Direction::Forward()
                                       : Acts::Direction::Backward();
          Acts::CombinatorialKalmanFilterOptions<SeqTrackContainer> extOpts(
              gctx, m_mctx, std::cref(m_cctx), ext, extPOpts, true, true);

          auto res = extCkf.findTracks(*start, extOpts, tracks);
          if (!res.ok()) continue;
          std::optional<Acts::TrackIndexType> bestIdx;
          std::size_t bestN = 0; double bestC2 = std::numeric_limits<double>::max();
          for (const auto& t : *res) {
            const std::size_t n = fingerprintOf(t).size();
            if (static_cast<int>(n) < minExt) continue;
            const double c2 = chi2NdfOf(t);
            if (n > bestN || (n == bestN && c2 < bestC2)) {
              bestN = n; bestC2 = c2; bestIdx = t.index();
            }
          }
          if (!bestIdx) {
            debug() << "[ACTSSequentialTracker] evt=" << evtNum
                    << " extend det" << src.det << " -> det" << nb
                    << (forward ? " fwd" : " bwd") << " from ("
                    << from.pos.z() << "," << from.pos.y() << ") mm: nothing picked up"
                    << endmsg;
            continue;
          }
          auto extTrack = tracks.getTrack(*bestIdx);
          auto fp = fingerprintOf(extTrack);
          auto ends = endStatesOf(extTrack);
          if (!ends) continue;

          SeqSegment e;
          e.det       = nb;
          e.iteration = 2;
          e.idx       = *bestIdx;
          e.fp        = fp;
          e.chi2      = extTrack.chi2();
          e.ndf       = std::max(1, static_cast<int>(extTrack.nDoF()) - kHelixParams);
          e.front     = ends->first;
          e.back      = ends->second;
          e.seedX     = src.seedX;
          e.seedY     = src.seedY;
          e.extension = true;
          owned.insert(fp.begin(), fp.end());
          debug() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " extend det" << src.det << "(" << src.fp.size() << "m) -> det"
                  << nb << (forward ? " fwd" : " bwd") << ": picked up "
                  << fp.size() << " measurements, chi2/ndf=" << bestC2 << endmsg;
          if (m_dumpSegmentCellIDs.value()) {
            std::ostringstream cids;
            for (std::size_t mIdx : fp) cids << ' ' << measurements[mIdx].cellID;
            info() << "[SegDump] evt=" << evtNum << " det=" << nb
                   << " completeness=1 nMeas=" << fp.size()
                   << " cellIDs:" << cids.str() << " extension" << endmsg;
          }
          extensions.push_back(std::move(e));
          ++nExtended;
        }
      }
      for (auto& e : extensions) segments.push_back(std::move(e));
    }

    // =====================================================================
    // PHASE 3 — splice adjacent segments and refit the chains
    // =====================================================================
    // Pairs are judged by matchAtBoundary (propagation through the field),
    // defined above PHASE 2.4.
    const double maxDist  = m_spliceMaxDist.value();
    const double maxAngle = m_spliceMaxAngle.value();
    const double angleW   = m_spliceAngleWeight.value();

    std::vector<SeqLink> links;
    for (int detUp = 0; detUp < 2; ++detUp) {         // 0->1 and 1->2 only
      const int detDn = detUp + 1;
      for (std::size_t i = 0; i < segments.size(); ++i) {
        if (segments[i].det != detUp) continue;
        for (std::size_t j = 0; j < segments.size(); ++j) {
          if (segments[j].det != detDn) continue;
          auto m = matchAtBoundary(segments[i].back, segments[j].front);
          if (!m) continue;
          const double dist  = m->first;
          const double angle = m->second;
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
    // On top of those, two states carry the FITTED MOMENTUM, which the AtOther
    // states do not: location=AtFirstHit and location=AtLastHit, with
    //   omega     = q/p [1/GeV] of the smoothed parameters (NOT a curvature:
    //               same local convention as omega=stereo tilt above),
    //   phi, tanLambda, referencePoint = direction and position there,
    //   covMatrix(phi,phi), (omega,omega), (tanLambda,tanLambda) from the
    //               smoothed covariance.
    // q/p is what the KF measures in the MTC iron field, and SND-PFA's muon
    // momentum comes from it; without these two states it was fitted and then
    // dropped. Consumers that decode AtIP/AtOther (EDM4HEP2RNTuple, the event
    // display) filter on `location` and are unaffected.
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
        // Most upstream / most downstream smoothed state along the beam, for
        // the AtFirstHit/AtLastHit momentum states written after the loop.
        struct MomState {
          bool  valid = false;
          float beamZ = 0.f, ddX = 0.f, ddY = 0.f;
          float phi = 0.f, tanLambda = 0.f, qOverP = 0.f;
          float varPhi = 0.f, varQoP = 0.f, varTanL = 0.f;
        };
        MomState firstHit, lastHit;
        while (true) {
          auto ts = mutableTraj.getTrackState(tipIdx);
          if (ts.hasCalibrated()) {
            // Point the Track at the TrackerHit3D this state was fitted to.
            // Without it a written track carries only geometry, and nothing
            // downstream can say which particle built it: MC-truth studies had
            // to parse the [SegDump] log lines, and SND-PFA's muon corridor
            // needs the same join at run time. The hit objects live in the
            // input measurement collections, which job4 also writes out, so
            // the relation resolves when the file is read back.
            try {
              const std::size_t mIdx =
                  ts.getUncalibratedSourceLink().template get<SeqSourceLink>().index;
              if (mIdx < measurements.size()) {
                const auto& m = measurements[mIdx];
                const edm4hep::TrackerHit3DCollection* srcColl =
                    (m.detectorID == 0) ? stHits :
                    (m.detectorID == 1) ? spHits : mtcHits;
                if (srcColl && m.srcIndex >= 0 &&
                    static_cast<std::size_t>(m.srcIndex) < srcColl->size()) {
                  track.addToTrackerHits((*srcColl)[m.srcIndex]);
                }
              }
            } catch (const std::exception&) {
              // A state without a usable source link is not fatal here: the
              // geometry of the track is still written.
            }

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
            edm4ts.omega = static_cast<float>(stereoTiltZ);

            // Local -> DD4hep through the surface's own transform. ACTS global
            // is (beam, ddY, ddX), so ddX = g.z(), ddY = g.y(). This used to be
            // hand-written trig with ddX = +loc0*cos(alpha), the mirror image
            // of the real frame (loc0 = -ddX, see SeqCalibrator): every track
            // was drawn reflected in x, invisible while the fit only ever
            // found near-axis tracks and glaring once it found the GENIE
            // event-1 muon at x = -138 mm.
            auto toDD = [&](double l0, double l1) {
              const Acts::Vector3 g = surf.localToGlobal(
                  gctx, Acts::Vector2(l0, l1), Acts::Vector3::UnitX());
              return std::pair<float, float>{static_cast<float>(g.z()),
                                             static_cast<float>(g.y())};
            };

            if (ts.hasSmoothed() || ts.hasFiltered()) {
              const Acts::BoundVector bv =
                  ts.hasSmoothed() ? ts.smoothed() : ts.filtered();
              if (ts.hasSmoothed()) {
                const double theta = bv[Acts::eBoundTheta];
                const auto cov = ts.smoothedCovariance();
                // var(tanLambda) from var(theta): tanLambda = cot(theta), so
                // d(tanLambda)/d(theta) = -1/sin^2(theta).
                const double s2 = std::sin(theta) * std::sin(theta);
                const auto g = surf.localToGlobal(
                    gctx, Acts::Vector2(bv[Acts::eBoundLoc0], bv[Acts::eBoundLoc1]),
                    Acts::Vector3::UnitX());
                MomState m;
                m.valid  = true;
                m.beamZ  = beamZ;
                m.ddX    = static_cast<float>(g.z());
                m.ddY    = static_cast<float>(g.y());
                m.phi    = static_cast<float>(bv[Acts::eBoundPhi]);
                m.tanLambda =
                    static_cast<float>(std::tan(M_PI / 2.0 - theta));
                // q/p in 1/GeV (ACTS native energy unit is already GeV; the
                // factor documents the unit rather than changing the number).
                //
                // The SIGN is flipped on purpose. This tracking geometry maps
                // DD4hep (x, y, z) onto ACTS (z, y, x) — beam along ACTS x —
                // which is an odd permutation, i.e. a reflection. A magnetic
                // field is a pseudovector, so the physical By = -1.9 T of the
                // MTC iron corresponds to +1.9 T in this frame; the fit is
                // given -1.9 T and is self-consistent with the equally
                // reflected measurements (the calibrator's loc0 = -ddX), so
                // the trajectory and |q/p| come out right and nothing about
                // the fit needs changing — but the CHARGE it infers is the
                // mirror image of the real one. Measured on the 2000-event
                // GENIE sample before this flip: q/p was negative for 99 % of
                // true mu+ and for only 13 % of true mu-. Flipping it here
                // keeps the convention where a consumer expects it: q/p < 0
                // means a negative particle.
                m.qOverP = static_cast<float>(-bv[Acts::eBoundQOverP] *
                                              Acts::UnitConstants::GeV);
                m.varPhi = static_cast<float>(
                    cov(Acts::eBoundPhi, Acts::eBoundPhi));
                m.varQoP = static_cast<float>(
                    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) *
                    Acts::UnitConstants::GeV * Acts::UnitConstants::GeV);
                m.varTanL = (s2 > 1e-12)
                    ? static_cast<float>(
                          cov(Acts::eBoundTheta, Acts::eBoundTheta) / (s2 * s2))
                    : 0.f;
                if (!firstHit.valid || m.beamZ < firstHit.beamZ) firstHit = m;
                if (!lastHit.valid  || m.beamZ > lastHit.beamZ)  lastHit  = m;
              }
              edm4ts.phi = static_cast<float>(bv[Acts::eBoundPhi]);
              edm4ts.tanLambda = ts.hasSmoothed()
                  ? static_cast<float>(std::tan(M_PI / 2.0 - bv[Acts::eBoundTheta]))
                  : 0.0f;
              const double l0 = bv[Acts::eBoundLoc0];
              const double l1 = bv[Acts::eBoundLoc1];
              edm4ts.D0 = static_cast<float>(l0);
              const auto [ddX, ddY] = toDD(l0, l1);
              edm4ts.referencePoint = edm4hep::Vector3f{ddX, ddY, beamZ};
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

        auto addMomState = [&track](const MomState& m, int location) {
          if (!m.valid) return;
          edm4hep::TrackState ms{};
          ms.location  = location;
          ms.D0        = 0.0f;
          ms.Z0        = 0.0f;
          ms.phi       = m.phi;
          ms.tanLambda = m.tanLambda;
          ms.omega     = m.qOverP;   // q/p [1/GeV], see the note above
          ms.referencePoint = edm4hep::Vector3f{m.ddX, m.ddY, m.beamZ};
          ms.setCovMatrix(m.varPhi,  edm4hep::TrackParams::phi,
                          edm4hep::TrackParams::phi);
          ms.setCovMatrix(m.varQoP,  edm4hep::TrackParams::omega,
                          edm4hep::TrackParams::omega);
          ms.setCovMatrix(m.varTanL, edm4hep::TrackParams::tanLambda,
                          edm4hep::TrackParams::tanLambda);
          track.addToTrackStates(ms);
        };
        addMomState(firstHit, edm4hep::TrackState::AtFirstHit);
        addMomState(lastHit,  edm4hep::TrackState::AtLastHit);
      } catch (const std::exception& e) {
        warning() << "[ACTSSequentialTracker] evt=" << evtNum
                  << " trackStates iteration failed: " << e.what() << endmsg;
      }
    };

    for (const auto& g : globals)
      writeTrack(globalOut, g.idx, g.chi2, g.ndf, g.seedX, g.seedY, g.typeMask);

    // PromoteUnspliced: a segment is eligible only if EVERY detector it could
    // have spliced through found zero segments this event -- "nothing to
    // match against", not "a match was attempted and lost". SiTarget(0) and
    // MTC(2) each have one neighbour (SiPad); SiPad(1) has two, and needs
    // both empty (a segment that could have linked upstream but failed while
    // its downstream neighbour is genuinely empty is still a failed match,
    // not an absent one).
    std::array<std::size_t, 3> nSegByDet{0, 0, 0};
    for (const auto& s : segments) if (!s.extension) ++nSegByDet[s.det];
    static const std::array<std::vector<int>, 3> kNeighbours{
        std::vector<int>{1}, std::vector<int>{0, 2}, std::vector<int>{1}};
    auto hasNoSpliceCandidate = [&](int det) {
      for (int n : kNeighbours[det])
        if (nSegByDet[n] != 0) return false;
      return true;
    };

    std::array<int, 3> nLocal{0, 0, 0};
    int nPromoted = 0;
    for (const auto& s : segments) {
      if (s.spliced) continue;
      if (s.extension) continue;   // only meaningful as part of a chain
      if (m_promoteUnspliced.value() && hasNoSpliceCandidate(s.det)) {
        writeTrack(globalOut, s.idx, s.chi2, s.ndf, s.seedX, s.seedY,
                   1 << s.det);
        ++nPromoted;
        continue;
      }
      writeTrack(localOut[s.det], s.idx, s.chi2, s.ndf, s.seedX, s.seedY,
                 1 << s.det);
      ++nLocal[s.det];
    }

    info() << "[ACTSSequentialTracker] evt=" << evtNum
           << " measurements=" << measurements.size()
           << " segments=" << segments.size()
           << " merged=" << nMerged
           << " extended=" << nExtended
           << " proxdup=" << nProxDup
           << " links=" << takenLinks.size() << "/" << links.size()
           << " global=" << globals.size()
           << " promoted=" << nPromoted
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
