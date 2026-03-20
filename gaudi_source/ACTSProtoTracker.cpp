// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/DataHandle.h"

// edm4hep input
#include "edm4hep/TrackerHit3DCollection.h"

// edm4hep output
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MutableTrack.h"

// ACTS
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"

// SND geometry service
#include "ISNDGeoSvc.h"

// Standard
#include <algorithm>
#include <atomic>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDMeasurement: holds all information needed from a single measurement
// for proto-tracking and eventually for the CKF.
// ---------------------------------------------------------------------------

struct SNDMeasurement {
  const Acts::Surface* surface     = nullptr; // ACTS surface for this measurement
  double               localCoord  = 0.0;     // primary local coordinate [mm]
  double               localCoord2 = 0.0;     // secondary local coord [mm] (2D only)
  double               variance    = 0.0;     // variance of localCoord [mm^2]
  double               variance2   = 0.0;     // variance of localCoord2 (2D only)
  bool                 is2D        = false;   // true for SiPad pixel hits
  int                  detectorID  = -1;      // 0=SiTarget, 1=SiPad
  int                  plane       = -1;      // 0=StripY, 1=StripZ, -1=pixel
  float                time        = 0.0f;    // hit time [ns]
  float                eDep        = 0.0f;    // energy deposit [GeV]

  SNDMeasurement(const Acts::Surface* sf,
                 double lc, double lc2,
                 double var, double var2,
                 bool twod,
                 int detID, int pl,
                 float t, float e)
      : surface(sf), localCoord(lc), localCoord2(lc2),
        variance(var), variance2(var2), is2D(twod),
        detectorID(detID), plane(pl), time(t), eDep(e) {}
};

// ---------------------------------------------------------------------------
// ACTSProtoTracker
// ---------------------------------------------------------------------------

class ACTSProtoTracker : public Gaudi::Algorithm {
public:
  ACTSProtoTracker(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  // ---- Properties ----
  Gaudi::Property<std::string> m_inputSiTarget{
      this, "InputSiTarget", "SiTargetMeasurements",
      "SiTarget TrackerHit3DCollection from SiTargetMeasConverter"};

  Gaudi::Property<std::string> m_inputSiPad{
      this, "InputSiPad", "SiPadMeasurements",
      "SiPad TrackerHit3DCollection from SiPadMeasConverter"};

  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "ACTSTracks",
      "Output edm4hep::TrackCollection name"};

  Gaudi::Property<double> m_bFieldX{
      this, "BFieldX", 0.0, "Magnetic field X [T]"};
  Gaudi::Property<double> m_bFieldY{
      this, "BFieldY", 0.0, "Magnetic field Y [T]"};
  Gaudi::Property<double> m_bFieldZ{
      this, "BFieldZ", 0.0, "Magnetic field Z [T]"};

  // ---- DataHandles ----
  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siTargetHandle;

  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siPadHandle;

  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::TrackCollection>>
      m_outputHandle;

  // ---- Services ----
  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  // ---- Event counter ----
  mutable std::atomic<long long> m_eventCount{0};

  // ---- Cached surface offset for SiPad ----
  // SiPad surfaces start at index 20 in allSurfaces() after 20 SiTarget surfaces.
  static constexpr std::size_t kSiPadOffset = 20;
};

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    // Build DataHandles
    m_siTargetHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiTarget.value(), Gaudi::DataHandle::Reader, this);

    m_siPadHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiPad.value(), Gaudi::DataHandle::Reader, this);

    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    // Retrieve geometry service
    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[ACTSProtoTracker] Failed to retrieve ACTSGeoSvc." << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "[ACTSProtoTracker] Initialized. GeoSvc has "
           << m_geoSvc->allSurfaces().size() << " surfaces." << endmsg;

    return sc;
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in initialize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    // Always create output collection
    auto* output = m_outputHandle->createAndPut();

    // Get pre-sorted surface list from geometry service
    const auto& allSurfaces = m_geoSvc->allSurfaces();
    const auto& gctx        = m_geoSvc->geometryContext();

    // =========================================================================
    // STEP 1: Collect measurements from both detectors into SNDMeasurement list
    // =========================================================================
    std::vector<SNDMeasurement> measurements;

    // ---- SiTarget measurements (1D strips) ----------------------------------
    const auto* stHits = m_siTargetHandle->get();
    if (stHits) {
      for (std::size_t i = 0; i < stHits->size(); ++i) {
        const auto& hit   = (*stHits)[i];
        const int   plane = hit.getQuality();  // 0=StripY, 1=StripZ

        // Surface index: layer * 2 + plane
        // Layer is not stored directly in TrackerHit3D — recover from
        // surface position: find the surface whose X center is closest
        // to hit position X within the SiTarget range (first 20 surfaces).
        const auto& pos = hit.getPosition();  // global position [mm]

        // Find closest SiTarget surface by X position
        // (more robust than computing layer index from X)
        const Acts::Surface* bestSurface = nullptr;
        double bestDist = 1e9;
        for (std::size_t j = 0; j < kSiPadOffset && j < allSurfaces.size(); ++j) {
          double dx = std::abs(allSurfaces[j]->center(gctx).x() - pos.x);
          if (dx < bestDist) {
            bestDist    = dx;
            bestSurface = allSurfaces[j];
          }
        }

        if (!bestSurface) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " SiTarget hit " << i << " has no matching surface."
                    << endmsg;
          continue;
        }

        // Local coordinate from position:
        //   plane=0 (StripY): localCoord = position.y
        //   plane=1 (StripZ): localCoord = position.z
        double localCoord = (plane == 0) ? pos.y : pos.z;

        // Variance from covMatrix:
        //   plane=0: cov[3] (yy element)
        //   plane=1: cov[5] (zz element)
        const auto& cov = hit.getCovMatrix();
        double var = (plane == 0) ? cov[3] : cov[5];

        measurements.emplace_back(
            bestSurface,
            localCoord, 0.0,   // localCoord, localCoord2
            var, 0.0,          // variance, variance2
            false,             // is2D
            0, plane,          // detectorID=0 (SiTarget), plane
            hit.getTime(), hit.getEDep());
      }
    }

    // ---- SiPad measurements (2D pixels) -------------------------------------
    const auto* spHits = m_siPadHandle->get();
    if (spHits) {
      for (std::size_t i = 0; i < spHits->size(); ++i) {
        const auto& hit = (*spHits)[i];
        const auto& pos = hit.getPosition();

        // Find closest SiPad surface by X position
        // (SiPad surfaces are at indices kSiPadOffset..end)
        const Acts::Surface* bestSurface = nullptr;
        double bestDist = 1e9;
        for (std::size_t j = kSiPadOffset; j < allSurfaces.size(); ++j) {
          double dx = std::abs(allSurfaces[j]->center(gctx).x() - pos.x);
          if (dx < bestDist) {
            bestDist    = dx;
            bestSurface = allSurfaces[j];
          }
        }

        if (!bestSurface) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " SiPad hit " << i << " has no matching surface."
                    << endmsg;
          continue;
        }

        // SiPad: 2D measurement, Y and Z
        const auto& cov = hit.getCovMatrix();
        double varY = cov[3];  // yy element
        double varZ = cov[5];  // zz element

        measurements.emplace_back(
            bestSurface,
            pos.y, pos.z,      // localCoord=Y, localCoord2=Z
            varY, varZ,        // variance in Y and Z
            true,              // is2D
            1, -1,             // detectorID=1 (SiPad), plane=-1 (pixel)
            hit.getTime(), hit.getEDep());
      }
    }

    if (measurements.empty()) {
      debug() << "[ACTSProtoTracker] evt=" << evtNum
              << " no measurements." << endmsg;
      return StatusCode::SUCCESS;
    }

    // Sort measurements by surface X position (beam axis order)
    std::sort(measurements.begin(), measurements.end(),
              [&](const SNDMeasurement& a, const SNDMeasurement& b) {
                return a.surface->center(gctx).x() <
                       b.surface->center(gctx).x();
              });

    debug() << "[ACTSProtoTracker] evt=" << evtNum
            << " SiTarget=" << (stHits ? stHits->size() : 0)
            << " SiPad=" << (spHits ? spHits->size() : 0)
            << " total measurements=" << measurements.size() << endmsg;

    // =========================================================================
    // STEP 2: CKF PLACEHOLDER
    // =========================================================================
    // TODO (Step 5): Replace the proto-track below with a real CKF.
    //
    // The CKF will need:
    //   - m_geoSvc->trackingGeometry()  : Acts::TrackingGeometry
    //   - m_geoSvc->geometryContext()   : Acts::GeometryContext
    //   - Acts::MagneticFieldContext    : from BField properties
    //   - Acts::CalibrationContext      : default empty context
    //   - measurements vector           : built above, contains source links
    //                                     and measurement covariances
    //   - seed parameters               : estimated from first/last measurement
    //
    // Pseudocode for Step 5:
    //
    //   Acts::ConstantBField bField(Acts::Vector3(BFieldX, BFieldY, BFieldZ));
    //   Stepper stepper(bField);
    //   Navigator navigator({&m_geoSvc->trackingGeometry()});
    //   Propagator propagator(stepper, navigator);
    //   CKF ckf(propagator, logger);
    //
    //   auto seed = estimateSeed(measurements, gctx);
    //   auto result = ckf.findTracks(seed, sourceLinkAccessor, ckfOptions, tracks);
    //
    //   if (result.ok()) { writeTracks(result.value(), output); }
    //
    // END CKF PLACEHOLDER
    // =========================================================================

    // =========================================================================
    // STEP 3: Proto-track — straight line estimate from all measurements
    // =========================================================================

    // Estimate direction from first to last measurement surface center
    Acts::Vector3 dir(1.0, 0.0, 0.0);  // default: along beam axis
    if (measurements.size() >= 2) {
      const Acts::Vector3& fc = measurements.front().surface->center(gctx);
      const Acts::Vector3& lc = measurements.back().surface->center(gctx);
      Acts::Vector3 d = lc - fc;
      if (d.norm() > 1e-6) dir = d.normalized();
    }

    // Write one proto-track per event window
    auto track = output->create();
    track.setType(1);   // muon hypothesis
    track.setChi2(0.0);
    track.setNdf(static_cast<int>(measurements.size()) - 4);

    const float phi = static_cast<float>(std::atan2(dir.y(), dir.x()));
    for (const auto& meas : measurements) {
      edm4hep::TrackState ts;
      ts.location  = edm4hep::TrackState::AtOther;
      ts.D0        = 0.0f;
      ts.Z0        = 0.0f;
      ts.phi       = phi;
      ts.tanLambda = 0.0f;
      ts.omega     = 0.0f;
      track.addToTrackStates(ts);
    }

    if (evtNum % 100 == 0 || evtNum < 5) {
      info() << "[ACTSProtoTracker] evt=" << evtNum
             << " measurements=" << measurements.size()
             << " ndf=" << track.getNdf()
             << " phi=" << phi
             << " tracks written=1 (proto-track)" << endmsg;
    }

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in execute()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::finalize() {
  try {
    m_siTargetHandle.reset();
    m_siPadHandle.reset();
    m_outputHandle.reset();
    info() << "[ACTSProtoTracker] Done. Total events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in finalize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(ACTSProtoTracker)
