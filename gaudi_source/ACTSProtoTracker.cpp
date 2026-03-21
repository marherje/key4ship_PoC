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
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"

// ACTS propagator
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"

// ACTS track finding
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"

// ACTS track fitting
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

// ACTS track containers
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"

// ACTS magnetic field
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

// ACTS calibration context
#include "Acts/Utilities/CalibrationContext.hpp"

// ACTS measurement helpers
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"

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
// SNDMeasurement
// ---------------------------------------------------------------------------

struct SNDMeasurement {
  const Acts::Surface* surface     = nullptr;
  double               localCoord  = 0.0;
  double               localCoord2 = 0.0;
  double               variance    = 0.0;
  double               variance2   = 0.0;
  bool                 is2D        = false;
  int                  detectorID  = -1;
  int                  plane       = -1;
  float                time        = 0.0f;
  float                eDep        = 0.0f;

  SNDMeasurement(const Acts::Surface* sf,
                 double lc, double lc2,
                 double var, double var2,
                 bool twod, int detID, int pl,
                 float t, float e)
      : surface(sf), localCoord(lc), localCoord2(lc2),
        variance(var), variance2(var2), is2D(twod),
        detectorID(detID), plane(pl), time(t), eDep(e) {}
};

// ---------------------------------------------------------------------------
// SNDSourceLink
// ---------------------------------------------------------------------------

struct SNDSourceLink {
  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
  std::size_t index = 0;
  void setGeometryId(Acts::GeometryIdentifier gid) { m_geometryId = gid; }
private:
  Acts::GeometryIdentifier m_geometryId;
};

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

using SNDStepper    = Acts::EigenStepper<>;
using SNDNavigator  = Acts::Navigator;
using SNDPropagator = Acts::Propagator<SNDStepper, SNDNavigator>;

using SNDTrackContainer = Acts::TrackContainer<
    Acts::VectorTrackContainer,
    Acts::VectorMultiTrajectory,
    std::shared_ptr>;

using SNDCKF = Acts::CombinatorialKalmanFilter<SNDPropagator, SNDTrackContainer>;

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
  Gaudi::Property<std::string> m_inputSiTarget{
      this, "InputSiTarget", "SiTargetMeasurements",
      "SiTarget TrackerHit3DCollection from SiTargetMeasConverter"};
  Gaudi::Property<std::string> m_inputSiPad{
      this, "InputSiPad", "SiPadMeasurements",
      "SiPad TrackerHit3DCollection from SiPadMeasConverter"};
  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "ACTSTracks",
      "Output edm4hep::TrackCollection name"};
  Gaudi::Property<double> m_bFieldX{this, "BFieldX", 0.0, "BField X [T]"};
  Gaudi::Property<double> m_bFieldY{this, "BFieldY", 0.0, "BField Y [T]"};
  Gaudi::Property<double> m_bFieldZ{this, "BFieldZ", 0.0, "BField Z [T]"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siTargetHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siPadHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackCollection>>
      m_outputHandle;

  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  mutable std::atomic<long long> m_eventCount{0};
  mutable Acts::MagneticFieldContext m_mctx;
  mutable Acts::CalibrationContext   m_cctx;

  // Dynamically computed in initialize(): number of SiTarget surfaces.
  // SiPad surfaces start at this index in allSurfaces().
  mutable std::size_t m_nSiTargetSurfaces{0};

  // GainMatrixUpdater stored as member to ensure valid lifetime
  // when the CKF lambda delegate captures a pointer to it.
  mutable Acts::GainMatrixUpdater m_gainMatrixUpdater;
};

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    m_siTargetHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiTarget.value(), Gaudi::DataHandle::Reader, this);
    m_siPadHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiPad.value(), Gaudi::DataHandle::Reader, this);
    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[ACTSProtoTracker] Failed to retrieve ACTSGeoSvc." << endmsg;
      return StatusCode::FAILURE;
    }

    // Compute the number of SiTarget surfaces dynamically.
    // Find the largest gap in X between consecutive surfaces —
    // that gap separates SiTarget from SiPad.
    const auto& allSurf = m_geoSvc->allSurfaces();
    const auto& gctx    = m_geoSvc->geometryContext();

    if (allSurf.size() < 2) {
      error() << "[ACTSProtoTracker] Too few surfaces in geometry." << endmsg;
      return StatusCode::FAILURE;
    }

    double maxGap      = 0.0;
    std::size_t gapIdx = 0;
    for (std::size_t i = 1; i < allSurf.size(); ++i) {
      double gap = allSurf[i]->center(gctx).x() -
                   allSurf[i-1]->center(gctx).x();
      if (gap > maxGap) {
        maxGap  = gap;
        gapIdx  = i;
      }
    }
    m_nSiTargetSurfaces = gapIdx;

    info() << "[ACTSProtoTracker] Initialized. GeoSvc has "
           << allSurf.size() << " surfaces. "
           << "SiTarget surfaces: " << m_nSiTargetSurfaces
           << " SiPad surfaces: " << (allSurf.size() - m_nSiTargetSurfaces)
           << " (gap=" << maxGap << " mm)" << endmsg;

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

    auto* output = m_outputHandle->createAndPut();

    const auto& allSurfaces = m_geoSvc->allSurfaces();
    const auto& gctx        = m_geoSvc->geometryContext();

    // =========================================================================
    // STEP 1: Collect measurements
    // =========================================================================
    std::vector<SNDMeasurement> measurements;

    // ---- SiTarget (1D strips) -----------------------------------------------
    const auto* stHits = m_siTargetHandle->get();
    if (stHits) {
      for (std::size_t i = 0; i < stHits->size(); ++i) {
        const auto& hit   = (*stHits)[i];
        const int   plane = hit.getQuality();
        const auto& pos   = hit.getPosition();

        const Acts::Surface* bestSurface = nullptr;
        double bestDist = 1e9;
        for (std::size_t j = 0; j < m_nSiTargetSurfaces; ++j) {
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

        const auto& cov   = hit.getCovMatrix();
        double localCoord = (plane == 0) ? pos.y : pos.z;
        double var        = (plane == 0) ? cov[3] : cov[5];

        measurements.emplace_back(bestSurface, localCoord, 0.0,
                                   var, 0.0, false, 0, plane,
                                   hit.getTime(), hit.getEDep());
      }
    }

    // ---- SiPad (2D pixels) --------------------------------------------------
    const auto* spHits = m_siPadHandle->get();
    if (spHits) {
      for (std::size_t i = 0; i < spHits->size(); ++i) {
        const auto& hit = (*spHits)[i];
        const auto& pos = hit.getPosition();

        const Acts::Surface* bestSurface = nullptr;
        double bestDist = 1e9;
        for (std::size_t j = m_nSiTargetSurfaces; j < allSurfaces.size(); ++j) {
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

        const auto& cov = hit.getCovMatrix();
        measurements.emplace_back(bestSurface, pos.y, pos.z,
                                   cov[3], cov[5], true, 1, -1,
                                   hit.getTime(), hit.getEDep());
      }
    }

    if (measurements.empty()) {
      debug() << "[ACTSProtoTracker] evt=" << evtNum
              << " no measurements." << endmsg;
      return StatusCode::SUCCESS;
    }

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
    // STEP 2: Build source link map
    // =========================================================================
    std::unordered_map<Acts::GeometryIdentifier,
                       std::vector<Acts::SourceLink>> sourceLinkMap;

    for (std::size_t i = 0; i < measurements.size(); ++i) {
      SNDSourceLink ssl;
      ssl.index = i;
      ssl.setGeometryId(measurements[i].surface->geometryId());
      sourceLinkMap[measurements[i].surface->geometryId()].push_back(
          Acts::SourceLink(ssl));
    }

    if (sourceLinkMap.empty()) {
      debug() << "[ACTSProtoTracker] evt=" << evtNum
              << " no source links." << endmsg;
      return StatusCode::SUCCESS;
    }

    // =========================================================================
    // STEP 3: Build CKF components
    // =========================================================================

    auto bField = std::make_shared<Acts::ConstantBField>(
        Acts::Vector3(m_bFieldX.value() * Acts::UnitConstants::T,
                      m_bFieldY.value() * Acts::UnitConstants::T,
                      m_bFieldZ.value() * Acts::UnitConstants::T));

    SNDStepper stepper(bField);

    auto tgPtr = std::shared_ptr<const Acts::TrackingGeometry>(
        &m_geoSvc->trackingGeometry(),
        [](const Acts::TrackingGeometry*) {});

    Acts::Navigator::Config navCfg{
        std::const_pointer_cast<Acts::TrackingGeometry>(tgPtr)};
    navCfg.resolvePassive   = false;
    navCfg.resolveMaterial  = false;
    navCfg.resolveSensitive = true;
    SNDNavigator navigator(navCfg);
    SNDPropagator propagator(std::move(stepper), std::move(navigator));

    // ---- Source link accessor -----------------------------------------------
    struct SourceLinkAccessor {
      const std::unordered_map<Acts::GeometryIdentifier,
                               std::vector<Acts::SourceLink>>* slMap = nullptr;
      using Iterator = std::vector<Acts::SourceLink>::const_iterator;  // ← AÑADIR ESTA LÍNEA
      std::pair<Iterator, Iterator> operator()(
        const Acts::Surface& surface) const {
        // TEMPORARY DIAGNOSTIC (no atomic, simpler)
        static int callCount = 0;
        if (callCount++ < 5) {
          std::size_t nLinks = slMap->count(surface.geometryId()) ?
          slMap->at(surface.geometryId()).size() : 0;
          std::cerr << "[SLAccessor] geoID=" << surface.geometryId().value()
              << " nLinks=" << nLinks << "\n";
        }
        // END DIAGNOSTIC
        auto it = slMap->find(surface.geometryId());
        if (it == slMap->end()) {
          static const std::vector<Acts::SourceLink> empty{};
          return {empty.begin(), empty.end()};
        }
        return {it->second.begin(), it->second.end()};
      } 
    };

    // ---- Calibrator ---------------------------------------------------------
    struct SNDCalibrator {
      const std::vector<SNDMeasurement>* meas = nullptr;

      void operator()(const Acts::GeometryContext& /*gctx*/,
                      const Acts::CalibrationContext& /*cctx*/,
                      const Acts::SourceLink& sl,
                      Acts::VectorMultiTrajectory::TrackStateProxy ts) const {
        const auto& ssl = sl.get<SNDSourceLink>();
        const auto& m   = (*meas)[ssl.index];

        if (m.is2D) {
          constexpr std::array<Acts::BoundIndices, 2> indices = {
              Acts::eBoundLoc0, Acts::eBoundLoc1};
          ts.allocateCalibrated(2);
          ts.template calibrated<2>() =
              Acts::ActsVector<2>(m.localCoord, m.localCoord2);
          ts.template calibratedCovariance<2>() =
              Acts::ActsSquareMatrix<2>{{m.variance, 0.0},
                                        {0.0, m.variance2}};
          ts.setProjectorSubspaceIndices(indices);
        } else {
          constexpr std::array<Acts::BoundIndices, 1> indices = {
              Acts::eBoundLoc0};
          ts.allocateCalibrated(1);
          ts.template calibrated<1>() =
              Acts::ActsVector<1>(m.localCoord);
          ts.template calibratedCovariance<1>() =
              Acts::ActsSquareMatrix<1>{{m.variance}};
          ts.setProjectorSubspaceIndices(indices);
        }
      }
    };

    // ---- Measurement selector -----------------------------------------------
    Acts::MeasurementSelectorCuts selCuts;
    selCuts.chi2CutOff.clear();
    selCuts.chi2CutOff.push_back(15.0);
    selCuts.numMeasurementsCutOff.clear();
    selCuts.numMeasurementsCutOff.push_back(std::size_t{10});
    Acts::MeasurementSelector measSelector(selCuts);

    // ---- TrackStateCreator --------------------------------------------------
    using SLIter    = std::vector<Acts::SourceLink>::const_iterator;
    using TSCreator = Acts::TrackStateCreator<SLIter, SNDTrackContainer>;

    TSCreator tsCreator;

    SourceLinkAccessor slAccessor;
    slAccessor.slMap = &sourceLinkMap;
    tsCreator.sourceLinkAccessor
        .template connect<&SourceLinkAccessor::operator()>(&slAccessor);

    SNDCalibrator calibrator;
    calibrator.meas = &measurements;
    tsCreator.calibrator
        .template connect<&SNDCalibrator::operator()>(&calibrator);

    tsCreator.measurementSelector
        .template connect<&Acts::MeasurementSelector::select<
            Acts::VectorMultiTrajectory>>(&measSelector);

    // ---- CKF extensions -----------------------------------------------------
    Acts::CombinatorialKalmanFilterExtensions<SNDTrackContainer> extensions;

    extensions.updater
        .connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
            &m_gainMatrixUpdater);

    extensions.createTrackStates
        .template connect<&TSCreator::createTrackStates>(&tsCreator);

    // =========================================================================
    // STEP 4: Seed parameters
    // =========================================================================
    const auto& firstM  = measurements.front();
    const auto& sfFirst = *firstM.surface;

    Acts::Vector4 seedPos4;
    seedPos4[Acts::ePos0] = sfFirst.center(gctx).x();
    seedPos4[Acts::ePos1] = firstM.localCoord;
    seedPos4[Acts::ePos2] = 0.0;
    seedPos4[Acts::eTime] = static_cast<double>(firstM.time);

    Acts::Vector3 seedDir(1.0, 0.0, 0.0);
    if (measurements.size() >= 2) {
      const auto& lastM = measurements.back();
      Acts::Vector3 d(
          lastM.surface->center(gctx).x() - sfFirst.center(gctx).x(),
          lastM.localCoord - firstM.localCoord,
          0.0);
      if (d.norm() > 1e-6) seedDir = d.normalized();
    }

    const double seedMomentum = 10.0 * Acts::UnitConstants::GeV;
    const double seedQoverP   = 1.0 / seedMomentum;

    Acts::BoundSquareMatrix seedCov = Acts::BoundSquareMatrix::Zero();
    seedCov(Acts::eBoundLoc0,   Acts::eBoundLoc0)   = 100.0;
    seedCov(Acts::eBoundLoc1,   Acts::eBoundLoc1)   = 100.0;
    seedCov(Acts::eBoundPhi,    Acts::eBoundPhi)     = 0.1;
    seedCov(Acts::eBoundTheta,  Acts::eBoundTheta)   = 0.1;
    seedCov(Acts::eBoundQOverP, Acts::eBoundQOverP)  = 1.0;
    seedCov(Acts::eBoundTime,   Acts::eBoundTime)    = 1e6;

    auto seedParamsResult = Acts::BoundTrackParameters::create(
        gctx, sfFirst.getSharedPtr(), seedPos4, seedDir, seedQoverP,
        seedCov, Acts::ParticleHypothesis::muon());

    if (!seedParamsResult.ok()) {
      warning() << "[ACTSProtoTracker] evt=" << evtNum
                << " seed creation failed: "
                << seedParamsResult.error() << endmsg;
      return StatusCode::SUCCESS;
    }
    const auto& seedParams = *seedParamsResult;

    // =========================================================================
    // STEP 5: Run CKF
    // =========================================================================
    Acts::PropagatorPlainOptions pOptions(gctx, m_mctx);
    pOptions.direction = Acts::Direction::Forward();
    pOptions.maxSteps  = 1000;

    Acts::CombinatorialKalmanFilterOptions<SNDTrackContainer> ckfOptions(
        gctx, m_mctx, std::cref(m_cctx), extensions, pOptions,
        false, false);

    auto trackBackend = std::make_shared<Acts::VectorTrackContainer>();
    auto trajBackend  = std::make_shared<Acts::VectorMultiTrajectory>();
    SNDTrackContainer ckfTracks(trackBackend, trajBackend);

    SNDCKF ckf(std::move(propagator),
               Acts::getDefaultLogger("CKF", Acts::Logging::WARNING));

    auto ckfResult = ckf.findTracks(seedParams, ckfOptions, ckfTracks);

    debug() << "[ACTSProtoTracker] evt=" << evtNum
            << " CKF ok=" << ckfResult.ok()
            << " nTracksInContainer=" << ckfTracks.size() << endmsg;
    if (ckfResult.ok()) {
      for (const auto& t : ckfTracks) {
        debug() << "  track nMeas=" << t.nMeasurements()
                << " nHoles=" << t.nHoles()
                << " chi2=" << t.chi2() << endmsg;
      }
    }

    // =========================================================================
    // STEP 6: Write output
    // =========================================================================
    if (!ckfResult.ok()) {
      warning() << "[ACTSProtoTracker] evt=" << evtNum
                << " CKF failed: " << ckfResult.error() << endmsg;
      return StatusCode::SUCCESS;
    }

    std::size_t nTracks = 0;
    for (const auto& ckfTrack : ckfTracks) {
      if (ckfTrack.nMeasurements() < 3) continue;

      auto track = output->create();
      track.setType(1);
      track.setChi2(static_cast<float>(ckfTrack.chi2()));
      track.setNdf(static_cast<int>(ckfTrack.nDoF()));

      for (const auto& ts : ckfTrack.trackStates()) {
        if (!ts.hasCalibrated()) continue;
        edm4hep::TrackState edm4ts;
        edm4ts.location  = edm4hep::TrackState::AtOther;
        edm4ts.D0        = 0.0f;
        edm4ts.Z0        = 0.0f;
        if (ts.hasSmoothed()) {
          edm4ts.phi       = static_cast<float>(
              ts.smoothed()[Acts::eBoundPhi]);
          edm4ts.tanLambda = static_cast<float>(
              std::tan(M_PI / 2.0 - ts.smoothed()[Acts::eBoundTheta]));
          edm4ts.omega     = static_cast<float>(
              ts.smoothed()[Acts::eBoundQOverP]);
        } else if (ts.hasFiltered()) {
          edm4ts.phi       = static_cast<float>(
              ts.filtered()[Acts::eBoundPhi]);
          edm4ts.tanLambda = 0.0f;
          edm4ts.omega     = 0.0f;
        } else {
          edm4ts.phi       = 0.0f;
          edm4ts.tanLambda = 0.0f;
          edm4ts.omega     = 0.0f;
        }
        track.addToTrackStates(edm4ts);
      }
      ++nTracks;
    }

    if (evtNum < 5 || evtNum % 100 == 0) {
      info() << "[ACTSProtoTracker] evt=" << evtNum
             << " measurements=" << measurements.size()
             << " CKF tracks=" << nTracks << endmsg;
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