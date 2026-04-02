// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/DataHandle.h"

// edm4hep input
#include "edm4hep/SimCalorimeterHitCollection.h"

// edm4hep output
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/MutableTrackerHit3D.h"
#include "edm4hep/CovMatrix3f.h"
#include "edm4hep/Vector3d.h"

// DD4hep segmentation
#include "DDSegmentation/BitFieldCoder.h"

// ACTS
#include "Acts/Definitions/Units.hpp"

// SND geometry service
#include "ISNDGeoSvc.h"

// Standard
#include <atomic>
#include <cmath>
#include <memory>
#include <string>

// ---------------------------------------------------------------------------
// SiTargetMeasConverter
// ---------------------------------------------------------------------------

class SiTargetMeasConverter : public Gaudi::Algorithm {
public:
  SiTargetMeasConverter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  // ---- Properties ----
  Gaudi::Property<std::string> m_inputCollection{
      this, "InputCollection", "SiTargetHitsWindowed",
      "Input SimCalorimeterHitCollection name"};

  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "SiTargetMeasurements",
      "Output TrackerHit3DCollection name"};

  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,layer:8,slice:4,plane:1,strip:14",
      "DD4hep BitField string for SiTarget cellID decoding"};

  Gaudi::Property<double> m_stripPitch{
      this, "StripPitch", 0.0755,
      "Strip pitch in mm — used to compute measurement variance"};

  // ---- DataHandles ----
  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>
      m_inputHandle;

  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_outputHandle;

  // ---- Services ----
  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  // ---- Internal state ----
  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::atomic<long long> m_eventCount{0};
};

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode SiTargetMeasConverter::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    // Build DataHandles
    m_inputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inputCollection.value(), Gaudi::DataHandle::Reader, this);

    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    // Build BitField decoder
    m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
        m_bitField.value());

    // Retrieve geometry service
    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[SiTargetMeasConverter] Failed to retrieve ACTSGeoSvc."
              << endmsg;
      return StatusCode::FAILURE;
    }

    const std::size_t nSurfaces = m_geoSvc->allSurfaces().size();
    info() << "[SiTargetMeasConverter] Initialized. GeoSvc has "
           << nSurfaces << " surfaces." << endmsg;

    return sc;
  } catch (const std::exception& e) {
    error() << "[SiTargetMeasConverter] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[SiTargetMeasConverter] Unknown exception in initialize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode SiTargetMeasConverter::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    // Always create output collection
    auto* output = m_outputHandle->createAndPut();

    // Read input hits
    const auto* hits = m_inputHandle->get();
    if (!hits || hits->size() == 0) {
      debug() << "[SiTargetMeasConverter] evt=" << evtNum
              << " no SiTarget hits." << endmsg;
      return StatusCode::SUCCESS;
    }

    // Get ordered surfaces from geometry service.
    // SiTarget surfaces are the first 20 in the sorted list (more negative X).
    // Index: layer * 2 + plane  (plane=0 StripY, plane=1 StripZ)
    const auto& allSurfaces = m_geoSvc->allSurfaces();

    // Precompute variance from strip pitch
    // variance = (pitch / sqrt(12))^2 = pitch^2 / 12
    const double pitch    = m_stripPitch.value();  // mm
    const double variance = (pitch * pitch) / 12.0;  // mm^2

    std::size_t nConverted = 0;
    std::size_t nSkipped   = 0;

    for (std::size_t i = 0; i < hits->size(); ++i) {
      const auto&    hit = (*hits)[i];
      const uint64_t cid = hit.getCellID();

      // Decode layer and plane from cellID
      int layer = static_cast<int>(m_decoder->get(cid, "layer"));
      int plane = static_cast<int>(m_decoder->get(cid, "plane"));

      // SiTarget: 2 surfaces per station, ordered as StripY (plane=0) then
      // StripZ (plane=1) within each station.
      // Global surface index = layer * 2 + plane
      const std::size_t surfIdx =
          static_cast<std::size_t>(layer * 2 + plane);

      if (surfIdx >= allSurfaces.size()) {
        warning() << "[SiTargetMeasConverter] evt=" << evtNum
                  << " hit=" << i << " layer=" << layer
                  << " plane=" << plane
                  << " surfIdx=" << surfIdx << " out of range, skipping."
                  << endmsg;
        ++nSkipped;
        continue;
      }

      // Get hit position in global coordinates [mm]
      // edm4hep stores positions in mm
      const auto& pos = hit.getPosition();

      // Local measurement coordinate (Z=beam convention):
      //   plane=0 (StripX): strips parallel to Y, measures X
      //   plane=1 (StripY): strips parallel to X, measures Y
      // Stored in covMatrix as variance; localCoord implicit in position.

      // Get hit time: use time of first contribution if available
      float hitTime = 0.0f;
      if (hit.contributions_size() > 0) {
        hitTime = static_cast<float>(hit.getContributions(0).getTime());
      }

      // Build CovMatrix3f: only the diagonal element for the measured
      // coordinate carries the strip variance.
      // The matrix is symmetric 3x3 stored as upper triangle:
      // [xx, xy, xz, yy, yz, zz] — 6 elements.
      //   plane=0 → xx = variance (position.x is the measurement)
      //   plane=1 → yy = variance (position.y is the measurement)
      edm4hep::CovMatrix3f cov{};  // zero-initialized
      // Z=beam: plane=0 (StripX) measures X → xx=cov[0]
      //         plane=1 (StripY) measures Y → yy=cov[3]
      if (plane == 0) {
        cov[0] = static_cast<float>(variance);  // xx element
      } else {
        cov[3] = static_cast<float>(variance);  // yy element
      }

      // Create TrackerHit3D
      auto mhit = output->create();
      mhit.setCellID(cid);
      mhit.setType(0);          // detector ID: 0 = SiTarget
      mhit.setQuality(plane);   // plane: 0=StripY, 1=StripZ
      mhit.setTime(hitTime);
      mhit.setEDep(hit.getEnergy());
      mhit.setEDepError(0.0f);
      mhit.setPosition(edm4hep::Vector3d{pos.x, pos.y, pos.z});
      mhit.setCovMatrix(cov);

      ++nConverted;
    }

    debug() << "[SiTargetMeasConverter] evt=" << evtNum
            << " input=" << hits->size()
            << " converted=" << nConverted
            << " skipped=" << nSkipped << endmsg;

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[SiTargetMeasConverter] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[SiTargetMeasConverter] Unknown exception in execute()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode SiTargetMeasConverter::finalize() {
  try {
    m_inputHandle.reset();
    m_outputHandle.reset();
    info() << "[SiTargetMeasConverter] Done. Events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[SiTargetMeasConverter] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[SiTargetMeasConverter] Unknown exception in finalize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(SiTargetMeasConverter)
