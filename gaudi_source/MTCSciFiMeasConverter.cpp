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
// MTCSciFiMeasConverter
//
// Converts MTCSciFiHitsWindowed (SimCalorimeterHitCollection of raw windowed
// hits) into TrackerHit3D measurements for ACTS tracking.
//
// Key unit note: position.x from SND_SciFiAction is in ROOT cm (not mm).
// It must be multiplied by 10 to get mm.  position.y (truth avg-Y for light
// propagation) is NOT stored here — it must never be used in tracking.
// position.z is in Geant4 mm and is stored as-is.
//
// Surface matching uses ACTSGeoSvc::surfaceByAddress(detID=2, station, layer,
// plane) decoded from the hit cellID — no Z-proximity arithmetic.
// ---------------------------------------------------------------------------

class MTCSciFiMeasConverter : public Gaudi::Algorithm {
public:
  MTCSciFiMeasConverter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  // ---- Properties ----
  Gaudi::Property<std::string> m_inputCollection{
      this, "InputCollection", "MTCSciFiHitsWindowed",
      "Input SimCalorimeterHitCollection name (windowed raw hits)"};

  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "MTCSciFiMeasurements",
      "Output TrackerHit3DCollection name"};

  Gaudi::Property<std::string> m_bitField{
      this, "BitField",
      "system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9",
      "DD4hep BitField string for MTC cellID decoding"};

  Gaudi::Property<double> m_stripPitch{
      this, "StripPitch", 1.0,
      "SciFi strip pitch [mm] — determines measurement variance"};

  Gaudi::Property<double> m_stereoAngleDeg{
      this, "StereoAngleDeg", 5.0,
      "Fiber stereo angle magnitude [degrees] (unused in conversion, stored for reference)"};

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

StatusCode MTCSciFiMeasConverter::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    m_inputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inputCollection.value(), Gaudi::DataHandle::Reader, this);

    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
        m_bitField.value());

    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[MTCSciFiMeasConverter] Failed to retrieve ACTSGeoSvc."
              << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "[MTCSciFiMeasConverter] Initialized. GeoSvc has "
           << m_geoSvc->allSurfaces().size() << " surfaces." << endmsg;

    return sc;
  } catch (const std::exception& e) {
    error() << "[MTCSciFiMeasConverter] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[MTCSciFiMeasConverter] Unknown exception in initialize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode MTCSciFiMeasConverter::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    auto* output = m_outputHandle->createAndPut();

    const auto* hits = m_inputHandle->get();
    if (!hits || hits->size() == 0) {
      debug() << "[MTCSciFiMeasConverter] evt=" << evtNum
              << " no MTC SciFi hits." << endmsg;
      return StatusCode::SUCCESS;
    }

    const double pitch    = m_stripPitch.value();      // mm
    const double variance = (pitch * pitch) / 12.0;   // mm^2

    std::size_t nConverted = 0;
    std::size_t nSkipped   = 0;

    for (std::size_t i = 0; i < hits->size(); ++i) {
      const auto&    hit = (*hits)[i];
      const uint64_t cid = hit.getCellID();

      // Decode detector address from cellID
      int plane   = static_cast<int>(m_decoder->get(cid, "plane"));
      int station = static_cast<int>(m_decoder->get(cid, "station"));
      int layer   = static_cast<int>(m_decoder->get(cid, "layer"));

      // Skip scintillator hits (plane 2); only track SciFi U (0) and V (1).
      if (plane >= 2) {
        ++nSkipped;
        continue;
      }

      // Decode position from hit.
      // IMPORTANT: position.x from SND_SciFiAction is in ROOT cm — multiply
      // by 10 to get mm.  position.y is the truth energy-weighted step-Y used
      // only by SciFiDigitizer for light propagation; must NOT be used here.
      const auto& pos = hit.getPosition();
      const double localCoord = pos.x * 10.0;  // ROOT cm → mm

      // Find matching ACTS surface via detector address lookup.
      const Acts::Surface* surf = m_geoSvc->surfaceByAddress(2, station, layer, plane);
      if (!surf) {
        warning() << "[MTCSciFiMeasConverter] evt=" << evtNum
                  << " hit=" << i
                  << " station=" << station
                  << " layer=" << layer
                  << " plane=" << plane
                  << " has no matching ACTS surface — skipping." << endmsg;
        ++nSkipped;
        continue;
      }

      float hitTime = 0.0f;
      if (hit.contributions_size() > 0)
        hitTime = static_cast<float>(hit.getContributions(0).getTime());

      // 1D covariance: only the xx element (eBoundLoc0) carries the strip
      // precision.  The surface is pre-rotated for stereo in ACTSGeoSvc, so
      // eBoundLoc0 on this surface corresponds to strip_centre_at_y0.
      edm4hep::CovMatrix3f cov{};
      cov[0] = static_cast<float>(variance);  // xx = σ²_strip

      auto mhit = output->create();
      mhit.setCellID(cid);
      mhit.setType(2);          // detector ID: 2 = MTC
      mhit.setQuality(plane);   // 0 = U, 1 = V
      mhit.setTime(hitTime);
      mhit.setEDep(hit.getEnergy());
      mhit.setEDepError(0.0f);
      // Store strip coordinate in position.x (mm), Y=0 (1D strip), Z from hit.
      mhit.setPosition(edm4hep::Vector3d{localCoord, 0.0, pos.z});
      mhit.setCovMatrix(cov);

      ++nConverted;
    }

    debug() << "[MTCSciFiMeasConverter] evt=" << evtNum
            << " input=" << hits->size()
            << " converted=" << nConverted
            << " skipped=" << nSkipped << endmsg;

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[MTCSciFiMeasConverter] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[MTCSciFiMeasConverter] Unknown exception in execute()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode MTCSciFiMeasConverter::finalize() {
  try {
    m_inputHandle.reset();
    m_outputHandle.reset();
    info() << "[MTCSciFiMeasConverter] Done. Events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[MTCSciFiMeasConverter] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[MTCSciFiMeasConverter] Unknown exception in finalize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(MTCSciFiMeasConverter)
