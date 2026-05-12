#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <atomic>
#include <cmath>
#include <memory>
#include <random>

// Simulates MTC SciFi light attenuation, photon statistics, SiPM
// saturation, and QDC conversion following the FairRoot MTCDetHit model.
//
// Requires SND_SciFiAction DDG4 plugin active during simulation so that
// hit.position.y carries the energy-weighted average step Y (mm) rather
// than the CartesianStripX centroid (always 0).
//
// Bit layout for MTCDetHits cellID:
//   system:8, station:2, layer:8, slice:4, plane:2, ...
//   station = (cellID >> 8) & 0x3   → index into EnvHeightHalf
//   plane   = (cellID >> 22) & 0x3  → 0=SciFi U, 1=SciFi V

class SciFiDigitizer : public Gaudi::Algorithm {
public:
  SciFiDigitizer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(),  Gaudi::DataHandle::Reader, this);
      m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_outputName.value(), Gaudi::DataHandle::Writer, this);

      const auto& hh = m_envHeightHalf.value();
      if (hh.size() < 3) {
        error() << "[SciFiDigitizer] EnvHeightHalf must have at least 3 entries (MTC40/50/60)" << endmsg;
        return StatusCode::FAILURE;
      }
      return Gaudi::Algorithm::initialize();
    } catch (const std::exception& e) {
      error() << "[SciFiDigitizer] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* input  = m_inputHandle->get();
      auto*       output = m_outputHandle->createAndPut();

      const double att_len    = m_attLen.value();     // cm
      const double att_off    = m_attOff.value();     // cm
      const double ph_per_gev = m_photonsPerGeV.value();
      const double R_mirror   = m_mirrorR.value();
      const double N_pix      = m_sipmNpix.value();
      const double n_thresh   = m_photonThresh.value();
      const double qdc_A      = m_qdcA.value();
      const double qdc_B      = m_qdcB.value();
      const double qdc_sA     = m_qdcSigA.value();
      const double qdc_sB     = m_qdcSigB.value();
      const int    sipm_side  = m_sipmSide.value();   // +1 or -1
      const double angle_rad  = m_fiberAngle.value() * M_PI / 180.0;
      const double tan_a      = std::tan(angle_rad);
      const double cos_a      = std::cos(angle_rad);
      const auto&  env_hh     = m_envHeightHalf.value();  // mm

      const long long evtNum = m_eventCount.fetch_add(1);
      const bool doPrint = (evtNum % m_debugFreq == 0);

      int n_pass = 0;
      for (const auto& hit : *input) {
        const uint64_t cellID   = hit.getCellID();
        const int station_id    = static_cast<int>((cellID >> 8)  & 0x3);
        const int plane         = static_cast<int>((cellID >> 22) & 0x3);

        // Only SciFi planes (0=U, 1=V). Skip scintillator (plane=2).
        if (plane >= 2) continue;
        if (station_id >= static_cast<int>(env_hh.size())) continue;

        const double h_half = env_hh[station_id];          // mm

        // Hit position: x = strip center (mm), y = energy-weighted avg Y (mm)
        const auto&  pos    = hit.getPosition();
        const double x_hit  = pos.x;
        const double y_hit  = pos.y;

        // SiPM end of this specific fiber.
        // U plane (plane=0) tilts +angle; V plane (plane=1) tilts -angle.
        const double sign    = (plane == 0) ? -1.0 : +1.0;
        const double y_sipm  = sipm_side * h_half;                        // mm
        const double x_sipm  = x_hit + sign * y_sipm * tan_a;            // mm

        // 2D distance hit → SiPM in the XY plane (Z same layer, ignored)
        const double dx    = x_hit - x_sipm;
        const double dy    = y_hit - y_sipm;
        const double d_mm  = std::sqrt(dx * dx + dy * dy);
        const double d_cm  = d_mm * 0.1;

        // Fiber total length (SiPM end → mirror end) in cm
        const double fiber_len_cm = (2.0 * h_half / cos_a) * 0.1;

        // Light yield: direct + mirror reflection (one end + mirror geometry)
        const double att_direct   = std::exp(-(d_cm              - att_off) / att_len);
        const double att_reflect  = R_mirror * std::exp(-((fiber_len_cm - d_cm) - att_off) / att_len);
        const double ly_raw = hit.getEnergy() * ph_per_gev * (att_direct + att_reflect);

        // Poisson fluctuation (clamp to avoid undefined behaviour for tiny negatives)
        const int smeared = std::poisson_distribution<int>{std::max(0.0, ly_raw)}(m_rng);

        // Photon threshold
        if (smeared < static_cast<int>(n_thresh)) continue;

        // SiPM saturation
        const double n_pixels = N_pix * (1.0 - std::exp(-static_cast<double>(smeared) / N_pix));

        // QDC conversion with Gaussian smearing
        const double A      = std::normal_distribution<double>{qdc_A, qdc_sA}(m_rng);
        const double B      = std::normal_distribution<double>{qdc_B, qdc_sB}(m_rng);
        const double signal = A * n_pixels + B;

        auto nh = output->create();
        nh.setCellID(cellID);
        nh.setEnergy(signal);    // QDC counts
        nh.setPosition(hit.getPosition());

        ++n_pass;

        if (doPrint) {
          debug() << "  SciFiHit st=" << station_id << " pl=" << plane
                  << " y=" << y_hit << "mm d=" << d_cm << "cm"
                  << " ly=" << ly_raw << " N=" << smeared
                  << " QDC=" << signal
                  << endmsg;
        }
      }

      debug() << "SciFiDigitizer: " << input->size() << " in, "
              << n_pass << " passing" << endmsg;
      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[SciFiDigitizer] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      m_inputHandle.reset();
      m_outputHandle.reset();
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[SciFiDigitizer] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  // --- I/O ---
  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection",  "MTCSciFiHitsWindowed", "Input SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "MTCSciFiHitsDigi",     "Output SimCalorimeterHit collection"};

  // --- Attenuation ---
  Gaudi::Property<double> m_attLen{
      this, "AttenuationLength", 300.0, "Fiber attenuation length (cm)"};
  Gaudi::Property<double> m_attOff{
      this, "AttenuationOffset", 20.0,  "Attenuation formula offset x0 (cm)"};
  Gaudi::Property<double> m_photonsPerGeV{
      this, "PhotonsPerGeV", 1.6e5,  "Photon yield (photons per GeV deposited)"};

  // --- Fiber geometry (configurable from steering file) ---
  Gaudi::Property<double> m_fiberAngle{
      this, "FiberAngleDeg", 5.0,
      "SciFi stereo angle from Y axis (degrees); from MTC_fiber_angle_deg in compact XML"};
  Gaudi::Property<int> m_sipmSide{
      this, "SiPMSide", +1,
      "+1 = SiPM at +y end, -1 = SiPM at -y end"};
  Gaudi::Property<std::vector<double>> m_envHeightHalf{
      this, "EnvHeightHalf", {200.0, 250.0, 300.0},
      "Per-station envelope half-heights in mm (MTC40, MTC50, MTC60); "
      "from MTC40_env_height/2 etc. in compact XML"};

  // --- SiPM + readout model (FairRoot MTCDetHit defaults) ---
  Gaudi::Property<double> m_mirrorR{
      this, "MirrorReflectivity", 0.9,  "Mirror reflectance at far fiber end"};
  Gaudi::Property<double> m_sipmNpix{
      this, "SiPMNpixels", 104.0, "SiPM pixel count (saturation parameter)"};
  Gaudi::Property<double> m_photonThresh{
      this, "PhotonThreshold", 3.5,  "Min detected photons to register a hit"};
  Gaudi::Property<double> m_sigSpeed{
      this, "SignalSpeedCmNs", 15.0, "Signal propagation speed in fiber (cm/ns)"};
  Gaudi::Property<double> m_qdcA{
      this, "QDC_A",      0.172,  "QDC linear factor A"};
  Gaudi::Property<double> m_qdcB{
      this, "QDC_B",     -1.31,   "QDC offset B"};
  Gaudi::Property<double> m_qdcSigA{
      this, "QDC_sigmaA", 0.006,  "Gaussian smearing on QDC factor A"};
  Gaudi::Property<double> m_qdcSigB{
      this, "QDC_sigmaB", 0.33,   "Gaussian smearing on QDC offset B"};

  // --- Debug ---
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print per-hit debug info every N events"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
  mutable std::atomic<long long> m_eventCount{0};
  mutable std::mt19937 m_rng{42};
};

DECLARE_COMPONENT(SciFiDigitizer)
