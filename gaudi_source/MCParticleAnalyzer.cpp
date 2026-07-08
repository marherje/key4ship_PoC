#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"

#include "edm4hep/MCParticleCollection.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include <array>
#include <atomic>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// MCParticleAnalyzer — basic MC-truth analysis for (multi-)Particle-Gun runs.
//
// Reads the MCParticles collection written by ddsim and fills histograms:
//   * energy of the generator primaries (generatorStatus == 1)
//   * angular distributions (theta, phi) and the opening angle between the
//     first two primaries
//   * primary vertex position
//   * number of MCParticles / primaries / final-state particles per event
//   * primaries vs final-state comparison (multiplicity and energy sums)
//
// Extend here later for calorimeter studies (deposited energy, hits, shower
// shape, PID): add DataHandles for the SimCalorimeterHit collections and fill
// additional histograms in execute().
// ---------------------------------------------------------------------------

class MCParticleAnalyzer : public Gaudi::Algorithm {
public:
  MCParticleAnalyzer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;

      m_mcHandle = std::make_unique<
          k4FWCore::DataHandle<edm4hep::MCParticleCollection>>(
          m_inputCollection.value(), Gaudi::DataHandle::Reader, this);

      // -- book histograms ---------------------------------------------------
      const double eMax = m_maxEnergy.value();
      h_nMC      = new TH1D("h_nMCParticles", "MCParticles per event;N_{MC};events",
                            200, 0.0, m_maxNMC.value());
      h_nPrim    = new TH1D("h_nPrimaries", "Generator primaries per event;N_{prim};events",
                            10, -0.5, 9.5);
      h_nFinal   = new TH1D("h_nFinal", "Final-state particles per event (no daughters);N_{final};events",
                            200, 0.0, m_maxNMC.value());
      h_pdgPrim  = new TH1D("h_pdg_primary", "PDG ID of primaries;PDG;entries",
                            501, -250.5, 250.5);
      h_EPrim    = new TH1D("h_E_primary", "Primary energy;E [GeV];entries",
                            200, 0.0, eMax);
      h_thetaPrim = new TH1D("h_theta_primary", "Primary polar angle;#theta [deg];entries",
                             300, 0.0, 15.0);
      h_phiPrim   = new TH1D("h_phi_primary", "Primary azimuth;#phi [deg];entries",
                             72, -180.0, 180.0);
      h_openAngle = new TH1D("h_openingAngle", "Opening angle between first two primaries;#alpha [deg];entries",
                             500, 0.0, 10.0);
      h_vtxX = new TH1D("h_vtx_x", "Primary vertex x;x [mm];entries", 100, -50.0, 50.0);
      h_vtxY = new TH1D("h_vtx_y", "Primary vertex y;y [mm];entries", 100, -50.0, 50.0);
      h_vtxZ = new TH1D("h_vtx_z", "Primary vertex z;z [mm];entries", 200, -1100.0, 900.0);
      h_EFinal = new TH1D("h_E_final", "Final-state particle energy;E [GeV];entries",
                          200, 0.0, eMax);
      h_ESumPrim  = new TH1D("h_Esum_primary", "#Sigma E of primaries per event;E [GeV];events",
                             200, 0.0, 2.0 * eMax);
      h_ESumFinal = new TH1D("h_Esum_final", "#Sigma E of final-state particles per event;E [GeV];events",
                             200, 0.0, 2.0 * eMax);
      h_nPrimVsFinal = new TH2D("h_nPrim_vs_nFinal",
                                "Primaries vs final-state;N_{prim};N_{final}",
                                10, -0.5, 9.5, 100, 0.0, m_maxNMC.value());

      for (TH1* h : allHistos()) h->SetDirectory(nullptr);

      info() << "[MCParticleAnalyzer] Reading '" << m_inputCollection.value()
             << "', writing histograms to '" << m_outputFile.value() << "'"
             << endmsg;
      return sc;
    } catch (const std::exception& e) {
      error() << "[MCParticleAnalyzer] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
      const auto* mc = m_mcHandle->get();
      if (mc == nullptr) {
        warning() << "[MCParticleAnalyzer] No MCParticles collection" << endmsg;
        return StatusCode::SUCCESS;
      }
      const long long evt = m_eventCount.fetch_add(1);

      int    nPrim = 0, nFinal = 0;
      double eSumPrim = 0.0, eSumFinal = 0.0;
      std::vector<std::array<double, 3>> primDirs;

      for (const auto& p : *mc) {
        const auto  mom = p.getMomentum();
        const double pMag = std::sqrt(mom.x * mom.x + mom.y * mom.y + mom.z * mom.z);
        const double e    = std::sqrt(pMag * pMag + p.getMass() * p.getMass());

        if (p.getGeneratorStatus() == 1) {
          ++nPrim;
          eSumPrim += e;
          h_pdgPrim->Fill(static_cast<double>(p.getPDG()));
          h_EPrim->Fill(e);
          const auto vtx = p.getVertex();
          h_vtxX->Fill(vtx.x);
          h_vtxY->Fill(vtx.y);
          h_vtxZ->Fill(vtx.z);
          if (pMag > 0.0) {
            const double theta = std::acos(mom.z / pMag) * 180.0 / M_PI;
            const double phi   = std::atan2(mom.y, mom.x) * 180.0 / M_PI;
            h_thetaPrim->Fill(theta);
            h_phiPrim->Fill(phi);
            primDirs.push_back({mom.x / pMag, mom.y / pMag, mom.z / pMag});
          }
          if (evt < 3) {
            info() << "[MCParticleAnalyzer] evt=" << evt
                   << " primary PDG=" << p.getPDG()
                   << " E=" << e << " GeV"
                   << " p=(" << mom.x << ", " << mom.y << ", " << mom.z << ")"
                   << " vtx=(" << p.getVertex().x << ", " << p.getVertex().y
                   << ", " << p.getVertex().z << ") mm" << endmsg;
          }
        }
        if (p.getDaughters().empty()) {
          ++nFinal;
          eSumFinal += e;
          h_EFinal->Fill(e);
        }
      }

      if (primDirs.size() >= 2) {
        const auto& a = primDirs[0];
        const auto& b = primDirs[1];
        double cosA = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        cosA = std::max(-1.0, std::min(1.0, cosA));
        h_openAngle->Fill(std::acos(cosA) * 180.0 / M_PI);
      }

      h_nMC->Fill(static_cast<double>(mc->size()));
      h_nPrim->Fill(static_cast<double>(nPrim));
      h_nFinal->Fill(static_cast<double>(nFinal));
      h_ESumPrim->Fill(eSumPrim);
      h_ESumFinal->Fill(eSumFinal);
      h_nPrimVsFinal->Fill(static_cast<double>(nPrim), static_cast<double>(nFinal));

      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[MCParticleAnalyzer] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      TFile fout(m_outputFile.value().c_str(), "RECREATE");
      for (TH1* h : allHistos()) h->Write();
      fout.Close();

      info() << "[MCParticleAnalyzer] Done. Events: " << m_eventCount.load()
             << "  <N primaries>=" << h_nPrim->GetMean()
             << "  <E primary>=" << h_EPrim->GetMean() << " GeV"
             << "  <opening angle>=" << h_openAngle->GetMean() << " deg"
             << "  histograms: " << m_outputFile.value() << endmsg;

      for (TH1* h : allHistos()) delete h;
      m_mcHandle.reset();
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[MCParticleAnalyzer] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  std::vector<TH1*> allHistos() const {
    return {h_nMC, h_nPrim, h_nFinal, h_pdgPrim, h_EPrim, h_thetaPrim,
            h_phiPrim, h_openAngle, h_vtxX, h_vtxY, h_vtxZ, h_EFinal,
            h_ESumPrim, h_ESumFinal, h_nPrimVsFinal};
  }

  Gaudi::Property<std::string> m_inputCollection{
      this, "InputCollection", "MCParticles",
      "Input edm4hep::MCParticleCollection name"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "mcparticle_histos.root",
      "Output ROOT file for the histograms"};
  Gaudi::Property<double> m_maxEnergy{
      this, "MaxEnergy", 20.0,
      "Upper edge [GeV] of the energy histograms"};
  Gaudi::Property<double> m_maxNMC{
      this, "MaxNMCParticles", 2000.0,
      "Upper edge of the multiplicity histograms"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::MCParticleCollection>>
      m_mcHandle;
  mutable std::atomic<long long> m_eventCount{0};

  TH1D* h_nMC = nullptr;
  TH1D* h_nPrim = nullptr;
  TH1D* h_nFinal = nullptr;
  TH1D* h_pdgPrim = nullptr;
  TH1D* h_EPrim = nullptr;
  TH1D* h_thetaPrim = nullptr;
  TH1D* h_phiPrim = nullptr;
  TH1D* h_openAngle = nullptr;
  TH1D* h_vtxX = nullptr;
  TH1D* h_vtxY = nullptr;
  TH1D* h_vtxZ = nullptr;
  TH1D* h_EFinal = nullptr;
  TH1D* h_ESumPrim = nullptr;
  TH1D* h_ESumFinal = nullptr;
  TH2D* h_nPrimVsFinal = nullptr;
};

DECLARE_COMPONENT(MCParticleAnalyzer)
