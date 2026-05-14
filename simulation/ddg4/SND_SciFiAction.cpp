//====================================================================
// SND_SciFiAction: custom DDG4 sensitive action for MTC SciFi planes
//
// Implements the standard calorimeter hit accumulation (one hit per cellID,
// N contributions per hit) and additionally encodes the energy-weighted
// average step-Y into hit.position.y at end-of-event.
//
// Why position.y?  CartesianStripX segmentation always leaves position.y = 0
// (no Y segmentation). The DDG4→EDM4HEP mapper does NOT transfer
// MonteCarloContrib::y (the contribution step position) to
// CaloHitContribution.stepPosition in key4hep-2026-02-01. Repurposing
// position.y is the only field the mapper reliably transfers that is
// otherwise unused for this detector type.  SciFiDigitizer reads it to
// compute the propagation distance d(hit, SiPM) for light attenuation.
//====================================================================

#define GEANT4_HIT_COLLECTION_WRAPPER_DEFINED 1

#include <DDG4/Geant4SensDetAction.h>
#include <DDG4/Geant4Data.h>
#include <DDG4/Geant4HitCollection.h>
#include <DDG4/Factories.h>
#include <G4Step.hh>

#include <unordered_map>

namespace dd4hep {
namespace sim {

class SND_SciFiAction : public Geant4Sensitive {
public:
  // Factory signature: DECLARE_GEANT4SENSITIVE dereferences DetElement* and Detector*
  // before calling the constructor, so we take them by value / reference.
  SND_SciFiAction(Geant4Context* ctxt, const std::string& nam,
                  DetElement det, Detector& desc)
      : Geant4Sensitive(ctxt, nam, det, desc) {}

  ~SND_SciFiAction() override = default;

  // Called once during detector initialisation to declare the hit collection.
  void defineCollections() override {
    defineCollection<Geant4Calorimeter::Hit>(m_readout.name());
  }

  // Reset per-event accumulators at the start of every G4 event.
  void begin(G4HCofThisEvent* hce) override {
    Geant4Sensitive::begin(hce);
    m_sumE.clear();
    m_sumEY.clear();
  }

  // Called for every G4 step inside a sensitive MTC SciFi volume.
  bool process(const G4Step* step, G4TouchableHistory* /*hist*/) override {
    const double edep = step->GetTotalEnergyDeposit();
    if (!accept(step)) return false;
    if (edep <= 0.0)   return false;

    // Standard contribution: trackID, PDG, deposit, time, step-length, step-pos, momentum
    auto contrib = Geant4HitData::extractContribution(step);
    const long long cid = cellID(step);

    Geant4HitCollection* coll = collection(0);
    auto* hit = coll->findByKey<Geant4Calorimeter::Hit>(cid);

    if (!hit) {
      // Initial position: strip-center X from segmentation (TGeo cm → mm), step-midpoint Y/Z
      const G4ThreeVector mid =
          0.5 * (step->GetPreStepPoint()->GetPosition() +
                 step->GetPostStepPoint()->GetPosition());  // global mm
      const Position seg_p = m_segmentation.position(cid); // local→approx global for centred placement
      hit = new Geant4Calorimeter::Hit();
      hit->cellID   = cid;
      hit->position = Position(seg_p.x() * 10.0, mid.y(), mid.z());
      coll->add(cid, hit);
    }

    hit->truth.push_back(contrib);
    hit->energyDeposit += edep;

    // Accumulate energy-weighted Y for this cell
    const double y_mid =
        0.5 * (step->GetPreStepPoint()->GetPosition().y() +
               step->GetPostStepPoint()->GetPosition().y());
    m_sumE [cid] += edep;
    m_sumEY[cid] += edep * y_mid;

    mark(step);
    return true;
  }

  // After all steps: overwrite hit.position.y with energy-weighted avg step-Y
  // so that SciFiDigitizer can compute the SiPM propagation distance.
  void end(G4HCofThisEvent* hce) override {
    Geant4HitCollection* coll = collection(0);
    if (coll) {
      for (const auto& [cid, sumE] : m_sumE) {
        if (sumE <= 0.0) continue;
        auto* hit = coll->findByKey<Geant4Calorimeter::Hit>(cid);
        if (hit) {
          const double y_avg = m_sumEY.at(cid) / sumE;
          hit->position = Position(hit->position.x(), y_avg, hit->position.z());
        }
      }
    }
    m_sumE.clear();
    m_sumEY.clear();
    Geant4Sensitive::end(hce);
  }

private:
  std::unordered_map<long long, double> m_sumE;
  std::unordered_map<long long, double> m_sumEY;
};

}  // namespace sim
}  // namespace dd4hep

DECLARE_GEANT4SENSITIVE_NS(dd4hep::sim, SND_SciFiAction)
