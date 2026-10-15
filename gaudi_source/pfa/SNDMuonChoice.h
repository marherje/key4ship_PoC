#pragma once

#include "edm4hep/ParticleIDCollection.h"

// ---------------------------------------------------------------------------
// Which muon candidate is THE muon of the event.
//
// SNDMuonFinder can return more than one corridor; phase 3.3 showed the extra
// ones are pions crossing the MTC. Every stage after it has to agree on which
// one is kept, or one stage subtracts a corridor that another counts, so the
// rule lives here: the candidate whose corridor hit the most MTC layers
// (MuonPID parameter 0).
// ---------------------------------------------------------------------------

namespace snd {

// Returns false when there is no candidate. *nCand, if given, receives the
// number of candidates.
inline bool chooseMuon(const edm4hep::ParticleIDCollection* pids,
                       edm4hep::ParticleID* best, int* nCand = nullptr) {
  bool found = false;
  float bestLayers = -1.f;
  int n = 0;
  if (pids) {
    for (const auto& pid : *pids) {
      ++n;
      const auto params = pid.getParameters();
      const float nLayers = params.empty() ? 0.f : params[0];
      if (nLayers > bestLayers) {
        bestLayers = nLayers;
        *best = pid;
        found = true;
      }
    }
  }
  if (nCand) *nCand = n;
  return found;
}

}  // namespace snd
