#pragma once

#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/Frame.h"

#include <map>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDTruthLink — which MC particle put the energy in a reconstructed object.
//
// The chain keeps that answer in two halves and this header is the only place
// that joins them:
//
//   * the *_Contributions of a SimCalorimeterHit hold the energies;
//   * the frame parameters *ContribParticleIDs and *ContribPDGs hold, one
//     entry per contribution IN COLLECTION ORDER, the MCParticle index and
//     PDG behind it -- the index being local to the source file it was read
//     from, which MCParticleSourceOffsets (paired with MCParticleSourceIDs)
//     shifts into the merged MCParticles collection. EventOverlay stores a
//     contribution's source id in CaloHitContribution::PDG.
//
// A reconstructed hit reaches its sim hit through the CaloHitSimCaloHitLink
// that SNDCaloHitBuilder writes; going by index instead would break the
// moment any stage filters or reorders hits, which the threshold already
// does. Everything here is keyed by cellID, which survives both.
// ---------------------------------------------------------------------------

namespace snd {

class TruthLink {
public:
  // `collections` lists the sim hit collections to index, each paired with the
  // name of its ContribParticleIDs frame parameter.
  TruthLink(const podio::Frame& frame,
            const std::vector<std::pair<std::string, std::string>>& collections,
            const std::string& mcCollection = "MCParticles") {
    m_mc = frame.get<edm4hep::MCParticleCollection>(mcCollection).size() > 0
               ? &frame.get<edm4hep::MCParticleCollection>(mcCollection)
               : nullptr;
    m_offsets = intParam(frame, "MCParticleSourceOffsets");
    m_sourceIDs = intParam(frame, "MCParticleSourceIDs");

    for (const auto& [collName, paramName] : collections) {
      const auto* coll = dynamic_cast<const edm4hep::SimCalorimeterHitCollection*>(
          frame.get(collName));
      if (!coll) continue;
      const auto ids = intParam(frame, paramName);
      if (ids.empty()) continue;

      std::size_t k = 0;
      for (const auto& hit : *coll) {
        auto& acc = m_byCell[hit.getCellID()];
        for (const auto& c : hit.getContributions()) {
          if (k >= ids.size()) break;
          const int idx = globalIndex(ids[k], c.getPDG());
          if (idx >= 0) acc[idx] += c.getEnergy();
          ++k;
        }
      }
    }
  }

  bool hasTruth() const { return m_mc != nullptr && !m_byCell.empty(); }

  // Energy each MC particle left in one cell.
  const std::map<int, double>& cellTruth(uint64_t cellID) const {
    static const std::map<int, double> kEmpty;
    auto it = m_byCell.find(cellID);
    return (it == m_byCell.end()) ? kEmpty : it->second;
  }

  // (MC index, purity) of a set of cells: the particle holding most of the
  // energy, and the fraction it holds.
  std::pair<int, double> dominantOf(const std::map<int, double>& weights) const {
    if (weights.empty()) return {-1, 0.0};
    double total = 0.0;
    int best = -1;
    double bestE = -1.0;
    for (const auto& [idx, e] : weights) {
      total += e;
      if (e > bestE) { bestE = e; best = idx; }
    }
    return {best, (total > 0.0) ? bestE / total : 0.0};
  }

  std::pair<int, double> hitTruth(const edm4hep::CalorimeterHit& hit) const {
    return dominantOf(cellTruth(hit.getCellID()));
  }

  std::pair<int, double> clusterTruth(const edm4hep::Cluster& cluster) const {
    std::map<int, double> acc;
    for (const auto& h : cluster.getHits())
      for (const auto& [idx, e] : cellTruth(h.getCellID())) acc[idx] += e;
    return dominantOf(acc);
  }

  const edm4hep::MCParticleCollection* mcParticles() const { return m_mc; }

private:
  static std::vector<int> intParam(const podio::Frame& frame,
                                   const std::string& name) {
    auto v = frame.getParameter<std::vector<int>>(name);
    return v.has_value() ? v.value() : std::vector<int>{};
  }

  int globalIndex(int localIndex, int sourceID) const {
    if (localIndex < 0 || !m_mc) return -1;
    int off = 0;
    if (!m_offsets.empty()) {
      auto it = std::find(m_sourceIDs.begin(), m_sourceIDs.end(), sourceID);
      if (it != m_sourceIDs.end())
        off = m_offsets[std::distance(m_sourceIDs.begin(), it)];
      else if (m_offsets.size() == 1)
        off = m_offsets[0];
    }
    const int idx = localIndex + off;
    return (idx < static_cast<int>(m_mc->size())) ? idx : -1;
  }

  const edm4hep::MCParticleCollection* m_mc = nullptr;
  std::vector<int> m_offsets;
  std::vector<int> m_sourceIDs;
  std::map<uint64_t, std::map<int, double>> m_byCell;
};

}  // namespace snd
