#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"

#include "SNDCellDecoder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDLayerClusterer — groups a layer's hits into clusters, one view at a time.
//
// This is CyberPFA's Calo1DCluster stage adapted to the SND: there, long
// crystal bars are read out in two orthogonal views and the reconstruction
// builds per-layer, per-view clusters before ever attempting three
// dimensions. Here the same shape appears three times over -- SiTarget X/Y
// strips, MTC SciFi U/V strips, and the 2D pads and tiles of the SiPad and
// the MTC scintillator -- so one algorithm, instantiated per collection,
// covers all of them.
//
// Two steps, both taken inside one (layer, view):
//
//   1. connectivity. Hits whose channels are neighbours join the same proto
//      cluster: adjacent strips in 1D (a gap of up to MaxGapChannels dead or
//      sub-threshold channels is bridged), 8-connected cells in 2D.
//
//   2. splitting on local maxima. A proto cluster holding more than one
//      energy maximum is usually more than one particle, so each maximum
//      above SeedMIP that is at least MinSeedSeparation channels away from a
//      stronger one seeds a cluster, and every hit goes to the nearest seed.
//      A single maximum is left alone.
//
// The splitting is deliberately the simple nearest-seed rule rather than
// CyberPFA's energy sharing between overlapping showers: at this stage the
// clusters feed the muon corridor, which asks whether a layer holds an
// isolated MIP, and that question is answered by separating the maxima, not
// by splitting the energy of the hits between them. Energy sharing belongs
// with the EM shower work of phase 6, if it is ever needed.
// ---------------------------------------------------------------------------

class SNDLayerClusterer : public Gaudi::Algorithm {
public:
  SNDLayerClusterer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    if (m_inputCollection.value().empty() ||
        m_outputCollection.value().empty() ||
        m_readoutDescriptor.value().empty()) {
      error() << "[SNDLayerClusterer] InputCollection, OutputCollection and "
                 "ReadoutDescriptor must all be set." << endmsg;
      return StatusCode::FAILURE;
    }

    m_decoder = std::make_unique<snd::SNDCellDecoder>(
        std::map<std::string, std::string>{
            {m_inputCollection.value(), m_readoutDescriptor.value()}},
        std::map<std::string, snd::Det>{
            {m_inputCollection.value(),
             static_cast<snd::Det>(m_detectorID.value())}},
        m_mtcLayersPerStation.value());

    m_in = std::make_unique<DataHandle<edm4hep::CalorimeterHitCollection>>(
        m_inputCollection.value(), Gaudi::DataHandle::Reader, this);
    m_out = std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    info() << "[SNDLayerClusterer] " << m_inputCollection.value() << " -> "
           << m_outputCollection.value()
           << "  maxGap=" << m_maxGapChannels.value()
           << " seedMIP=" << m_seedMip.value()
           << " minSeedSep=" << m_minSeedSeparation.value() << endmsg;
    return sc;
  }

  StatusCode execute(const EventContext&) const override {
    auto clusters = std::make_unique<edm4hep::ClusterCollection>();
    const auto* hits = m_in->get();
    if (!hits) {
      m_out->put(std::move(clusters));
      return StatusCode::SUCCESS;
    }

    // Group by (layer, view). Both come from the cellID, so a collection that
    // mixes views -- the MTC SciFi one holds U and V -- is handled here and
    // not by the caller.
    std::map<std::pair<int, int>, std::vector<HitRef>> byLayerView;

    for (const auto& hit : *hits) {
      const auto a = m_decoder->decode(m_inputCollection.value(),
                                       hit.getCellID());
      if (!a.valid) continue;
      HitRef r;
      r.hit = hit;
      m_decoder->channel(m_inputCollection.value(), hit.getCellID(),
                         &r.cx, &r.cy);
      r.twoD = (a.view == snd::View::TwoD);
      byLayerView[{a.layerOrdinal, static_cast<int>(a.view)}].push_back(r);
    }

    for (auto& [key, group] : byLayerView) {
      const int layerOrdinal = key.first;
      const int view = key.second;
      const bool twoD = !group.empty() && group.front().twoD;

      for (auto& proto : connectedGroups(group, twoD)) {
        for (auto& part : splitOnMaxima(proto, twoD)) {
          writeCluster(*clusters, part, layerOrdinal, view);
        }
      }
    }

    m_nClusters += clusters->size();
    ++m_nEvents;
    m_out->put(std::move(clusters));
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    if (m_nEvents > 0) {
      info() << "[SNDLayerClusterer] " << m_outputCollection.value() << ": "
             << m_nClusters << " clusters in " << m_nEvents << " events ("
             << static_cast<double>(m_nClusters) / m_nEvents
             << " per event)" << endmsg;
    }
    return Gaudi::Algorithm::finalize();
  }

private:
  // The hit is held BY VALUE. A podio object handed out by a collection's
  // iterator is a temporary handle, so keeping its address gives a pointer
  // into a buffer that the next iteration reuses -- every stored hit then
  // reads back as the same one (found 2026-09-23: clusters whose hits all had
  // one cellID). Copying the handle is cheap and refers to the real object.
  struct HitRef {
    edm4hep::CalorimeterHit hit;
    int cx = 0, cy = 0;
    bool twoD = false;
  };

  // Hits of one (layer, view), split into connected groups.
  std::vector<std::vector<HitRef>> connectedGroups(
      std::vector<HitRef>& group, bool twoD) const {
    std::vector<std::vector<HitRef>> out;
    if (group.empty()) return out;

    if (!twoD) {
      std::sort(group.begin(), group.end(),
                [](const HitRef& a, const HitRef& b) { return a.cx < b.cx; });
      std::vector<HitRef> cur{group.front()};
      for (std::size_t i = 1; i < group.size(); ++i) {
        if (group[i].cx - cur.back().cx <= m_maxGapChannels.value() + 1) {
          cur.push_back(group[i]);
        } else {
          out.push_back(std::move(cur));
          cur = {group[i]};
        }
      }
      out.push_back(std::move(cur));
      return out;
    }

    // 2D: flood fill over 8-connected neighbours, the gap allowance widening
    // the neighbourhood the same way it bridges dead strips in 1D.
    const int reach = m_maxGapChannels.value() + 1;
    std::vector<bool> used(group.size(), false);
    for (std::size_t i = 0; i < group.size(); ++i) {
      if (used[i]) continue;
      std::vector<HitRef> cur;
      std::vector<std::size_t> stack{i};
      used[i] = true;
      while (!stack.empty()) {
        const std::size_t k = stack.back();
        stack.pop_back();
        cur.push_back(group[k]);
        for (std::size_t j = 0; j < group.size(); ++j) {
          if (used[j]) continue;
          if (std::abs(group[j].cx - group[k].cx) <= reach &&
              std::abs(group[j].cy - group[k].cy) <= reach) {
            used[j] = true;
            stack.push_back(j);
          }
        }
      }
      out.push_back(std::move(cur));
    }
    return out;
  }

  // Splits a proto cluster at its local maxima; returns it unchanged when it
  // holds at most one.
  std::vector<std::vector<HitRef>> splitOnMaxima(
      const std::vector<HitRef>& proto, bool twoD) const {
    std::vector<std::vector<HitRef>> out;
    if (proto.size() < 3) {
      out.push_back(proto);
      return out;
    }

    // Candidate seeds: above SeedMIP and no neighbour strictly higher.
    const int reach = m_maxGapChannels.value() + 1;
    std::vector<std::size_t> seeds;
    for (std::size_t i = 0; i < proto.size(); ++i) {
      const float e = proto[i].hit.getEnergyError();   // MIPs
      if (e < m_seedMip.value()) continue;
      bool isMax = true;
      for (std::size_t j = 0; j < proto.size() && isMax; ++j) {
        if (i == j) continue;
        const int dx = std::abs(proto[j].cx - proto[i].cx);
        const int dy = twoD ? std::abs(proto[j].cy - proto[i].cy) : 0;
        if (dx <= reach && dy <= reach &&
            proto[j].hit.getEnergyError() > e) isMax = false;
      }
      if (isMax) seeds.push_back(i);
    }

    // Keep the strongest seeds first and drop any that sits too close to one
    // already kept: two maxima a channel apart are one shower fluctuating,
    // not two particles.
    std::sort(seeds.begin(), seeds.end(),
              [&](std::size_t a, std::size_t b) {
                return proto[a].hit.getEnergyError() >
                       proto[b].hit.getEnergyError();
              });
    std::vector<std::size_t> kept;
    for (std::size_t s : seeds) {
      bool ok = true;
      for (std::size_t k : kept) {
        const int dx = std::abs(proto[s].cx - proto[k].cx);
        const int dy = twoD ? std::abs(proto[s].cy - proto[k].cy) : 0;
        if (std::max(dx, dy) < m_minSeedSeparation.value()) { ok = false; break; }
      }
      if (ok) kept.push_back(s);
    }

    if (kept.size() < 2) {
      out.push_back(proto);
      return out;
    }

    out.resize(kept.size());
    for (const auto& h : proto) {
      std::size_t best = 0;
      int bestD = std::numeric_limits<int>::max();
      for (std::size_t k = 0; k < kept.size(); ++k) {
        const int dx = std::abs(h.cx - proto[kept[k]].cx);
        const int dy = twoD ? std::abs(h.cy - proto[kept[k]].cy) : 0;
        const int d = std::max(dx, dy);
        if (d < bestD) { bestD = d; best = k; }
      }
      out[best].push_back(h);
    }
    out.erase(std::remove_if(out.begin(), out.end(),
                             [](const std::vector<HitRef>& v) { return v.empty(); }),
              out.end());
    return out;
  }

  void writeCluster(edm4hep::ClusterCollection& coll,
                    const std::vector<HitRef>& hits,
                    int layerOrdinal, int view) const {
    if (hits.empty()) return;
    double eSum = 0.0, mipSum = 0.0;
    double sx = 0.0, sy = 0.0, sz = 0.0;
    for (const auto& h : hits) {
      const double e = h.hit.getEnergy();
      eSum += e;
      mipSum += h.hit.getEnergyError();
      const auto& p = h.hit.getPosition();
      sx += e * p.x; sy += e * p.y; sz += e * p.z;
    }
    if (mipSum < m_minClusterMip.value()) return;
    const double norm = (eSum > 0.0) ? eSum : static_cast<double>(hits.size());
    if (eSum <= 0.0) {
      for (const auto& h : hits) {
        const auto& p = h.hit.getPosition();
        sx += p.x; sy += p.y; sz += p.z;
      }
    }

    auto cl = coll.create();
    cl.setEnergy(static_cast<float>(eSum));
    cl.setPosition({static_cast<float>(sx / norm),
                    static_cast<float>(sy / norm),
                    static_cast<float>(sz / norm)});
    cl.setType(m_detectorID.value() * 10 + view);

    // Extent of the cluster in the coordinates its view measures. The MTC's
    // stereo strips report their centre at y = 0, so for them -- as for the
    // SiTarget X view -- that coordinate is x, and only the SiTarget Y view
    // measures y. A 2D readout is measured in both: tiles that are
    // contiguous in y share one x, so an x-only extent reads back as zero.
    const bool measuresY = (m_detectorID.value() == 0 && view == 1);
    const bool twoDView = (view == 2);
    double loX = 1e30, hiX = -1e30, loY = 1e30, hiY = -1e30;
    for (const auto& h : hits) {
      const auto& p = h.hit.getPosition();
      loX = std::min(loX, static_cast<double>(p.x));
      hiX = std::max(hiX, static_cast<double>(p.x));
      loY = std::min(loY, static_cast<double>(p.y));
      hiY = std::max(hiY, static_cast<double>(p.y));
    }
    const double extX = hiX - loX;
    const double extY = hiY - loY;
    const double width = twoDView ? std::hypot(extX, extY)
                                  : (measuresY ? extY : extX);

    cl.addToShapeParameters(static_cast<float>(layerOrdinal));
    cl.addToShapeParameters(static_cast<float>(hits.size()));
    cl.addToShapeParameters(static_cast<float>(width));
    cl.addToShapeParameters(static_cast<float>(sz / norm));
    cl.addToShapeParameters(static_cast<float>(mipSum));
    for (const auto& h : hits) cl.addToHits(h.hit);
  }

  Gaudi::Property<std::string> m_inputCollection{
      this, "InputCollection", "", "CalorimeterHit collection to cluster"};
  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "", "Cluster collection to write"};
  Gaudi::Property<std::string> m_readoutDescriptor{
      this, "ReadoutDescriptor", "", "cellID descriptor of the readout"};
  Gaudi::Property<int> m_detectorID{
      this, "DetectorID", -1, "0=SiTarget, 1=SiPad, 2=MTC"};
  Gaudi::Property<int> m_mtcLayersPerStation{
      this, "MTCLayersPerStation", 14, "DD4hep constant MTC_NLayers"};
  Gaudi::Property<int> m_maxGapChannels{
      this, "MaxGapChannels", 1,
      "Dead or sub-threshold channels bridged when joining neighbours. 1 for "
      "strips, 0 for pads and tiles, whose cells are large enough that a gap "
      "is a real separation."};
  Gaudi::Property<double> m_seedMip{
      this, "SeedMIP", 2.0,
      "A local maximum must reach this many MIPs to seed a split"};
  Gaudi::Property<int> m_minSeedSeparation{
      this, "MinSeedSeparation", 3,
      "Channels two seeds must be apart; closer maxima are one shower "
      "fluctuating rather than two particles"};
  Gaudi::Property<double> m_minClusterMip{
      this, "MinClusterMIP", 0.5, "Drop clusters below this many MIPs"};

  std::unique_ptr<snd::SNDCellDecoder> m_decoder;
  std::unique_ptr<DataHandle<edm4hep::CalorimeterHitCollection>> m_in;
  std::unique_ptr<DataHandle<edm4hep::ClusterCollection>> m_out;

  mutable long m_nClusters = 0;
  mutable long m_nEvents = 0;
};

DECLARE_COMPONENT(SNDLayerClusterer)
