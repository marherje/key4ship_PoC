#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "podio/detail/OrderKey.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ParticleIDCollection.h"

#include "SNDClue3D.h"
#include "SNDMuonChoice.h"
#include "SNDSamplingWeights.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDHadronicClusterer — 3D clusters of the hadronic system, with CLUE3D.
//
// Phase 4.3 of the plan asked for CyberPFA's cone linking; this uses CLUE3D
// instead (decision of 2026-10-13), the density-based clustering CMS runs on
// HGCAL layer clusters. Two steps:
//
//   1. LAYER POINTS. CLUE3D wants (x, y) per point, and most SND layer
//      clusters measure one coordinate: the SiTarget reads X and Y on
//      separate planes, the MTC SciFi reads two stereo views. Within each
//      layer the views are paired, strongest first, when their energies agree
//      (|E1 - E2| / (E1 + E2) < MaxPairAsymmetry -- the two views sample the
//      same step of the shower, so a real pair carries the same energy). A U/V
//      pair gives (x, y) through the stereo angle exactly as
//      MTCSciFiMeasConverter combines them. SiPad pads and MTC tiles are 2D
//      already. A view left unpaired borrows the missing coordinate from the
//      energy barycentre of the other view in that layer, or of the event.
//
//   2. CLUE3D over the points (SNDClue3D.h), with per-detector critical
//      distances and densities.
//
// The input is every layer cluster except the event muon's corridor, chosen
// by the same rule the PFO builder uses (SNDMuonChoice.h). Corridors of the
// other candidates are pions crossing the MTC and are clustered like any
// other hadron. With ExcludeMuonCorridor = False the corridor goes in too:
// that is the setting for SNDPFOBuilder's MuonSubtraction = "trajectory",
// which takes the muon's expected deposit out of whichever cluster holds it
// instead of removing whole corridor clusters beforehand.
//
// Output: HadronicClusters, type 100 + k, with the layer clusters attached
// (addToClusters), energy, energy barycentre, and in shapeParameters
//   [0] number of layer points   [1] layers spanned
//   [2] first layer index        [3] last layer index
//   [4..6] axis direction (unit, DD4hep axes; (0,0,1) if one layer)
//   [7] density of the seed      [8] energy in the MTC [GeV]
// The axis is a straight energy-weighted fit x(z), y(z) through the points.
// ---------------------------------------------------------------------------

namespace {

// Identity of a podio object, for sets and maps. NOT (collectionID, index):
// a collection made earlier in the same job has no collectionID until it is
// written (every one reads 0xffffffff here), so that pair collides between
// collections -- index 3 of the SiTarget clusters is index 3 of the MTC ones.
// The address of the object's data is what podio's own operator< compares.
using Key = podio::detail::OrderKey;

template <typename T>
Key keyOf(const T& obj) {
  return podio::detail::getOrderKey(obj);
}

}  // namespace

class SNDHadronicClusterer : public Gaudi::Algorithm {
public:
  SNDHadronicClusterer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    for (const auto& c : m_clusterCollections.value())
      m_clusterHandles.push_back(
          std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
              c, Gaudi::DataHandle::Reader, this));
    if (m_clusterHandles.empty()) {
      error() << "[SNDHadronicClusterer] ClusterCollections must be set."
              << endmsg;
      return StatusCode::FAILURE;
    }
    m_inMuonPID = std::make_unique<DataHandle<edm4hep::ParticleIDCollection>>(
        m_muonPID.value(), Gaudi::DataHandle::Reader, this);
    m_out = std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    if (m_dc.value().size() != 3 || m_rhoc.value().size() != 3) {
      error() << "[SNDHadronicClusterer] CriticalDistance and CriticalDensity "
                 "need one value per detector (SiTarget, SiPad, MTC)." << endmsg;
      return StatusCode::FAILURE;
    }
    for (int d = 0; d < 3; ++d) {
      m_par.dc[d] = m_dc.value()[d];
      m_par.rhoc[d] = m_rhoc.value()[d];
    }
    m_par.densityLayers = m_densityLayers.value();
    m_par.nearestHigherLayers = m_nearestHigherLayers.value();
    m_par.outlierDeltaFactor = m_outlierDeltaFactor.value();
    m_tanStereo = std::tan(m_mtcStereoDeg.value() * M_PI / 180.0);

    if (StatusCode gsc = loadLayerSequence(); gsc.isFailure()) return gsc;

    info() << "[SNDHadronicClusterer] CLUE3D dc=(" << m_par.dc[0] << ","
           << m_par.dc[1] << "," << m_par.dc[2] << ") mm rhoc=("
           << m_par.rhoc[0] << "," << m_par.rhoc[1] << "," << m_par.rhoc[2]
           << ") GeV densityLayers=" << m_par.densityLayers
           << " nearestHigherLayers=" << m_par.nearestHigherLayers
           << "; " << m_layerIndex.size() << " layers" << endmsg;
    return sc;
  }

  // Same beam-ordered sequence as SNDPFOBuilder: a layer window has to count
  // the empty layers too.
  StatusCode loadLayerSequence() {
    auto& detector = dd4hep::Detector::getInstance();
    if (detector.state() != dd4hep::Detector::READY) {
      if (m_compactFile.value().empty()) {
        error() << "[SNDHadronicClusterer] no geometry loaded and CompactFile "
                   "empty." << endmsg;
        return StatusCode::FAILURE;
      }
      detector.fromXML(m_compactFile.value());
    }
    std::map<std::pair<int, int>, double> zOfLayer;
    for (const auto& p : snd::collectSensitivePlanes()) {
      const std::pair<int, int> key{p.det, p.layerOrdinal};
      auto it = zOfLayer.find(key);
      if (it == zOfLayer.end() || p.z < it->second) zOfLayer[key] = p.z;
    }
    std::vector<std::pair<double, std::pair<int, int>>> seq;
    for (const auto& [key, z] : zOfLayer) seq.push_back({z, key});
    std::sort(seq.begin(), seq.end());
    for (std::size_t i = 0; i < seq.size(); ++i)
      m_layerIndex[seq[i].second] = static_cast<int>(i);
    return m_layerIndex.empty() ? StatusCode::FAILURE : StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    auto out = std::make_unique<edm4hep::ClusterCollection>();

    // The muon's corridor stays out, unless the PFO builder subtracts the
    // muon along its trajectory instead.
    std::set<Key> corridor;
    edm4hep::ParticleID mpid;
    if (m_excludeCorridor.value() &&
        snd::chooseMuon(m_inMuonPID->get(), &mpid))
      for (const auto& cl : mpid.getParticle().getClusters())
        corridor.insert(keyOf(cl));

    struct Clu {
      edm4hep::Cluster cluster;
      int det, view, layer;
      double z;
    };
    std::map<std::pair<int, int>, std::vector<Clu>> byLayer;
    double ex = 0.0, wx = 0.0, ey = 0.0, wy = 0.0;
    for (const auto& handle : m_clusterHandles) {
      const auto* coll = handle->get();
      if (!coll) continue;
      for (const auto& cl : *coll) {
        if (corridor.count(keyOf(cl))) continue;
        const auto sp = cl.getShapeParameters();
        if (sp.size() < 4) continue;
        Clu c{cl, cl.getType() / 10, cl.getType() % 10,
              static_cast<int>(sp[0]), sp[3]};
        byLayer[{c.det, c.layer}].push_back(c);
        const double e = cl.getEnergy();
        const auto& p = cl.getPosition();
        if (!(c.det == 0 && c.view == 1)) { ex += e * p.x; wx += e; }
        if ((c.det == 0 && c.view == 1) || c.view == 2) { ey += e * p.y; wy += e; }
      }
    }
    const double evX = wx > 0.0 ? ex / wx : 0.0;
    const double evY = wy > 0.0 ? ey / wy : 0.0;

    // --- 1. layer points ----------------------------------------------------
    std::vector<snd::clue::Point> pts;
    std::vector<std::vector<edm4hep::Cluster>> members;
    std::vector<double> ptZ;
    for (auto& [key, clus] : byLayer) {
      const auto [det, layer] = key;
      auto lit = m_layerIndex.find(key);
      if (lit == m_layerIndex.end()) continue;
      const int li = lit->second;

      std::vector<const Clu*> v0, v1;
      for (const auto& c : clus) {
        if (c.view == 2) {
          const auto& p = c.cluster.getPosition();
          pts.push_back({p.x, p.y, li, det, c.cluster.getEnergy()});
          members.push_back({c.cluster});
          ptZ.push_back(c.z);
        } else {
          (c.view == 0 ? v0 : v1).push_back(&c);
        }
      }
      if (v0.empty() && v1.empty()) continue;

      auto byE = [](const Clu* a, const Clu* b) {
        return a->cluster.getEnergy() > b->cluster.getEnergy();
      };
      std::sort(v0.begin(), v0.end(), byE);
      std::sort(v1.begin(), v1.end(), byE);

      // Barycentre of each view in this layer, for the unpaired ones.
      auto bary = [](const std::vector<const Clu*>& v) {
        double s = 0.0, w = 0.0;
        for (const auto* c : v) {
          s += c->cluster.getEnergy() * c->cluster.getPosition().x;
          w += c->cluster.getEnergy();
        }
        return w > 0.0 ? s / w : 0.0;
      };
      auto baryY = [](const std::vector<const Clu*>& v) {
        double s = 0.0, w = 0.0;
        for (const auto* c : v) {
          s += c->cluster.getEnergy() * c->cluster.getPosition().y;
          w += c->cluster.getEnergy();
        }
        return w > 0.0 ? s / w : 0.0;
      };

      // Coordinate a view measures: SiTarget X -> x, SiTarget Y -> y, MTC
      // U/V -> the stereo coordinate s, stored in the position's x.
      std::vector<bool> used1(v1.size(), false);
      for (const auto* a : v0) {
        const double ea = a->cluster.getEnergy();
        int bestJ = -1;
        double bestAsym = m_maxPairAsym.value();
        for (std::size_t j = 0; j < v1.size(); ++j) {
          if (used1[j]) continue;
          const double eb = v1[j]->cluster.getEnergy();
          const double asym = std::abs(ea - eb) / (ea + eb);
          if (asym < bestAsym) { bestAsym = asym; bestJ = static_cast<int>(j); }
        }
        double x, y;
        std::vector<edm4hep::Cluster> mem{a->cluster};
        double e = ea, z = a->z;
        if (bestJ >= 0) {
          const auto* b = v1[bestJ];
          used1[bestJ] = true;
          mem.push_back(b->cluster);
          e += b->cluster.getEnergy();
          z = 0.5 * (a->z + b->z);
          xyOf(det, a, b, &x, &y);
        } else {
          xyOf<Clu>(det, a, nullptr, &x, &y);
          if (det == 0) y = v1.empty() ? evY : baryY(v1);
          else y = evY;
        }
        pts.push_back({x, y, li, det, e});
        members.push_back(mem);
        ptZ.push_back(z);
      }
      for (std::size_t j = 0; j < v1.size(); ++j) {
        if (used1[j]) continue;
        const auto* b = v1[j];
        double x, y;
        if (det == 0) {
          y = b->cluster.getPosition().y;
          x = v0.empty() ? evX : bary(v0);
        } else {
          x = b->cluster.getPosition().x;
          y = evY;
        }
        pts.push_back({x, y, li, det, b->cluster.getEnergy()});
        members.push_back({b->cluster});
        ptZ.push_back(b->z);
      }
    }

    // --- 2. CLUE3D ----------------------------------------------------------
    const auto res = snd::clue::run(pts, m_par);

    std::vector<std::vector<int>> byCluster(res.nClusters);
    std::vector<double> seedRho(res.nClusters, 0.0);
    for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
      if (res.cluster[i] < 0) continue;
      byCluster[res.cluster[i]].push_back(i);
      if (res.isSeed[i]) seedRho[res.cluster[i]] = res.rho[i];
    }

    for (int k = 0; k < res.nClusters; ++k) {
      const auto& idx = byCluster[k];
      if (idx.empty()) continue;
      double e = 0.0, sx = 0.0, sy = 0.0, sz = 0.0, eMtc = 0.0;
      int lo = 1 << 30, hi = -1;
      std::set<int> layers;
      for (int i : idx) {
        const double w = pts[i].energy;
        e += w; sx += w * pts[i].x; sy += w * pts[i].y; sz += w * ptZ[i];
        if (pts[i].det == 2) eMtc += w;
        lo = std::min(lo, pts[i].layer);
        hi = std::max(hi, pts[i].layer);
        layers.insert(pts[i].layer);
      }
      if (e <= 0.0) continue;
      const double cx = sx / e, cy = sy / e, cz = sz / e;

      // Axis: energy-weighted straight-line fit x(z), y(z).
      double dx = 0.0, dy = 0.0, dz = 1.0;
      if (layers.size() >= 2) {
        double szz = 0.0, szx = 0.0, szy = 0.0;
        for (int i : idx) {
          const double w = pts[i].energy, u = ptZ[i] - cz;
          szz += w * u * u;
          szx += w * u * (pts[i].x - cx);
          szy += w * u * (pts[i].y - cy);
        }
        if (szz > 0.0) { dx = szx / szz; dy = szy / szz; }
        const double nrm = std::sqrt(dx * dx + dy * dy + 1.0);
        dx /= nrm; dy /= nrm; dz = 1.0 / nrm;
      }

      auto c = out->create();
      c.setType(100 + k);
      c.setEnergy(static_cast<float>(e));
      c.setPosition({static_cast<float>(cx), static_cast<float>(cy),
                     static_cast<float>(cz)});
      c.addToShapeParameters(static_cast<float>(idx.size()));
      c.addToShapeParameters(static_cast<float>(layers.size()));
      c.addToShapeParameters(static_cast<float>(lo));
      c.addToShapeParameters(static_cast<float>(hi));
      c.addToShapeParameters(static_cast<float>(dx));
      c.addToShapeParameters(static_cast<float>(dy));
      c.addToShapeParameters(static_cast<float>(dz));
      c.addToShapeParameters(static_cast<float>(seedRho[k]));
      c.addToShapeParameters(static_cast<float>(eMtc));
      for (int i : idx)
        for (const auto& m : members[i]) c.addToClusters(m);
      ++m_nClusters;
    }

    long nOut = 0;
    for (int c : res.cluster) if (c < 0) ++nOut;
    m_nPoints += pts.size();
    m_nOutliers += nOut;
    ++m_nEvents;
    m_out->put(std::move(out));
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    if (m_nEvents > 0)
      info() << "[SNDHadronicClusterer] " << m_nEvents << " events: "
             << double(m_nPoints) / m_nEvents << " layer points and "
             << double(m_nClusters) / m_nEvents
             << " 3D clusters per event; "
             << (m_nPoints ? 100.0 * m_nOutliers / m_nPoints : 0.0)
             << " % of the points left unclustered" << endmsg;
    return Gaudi::Algorithm::finalize();
  }

private:
  // (x, y) of a view pair, or of a single view's measured coordinate.
  template <typename C>
  void xyOf(int det, const C* a, const C* b, double* x, double* y) const {
    if (det == 0) {           // SiTarget: a = X view, b = Y view
      *x = a->cluster.getPosition().x;
      *y = b ? b->cluster.getPosition().y : 0.0;
      return;
    }
    // MTC SciFi: a = U (+alpha), b = V (-alpha); s_raw = x -+ y tan(alpha)
    const double sU = a->cluster.getPosition().x;
    if (!b) { *x = sU; *y = 0.0; return; }
    const double sV = b->cluster.getPosition().x;
    *x = 0.5 * (sU + sV);
    *y = (m_tanStereo != 0.0) ? (sV - sU) / (2.0 * m_tanStereo) : 0.0;
  }

  Gaudi::Property<std::vector<std::string>> m_clusterCollections{
      this, "ClusterCollections", {}, "Layer cluster collections"};
  Gaudi::Property<std::string> m_muonPID{
      this, "MuonPID", "MuonPID", "SNDMuonFinder's ParticleIDs"};
  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "HadronicClusters", "3D clusters to write"};
  Gaudi::Property<std::string> m_compactFile{
      this, "CompactFile", "",
      "Compact for the layer sequence; only read if no geometry is loaded"};
  Gaudi::Property<double> m_mtcStereoDeg{
      this, "MTCStereoAngleDeg", 5.0, "SciFi stereo angle [deg]"};
  Gaudi::Property<bool> m_excludeCorridor{
      this, "ExcludeMuonCorridor", true,
      "Leave the event muon's corridor clusters out of CLUE3D (False when "
      "SNDPFOBuilder subtracts the muon along its trajectory)"};
  Gaudi::Property<double> m_maxPairAsym{
      this, "MaxPairAsymmetry", 0.5,
      "Two views of a layer pair up only if |E1-E2|/(E1+E2) is below this"};
  Gaudi::Property<std::vector<double>> m_dc{
      this, "CriticalDistance", {10.0, 15.0, 40.0},
      "CLUE critical distance per detector (SiTarget, SiPad, MTC) [mm]"};
  Gaudi::Property<std::vector<double>> m_rhoc{
      this, "CriticalDensity", {0.02, 0.04, 0.2},
      "CLUE seed density per detector (SiTarget, SiPad, MTC) [GeV]"};
  Gaudi::Property<int> m_densityLayers{
      this, "DensitySiblingLayers", 2, "+- layers summed into the density"};
  Gaudi::Property<int> m_nearestHigherLayers{
      this, "NearestHigherLayers", 3,
      "+- layers searched for the nearest point of higher density"};
  Gaudi::Property<double> m_outlierDeltaFactor{
      this, "OutlierDeltaFactor", 2.0,
      "A point below the seed density further than this x dc from any denser "
      "one is an outlier; also the transverse window of the nearest-higher"};

  std::vector<std::unique_ptr<DataHandle<edm4hep::ClusterCollection>>>
      m_clusterHandles;
  std::unique_ptr<DataHandle<edm4hep::ParticleIDCollection>> m_inMuonPID;
  std::unique_ptr<DataHandle<edm4hep::ClusterCollection>> m_out;

  snd::clue::Params m_par;
  double m_tanStereo = 0.0;
  std::map<std::pair<int, int>, int> m_layerIndex;

  mutable long m_nEvents = 0;
  mutable long m_nPoints = 0;
  mutable long m_nOutliers = 0;
  mutable long m_nClusters = 0;
};

DECLARE_COMPONENT(SNDHadronicClusterer)
