#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "podio/detail/OrderKey.h"
#include "k4FWCore/MetadataUtils.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/TrackCollection.h"

#include "SNDMuonChoice.h"
#include "SNDSamplingWeights.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

// ---------------------------------------------------------------------------
// SNDPFOBuilder — turns the muon and the calorimeter into an event.
//
// CyberPFA's PFOCreatingAlg reduced to what a neutrino interaction in the SND
// needs: at most one muon, and a hadronic system that is "everything else".
// Two ways of building the hadronic energy:
//   * without HadronicClusters (phase 4.1): every layer cluster that the
//     muon's corridor did not claim counts towards it;
//   * with HadronicClusters (phase 4.3, CLUE3D): one hadronic PFO per 3D
//     cluster, and E_had sums only those whose axis passes within
//     VertexLinkMM of the vertex -- a cluster pointing elsewhere is not from
//     this interaction. The 4.1 sum is still written, as pfa_E_had_layers.
//
// THE MUON'S OWN DEPOSIT (MuonSubtraction). A muon crossing the whole
// detector leaves 1-4 GeV in it once its MIPs are weighted like a shower,
// and none of that is hadronic energy. Two ways of taking it out:
//   * "corridor" (phase 4.1): drop the corridor clusters SNDMuonFinder
//     collected, keeping what each holds above one MIP. It leaves +16 % of
//     E_had in nu_mu CC (phase 4.2) for two reasons the 10 GeV muon gun
//     shows: one "MIP" is the Landau PEAK (MIPEnergyOverride), while the
//     muon deposits on average 1.35-1.8 of those per plane once delta rays
//     and radiative losses in the tungsten and iron are counted; and the
//     corridor refuses clusters above CorridorMaxMIP, which is exactly that
//     tail, so those stay whole in E_had;
//   * "trajectory" (option (a) of the 4.2 review, 2026-10-14): follow the
//     muon's fitted track through every sensitive plane from the vertex on,
//     and take from the cluster nearest to it in each plane the MEAN deposit
//     the guns measured for that kind of plane at the muon's reconstructed
//     momentum (MuonMeanMIP, tabulated at 3/10/30/100 GeV and interpolated in
//     log p, in units of the plane's MIP, scaled by the path length
//     1/cos(theta)), never more than the cluster holds. What a cluster
//     cannot give (the deposit is a Landau: the typical plane holds less than
//     the mean, a few far more) is taken from the event total instead, so
//     the subtraction is unbiased on E_had; "trajectory-capped" drops that
//     remainder (the first version, which took 85 % of the muon). Every cluster then counts towards E_had with what is
//     left, so a hadron sharing a cluster with the muon keeps its energy and
//     a pure muon cluster drops out on its own. No cut on the corridor
//     matters any more, and CLUE3D can cluster the muon's hits with the rest
//     (SNDHadronicClusterer.ExcludeMuonCorridor = False): the subtraction
//     follows each layer cluster into whichever 3D cluster took it.
//
// ENERGY SCALES. The hits arrive on the electromagnetic scale in the silicon
// (GlobalScaleECAL, fitted on electrons) and on the HADRONIC scale in the MTC
// (GlobalScaleMTC was fitted on pions fired from inside MTC40 -- phase 1.4).
// So the e/h correction belongs to the silicon only: EoverHECAL = 1.075 is
// slope(e) / slope(pi) measured on the same guns with the scales in place,
// and EoverHMTC stays 1. The plan's single EoverH = 1.4 predates that
// calibration; applied to both it would have pushed the MTC 40 % high.
//
// What this step deliberately does NOT correct: the longitudinal leakage of
// phase 1.5, which depends on where the neutrino interacted (0.95 response in
// the SiTarget, 0.64 at the back of the MTC). Step 4.2 measures it against
// that bound before anything is fitted.
//
// Frame parameters cannot be written from an algorithm under IOSvc (only the
// run metadata can), so the per-event summary travels as a one-element
// ParticleIDCollection, PFASummary, whose parameters follow the order of the
// "PFASummaryNames" metadata entry -- the reader never has to hard-code it.
// ---------------------------------------------------------------------------

namespace {

// Order of PFASummary's parameters. Written as metadata in initialize().
const std::vector<std::string> kSummaryNames = {
    "pfa_E_nu", "pfa_E_had", "pfa_E_had_ecal", "pfa_E_had_mtc",
    "pfa_p_mu", "pfa_q_mu", "pfa_p_mu_source", "pfa_flavour",
    "pfa_vtx_x", "pfa_vtx_y", "pfa_vtx_z", "pfa_vtx_source",
    "pfa_vtx_clu_x", "pfa_vtx_clu_y", "pfa_vtx_clu_z",
    "pfa_vtx_trk_x", "pfa_vtx_trk_y", "pfa_vtx_trk_z",
    "pfa_n_muon_cand", "pfa_E_corridor_extra", "pfa_n_had_clusters",
    // phase 4.3: E_had from the CLUE3D clusters linked to the vertex. The
    // layer-cluster sum of 4.1 stays alongside so the two can be compared.
    "pfa_E_had_layers", "pfa_n_3d", "pfa_n_3d_linked", "pfa_E_had_unlinked",
    "pfa_E_had_unclustered",
    // MuonSubtraction "trajectory": the muon energy taken out (e/h applied),
    // the planes its path crossed from the vertex on, and in how many of
    // those a cluster lay close enough to take it from.
    "pfa_E_mu_sub", "pfa_n_mu_planes", "pfa_n_mu_matched"};

constexpr double kMuonMass = 0.10566;  // GeV

// Vertex source codes in PFASummary: 0 = none, 1 = clusters, 2 = track.
enum VtxSource { kNoVtx = 0, kVtxClusters = 1, kVtxTrack = 2 };

struct Vtx {
  bool   ok = false;
  double x = 0.0, y = 0.0, z = 0.0;
};

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

class SNDPFOBuilder : public Gaudi::Algorithm {
public:
  SNDPFOBuilder(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    if (m_vertexMode.value() != "clusters" && m_vertexMode.value() != "track") {
      error() << "[SNDPFOBuilder] VertexMode must be 'clusters' or 'track', "
                 "not '" << m_vertexMode.value() << "'" << endmsg;
      return StatusCode::FAILURE;
    }

    if (m_muonSubtraction.value() != "corridor" &&
        m_muonSubtraction.value() != "trajectory" &&
        m_muonSubtraction.value() != "trajectory-capped") {
      error() << "[SNDPFOBuilder] MuonSubtraction must be 'corridor', "
                 "'trajectory' or 'trajectory-capped', not '" << m_muonSubtraction.value() << "'"
              << endmsg;
      return StatusCode::FAILURE;
    }
    const auto& mipP = m_muonMeanMipP.value();
    bool rising = !mipP.empty() && mipP.front() > 0.0;
    for (std::size_t i = 1; i < mipP.size(); ++i)
      rising = rising && mipP[i] > mipP[i - 1];
    if (!rising || m_muonMeanMip.value().size() != 6 * mipP.size() ||
        m_muonHalfWidth.value().size() != 4) {
      error() << "[SNDPFOBuilder] MuonMeanMIPMomenta must rise and be > 0, "
                 "MuonMeanMIP needs 6 values per momentum (SiTarget X, Y, "
                 "SiPad, MTC U, V, scintillator) and MuonTrajectoryHalfWidth "
                 "4 (SiTarget, SiPad, MTC SciFi, MTC scintillator)." << endmsg;
      return StatusCode::FAILURE;
    }
    m_tanStereo = std::tan(m_mtcStereoDeg.value() * M_PI / 180.0);

    for (const auto& c : m_clusterCollections.value())
      m_clusterHandles.push_back(
          std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
              c, Gaudi::DataHandle::Reader, this));
    for (const auto& c : m_vertexTrackCollections.value())
      m_trackHandles.push_back(
          std::make_unique<DataHandle<edm4hep::TrackCollection>>(
              c, Gaudi::DataHandle::Reader, this));
    if (m_clusterHandles.empty()) {
      error() << "[SNDPFOBuilder] ClusterCollections must be set." << endmsg;
      return StatusCode::FAILURE;
    }
    m_inMuonPID = std::make_unique<DataHandle<edm4hep::ParticleIDCollection>>(
        m_muonPID.value(), Gaudi::DataHandle::Reader, this);
    // Truth is optional: the links are for validation, not for the
    // reconstruction, and a data-like run has neither collection.
    if (!m_truthClusters.value().empty())
      m_inTruth = std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
          m_truthClusters.value(), Gaudi::DataHandle::Reader, this);
    if (!m_hadronicClusters.value().empty())
      m_inHad = std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
          m_hadronicClusters.value(), Gaudi::DataHandle::Reader, this);
    if (!m_mcParticles.value().empty())
      m_inMC = std::make_unique<DataHandle<edm4hep::MCParticleCollection>>(
          m_mcParticles.value(), Gaudi::DataHandle::Reader, this);

    m_outPFOs = std::make_unique<
        DataHandle<edm4hep::ReconstructedParticleCollection>>(
        m_outputPFOs.value(), Gaudi::DataHandle::Writer, this);
    m_outLinks = std::make_unique<
        DataHandle<edm4hep::RecoMCParticleLinkCollection>>(
        m_outputLinks.value(), Gaudi::DataHandle::Writer, this);
    m_outSummary = std::make_unique<DataHandle<edm4hep::ParticleIDCollection>>(
        m_outputSummary.value(), Gaudi::DataHandle::Writer, this);

    if (StatusCode gsc = loadLayerSequence(); gsc.isFailure()) return gsc;

    k4FWCore::putParameter("PFASummaryNames", kSummaryNames, this);

    info() << "[SNDPFOBuilder] VertexMode=" << m_vertexMode.value()
           << " EoverHECAL=" << m_eOverHECAL.value()
           << " EoverHMTC=" << m_eOverHMTC.value()
           << " VertexSeedMIP=" << m_vertexSeedMip.value()
           << " MuonSubtraction=" << m_muonSubtraction.value()
           << "; " << m_layerSeq.size() << " layers, " << m_planes.size()
           << " planes along the beam" << endmsg;
    return sc;
  }

  // Every layer of the three detectors in beam order. "Activity in the next
  // three layers" has to count empty layers too, so the sequence comes from
  // the geometry, not from the layers that happen to hold clusters.
  StatusCode loadLayerSequence() {
    auto& detector = dd4hep::Detector::getInstance();
    if (detector.state() != dd4hep::Detector::READY) {
      if (m_compactFile.value().empty()) {
        error() << "[SNDPFOBuilder] no geometry loaded and CompactFile empty."
                << endmsg;
        return StatusCode::FAILURE;
      }
      detector.fromXML(m_compactFile.value());
    }
    std::map<std::pair<int, int>, double> zOfLayer;
    std::map<std::tuple<int, int, int>, double> zOfPlane;
    for (const auto& p : snd::collectSensitivePlanes()) {
      const std::pair<int, int> key{p.det, p.layerOrdinal};
      auto it = zOfLayer.find(key);
      if (it == zOfLayer.end() || p.z < it->second) zOfLayer[key] = p.z;
      const std::tuple<int, int, int> pk{p.det, p.view, p.layerOrdinal};
      auto ip = zOfPlane.find(pk);
      if (ip == zOfPlane.end() || p.z < ip->second) zOfPlane[pk] = p.z;
    }
    m_planes.clear();
    for (const auto& [pk, z] : zOfPlane)
      m_planes.push_back({z, std::get<0>(pk), std::get<1>(pk),
                          std::get<2>(pk)});
    std::sort(m_planes.begin(), m_planes.end(),
              [](const Plane& a, const Plane& b) { return a.z < b.z; });
    std::vector<std::pair<double, std::pair<int, int>>> seq;
    for (const auto& [key, z] : zOfLayer) seq.push_back({z, key});
    std::sort(seq.begin(), seq.end());
    for (std::size_t i = 0; i < seq.size(); ++i)
      m_layerIndex[seq[i].second] = static_cast<int>(i);
    m_layerSeq = seq;
    if (m_layerSeq.empty()) {
      error() << "[SNDPFOBuilder] no sensitive planes found in the geometry."
              << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    auto pfos = std::make_unique<edm4hep::ReconstructedParticleCollection>();
    auto links = std::make_unique<edm4hep::RecoMCParticleLinkCollection>();
    auto summary = std::make_unique<edm4hep::ParticleIDCollection>();

    // --- the muon ---------------------------------------------------------
    // The candidate with the most MTC layers in its corridor. Any other
    // candidate is, by this choice, not the event's muon -- in phase 3.3 the
    // extra corridors were pions crossing the MTC -- so its clusters return
    // to the hadronic system instead of vanishing from the energy sum.
    edm4hep::ParticleID bestPID;
    int nMuonCand = 0;
    const bool haveMuon = snd::chooseMuon(m_inMuonPID->get(), &bestPID,
                                          &nMuonCand);
    std::set<Key> corridor;
    edm4hep::ReconstructedParticle muon;
    if (haveMuon) {
      muon = bestPID.getParticle();
      for (const auto& cl : muon.getClusters()) corridor.insert(keyOf(cl));
    }

    // --- the clusters -----------------------------------------------------
    struct Clu {
      edm4hep::Cluster cluster;
      int det = -1, view = -1, layer = -1;
      double z = 0.0, mips = 0.0;
      bool inCorridor = false;
    };
    std::vector<Clu> clus;
    for (const auto& handle : m_clusterHandles) {
      const auto* coll = handle->get();
      if (!coll) continue;
      for (const auto& cl : *coll) {
        const auto sp = cl.getShapeParameters();
        if (sp.size() < 5) continue;
        Clu c;
        c.cluster = cl;
        c.det = cl.getType() / 10;
        c.view = cl.getType() % 10;
        c.layer = static_cast<int>(sp[0]);
        c.z = sp[3];
        c.mips = sp[4];
        c.inCorridor = corridor.count(keyOf(cl)) > 0;
        clus.push_back(c);
      }
    }

    // --- vertex -----------------------------------------------------------
    const Vtx vClu = vertexFromClusters(clus);
    const Vtx vTrk = vertexFromTracks();
    Vtx vtx;
    int vtxSource = kNoVtx;
    if (m_vertexMode.value() == "track" && vTrk.ok) {
      vtx = vTrk; vtxSource = kVtxTrack;
    } else if (vClu.ok) {
      vtx = vClu; vtxSource = kVtxClusters;
    } else if (vTrk.ok) {
      vtx = vTrk; vtxSource = kVtxTrack;
    }

    // --- the muon's expected deposit, plane by plane ----------------------
    const bool alongPath = m_muonSubtraction.value() != "corridor";
    std::map<Key, double> muSub;   // layer cluster -> energy taken (no e/h)
    int nMuPlanes = 0, nMuMatched = 0;
    // Expected deposit no cluster could hold: (silicon, MTC), raw energy.
    std::pair<double, double> muResid{0.0, 0.0};
    if (alongPath && haveMuon)
      muSub = subtractAlongPath(muon, clus, vtx, &nMuPlanes, &nMuMatched,
                                &muResid);
    if (m_muonSubtraction.value() == "trajectory-capped") muResid = {0.0, 0.0};
    auto subOf = [&muSub](const edm4hep::Cluster& cl) {
      auto it = muSub.find(keyOf(cl));
      return it == muSub.end() ? 0.0 : it->second;
    };
    double eMuSub = m_eOverHECAL.value() * muResid.first +
                    m_eOverHMTC.value() * muResid.second;
    for (const auto& c : clus)
      eMuSub += subOf(c.cluster) *
          (c.det == 2 ? m_eOverHMTC.value() : m_eOverHECAL.value());

    // --- hadronic energy --------------------------------------------------
    // "trajectory": every cluster, less what the muon is expected to have
    // left in it. "corridor": everything the corridor did not claim, plus
    // what the corridor holds
    // above one MIP per cluster: a hadron that overlaps the muon in a layer
    // shares its cluster, and that excess is hadronic. A cluster's energy is
    // linear in its MIP count (same weight, same scale for all its hits), so
    // E x (MIP - 1) / MIP is the plan's (E_MIP - 1) x E_active x w with the
    // global scale already inside.
    double eEcal = 0.0, eMtc = 0.0, eCorridorExtra = 0.0;
    double eCorrEcal = 0.0, eCorrMtc = 0.0;
    double bx = 0.0, wx = 0.0, by = 0.0, wy = 0.0, bz = 0.0, wz = 0.0;
    std::vector<edm4hep::Cluster> hadClusters;
    for (const auto& c : clus) {
      double e = c.cluster.getEnergy();
      if (alongPath) {
        e -= subOf(c.cluster);
        if (e <= 0.0) continue;
        hadClusters.push_back(c.cluster);
      } else if (c.inCorridor) {
        if (c.mips <= 1.0) continue;
        e *= (c.mips - 1.0) / c.mips;
        eCorridorExtra += e;
        (c.det == 2 ? eCorrMtc : eCorrEcal) += e;
      } else {
        hadClusters.push_back(c.cluster);
      }
      (c.det == 2 ? eMtc : eEcal) += e;

      // Barycentre, each coordinate from the clusters that measure it: a
      // stereo SciFi strip or a SiTarget X strip knows nothing about y.
      const auto& p = c.cluster.getPosition();
      if (measuresX(c.det, c.view)) { bx += e * p.x; wx += e; }
      if (measuresY(c.det, c.view)) { by += e * p.y; wy += e; }
      bz += e * c.z; wz += e;
    }
    const double eHadLayers =
        std::max(0.0, m_eOverHECAL.value() * (eEcal - muResid.first)) +
        std::max(0.0, m_eOverHMTC.value() * (eMtc - muResid.second));
    double eHadEcal = m_eOverHECAL.value() * eEcal;
    double eHadMtc = m_eOverHMTC.value() * eMtc;

    // --- 3D clusters (phase 4.3) --------------------------------------------
    struct Had3D {
      edm4hep::Cluster cluster;
      double eEcal = 0.0, eMtc = 0.0;   // e/h applied
      double dir[3] = {0.0, 0.0, 1.0};
      bool linked = true;
    };
    std::vector<Had3D> had3d;
    double eUnlinked = 0.0, eUnclustered = 0.0;
    int nLinked = 0;
    const auto* hadColl = m_inHad ? m_inHad->get() : nullptr;
    if (hadColl) {
      std::set<Key> clustered;
      for (const auto& hc : *hadColl) {
        const auto sp = hc.getShapeParameters();
        Had3D h;
        h.cluster = hc;
        const double eMtcRaw = sp.size() > 8 ? sp[8] : 0.0;
        // The muon's share of its layer clusters, per detector.
        double subMtc = 0.0, subEcal = 0.0;
        for (const auto& cl : hc.getClusters())
          (cl.getType() / 10 == 2 ? subMtc : subEcal) += subOf(cl);
        h.eMtc = m_eOverHMTC.value() * std::max(0.0, eMtcRaw - subMtc);
        h.eEcal = m_eOverHECAL.value() *
            std::max(0.0, hc.getEnergy() - eMtcRaw - subEcal);
        if (sp.size() > 6) { h.dir[0] = sp[4]; h.dir[1] = sp[5]; h.dir[2] = sp[6]; }
        // Distance from the vertex to the cluster axis. A cluster only
        // belongs to this interaction if it points back at the vertex.
        if (vtxSource != kNoVtx) {
          const auto& p = hc.getPosition();
          const double rx = vtx.x - p.x, ry = vtx.y - p.y, rz = vtx.z - p.z;
          const double t = rx * h.dir[0] + ry * h.dir[1] + rz * h.dir[2];
          const double d = std::sqrt(std::max(
              0.0, rx * rx + ry * ry + rz * rz - t * t));
          h.linked = d < m_vertexLinkMm.value();
        }
        if (h.linked) ++nLinked;
        else eUnlinked += h.eEcal + h.eMtc;
        for (const auto& cl : hc.getClusters()) clustered.insert(keyOf(cl));
        had3d.push_back(h);
      }
      for (const auto& c : clus)
        if ((alongPath || !c.inCorridor) && !clustered.count(keyOf(c.cluster)))
          eUnclustered += std::max(0.0, c.cluster.getEnergy() -
                                            subOf(c.cluster)) *
              (c.det == 2 ? m_eOverHMTC.value() : m_eOverHECAL.value());

      // With the corridor rule its clusters never reach CLUE3D, so their
      // excess over one MIP is added back here; along the path they are
      // inside the 3D clusters already.
      eHadEcal = alongPath ? 0.0 : m_eOverHECAL.value() * eCorrEcal;
      eHadMtc = alongPath ? 0.0 : m_eOverHMTC.value() * eCorrMtc;
      for (const auto& h : had3d) {
        if (!h.linked) continue;
        eHadEcal += h.eEcal;
        eHadMtc += h.eMtc;
      }
    }
    eHadEcal = std::max(0.0, eHadEcal - m_eOverHECAL.value() * muResid.first);
    eHadMtc = std::max(0.0, eHadMtc - m_eOverHMTC.value() * muResid.second);
    const double eHad = eHadEcal + eHadMtc;

    // --- muon numbers -----------------------------------------------------
    double pMu = 0.0, qMu = 0.0, pMuSource = -1.0;
    if (haveMuon) {
      const auto& m = muon.getMomentum();
      pMu = std::sqrt(double(m.x) * m.x + double(m.y) * m.y +
                      double(m.z) * m.z);
      qMu = muon.getCharge();
      const auto params = bestPID.getParameters();
      if (params.size() > 7) pMuSource = params[7];
    }
    const double eMu = haveMuon ? std::sqrt(pMu * pMu + kMuonMass * kMuonMass)
                                : 0.0;
    const double eNu = eMu + eHad;
    const int flavour = haveMuon ? 14 : 0;

    // --- PFOs -------------------------------------------------------------
    edm4hep::MutableReconstructedParticle muPFO;
    if (haveMuon) {
      muPFO = muon.clone();
      pfos->push_back(muPFO);
    }

    // Phase 4.3: one PFO per 3D cluster. goodnessOfPID carries whether it is
    // linked to the vertex (1) or not (0) -- EDM4hep has no better field --
    // and the summary points at the most energetic linked one.
    std::vector<edm4hep::MutableReconstructedParticle> pfo3d;
    edm4hep::MutableReconstructedParticle lead3d;
    double leadE = -1.0;
    for (const auto& h : had3d) {
      auto q = pfos->create();
      const double e = h.eEcal + h.eMtc;
      q.setPDG(0);
      q.setEnergy(static_cast<float>(e));
      q.setMomentum({static_cast<float>(e * h.dir[0]),
                     static_cast<float>(e * h.dir[1]),
                     static_cast<float>(e * h.dir[2])});
      q.setReferencePoint(h.cluster.getPosition());
      q.setGoodnessOfPID(h.linked ? 1.f : 0.f);
      q.addToClusters(h.cluster);
      pfo3d.push_back(q);
      if (h.linked && e > leadE) { leadE = e; lead3d = q; }
    }

    // Without 3D clusters (phase 4.1), one PFO for the whole hadronic system.
    edm4hep::MutableReconstructedParticle had;
    if (!hadColl) {
    had = pfos->create();
    // EDM4hep 1.x dropped ReconstructedParticle::type, so PDG 0 is what tells
    // the hadronic system from the muon (+-13).
    had.setPDG(0);
    had.setEnergy(static_cast<float>(eHad));
    {
      // Direction vertex -> energy barycentre. Without a vertex, or with a
      // barycentre on top of it, the beam axis is the honest default.
      double dx = 0.0, dy = 0.0, dz = 1.0;
      if (vtxSource != kNoVtx && wz > 0.0) {
        dx = (wx > 0.0 ? bx / wx : vtx.x) - vtx.x;
        dy = (wy > 0.0 ? by / wy : vtx.y) - vtx.y;
        dz = bz / wz - vtx.z;
      }
      double n = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (n <= 0.0) { dx = 0.0; dy = 0.0; dz = 1.0; n = 1.0; }
      had.setMomentum({static_cast<float>(eHad * dx / n),
                       static_cast<float>(eHad * dy / n),
                       static_cast<float>(eHad * dz / n)});
      had.setReferencePoint({static_cast<float>(vtx.x),
                             static_cast<float>(vtx.y),
                             static_cast<float>(vtx.z)});
    }
    for (const auto& cl : hadClusters) had.addToClusters(cl);
    }

    // --- truth links --------------------------------------------------------
    const auto* mcs = m_inMC ? m_inMC->get() : nullptr;
    if (mcs && !mcs->empty()) {
      if (haveMuon && m_inTruth) {
        std::vector<edm4hep::Cluster> mcl(muon.getClusters().begin(),
                                          muon.getClusters().end());
        const auto [idx, purity] = dominantParticle(mcl);
        if (idx >= 0 && idx < static_cast<int>(mcs->size())) {
          auto l = links->create();
          l.setFrom(muPFO);
          l.setTo((*mcs)[idx]);
          l.setWeight(static_cast<float>(purity));
        }
      }
      // The hadronic system stands for the whole interaction: link it to the
      // neutrino (generator status 4, "beam"), the most energetic one if an
      // overlay window ever holds more than one.
      int nuIdx = -1;
      double nuE = -1.0;
      for (std::size_t i = 0; i < mcs->size(); ++i) {
        const auto& p = (*mcs)[i];
        if (p.getGeneratorStatus() == 4 && p.getEnergy() > nuE) {
          nuE = p.getEnergy();
          nuIdx = static_cast<int>(i);
        }
      }
      if (nuIdx >= 0 && !hadColl) {
        auto l = links->create();
        l.setFrom(had);
        l.setTo((*mcs)[nuIdx]);
        l.setWeight(1.f);
      }
      // A 3D cluster is a particle candidate in its own right: link it to the
      // MC particle that dominates its energy, with the purity as weight.
      if (m_inTruth) {
        for (std::size_t k = 0; k < had3d.size(); ++k) {
          std::vector<edm4hep::Cluster> cl(had3d[k].cluster.getClusters().begin(),
                                           had3d[k].cluster.getClusters().end());
          const auto [idx, purity] = dominantParticle(cl);
          if (idx < 0 || idx >= static_cast<int>(mcs->size())) continue;
          auto l = links->create();
          l.setFrom(pfo3d[k]);
          l.setTo((*mcs)[idx]);
          l.setWeight(static_cast<float>(purity));
        }
      }
    }

    // --- summary ------------------------------------------------------------
    auto s = summary->create();
    s.setType(0);
    s.setPDG(flavour);
    s.setAlgorithmType(0);
    const std::vector<double> values = {
        eNu, eHad, eHadEcal, eHadMtc,
        pMu, qMu, pMuSource, double(flavour),
        vtx.x, vtx.y, vtx.z, double(vtxSource),
        vClu.ok ? vClu.x : 0.0, vClu.ok ? vClu.y : 0.0, vClu.ok ? vClu.z : 0.0,
        vTrk.ok ? vTrk.x : 0.0, vTrk.ok ? vTrk.y : 0.0, vTrk.ok ? vTrk.z : 0.0,
        double(nMuonCand), eCorridorExtra, double(hadClusters.size()),
        eHadLayers, double(had3d.size()), double(nLinked), eUnlinked,
        eUnclustered, eMuSub, double(nMuPlanes), double(nMuMatched)};
    for (double v : values) s.addToParameters(static_cast<float>(v));
    if (!hadColl) s.setParticle(had);
    else if (leadE >= 0.0) s.setParticle(lead3d);

    ++m_nEvents;
    if (haveMuon) ++m_nWithMuon;
    if (vClu.ok) ++m_nVtxClu;
    if (vTrk.ok) ++m_nVtxTrk;

    m_outPFOs->put(std::move(pfos));
    m_outLinks->put(std::move(links));
    m_outSummary->put(std::move(summary));
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    info() << "[SNDPFOBuilder] " << m_nEvents << " events: " << m_nWithMuon
           << " with a muon, vertex from clusters in " << m_nVtxClu
           << ", from tracks in " << m_nVtxTrk << endmsg;
    if (m_muonSubtraction.value() == "trajectory") {
      static const char* kNames[6] = {"SiTarget X", "SiTarget Y", "SiPad",
                                      "MTC U", "MTC V", "MTC scint"};
      info() << "[SNDPFOBuilder] muon path: " << m_nMuNoPath
             << " muons without a usable track" << endmsg;
      for (int k = 0; k < 6; ++k)
        info() << "[SNDPFOBuilder]   " << kNames[k] << ": " << m_muMatched[k]
               << " / " << m_muPlanes[k] << " planes matched, mean miss "
               << (m_muMatched[k] ? m_muMiss[k] / m_muMatched[k] : 0.0)
               << " mm" << endmsg;
    }
    return Gaudi::Algorithm::finalize();
  }

private:
  // Which coordinates a layer cluster's position actually measures. The MTC
  // stereo strips report their centre at y = 0 and the SiTarget X strips sit
  // on the tile centre in y, so those values are placeholders, not data.
  static bool measuresX(int det, int view) {
    return !(det == 0 && view == 1);
  }
  static bool measuresY(int det, int view) {
    return (det == 0 && view == 1) || view == 2;
  }

  // VertexMode "clusters": the most upstream layer holding a cluster above
  // VertexSeedMIP, followed by activity in at least two of the three layers
  // after it -- a single hard scatter or a noise hit has no such follow-up,
  // an interaction does.
  template <typename CluVec>
  Vtx vertexFromClusters(const CluVec& clus) const {
    Vtx v;
    const int n = static_cast<int>(m_layerSeq.size());
    std::vector<double> maxMip(n, 0.0);
    std::vector<bool> active(n, false);
    for (const auto& c : clus) {
      auto it = m_layerIndex.find({c.det, c.layer});
      if (it == m_layerIndex.end()) continue;
      active[it->second] = true;
      maxMip[it->second] = std::max(maxMip[it->second], c.mips);
    }
    int seed = -1;
    for (int i = 0; i < n && seed < 0; ++i) {
      if (maxMip[i] < m_vertexSeedMip.value()) continue;
      int follow = 0;
      for (int j = i + 1; j <= std::min(n - 1, i + 3); ++j)
        if (active[j]) ++follow;
      if (follow >= 2) seed = i;
    }
    if (seed < 0) return v;
    const auto [det, layer] = m_layerSeq[seed].second;

    // Position in the seed layer: x and y each from the most energetic
    // cluster of a view that measures it. In the MTC, a U/V pair gives y
    // through the stereo angle, exactly as MTCSciFiMeasConverter combines
    // them: X = (s_U + s_V)/2, Y = (s_V - s_U)/(2 tan a).
    const edm4hep::Cluster* best[3] = {nullptr, nullptr, nullptr};
    double bestE[3] = {-1.0, -1.0, -1.0};
    double zSum = 0.0, zW = 0.0;
    for (const auto& c : clus) {
      if (c.det != det || c.layer != layer) continue;
      const double e = c.cluster.getEnergy();
      if (c.view >= 0 && c.view < 3 && e > bestE[c.view]) {
        bestE[c.view] = e;
        best[c.view] = &c.cluster;
      }
      zSum += c.z; zW += 1.0;
    }
    v.z = zW > 0.0 ? zSum / zW : m_layerSeq[seed].first;
    if (det == 0) {
      if (best[0]) v.x = best[0]->getPosition().x;
      if (best[1]) v.y = best[1]->getPosition().y;
    } else if (det == 1) {
      if (best[2]) { v.x = best[2]->getPosition().x; v.y = best[2]->getPosition().y; }
    } else {
      if (best[0] && best[1]) {
        const double sU = best[0]->getPosition().x;
        const double sV = best[1]->getPosition().x;
        const double tanA = std::tan(m_mtcStereoDeg.value() * M_PI / 180.0);
        v.x = 0.5 * (sU + sV);
        v.y = (tanA != 0.0) ? (sV - sU) / (2.0 * tanA) : 0.0;
      } else if (best[0] || best[1]) {
        v.x = (best[0] ? best[0] : best[1])->getPosition().x;
      }
      // The scintillator measures y directly, if coarsely (50 mm cells);
      // prefer it to a stereo y from a single unpaired view.
      if (best[2] && !(best[0] && best[1])) {
        v.y = best[2]->getPosition().y;
        if (!best[0] && !best[1]) v.x = best[2]->getPosition().x;
      }
    }
    v.ok = true;
    return v;
  }

  // MuonSubtraction "trajectory". The path is the muon track's fitted states
  // (one per surface it crossed, DD4hep axes), interpolated between them and
  // extrapolated straight along the end slopes outside them -- the same model
  // SNDMuonFinder builds its corridor on. It starts at the vertex: upstream
  // of it there was no muon. In every sensitive plane after that, the
  // cluster nearest to the path within MuonTrajectoryHalfWidth gives up the
  // gun's mean deposit for that kind of plane, in that plane's own MIPs
  // (E / MIPs of the cluster is the plane's calibrated energy per MIP), times
  // the path length through it. Returns the energy taken per layer cluster.
  template <typename CluVec>
  std::map<Key, double> subtractAlongPath(
      const edm4hep::ReconstructedParticle& muon, const CluVec& clus,
      const Vtx& vtx, int* nPlanes, int* nMatched,
      std::pair<double, double>* resid) const {
    std::map<Key, double> taken;
    std::vector<std::array<double, 3>> path;   // (z, x, y)
    for (const auto& track : muon.getTracks())
      for (const auto& st : track.getTrackStates()) {
        if (st.location != edm4hep::TrackState::AtOther) continue;
        const auto& rp = st.referencePoint;
        path.push_back({rp.z, rp.x, rp.y});
      }
    if (path.size() < 2) {
      ++m_nMuNoPath;
      return taken;
    }
    std::sort(path.begin(), path.end());

    // Position at z from the segment that contains it (or the end segment
    // beyond the states). The slope, which sets the path length, comes from
    // states at least MuonSlopeBaseline apart instead: the two views of a
    // layer sit a few mm apart in z, and the fit's kink between them would
    // otherwise read as a steep track.
    const double half = 0.5 * m_muonSlopeBaseline.value();
    auto at = [&path, half](double z, double* x, double* y, double* sx,
                            double* sy) {
      std::size_t i = 1;
      while (i + 1 < path.size() && z > path[i][0]) ++i;
      const auto& a = path[i - 1];
      const auto& b = path[i];
      const double dz = b[0] - a[0];
      const double t = dz != 0.0 ? (z - a[0]) / dz : 0.0;
      *x = a[1] + t * (b[1] - a[1]);
      *y = a[2] + t * (b[2] - a[2]);
      std::size_t lo = 0, hi = path.size() - 1;
      for (std::size_t k = 0; k < path.size(); ++k) {
        if (path[k][0] <= z - half) lo = k;
        if (path[k][0] >= z + half) { hi = k; break; }
      }
      const double bz = path[hi][0] - path[lo][0];
      *sx = bz > 0.0 ? (path[hi][1] - path[lo][1]) / bz : 0.0;
      *sy = bz > 0.0 ? (path[hi][2] - path[lo][2]) / bz : 0.0;
    };

    std::map<std::tuple<int, int, int>, std::vector<int>> onPlane;
    for (std::size_t i = 0; i < clus.size(); ++i)
      onPlane[{clus[i].det, clus[i].view, clus[i].layer}].push_back(
          static_cast<int>(i));

    const auto& mom = muon.getMomentum();
    const double pMu = std::sqrt(double(mom.x) * mom.x +
                                 double(mom.y) * mom.y +
                                 double(mom.z) * mom.z);
    std::array<double, 6> meanMip;
    for (int k = 0; k < 6; ++k) meanMip[k] = muonMeanMip(k, pMu);

    const double zStart = vtx.ok ? vtx.z - m_muonStartTolMm.value()
                                 : path.front()[0];
    for (const auto& pl : m_planes) {
      if (pl.z < zStart) continue;
      ++*nPlanes;
      const int kind = pl.det == 0 ? pl.view
                     : pl.det == 1 ? 2 : 3 + pl.view;   // MuonMeanMIP index
      const int wKind = pl.det < 2 ? pl.det : (pl.view == 2 ? 3 : 2);
      ++m_muPlanes[kind];
      auto it = onPlane.find({pl.det, pl.view, pl.layer});
      if (it == onPlane.end()) continue;

      double x, y, sx, sy;
      at(pl.z, &x, &y, &sx, &sy);
      // Only the coordinate the plane measures; the MTC stereo views read
      // s = x -+ y tan(a) (U, V), checked on truth muons on 2026-10-14.
      auto miss = [&](const edm4hep::Cluster& cl) {
        const auto& p = cl.getPosition();
        if (pl.det == 0) return std::abs(pl.view == 0 ? p.x - x : p.y - y);
        if (pl.view == 2) return std::hypot(p.x - x, p.y - y);
        const double s = x + (pl.view == 0 ? -1.0 : 1.0) * y * m_tanStereo;
        return std::abs(p.x - s);
      };
      int best = -1;
      double bestMiss = m_muonHalfWidth.value()[wKind];
      for (int i : it->second) {
        const double d = miss(clus[i].cluster);
        if (d <= bestMiss) { bestMiss = d; best = i; }
      }
      if (best < 0) continue;
      const auto& c = clus[best];
      if (c.mips <= 0.0) continue;
      const double eCl = c.cluster.getEnergy();
      const double expect = meanMip[kind] * eCl / c.mips *
                            std::sqrt(1.0 + sx * sx + sy * sy);
      double& t = taken[keyOf(c.cluster)];
      const double sub = std::max(0.0, std::min(expect, eCl - t));
      // What does not fit in the cluster is still the muon's on average: the
      // deposit is a Landau, so the typical plane holds less than the mean
      // and a few hold far more. Capping plane by plane takes only 56-83 %
      // of what the guns deposit; the remainder is taken from the event.
      (pl.det == 2 ? resid->second : resid->first) += expect - sub;
      ++*nMatched;
      ++m_muMatched[kind];
      m_muMiss[kind] += bestMiss;
      if (sub <= 0.0) continue;
      t += sub;
    }
    return taken;
  }

  // MuonMeanMIP for plane kind k at momentum p: linear in log p between the
  // tabulated momenta, flat beyond them (and at the first entry when the
  // muon has no momentum).
  double muonMeanMip(int k, double p) const {
    const auto& ps = m_muonMeanMipP.value();
    const auto& v = m_muonMeanMip.value();
    const std::size_t n = ps.size();
    if (n == 1 || p <= ps.front()) return v[k];
    if (p >= ps.back()) return v[6 * (n - 1) + k];
    std::size_t i = 1;
    while (ps[i] < p) ++i;
    const double t = std::log(p / ps[i - 1]) / std::log(ps[i] / ps[i - 1]);
    return v[6 * (i - 1) + k] + t * (v[6 * i + k] - v[6 * (i - 1) + k]);
  }

  // VertexMode "track": the most upstream fitted state of any track that
  // starts in the SiTarget (type bit 1).
  Vtx vertexFromTracks() const {
    Vtx v;
    double zMin = std::numeric_limits<double>::max();
    for (const auto& handle : m_trackHandles) {
      const auto* coll = handle->get();
      if (!coll) continue;
      for (const auto& t : *coll) {
        if (!(t.getType() & 1)) continue;
        for (const auto& st : t.getTrackStates()) {
          if (st.location != edm4hep::TrackState::AtOther) continue;
          const auto& rp = st.referencePoint;
          if (rp.z < zMin) {
            zMin = rp.z;
            v.x = rp.x; v.y = rp.y; v.z = rp.z;
            v.ok = true;
          }
        }
      }
    }
    return v;
  }

  // MC particle that deposited most of a set of layer clusters' energy, and
  // its share: the truth clusters already hold, for every reco hit, the particle
  // that dominates it (SNDTruthClusterer, shapeParameters[0]).
  std::pair<int, double> dominantParticle(
      const std::vector<edm4hep::Cluster>& clusters) const {
    const auto* truth = m_inTruth->get();
    if (!truth) return {-1, 0.0};
    std::map<Key, int> owner;
    for (const auto& tc : *truth) {
      const auto sp = tc.getShapeParameters();
      if (sp.empty()) continue;
      const int idx = static_cast<int>(sp[0]);
      for (const auto& h : tc.getHits()) owner[keyOf(h)] = idx;
    }
    std::map<int, double> eBy;
    double eTot = 0.0;
    for (const auto& cl : clusters) {
      for (const auto& h : cl.getHits()) {
        const double e = h.getEnergy();
        eTot += e;
        auto it = owner.find(keyOf(h));
        if (it != owner.end()) eBy[it->second] += e;
      }
    }
    int best = -1;
    double bestE = 0.0;
    for (const auto& [idx, e] : eBy)
      if (e > bestE) { bestE = e; best = idx; }
    return {best, eTot > 0.0 ? bestE / eTot : 0.0};
  }

  Gaudi::Property<std::vector<std::string>> m_clusterCollections{
      this, "ClusterCollections", {}, "Layer cluster collections"};
  Gaudi::Property<std::vector<std::string>> m_vertexTrackCollections{
      this, "VertexTrackCollections", {"SiTargetTracks", "GlobalTracks"},
      "Tracks searched for the most upstream SiTarget state (VertexMode "
      "'track'); spliced tracks live in GlobalTracks, so both are needed"};
  Gaudi::Property<std::string> m_muonPID{
      this, "MuonPID", "MuonPID",
      "ParticleIDs written by SNDMuonFinder, one per muon candidate"};
  Gaudi::Property<std::string> m_truthClusters{
      this, "TruthClusters", "TruthClusters",
      "Truth clusters for the muon's MC link; empty disables it"};
  Gaudi::Property<std::string> m_hadronicClusters{
      this, "HadronicClusters", "",
      "3D clusters from SNDHadronicClusterer; empty keeps the phase-4.1 sum"};
  Gaudi::Property<double> m_vertexLinkMm{
      this, "VertexLinkMM", 30.0,
      "A 3D cluster counts towards E_had if its axis passes this close to "
      "the vertex [mm]"};
  Gaudi::Property<std::string> m_mcParticles{
      this, "MCParticles", "MCParticles",
      "MC particles for the links; empty disables them"};
  Gaudi::Property<std::string> m_outputPFOs{
      this, "OutputPFOs", "ReconstructedParticles", "PFOs to write"};
  Gaudi::Property<std::string> m_outputLinks{
      this, "OutputLinks", "RecoMCParticleLinks", "PFO -> MC links"};
  Gaudi::Property<std::string> m_outputSummary{
      this, "OutputSummary", "PFASummary",
      "One ParticleID per event carrying the pfa_* numbers (see the "
      "PFASummaryNames metadata for their order)"};
  Gaudi::Property<std::string> m_compactFile{
      this, "CompactFile", "",
      "Compact for the layer sequence; only read if no geometry is loaded"};
  Gaudi::Property<std::string> m_vertexMode{
      this, "VertexMode", "clusters",
      "'clusters' (first layer with a seed and follow-up activity) or "
      "'track' (first SiTarget track state, falling back to clusters)"};
  Gaudi::Property<double> m_vertexSeedMip{
      this, "VertexSeedMIP", 2.0,
      "MIPs a cluster needs to seed the vertex layer"};
  Gaudi::Property<double> m_eOverHECAL{
      this, "EoverHECAL", 1.075,
      "e/h of SiTarget+SiPad: slope(e)/slope(pi) on the phase-1.4 guns with "
      "the current scales"};
  Gaudi::Property<double> m_eOverHMTC{
      this, "EoverHMTC", 1.0,
      "e/h of the MTC; 1 because GlobalScaleMTC was fitted on pions"};
  Gaudi::Property<double> m_mtcStereoDeg{
      this, "MTCStereoAngleDeg", 5.0,
      "SciFi stereo angle, for the y of an MTC vertex from a U/V pair"};
  Gaudi::Property<std::string> m_muonSubtraction{
      this, "MuonSubtraction", "trajectory",
      "'trajectory' (the guns' mean muon deposit, plane by plane along the "
      "track, what a cluster cannot hold taken from the event total), "
      "'trajectory-capped' (the same without that remainder) or 'corridor' "
      "(phase 4.1: drop corridor clusters, keep what they hold above one "
      "MIP)"};
  // Muon guns (pfa_calib_mu-_<p>GeV.root), 500 events each: mean of the MIPs
  // per crossed plane, where 1 MIP is the Landau peak of MIPEnergyOverride.
  // The medians are 1.08-1.16; the means sit well above because some planes
  // catch a delta ray or a radiative loss, and that tail grows with p (the
  // fraction of planes above 3 MIP goes from 6 % at 3 GeV to 15 % at 100).
  Gaudi::Property<std::vector<double>> m_muonMeanMipP{
      this, "MuonMeanMIPMomenta", {3.0, 10.0, 30.0, 100.0},
      "Muon momenta MuonMeanMIP is tabulated at [GeV], rising"};
  Gaudi::Property<std::vector<double>> m_muonMeanMip{
      this, "MuonMeanMIP",
      {1.529, 1.559, 1.397, 1.526, 1.563, 1.261,    //   3 GeV
       1.782, 1.810, 1.548, 1.665, 1.723, 1.350,    //  10 GeV
       1.956, 1.998, 1.916, 1.882, 1.929, 1.508,    //  30 GeV
       2.875, 2.919, 2.738, 2.490, 2.571, 1.928},   // 100 GeV
      "Mean MIPs a muon leaves per plane, 6 values per entry of "
      "MuonMeanMIPMomenta (SiTarget X, Y, SiPad, MTC U, V, MTC scintillator), "
      "interpolated linearly in log p and held flat beyond the table"};
  Gaudi::Property<std::vector<double>> m_muonHalfWidth{
      this, "MuonTrajectoryHalfWidth", {5.0, 15.0, 20.0, 50.0},
      "How far from the muon path a cluster may sit and still be the one "
      "that holds its deposit: SiTarget, SiPad, MTC SciFi, MTC scintillator "
      "[mm]"};
  Gaudi::Property<double> m_muonSlopeBaseline{
      this, "MuonSlopeBaselineMM", 50.0,
      "z span of the track states the muon's local slope is taken over [mm]"};
  Gaudi::Property<double> m_muonStartTolMm{
      this, "MuonStartToleranceMM", 5.0,
      "Planes this far upstream of the vertex still count as crossed by the "
      "muon [mm]"};

  std::vector<std::unique_ptr<DataHandle<edm4hep::ClusterCollection>>>
      m_clusterHandles;
  std::vector<std::unique_ptr<DataHandle<edm4hep::TrackCollection>>>
      m_trackHandles;
  std::unique_ptr<DataHandle<edm4hep::ParticleIDCollection>> m_inMuonPID;
  std::unique_ptr<DataHandle<edm4hep::ClusterCollection>> m_inTruth;
  std::unique_ptr<DataHandle<edm4hep::MCParticleCollection>> m_inMC;
  std::unique_ptr<DataHandle<edm4hep::ClusterCollection>> m_inHad;
  std::unique_ptr<DataHandle<edm4hep::ReconstructedParticleCollection>>
      m_outPFOs;
  std::unique_ptr<DataHandle<edm4hep::RecoMCParticleLinkCollection>>
      m_outLinks;
  std::unique_ptr<DataHandle<edm4hep::ParticleIDCollection>> m_outSummary;

  struct Plane {
    double z;
    int det, view, layer;
  };
  std::vector<Plane> m_planes;   // every sensitive plane, in beam order
  double m_tanStereo = 0.0;

  // Layers in beam order: (z of the first plane, (det, ordinal)).
  std::vector<std::pair<double, std::pair<int, int>>> m_layerSeq;
  std::map<std::pair<int, int>, int> m_layerIndex;

  mutable long m_nEvents = 0;
  mutable long m_nWithMuon = 0;
  mutable long m_nVtxClu = 0;
  mutable long m_nVtxTrk = 0;
  mutable long m_nMuNoPath = 0;
  mutable std::array<long, 6> m_muPlanes{};
  mutable std::array<long, 6> m_muMatched{};
  mutable std::array<double, 6> m_muMiss{};
};

DECLARE_COMPONENT(SNDPFOBuilder)
