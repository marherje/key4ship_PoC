#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "TFile.h"

#include "SNDTruthLink.h"

#include <cmath>
#include <cstdint>
#include <utility>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// PFASummaryWriter — the SND-PFA output as flat RNTuples (ShipPFA.root), for
// the MOBO metrics and for analysis without podio.
//
//   EventPFA      one row per event: the truth the metrics compare against
//                 (neutrino, primary lepton, E_had,true, vertex, whether a
//                 muon candidate IS the primary muon) next to the PFASummary
//                 numbers of SNDPFOBuilder;
//   PFOs          one row per ReconstructedParticle, with its MC link;
//   LayerClusters one row per layer cluster (only if WriteClusters), with
//                 what claimed it: 0 nothing, 1 the chosen muon, 2 hadronic.
//
// Like EDM4HEP2RNTuple and SNDTruthClusterer it reads the input file with its
// own podio::ROOTReader instead of the event store: the truth lives in frame
// parameters (*ContribParticleIDs) and the PFASummary column names in the
// run metadata, and an algorithm under IOSvc can read neither. The event loop
// still runs through IOSvc, which only paces it.
//
// The truth is computed here, not in the reconstruction, with the same
// definitions as the validation scripts (Optimization/Analysis/pfa): E_had,true
// is the kinetic energy of the primary final-state hadrons (phase 1.5), the
// vertex is where the primaries start (the neutrino itself sits at z = 0),
// and a muon counts as found when a MuonCandidate's clusters are dominated by
// the primary muon (muon_perf.py).
// ---------------------------------------------------------------------------

namespace {

bool isNeutrino(int pdg) {
  pdg = std::abs(pdg);
  return pdg == 12 || pdg == 14 || pdg == 16;
}
bool isLepton(int pdg) {
  pdg = std::abs(pdg);
  return pdg == 11 || pdg == 13 || pdg == 15;
}
double norm(const edm4hep::Vector3d& v) {
  return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}
double norm(const edm4hep::Vector3f& v) {
  return std::sqrt(double(v.x) * v.x + double(v.y) * v.y + double(v.z) * v.z);
}

// Objects read back from a file have their collectionID assigned, so
// (collectionID, index) identifies them -- unlike inside the job that made
// them (see SNDPFOBuilder).
using Key = std::pair<uint32_t, int>;
Key keyOf(const podio::ObjectID& id) { return {id.collectionID, id.index}; }

template <typename T>
const T* getColl(const podio::Frame& frame, const std::string& name) {
  return dynamic_cast<const T*>(frame.get(name));
}

}  // namespace

class PFASummaryWriter : public Gaudi::Algorithm {
public:
  PFASummaryWriter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    if (m_inputFile.value().empty() || m_outputFile.value().empty()) {
      error() << "[PFASummaryWriter] InputFile and OutputFile must be set."
              << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_simCollections.value().size() != m_particleIDParams.value().size()) {
      error() << "[PFASummaryWriter] SimCollections and ParticleIDParameters "
                 "must have the same size." << endmsg;
      return StatusCode::FAILURE;
    }
    m_reader = std::make_unique<podio::ROOTReader>();
    m_reader->openFile(m_inputFile.value());
    m_nEntries = m_reader->getEntries("events");

    // The PFASummary column names, written by SNDPFOBuilder as run metadata.
    if (m_reader->getEntries("metadata") > 0) {
      const podio::Frame md(m_reader->readEntry("metadata", 0));
      auto names =
          md.getParameter<std::vector<std::string>>(m_summaryNamesKey.value());
      if (names.has_value()) m_summaryNames = names.value();
    }
    if (m_summaryNames.empty()) {
      error() << "[PFASummaryWriter] no '" << m_summaryNamesKey.value()
              << "' metadata in " << m_inputFile.value()
              << ": was it written by SNDPFOBuilder?" << endmsg;
      return StatusCode::FAILURE;
    }
    for (std::size_t i = 0; i < m_summaryNames.size(); ++i)
      m_summaryIndex[m_summaryNames[i]] = i;

    m_file.reset(TFile::Open(m_outputFile.value().c_str(), "RECREATE"));
    if (!m_file || m_file->IsZombie()) {
      error() << "[PFASummaryWriter] could not open " << m_outputFile.value()
              << endmsg;
      return StatusCode::FAILURE;
    }

    // --- EventPFA --------------------------------------------------------
    {
      auto model = ROOT::RNTupleModel::Create();
      m_ev.window_id = model->MakeField<int>("window_id");
      for (const char* n : {"flavour_true", "is_cc", "n_muon_cand",
                            "muon_matched", "p_mu_source", "flavour_reco",
                            "vtx_source"})
        m_ev.ints[n] = model->MakeField<int>(n);
      for (const char* n :
           {"E_nu_true", "E_lep_true", "E_had_true", "vtx_true_x",
            "vtx_true_y", "vtx_true_z", "p_mu_true", "E_vis_nomu",
            "E_mu_dep_true", "p_mu_reco", "q_mu_reco", "E_had_reco_ecal",
            "E_had_reco_mtc", "E_had_reco", "E_nu_reco", "vtx_reco_x",
            "vtx_reco_y", "vtx_reco_z", "E_mu_sub", "E_had_layers"})
        m_ev.floats[n] = model->MakeField<float>(n);
      m_wEvent = ROOT::RNTupleWriter::Append(std::move(model), "EventPFA",
                                             *m_file);
    }
    // --- PFOs ------------------------------------------------------------
    {
      auto model = ROOT::RNTupleModel::Create();
      m_pfo.window_id = model->MakeField<int>("window_id");
      for (const char* n : {"type", "pdg", "n_clusters", "n_tracks",
                            "mc_index", "mc_pdg"})
        m_pfo.ints[n] = model->MakeField<int>(n);
      for (const char* n : {"E", "px", "py", "pz", "mc_E", "mc_weight"})
        m_pfo.floats[n] = model->MakeField<float>(n);
      m_wPFO = ROOT::RNTupleWriter::Append(std::move(model), "PFOs", *m_file);
    }
    // --- LayerClusters ---------------------------------------------------
    if (m_writeClusters.value()) {
      auto model = ROOT::RNTupleModel::Create();
      m_clu.window_id = model->MakeField<int>("window_id");
      for (const char* n : {"det", "view", "layer", "nhits", "used_by"})
        m_clu.ints[n] = model->MakeField<int>(n);
      for (const char* n : {"E", "x", "y", "z", "width", "mips"})
        m_clu.floats[n] = model->MakeField<float>(n);
      m_wClu = ROOT::RNTupleWriter::Append(std::move(model), "LayerClusters",
                                           *m_file);
    }

    info() << "[PFASummaryWriter] " << m_inputFile.value() << " ("
           << m_nEntries << " frames, " << m_summaryNames.size()
           << " PFASummary values) -> " << m_outputFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    if (m_iEvent >= m_nEntries) {
      warning() << "[PFASummaryWriter] the event loop went past the input "
                   "file's " << m_nEntries << " frames." << endmsg;
      ++m_iEvent;
      return StatusCode::SUCCESS;
    }
    const int window = static_cast<int>(m_iEvent);
    const podio::Frame frame(m_reader->readEntry("events", m_iEvent));
    ++m_iEvent;

    std::vector<std::pair<std::string, std::string>> spec;
    for (std::size_t i = 0; i < m_simCollections.value().size(); ++i)
      spec.emplace_back(m_simCollections.value()[i],
                        m_particleIDParams.value()[i]);
    const snd::TruthLink truth(frame, spec, m_mcCollection.value());
    const auto* mcs = truth.mcParticles();

    // --- truth -----------------------------------------------------------
    int nuPdg = 0, muIdx = -1;
    double eNu = NAN, eLep = 0.0, eHadTrue = 0.0, pMuTrue = NAN;
    double vtx[3] = {NAN, NAN, NAN};
    bool isCC = false;
    if (mcs) {
      int lepIdx = -1;
      double lepE = -1.0;
      for (std::size_t i = 0; i < mcs->size(); ++i) {
        const auto p = (*mcs)[i];
        if (nuPdg == 0 && isNeutrino(p.getPDG())) {
          nuPdg = p.getPDG();
          eNu = p.getEnergy();
        }
        if (p.getGeneratorStatus() != 1) continue;
        if (std::isnan(vtx[0])) {
          const auto& v = p.getVertex();
          vtx[0] = v.x; vtx[1] = v.y; vtx[2] = v.z;
        }
        if (isLepton(p.getPDG())) {
          if (p.getEnergy() > lepE) { lepE = p.getEnergy(); lepIdx = int(i); }
          continue;
        }
        if (isNeutrino(p.getPDG())) continue;
        eHadTrue += std::max(0.0, double(p.getEnergy()) - p.getMass());
      }
      if (lepIdx >= 0) {
        isCC = true;
        eLep = lepE;
        const auto lep = (*mcs)[lepIdx];
        if (std::abs(lep.getPDG()) == 13) {
          muIdx = lepIdx;
          pMuTrue = norm(lep.getMomentum());
        }
      }
    }

    // Ceiling of phase 1.5 and the muon's own calibrated deposit: every hit,
    // split by whether the primary muon dominates it, with the e/h of
    // SNDPFOBuilder.
    double eVisNoMu = 0.0, eMuDep = 0.0;
    for (const auto& name : m_hitCollections.value()) {
      const auto* hits = getColl<edm4hep::CalorimeterHitCollection>(frame, name);
      if (!hits) continue;
      for (const auto& h : *hits) {
        const double k = (h.getType() / 10 == 2) ? m_eOverHMTC.value()
                                                 : m_eOverHECAL.value();
        const int idx = truth.hitTruth(h).first;
        ((muIdx >= 0 && idx == muIdx) ? eMuDep : eVisNoMu) +=
            k * h.getEnergy();
      }
    }

    // --- reconstruction --------------------------------------------------
    std::vector<float> summary;
    if (const auto* s = getColl<edm4hep::ParticleIDCollection>(
            frame, m_summaryCollection.value());
        s && s->size() > 0) {
      const auto pars = (*s)[0].getParameters();
      summary.assign(pars.begin(), pars.end());
    } else {
      ++m_nNoSummary;
    }
    auto val = [&](const std::string& n) -> double {
      auto it = m_summaryIndex.find(n);
      if (it == m_summaryIndex.end() || it->second >= summary.size())
        return NAN;
      return summary[it->second];
    };

    bool matched = false;
    std::set<Key> muonClusters;
    const auto* cands = getColl<edm4hep::ReconstructedParticleCollection>(
        frame, m_muonCollection.value());
    const int nCand = cands ? static_cast<int>(cands->size()) : 0;
    if (cands) {
      for (const auto& c : *cands) {
        std::map<int, double> acc;
        for (const auto& cl : c.getClusters())
          for (const auto& h : cl.getHits())
            for (const auto& [idx, e] : truth.cellTruth(h.getCellID()))
              acc[idx] += e;
        if (muIdx >= 0 && truth.dominantOf(acc).first == muIdx)
          matched = true;
      }
    }

    // The PFOs; their clusters say what claimed each layer cluster.
    std::set<Key> hadClusters;
    const auto* pfos = getColl<edm4hep::ReconstructedParticleCollection>(
        frame, m_pfoCollection.value());
    std::map<Key, std::pair<int, float>> linkOf;  // pfo -> (mc, w)
    if (const auto* links = getColl<edm4hep::RecoMCParticleLinkCollection>(
            frame, m_linkCollection.value())) {
      for (const auto& l : *links) {
        auto& best = linkOf[keyOf(l.getFrom().getObjectID())];
        if (best.second <= 0.f || l.getWeight() > best.second)
          best = {l.getTo().getObjectID().index, l.getWeight()};
      }
    }
    if (pfos) {
      for (const auto& p : *pfos) {
        const bool isMu = std::abs(p.getPDG()) == 13;
        for (const auto& cl : p.getClusters())
          (isMu ? muonClusters : hadClusters).insert(keyOf(cl.getObjectID()));
        *m_pfo.window_id = window;
        *m_pfo.ints.at("type") = isMu ? 1 : 2;
        *m_pfo.ints.at("pdg") = p.getPDG();
        *m_pfo.ints.at("n_clusters") = static_cast<int>(p.clusters_size());
        *m_pfo.ints.at("n_tracks") = static_cast<int>(p.tracks_size());
        const auto& m = p.getMomentum();
        *m_pfo.floats.at("E") = p.getEnergy();
        *m_pfo.floats.at("px") = m.x;
        *m_pfo.floats.at("py") = m.y;
        *m_pfo.floats.at("pz") = m.z;
        int mcIdx = -1, mcPdg = 0;
        float mcE = NAN, w = 0.f;
        auto it = linkOf.find(keyOf(p.getObjectID()));
        if (it != linkOf.end() && mcs && it->second.first >= 0 &&
            it->second.first < static_cast<int>(mcs->size())) {
          mcIdx = it->second.first;
          w = it->second.second;
          mcPdg = (*mcs)[mcIdx].getPDG();
          mcE = (*mcs)[mcIdx].getEnergy();
        }
        *m_pfo.ints.at("mc_index") = mcIdx;
        *m_pfo.ints.at("mc_pdg") = mcPdg;
        *m_pfo.floats.at("mc_E") = mcE;
        *m_pfo.floats.at("mc_weight") = w;
        m_wPFO->Fill();
      }
    }

    if (m_wClu) {
      for (const auto& name : m_clusterCollections.value()) {
        const auto* clus = getColl<edm4hep::ClusterCollection>(frame, name);
        if (!clus) continue;
        for (const auto& cl : *clus) {
          const auto sp = cl.getShapeParameters();
          const auto& pos = cl.getPosition();
          const auto id = keyOf(cl.getObjectID());
          *m_clu.window_id = window;
          *m_clu.ints.at("det") = cl.getType() / 10;
          *m_clu.ints.at("view") = cl.getType() % 10;
          *m_clu.ints.at("layer") = sp.size() > 0 ? int(sp[0]) : -1;
          *m_clu.ints.at("nhits") = static_cast<int>(cl.hits_size());
          *m_clu.ints.at("used_by") = muonClusters.count(id) ? 1
                                      : hadClusters.count(id) ? 2 : 0;
          *m_clu.floats.at("E") = cl.getEnergy();
          *m_clu.floats.at("x") = pos.x;
          *m_clu.floats.at("y") = pos.y;
          *m_clu.floats.at("z") = pos.z;
          *m_clu.floats.at("width") = sp.size() > 2 ? sp[2] : NAN;
          *m_clu.floats.at("mips") = sp.size() > 4 ? sp[4] : NAN;
          m_wClu->Fill();
        }
      }
    }

    // --- the event row ---------------------------------------------------
    *m_ev.window_id = window;
    *m_ev.ints.at("flavour_true") = nuPdg;
    *m_ev.ints.at("is_cc") = isCC ? 1 : 0;
    *m_ev.ints.at("n_muon_cand") = nCand;
    *m_ev.ints.at("muon_matched") = matched ? 1 : 0;
    auto asInt = [&](const std::string& n, int dflt) {
      const double v = val(n);
      return std::isnan(v) ? dflt : static_cast<int>(std::lround(v));
    };
    *m_ev.ints.at("p_mu_source") = asInt("pfa_p_mu_source", -1);
    *m_ev.ints.at("flavour_reco") = asInt("pfa_flavour", 0);
    *m_ev.ints.at("vtx_source") = asInt("pfa_vtx_source", -1);
    auto& f = m_ev.floats;
    *f.at("E_nu_true") = eNu;
    *f.at("E_lep_true") = eLep;
    *f.at("E_had_true") = eHadTrue;
    *f.at("vtx_true_x") = vtx[0];
    *f.at("vtx_true_y") = vtx[1];
    *f.at("vtx_true_z") = vtx[2];
    *f.at("p_mu_true") = pMuTrue;
    *f.at("E_vis_nomu") = eVisNoMu;
    *f.at("E_mu_dep_true") = eMuDep;
    *f.at("p_mu_reco") = val("pfa_p_mu");
    *f.at("q_mu_reco") = val("pfa_q_mu");
    *f.at("E_had_reco_ecal") = val("pfa_E_had_ecal");
    *f.at("E_had_reco_mtc") = val("pfa_E_had_mtc");
    *f.at("E_had_reco") = val("pfa_E_had");
    *f.at("E_nu_reco") = val("pfa_E_nu");
    *f.at("vtx_reco_x") = val("pfa_vtx_x");
    *f.at("vtx_reco_y") = val("pfa_vtx_y");
    *f.at("vtx_reco_z") = val("pfa_vtx_z");
    *f.at("E_mu_sub") = val("pfa_E_mu_sub");
    *f.at("E_had_layers") = val("pfa_E_had_layers");
    m_wEvent->Fill();
    ++m_nWritten;
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    // The writers must be gone before the file they append to.
    m_wEvent.reset();
    m_wPFO.reset();
    m_wClu.reset();
    if (m_file) m_file->Close();
    info() << "[PFASummaryWriter] events written=" << m_nWritten
           << " without PFASummary=" << m_nNoSummary << endmsg;
    return Gaudi::Algorithm::finalize();
  }

private:
  struct Fields {
    std::shared_ptr<int> window_id;
    std::map<std::string, std::shared_ptr<int>> ints;
    std::map<std::string, std::shared_ptr<float>> floats;
  };

  Gaudi::Property<std::string> m_inputFile{
      this, "InputFile", "", "SND-PFA output (job6_pfa.py) to read"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "ShipPFA.root", "RNTuple file to write"};
  Gaudi::Property<bool> m_writeClusters{
      this, "WriteClusters", false, "Also write the LayerClusters RNTuple"};
  Gaudi::Property<std::string> m_summaryCollection{
      this, "SummaryCollection", "PFASummary", ""};
  Gaudi::Property<std::string> m_summaryNamesKey{
      this, "SummaryNamesKey", "PFASummaryNames",
      "Run metadata entry with the PFASummary column names"};
  Gaudi::Property<std::string> m_muonCollection{
      this, "MuonCandidates", "MuonCandidates", ""};
  Gaudi::Property<std::string> m_pfoCollection{
      this, "PFOs", "ReconstructedParticles", ""};
  Gaudi::Property<std::string> m_linkCollection{
      this, "Links", "RecoMCParticleLinks", ""};
  Gaudi::Property<std::string> m_mcCollection{
      this, "MCParticles", "MCParticles", ""};
  Gaudi::Property<std::vector<std::string>> m_hitCollections{
      this, "HitCollections",
      {"SiTargetCaloHits", "SiPadCaloHits", "MTCSciFiCaloHits",
       "MTCScintCaloHits"},
      "Calibrated hits, for the truth ceiling E_vis_nomu"};
  Gaudi::Property<std::vector<std::string>> m_clusterCollections{
      this, "ClusterCollections",
      {"SiTargetLayerClusters", "SiPadLayerClusters", "MTCSciFiLayerClusters",
       "MTCScintLayerClusters"},
      "Layer clusters"};
  Gaudi::Property<std::vector<std::string>> m_simCollections{
      this, "SimCollections",
      {"SiTargetHitsWindowed", "SiPadHitsWindowed", "MTCSciFiHitsWindowed",
       "MTCScintHitsWindowed"},
      "Sim hits whose contributions carry the truth"};
  Gaudi::Property<std::vector<std::string>> m_particleIDParams{
      this, "ParticleIDParameters",
      {"SiTargetContribParticleIDs", "SiPadContribParticleIDs",
       "MTCSciFiContribParticleIDs", "MTCScintContribParticleIDs"},
      "Frame parameters with each contribution's MCParticle index"};
  Gaudi::Property<double> m_eOverHECAL{
      this, "EoverHECAL", 1.07546, "Same as SNDPFOBuilder"};
  Gaudi::Property<double> m_eOverHMTC{
      this, "EoverHMTC", 1.0, "Same as SNDPFOBuilder"};

  std::unique_ptr<podio::ROOTReader> m_reader;
  std::size_t m_nEntries = 0;
  mutable std::size_t m_iEvent = 0;
  std::vector<std::string> m_summaryNames;
  std::map<std::string, std::size_t> m_summaryIndex;

  std::unique_ptr<TFile> m_file;
  Fields m_ev, m_pfo, m_clu;
  std::unique_ptr<ROOT::RNTupleWriter> m_wEvent, m_wPFO, m_wClu;
  mutable long m_nWritten = 0;
  mutable long m_nNoSummary = 0;
};

DECLARE_COMPONENT(PFASummaryWriter)
