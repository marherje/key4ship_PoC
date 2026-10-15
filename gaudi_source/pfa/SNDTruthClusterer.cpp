#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "SNDTruthLink.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDTruthClusterer — the clustering a perfect algorithm would produce.
//
// One cluster per MC particle, holding every reconstructed hit whose energy
// that particle dominates, across all four collections. It is CyberPFA's
// TruthClusteringAlg: not a reconstruction, but the ceiling the real one is
// measured against. When the muon finder of phase 3 identifies 70 % of the
// muons, this says how much of the remaining 30 % was ever findable.
//
// The MC particle index is written into shapeParameters[0]: EDM4hep 1.x has
// no Cluster-to-MCParticle link type, and inventing a parallel collection
// would be heavier than one number a consumer already knows how to read.
//
// The truth itself lives in frame PARAMETERS (*ContribParticleIDs), and with
// k4FWCore's IOSvc the event store is not a PodioDataSvc, so there is no
// in-memory frame to ask. The algorithm therefore opens the input file again
// with its own podio::ROOTReader and walks it in step with the event loop --
// the same thing EDM4HEP2RNTuple does for the same reason. It costs a second
// read of the file, which is acceptable for a stage that exists to measure
// the reconstruction rather than to run in production.
// ---------------------------------------------------------------------------

class SNDTruthClusterer : public Gaudi::Algorithm {
public:
  SNDTruthClusterer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    if (m_inputFile.value().empty()) {
      error() << "[SNDTruthClusterer] InputFile must be set: the MC truth is "
                 "in the frame parameters of the input file." << endmsg;
      return StatusCode::FAILURE;
    }
    m_reader = std::make_unique<podio::ROOTReader>();
    m_reader->openFile(m_inputFile.value());
    m_nEntries = m_reader->getEntries("events");
    info() << "[SNDTruthClusterer] truth read from " << m_inputFile.value()
           << " (" << m_nEntries << " frames)" << endmsg;

    const auto& hitColls = m_hitCollections.value();
    const auto& simColls = m_simCollections.value();
    const auto& params = m_particleIDParams.value();
    if (hitColls.empty() || hitColls.size() != simColls.size() ||
        hitColls.size() != params.size()) {
      error() << "[SNDTruthClusterer] HitCollections, SimCollections and "
                 "ParticleIDParameters must have the same (nonzero) size."
              << endmsg;
      return StatusCode::FAILURE;
    }
    for (const auto& c : hitColls) {
      m_hitHandles.push_back(
          std::make_unique<DataHandle<edm4hep::CalorimeterHitCollection>>(
              c, Gaudi::DataHandle::Reader, this));
    }
    m_out = std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);
    return sc;
  }

  StatusCode execute(const EventContext&) const override {
    auto clusters = std::make_unique<edm4hep::ClusterCollection>();
    if (m_iEvent >= m_nEntries) {
      warning() << "[SNDTruthClusterer] the event loop went past the input "
                   "file's " << m_nEntries << " frames; no truth for this event."
                << endmsg;
      ++m_iEvent;
      m_out->put(std::move(clusters));
      return StatusCode::SUCCESS;
    }
    const podio::Frame frame(m_reader->readEntry("events", m_iEvent));
    ++m_iEvent;

    std::vector<std::pair<std::string, std::string>> spec;
    for (std::size_t i = 0; i < m_simCollections.value().size(); ++i) {
      spec.emplace_back(m_simCollections.value()[i],
                        m_particleIDParams.value()[i]);
    }
    const snd::TruthLink truth(frame, spec, m_mcCollection.value());
    if (!truth.hasTruth()) {
      ++m_nEventsNoTruth;
      m_out->put(std::move(clusters));
      return StatusCode::SUCCESS;
    }

    // Collect the hits of each MC particle, then turn each into a cluster.
    // Hits are stored BY VALUE: a podio iterator hands out a temporary
    // handle, so a stored pointer would alias whatever the next iteration
    // puts in that buffer.
    std::map<int, std::vector<edm4hep::CalorimeterHit>> byParticle;
    for (const auto& handle : m_hitHandles) {
      const auto* hits = handle->get();
      if (!hits) continue;
      for (const auto& hit : *hits) {
        const auto [idx, purity] = truth.hitTruth(hit);
        if (idx >= 0 && purity >= m_minHitPurity.value())
          byParticle[idx].push_back(hit);
      }
    }

    for (const auto& [idx, hits] : byParticle) {
      double eSum = 0.0, sx = 0.0, sy = 0.0, sz = 0.0;
      for (const auto& h : hits) {
        const double e = h.getEnergy();
        eSum += e;
        const auto& p = h.getPosition();
        sx += e * p.x; sy += e * p.y; sz += e * p.z;
      }
      const double norm = (eSum > 0.0) ? eSum : static_cast<double>(hits.size());
      auto cl = clusters->create();
      cl.setEnergy(static_cast<float>(eSum));
      cl.setPosition({static_cast<float>(sx / norm),
                      static_cast<float>(sy / norm),
                      static_cast<float>(sz / norm)});
      cl.setType(200);
      cl.addToShapeParameters(static_cast<float>(idx));
      cl.addToShapeParameters(static_cast<float>(hits.size()));
      const auto* mcs = truth.mcParticles();
      cl.addToShapeParameters(
          (mcs && idx < static_cast<int>(mcs->size()))
              ? static_cast<float>((*mcs)[idx].getPDG()) : 0.f);
      for (const auto& h : hits) cl.addToHits(h);
      ++m_nClusters;
    }

    ++m_nEvents;
    m_out->put(std::move(clusters));
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    info() << "[SNDTruthClusterer] " << m_nClusters << " truth clusters in "
           << m_nEvents << " events";
    if (m_nEventsNoTruth > 0)
      info() << " (" << m_nEventsNoTruth << " events had no usable MC truth)";
    info() << endmsg;
    return Gaudi::Algorithm::finalize();
  }

private:
  Gaudi::Property<std::vector<std::string>> m_hitCollections{
      this, "HitCollections", {}, "Reconstructed CalorimeterHit collections"};
  Gaudi::Property<std::vector<std::string>> m_simCollections{
      this, "SimCollections", {},
      "SimCalorimeterHit collection behind each of them, same order"};
  Gaudi::Property<std::vector<std::string>> m_particleIDParams{
      this, "ParticleIDParameters", {},
      "*ContribParticleIDs frame parameter of each, same order"};
  Gaudi::Property<std::string> m_mcCollection{
      this, "MCParticleCollection", "MCParticles", "MCParticle collection"};
  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "TruthClusters", "Cluster collection to write"};
  Gaudi::Property<double> m_minHitPurity{
      this, "MinHitPurity", 0.0,
      "Assign a hit to its dominant particle only above this purity; 0 means "
      "always, which keeps the truth cluster a true upper bound"};

  Gaudi::Property<std::string> m_inputFile{
      this, "InputFile", "",
      "The file IOSvc is reading. Re-opened here because the MC truth lives "
      "in frame parameters and IOSvc exposes no frame."};

  std::unique_ptr<podio::ROOTReader> m_reader;
  mutable std::size_t m_iEvent = 0;
  std::size_t m_nEntries = 0;
  std::vector<std::unique_ptr<DataHandle<edm4hep::CalorimeterHitCollection>>>
      m_hitHandles;
  std::unique_ptr<DataHandle<edm4hep::ClusterCollection>> m_out;

  mutable long m_nClusters = 0;
  mutable long m_nEvents = 0;
  mutable long m_nEventsNoTruth = 0;
};

DECLARE_COMPONENT(SNDTruthClusterer)
