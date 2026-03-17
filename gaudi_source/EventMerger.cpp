#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "podio/ROOTWriter.h"
#include "podio/Frame.h"

#include <map>
#include <memory>
#include <numeric>
#include <vector>

// Plain struct to buffer contribution data between events.
struct ContribData {
  float    energy;
  float    time;
  int      pdg;
  // Exact Geant4 step position inside the sensitive volume
  float    stepX, stepY, stepZ;
  // Centre of the parent SimCalorimeterHit sensitive volume
  float    hitX, hitY, hitZ;
};

class EventMerger : public Gaudi::Algorithm {
public:
  EventMerger(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;
    m_inTargetHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inTargetName.value(), Gaudi::DataHandle::Reader, this);
    m_inPadHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inPadName.value(),  Gaudi::DataHandle::Reader, this);
    return sc;
  }

  StatusCode execute(const EventContext&) const override {
    const auto* inTarget = m_inTargetHandle->get();
    const auto* inPad  = m_inPadHandle->get();

    accumulate(*inTarget, m_bufferSiTarget);
    accumulate(*inPad,  m_bufferSiPad);

    ++m_eventCount;

    if (m_eventCount % 100 == 0) {
      info() << "[EventMerger] Accumulated evt=" << m_eventCount
             << "  SiTarget cellIDs=" << m_bufferSiTarget.size()
             << "  SiPad cellIDs="  << m_bufferSiPad.size()
             << endmsg;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    m_inTargetHandle.reset();
    m_inPadHandle.reset();

    // Write the single super-event directly via podio::ROOTWriter.
    podio::ROOTWriter writer(m_outputFile.value());
    {
      podio::Frame frame;

      auto [siTargetColl, siTargetContribs] = buildCollection(m_bufferSiTarget);
      auto [siPadColl,  siPadContribs]  = buildCollection(m_bufferSiPad);

      frame.put(std::move(siTargetColl),    m_outTargetName.value());
      frame.put(std::move(siTargetContribs), m_outTargetName.value() + "_Contributions");
      frame.put(std::move(siPadColl),     m_outPadName.value());
      frame.put(std::move(siPadContribs),  m_outPadName.value()  + "_Contributions");

      writer.writeFrame(frame, "events");
    }
    writer.finish();

    info() << "[EventMerger] Super-event written: SiTarget cellIDs=" << m_bufferSiTarget.size()
           << "  SiPad cellIDs=" << m_bufferSiPad.size() << endmsg;
    info() << "[EventMerger] Total input events merged: " << m_eventCount << endmsg;

    return Gaudi::Algorithm::finalize();
  }

private:
  static void accumulate(const edm4hep::SimCalorimeterHitCollection& input,
                         std::map<uint64_t, std::vector<ContribData>>& buffer) {
    for (const auto& hit : input) {
      auto& vec = buffer[hit.getCellID()];
      const auto& hp = hit.getPosition();
      for (const auto& contrib : hit.getContributions()) {
        const auto& sp = contrib.getStepPosition();
        vec.push_back({contrib.getEnergy(), contrib.getTime(), contrib.getPDG(),
                       sp.x, sp.y, sp.z,
                       hp.x, hp.y, hp.z});
      }
    }
  }

  static std::pair<edm4hep::SimCalorimeterHitCollection,
                   edm4hep::CaloHitContributionCollection>
  buildCollection(const std::map<uint64_t, std::vector<ContribData>>& buffer) {
    edm4hep::SimCalorimeterHitCollection  hits;
    edm4hep::CaloHitContributionCollection contribs;

    for (const auto& [cellID, contribVec] : buffer) {
      auto hit    = hits.create();
      hit.setCellID(cellID);
      // Use the hit position from the first contribution for this cellID
      hit.setPosition({contribVec[0].hitX, contribVec[0].hitY, contribVec[0].hitZ});

      float totalE = 0.0f;
      for (const auto& cd : contribVec) {
        auto contrib = contribs.create();
        contrib.setEnergy(cd.energy);
        contrib.setTime(cd.time);
        contrib.setPDG(cd.pdg);
        contrib.setStepPosition({cd.stepX, cd.stepY, cd.stepZ});
        hit.addToContributions(contrib);
        totalE += cd.energy;
      }
      hit.setEnergy(totalE);
    }
    return {std::move(hits), std::move(contribs)};
  }

  // Python-configurable properties
  Gaudi::Property<std::string> m_inTargetName{
      this, "InputCollectionSiTarget", "SiTargetHitsDelayed", "Input SiTarget collection"};
  Gaudi::Property<std::string> m_inPadName{
      this, "InputCollectionSiPad", "SiPadHitsDelayed", "Input SiPad collection"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "superevent.edm4hep.root", "Output ROOT file for the super-event"};
  Gaudi::Property<std::string> m_outTargetName{
      this, "OutputCollectionSiTarget", "SiTargetHitsMerged", "Output SiTarget collection name"};
  Gaudi::Property<std::string> m_outPadName{
      this, "OutputCollectionSiPad", "SiPadHitsMerged", "Output SiPad collection name"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inTargetHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inPadHandle;

  // Accumulation buffers: cellID -> list of contributions from all input events.
  mutable std::map<uint64_t, std::vector<ContribData>> m_bufferSiTarget;
  mutable std::map<uint64_t, std::vector<ContribData>> m_bufferSiPad;
  mutable long long m_eventCount{0};
};

DECLARE_COMPONENT(EventMerger)
