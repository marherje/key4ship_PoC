#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

#include <atomic>
#include <memory>

class DelayTagger : public Gaudi::Algorithm {
public:
  DelayTagger(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    m_inTargetHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inTargetName.value(),  Gaudi::DataHandle::Reader, this);
    m_inPadHandle   = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inPadName.value(),   Gaudi::DataHandle::Reader, this);
    m_outTargetHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_outTargetName.value(), Gaudi::DataHandle::Writer, this);
    m_outPadHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_outPadName.value(),  Gaudi::DataHandle::Writer, this);
    // Contributions collections must also be registered so they are
    // persisted alongside the hit collections in the output file.
    m_outTargetContribsHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::CaloHitContributionCollection>>(
        m_outTargetName.value() + "_Contributions", Gaudi::DataHandle::Writer, this);
    m_outPadContribsHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::CaloHitContributionCollection>>(
        m_outPadName.value()  + "_Contributions", Gaudi::DataHandle::Writer, this);
    return sc;
  }

  StatusCode execute(const EventContext&) const override {
    const long long evtNum  = m_eventCount.fetch_add(1);
    const bool      doPrint = (evtNum % m_debugFreq == 0);

    const auto* inTarget = m_inTargetHandle->get();
    const auto* inPad  = m_inPadHandle->get();
    auto* outTarget        = m_outTargetHandle->createAndPut();
    auto* outPad         = m_outPadHandle->createAndPut();
    auto* outTargetContribs = m_outTargetContribsHandle->createAndPut();
    auto* outPadContribs  = m_outPadContribsHandle->createAndPut();

    tagCollection(*inTarget, *outTarget, *outTargetContribs,
                  m_inTargetName.value(), evtNum, doPrint);
    tagCollection(*inPad,  *outPad,  *outPadContribs,
                  m_inPadName.value(),  evtNum, doPrint);

    if (doPrint) {
      info() << "[DelayTagger] evt=" << evtNum
             << ": SiTarget hits=" << inTarget->size()
             << "  SiPad hits=" << inPad->size()
             << "  delay=" << static_cast<double>(evtNum) * m_eventDelay << " ns"
             << endmsg;
    }
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    m_inTargetHandle.reset();
    m_inPadHandle.reset();
    m_outTargetHandle.reset();
    m_outPadHandle.reset();
    m_outTargetContribsHandle.reset();
    m_outPadContribsHandle.reset();
    return Gaudi::Algorithm::finalize();
  }

private:
  void tagCollection(const edm4hep::SimCalorimeterHitCollection& input,
                     edm4hep::SimCalorimeterHitCollection&        output,
                     edm4hep::CaloHitContributionCollection&      outContribs,
                     const std::string& collName,
                     long long evtNum, bool doPrint) const {
    for (const auto& hit : input) {
      auto outHit = output.create();
      outHit.setCellID(hit.getCellID());
      outHit.setPosition(hit.getPosition());

      float totalE = 0.0f;
      int   j      = 0;
      for (const auto& contrib : hit.getContributions()) {
        const float t_abs = static_cast<float>(
            static_cast<double>(evtNum) * m_eventDelay + contrib.getTime());

        auto newContrib = outContribs.create();
        newContrib.setEnergy(contrib.getEnergy());
        newContrib.setTime(t_abs);
        newContrib.setPDG(contrib.getPDG());
        newContrib.setStepPosition(contrib.getStepPosition());
        outHit.addToContributions(newContrib);

        totalE += contrib.getEnergy();

        if (doPrint) {
          debug() << "[DelayTagger] evt=" << evtNum
                  << "  cellID=" << hit.getCellID()
                  << "  contrib=" << j
                  << "  t_orig=" << contrib.getTime() << " ns"
                  << "  t_abs=" << t_abs << " ns"
                  << endmsg;
        }
        ++j;
      }
      outHit.setEnergy(totalE);
    }
  }

  // Python-configurable properties
  Gaudi::Property<std::string> m_inTargetName{
      this, "InputCollectionSiTarget", "SiTargetHits", "Input SiTarget SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_inPadName{
      this, "InputCollectionSiPad", "SiPadHits", "Input SiPad SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_outTargetName{
      this, "OutputCollectionSiTarget", "SiTargetHitsDelayed", "Output SiTarget SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_outPadName{
      this, "OutputCollectionSiPad", "SiPadHitsDelayed", "Output SiPad SimCalorimeterHit collection"};
  Gaudi::Property<double> m_eventDelay{
      this, "EventDelay", 25.0, "Time delay between consecutive input events (ns)"};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print debug/info lines every N events"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>   m_inTargetHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>   m_inPadHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>   m_outTargetHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>   m_outPadHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::CaloHitContributionCollection>> m_outTargetContribsHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::CaloHitContributionCollection>> m_outPadContribsHandle;
  mutable std::atomic<long long> m_eventCount{0};
};

DECLARE_COMPONENT(DelayTagger)
