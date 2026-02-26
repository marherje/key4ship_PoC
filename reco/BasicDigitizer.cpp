#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <memory>

class BasicDigitizer : public Gaudi::Algorithm {
public:
  BasicDigitizer(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inputName.value(),  Gaudi::DataHandle::Reader, this);
    m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_outputName.value(), Gaudi::DataHandle::Writer, this);
    return Gaudi::Algorithm::initialize();
  }

  StatusCode execute(const EventContext&) const override {
    const auto* input  = m_inputHandle->get();
    auto*       output = m_outputHandle->createAndPut();

    int n_pass = 0;
    for (const auto& hit : *input) {
      if (hit.getEnergy() > m_threshold) {
        auto nh = output->create();
        nh.setCellID(hit.getCellID());
        nh.setEnergy(hit.getEnergy());
        nh.setPosition(hit.getPosition());
        ++n_pass;
      }
    }

    info() << "BasicDigitizer: " << input->size() << " in, "
           << n_pass << " passing threshold" << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    m_inputHandle.reset();
    m_outputHandle.reset();
    return Gaudi::Algorithm::finalize();
  }

private:
  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiTargetHitsMIP", "Input SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "SiTargetDigiHits", "Output SimCalorimeterHit collection"};
  Gaudi::Property<double> m_threshold{
      this, "Threshold", 0.5, "Minimum energy to keep a hit"};

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
};

DECLARE_COMPONENT(BasicDigitizer)
