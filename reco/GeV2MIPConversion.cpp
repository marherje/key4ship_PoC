#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <memory>

class GeV2MIPConversion : public Gaudi::Algorithm {
public:
  GeV2MIPConversion(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    // DataHandles are built here, AFTER Gaudi has applied all Python-set
    // property values, so m_inputName / m_outputName hold the user values.
    m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inputName.value(),  Gaudi::DataHandle::Reader, this);
    m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_outputName.value(), Gaudi::DataHandle::Writer, this);
    return Gaudi::Algorithm::initialize();
  }

  StatusCode execute(const EventContext&) const override {
    const auto* input  = m_inputHandle->get();
    auto*       output = m_outputHandle->createAndPut();

    for (const auto& hit : *input) {
      auto nh = output->create();
      nh.setCellID(hit.getCellID());
      nh.setEnergy(static_cast<float>(hit.getEnergy() * (1/m_MIPValue)));
      nh.setPosition(hit.getPosition());
    }

    info() << "GeV2MIPConversion: " << input->size() << " hits processed" << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    m_inputHandle.reset();
    m_outputHandle.reset();
    return Gaudi::Algorithm::finalize();
  }

private:
  // Python-configurable properties
  Gaudi::Property<std::string> m_inputName{
      this, "InputCollection", "SiTargetHits", "Input SimCalorimeterHit collection"};
  Gaudi::Property<std::string> m_outputName{
      this, "OutputCollection", "SiTargetHitsMIP", "Output SimCalorimeterHit collection"};
  Gaudi::Property<double> m_MIPValue{
      this, "MIPValue", 0.000009, "MIP value for conversion (in GeV)"}; // Default is 9 keV, for 0.3mm silicon

  // DataHandles constructed in initialize() from the properties above
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inputHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_outputHandle;
};

DECLARE_COMPONENT(GeV2MIPConversion)
