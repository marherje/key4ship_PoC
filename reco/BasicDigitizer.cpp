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
    try {
      m_inputHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_inputName.value(),  Gaudi::DataHandle::Reader, this);
      m_outputHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_outputName.value(), Gaudi::DataHandle::Writer, this);
      return Gaudi::Algorithm::initialize();
    } catch (const std::exception& e) {
      error() << "[BasicDigitizer] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[BasicDigitizer] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode execute(const EventContext&) const override {
    try {
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
    } catch (const std::exception& e) {
      error() << "[BasicDigitizer] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[BasicDigitizer] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      m_inputHandle.reset();
      m_outputHandle.reset();
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[BasicDigitizer] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[BasicDigitizer] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
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
