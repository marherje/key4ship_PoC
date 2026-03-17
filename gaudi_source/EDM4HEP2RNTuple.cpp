#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "DD4hep/BitFieldCoder.h"
#include "podio/ROOTReader.h"
#include "podio/Frame.h"
#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "TFile.h"
#include "ValidationUtils.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class EDM4HEP2RNTuple : public Gaudi::Algorithm {
public:
  EDM4HEP2RNTuple(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;

      const auto& colls       = m_collections.value();
      const auto& bfs         = m_bitFields.value();
      const auto& params      = m_sourceIDParams.value();
      const auto& detIDs      = m_detectorIDs.value();
      const auto& ntupleNames = m_ntupleNames.value();

      // 1. Validate parallel vectors
      if (colls.size() != bfs.size()         ||
          colls.size() != params.size()      ||
          colls.size() != detIDs.size()      ||
          colls.size() != ntupleNames.size()) {
        error() << "[EDM4HEP2RNTuple] Collections, BitFields, SourceIDParams,"
                << " DetectorIDs and NTupleNames must all have the same size."
                << " Got: "
                << colls.size()       << " " << bfs.size()   << " "
                << params.size()      << " " << detIDs.size() << " "
                << ntupleNames.size()
                << endmsg;
        return StatusCode::FAILURE;
      }

      // 2. Build one DataHandle per collection
      for (const auto& cname : colls) {
        m_inputHandles.push_back(
            std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
                cname, Gaudi::DataHandle::Reader, this));
      }

      // 3. Build one BitFieldCoder per collection
      try {
        for (const auto& bf : bfs) {
          m_decoders.push_back(
              std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(bf));
        }
      } catch (const std::exception& e) {
        error() << "[EDM4HEP2RNTuple] Invalid BitField string: " << e.what() << endmsg;
        return StatusCode::FAILURE;
      }

      // 4. Inspect BitField field names per detector
      for (size_t i = 0; i < m_decoders.size(); ++i) {
        const auto& decoder = *m_decoders[i];
        std::vector<std::string> names;
        for (size_t f = 0; f < decoder.size(); ++f) {
          names.push_back(decoder[f].name());
        }
        m_detectorFieldNames.push_back(names);

        info() << "[EDM4HEP2RNTuple] detector_id=" << detIDs[i]
               << " ntuple='" << ntupleNames[i] << "' BitField fields: ";
        for (const auto& n : names) info() << n << " ";
        info() << endmsg;
      }

      // 5. Open podio::ROOTReader for frame parameter access
      m_reader = std::make_unique<podio::ROOTReader>();
      m_reader->openFile(m_inputFile.value());

      if (m_reader->getEntries("events") == 0) {
        error() << "[EDM4HEP2RNTuple] Input file has 0 events." << endmsg;
        return StatusCode::FAILURE;
      }

      // Validate first frame
      {
        auto firstFrame = podio::Frame(m_reader->readEntry("events", 0));
        for (const auto& cname : colls)
          validateCollection(firstFrame, cname, "EDM4HEP2RNTuple", msgStream());
        for (const auto& pname : params) {
          if (!firstFrame.getParameter<std::vector<int>>(pname).has_value()) {
            warning() << "[EDM4HEP2RNTuple] Frame parameter '" << pname
                      << "' not found in first window. source_id will be 0." << endmsg;
          }
        }
      }

      // 6. Open output TFile once; all per-detector writers share it
      m_outputTFile.reset(
          TFile::Open(m_outputFile.value().c_str(), "RECREATE"));
      if (!m_outputTFile || m_outputTFile->IsZombie()) {
        error() << "[EDM4HEP2RNTuple] Could not open output file: "
                << m_outputFile.value() << endmsg;
        return StatusCode::FAILURE;
      }

      // 7. Build one RNTuple model and writer per detector
      for (size_t i = 0; i < colls.size(); ++i) {
        auto model = ROOT::RNTupleModel::Create();
        DetectorFields df;

        // Fixed fields (no detector_id: each RNTuple is already detector-specific)
        auto fWindowID     = model->MakeField<int>("window_id");
        auto fCellID       = model->MakeField<uint64_t>("cellID");
        auto fX            = model->MakeField<float>("x");
        auto fY            = model->MakeField<float>("y");
        auto fZ            = model->MakeField<float>("z");
        auto fE            = model->MakeField<float>("E");
        auto fT            = model->MakeField<float>("t");
        auto fSourceID     = model->MakeField<int>("source_id");
        auto fTWindowStart = model->MakeField<float>("t_window_start");

        df.fWindowID     = fWindowID.get();
        df.fCellID       = fCellID.get();
        df.fX            = fX.get();
        df.fY            = fY.get();
        df.fZ            = fZ.get();
        df.fE            = fE.get();
        df.fT            = fT.get();
        df.fSourceID     = fSourceID.get();
        df.fTWindowStart = fTWindowStart.get();

        // Dynamic BitField columns for THIS detector only — no cross-detector padding
        info() << "[EDM4HEP2RNTuple] NTuple '" << ntupleNames[i]
               << "' dynamic BitField columns: ";
        for (const auto& fname : m_detectorFieldNames[i]) {
          const std::string colName = "bf_" + fname;
          auto fptr = model->MakeField<int>(colName);
          df.fieldPtrs[fname] = fptr.get();
          df.fieldSharedPtrs.push_back(std::move(fptr));
          info() << colName << " ";
        }
        info() << endmsg;

        // File was opened empty with "RECREATE"; Append works for all detectors.
        m_writers.push_back(
            ROOT::RNTupleWriter::Append(
                std::move(model), ntupleNames[i], *m_outputTFile));

        m_detFields.push_back(std::move(df));

        info() << "[EDM4HEP2RNTuple]   detector_id=" << detIDs[i]
               << "  collection=" << colls[i]
               << "  ntuple=" << ntupleNames[i]
               << "  bitfield=" << bfs[i]
               << endmsg;
      }

      info() << "[EDM4HEP2RNTuple] Initialized."
             << "  Output: " << m_outputFile.value()
             << "  Detectors: " << colls.size() << endmsg;

      return sc;
    } catch (const std::exception& e) {
      error() << "[EDM4HEP2RNTuple] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EDM4HEP2RNTuple] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Called once per time window by Gaudi (EvtMax=-1, IOSvc reads the file).
  StatusCode execute(const EventContext&) const override {
    try {
      const long long windowID = m_windowCount.fetch_add(1);
      const auto& colls  = m_collections.value();
      const auto& params = m_sourceIDParams.value();

      // Read the frame once per window to get all parameters
      std::optional<podio::Frame> maybeFrame;
      float t_window_start = 0.f;
      try {
        auto entry = m_reader->readEntry("events", static_cast<size_t>(windowID));
        maybeFrame.emplace(std::move(entry));
        auto ts = maybeFrame->getParameter<float>("t_window_start");
        if (ts.has_value()) t_window_start = ts.value();
      } catch (...) {
        warning() << "[EDM4HEP2RNTuple] Could not read frame parameters for"
                  << " window " << windowID << ". source_id and t_window_start will be 0."
                  << endmsg;
      }

      for (size_t i = 0; i < colls.size(); ++i) {
        const auto* hits = m_inputHandles[i]->get();
        auto& df = m_detFields[i];

        std::vector<int> sourceIDs;
        if (maybeFrame) {
          sourceIDs = maybeFrame->getParameter<std::vector<int>>(params[i])
                          .value_or(std::vector<int>{});
        }

        for (int hitIdx = 0; hitIdx < static_cast<int>(hits->size()); ++hitIdx) {
          const auto& hit       = (*hits)[hitIdx];
          const uint64_t cellID = hit.getCellID();
          const auto& pos       = hit.getPosition();

          float t = 0.f;
          if (hit.contributions_size() > 0)
            t = hit.getContributions(0).getTime();

          const int source_id = (hitIdx < static_cast<int>(sourceIDs.size()))
                                ? sourceIDs[hitIdx] : 0;

          // Decode BitField fields for this detector only (no zero-padding needed)
          for (auto& [fname, ptr] : df.fieldPtrs) {
            *ptr = static_cast<int>(m_decoders[i]->get(cellID, fname));
          }

          *df.fWindowID     = static_cast<int>(windowID);
          *df.fCellID       = cellID;
          *df.fX            = pos.x;
          *df.fY            = pos.y;
          *df.fZ            = pos.z;
          *df.fE            = hit.getEnergy();
          *df.fT            = t;
          *df.fSourceID     = source_id;
          *df.fTWindowStart = t_window_start;

          m_writers[i]->Fill();
        }
      }

      if (windowID % 100 == 0) {
        info() << "[EDM4HEP2RNTuple] Processed window " << windowID
               << "  t_window_start=" << t_window_start << " ns" << endmsg;
        for (size_t i = 0; i < colls.size(); ++i) {
          info() << "[EDM4HEP2RNTuple]   collection=" << colls[i]
                 << "  hits=" << m_inputHandles[i]->get()->size() << endmsg;
        }
      }

      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[EDM4HEP2RNTuple] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EDM4HEP2RNTuple] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      for (auto& w : m_writers) w.reset();  // flush all writers before closing file
      m_writers.clear();
      m_detFields.clear();
      m_outputTFile.reset();  // closes the TFile after all writers are flushed
      m_reader.reset();
      for (auto& h : m_inputHandles) h.reset();

      info() << "[EDM4HEP2RNTuple] Done. Total windows processed: "
             << m_windowCount.load() << endmsg;
      info() << "[EDM4HEP2RNTuple] RNTuples written to: " << m_outputFile.value() << endmsg;

      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[EDM4HEP2RNTuple] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EDM4HEP2RNTuple] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  // Python-configurable properties — parallel vectors, one entry per detector.
  Gaudi::Property<std::vector<std::string>> m_collections{
      this, "Collections",
      {"SiTargetHitsWindowed", "SiPadHitsWindowed"},
      "Input SimCalorimeterHit collection names, one per detector"};
  Gaudi::Property<std::vector<std::string>> m_bitFields{
      this, "BitFields",
      {"system:8,layer:8,slice:4,plane:1,strip:14",
       "system:8,layer:8,slice:4,x:9,y:9"},
      "DD4hep BitField encoding string for each collection"};
  Gaudi::Property<std::vector<std::string>> m_sourceIDParams{
      this, "SourceIDParams",
      {"SiTargetSourceIDs", "SiPadSourceIDs"},
      "Frame parameter name carrying source_id vector for each collection"};
  Gaudi::Property<std::vector<int>> m_detectorIDs{
      this, "DetectorIDs", {0, 1},
      "detector_id value (for logging) for each collection"};
  Gaudi::Property<std::vector<std::string>> m_ntupleNames{
      this, "NTupleNames", {"SiTarget", "SiPad"},
      "RNTuple name inside the output file for each collection"};
  Gaudi::Property<std::string> m_inputFile{
      this, "InputFile", "timewindows.root",
      "Input edm4hep file (same file read by IOSvc)"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "ShipHits.root",
      "Output ROOT file containing all per-detector RNTuples"};

  // One DataHandle per collection, built in initialize() from m_collections.
  mutable std::vector<
      std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>
  > m_inputHandles;

  // One BitFieldCoder per collection.
  mutable std::vector<
      std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>
  > m_decoders;

  // Per-detector field names, populated in initialize() from the BitFieldCoders.
  mutable std::vector<std::vector<std::string>> m_detectorFieldNames;

  // podio::ROOTReader opened in initialize(), used in execute() to read
  // frame parameters (source_id vectors, t_window_start) by window index.
  mutable std::unique_ptr<podio::ROOTReader> m_reader;

  // Output TFile shared by all per-detector writers.
  mutable std::unique_ptr<TFile> m_outputTFile;

  // Per-detector RNTuple writers (one per collection).
  mutable std::vector<std::unique_ptr<ROOT::RNTupleWriter>> m_writers;

  // Per-detector field pointer bundles.
  struct DetectorFields {
    int*      fWindowID     = nullptr;
    uint64_t* fCellID       = nullptr;
    float*    fX            = nullptr;
    float*    fY            = nullptr;
    float*    fZ            = nullptr;
    float*    fE            = nullptr;
    float*    fT            = nullptr;
    int*      fSourceID     = nullptr;
    float*    fTWindowStart = nullptr;
    // Dynamic BitField field pointers for this detector only.
    std::map<std::string, int*>        fieldPtrs;
    // Keeps shared_ptrs alive after model is moved into writer.
    std::vector<std::shared_ptr<int>>  fieldSharedPtrs;
  };
  mutable std::vector<DetectorFields> m_detFields;

  // Window counter: incremented in each execute() call.
  mutable std::atomic<long long> m_windowCount{0};
};

DECLARE_COMPONENT(EDM4HEP2RNTuple)
