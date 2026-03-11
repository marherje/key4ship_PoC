#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "podio/ROOTReader.h"
#include "podio/ROOTWriter.h"
#include "podio/Frame.h"
#include "ValidationUtils.h"

#include <algorithm>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

struct ContribData {
  float energy;
  float time;
  // source_id is encoded in the PDG field of CaloHitContribution by EventShuffler
  // because edm4hep::CaloHitContribution has no dedicated source field.
  // PDG is not used downstream in this pipeline.
  int   source_id;
  // Exact Geant4 step position inside the sensitive volume
  float stepX, stepY, stepZ;
  // Centre of the parent SimCalorimeterHit sensitive volume
  float hitX, hitY, hitZ;
};

struct TimedContrib {
  uint64_t   cellID;
  ContribData data;
};

class EventWindowSplitter : public Gaudi::Algorithm {
public:
  EventWindowSplitter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      return Gaudi::Algorithm::initialize();
    } catch (const std::exception& e) {
      error() << "[EventWindowSplitter] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventWindowSplitter] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // execute() runs exactly once (EvtMax=1). All work is done here by reading
  // the input file directly via podio::ROOTReader (bypassing IOSvc).
  StatusCode execute(const EventContext&) const override {
    try {
      // --- 0. Open input super-event file ---
      podio::ROOTReader reader;
      reader.openFile(m_inputFile.value());

      if (reader.getEntries("events") == 0) {
        error() << "[EventWindowSplitter] Input file has 0 events." << endmsg;
        return StatusCode::FAILURE;
      }

      auto frame = podio::Frame(reader.readEntry("events", 0));

      // Validate collections and source ID parameters before processing
      bool ok = true;
      ok &= validateCollection(frame, m_inTargetName.value(), "EventWindowSplitter", msgStream());
      ok &= validateCollection(frame, m_inPixelName.value(),  "EventWindowSplitter", msgStream());
      if (!ok) return StatusCode::FAILURE;

      const auto& inTarget = frame.get<edm4hep::SimCalorimeterHitCollection>(m_inTargetName.value());
      const auto& inPixel  = frame.get<edm4hep::SimCalorimeterHitCollection>(m_inPixelName.value());

      validateParameter(frame, "SiTargetSourceIDs", inTarget.size(),
                        "EventWindowSplitter", msgStream());
      validateParameter(frame, "SiPixelSourceIDs",  inPixel.size(),
                        "EventWindowSplitter", msgStream());

      // Read source_id parameters written by EventShuffler
      std::vector<int> sourceIDsSiTarget =
          frame.getParameter<std::vector<int>>("SiTargetSourceIDs")
               .value_or(std::vector<int>{});
      std::vector<int> sourceIDsSiPixel =
          frame.getParameter<std::vector<int>>("SiPixelSourceIDs")
               .value_or(std::vector<int>{});

      // Build cellID -> source_id maps for fast lookup during window construction
      std::unordered_map<uint64_t, int> siTargetSourceMap;
      std::unordered_map<uint64_t, int> siPixelSourceMap;
      {
        size_t idx = 0;
        for (const auto& hit : inTarget) {
          siTargetSourceMap[hit.getCellID()] =
              (idx < sourceIDsSiTarget.size()) ? sourceIDsSiTarget[idx] : 0;
          ++idx;
        }
      }
      {
        size_t idx = 0;
        for (const auto& hit : inPixel) {
          siPixelSourceMap[hit.getCellID()] =
              (idx < sourceIDsSiPixel.size()) ? sourceIDsSiPixel[idx] : 0;
          ++idx;
        }
      }

      // --- 1. Build flat, sorted lists of timed contributions ---
      auto buildList = [](const edm4hep::SimCalorimeterHitCollection& coll)
          -> std::vector<TimedContrib> {
        std::vector<TimedContrib> list;
        for (const auto& hit : coll) {
          const auto& hp = hit.getPosition();
          for (const auto& contrib : hit.getContributions()) {
            const auto& sp = contrib.getStepPosition();
            list.push_back({hit.getCellID(),
                            {contrib.getEnergy(), contrib.getTime(),
                             contrib.getPDG(),  // source_id encoded in PDG field by EventShuffler
                             sp.x, sp.y, sp.z,
                             hp.x, hp.y, hp.z}});
          }
        }
        std::sort(list.begin(), list.end(),
                  [](const TimedContrib& a, const TimedContrib& b) {
                    return a.data.time < b.data.time;
                  });
        return list;
      };

      const auto siTargetList = buildList(inTarget);
      const auto siPixelList  = buildList(inPixel);

      if (siTargetList.empty() && siPixelList.empty()) {
        info() << "[EventWindowSplitter] No contributions found; nothing to write." << endmsg;
        return StatusCode::SUCCESS;
      }

      // --- 2. Determine global t_start ---
      float t_start = std::numeric_limits<float>::max();
      if (!siTargetList.empty()) t_start = std::min(t_start, siTargetList.front().data.time);
      if (!siPixelList.empty())  t_start = std::min(t_start, siPixelList.front().data.time);

      const float t_end_global = [](const auto& a, const auto& b) -> float {
        float t = std::numeric_limits<float>::lowest();
        if (!a.empty()) t = std::max(t, a.back().data.time);
        if (!b.empty()) t = std::max(t, b.back().data.time);
        return t;
      }(siTargetList, siPixelList);

      // --- 3. Assign each contribution to a window index ---
      // Window k covers [t_start + k*W, t_start + (k+1)*W).
      const float W = static_cast<float>(m_windowSize);

      auto distribute =
          [&](const std::vector<TimedContrib>& list)
          -> std::map<int, std::map<uint64_t, std::vector<ContribData>>> {
        std::map<int, std::map<uint64_t, std::vector<ContribData>>> windows;
        for (const auto& tc : list) {
          const int wIdx = std::max(0, static_cast<int>((tc.data.time - t_start) / W));
          windows[wIdx][tc.cellID].push_back(tc.data);
        }
        return windows;
      };

      const auto windowsSiTarget = distribute(siTargetList);
      const auto windowsSiPixel  = distribute(siPixelList);

      int numWindows = 0;
      if (!windowsSiTarget.empty()) numWindows = std::max(numWindows, windowsSiTarget.rbegin()->first + 1);
      if (!windowsSiPixel.empty())  numWindows = std::max(numWindows, windowsSiPixel.rbegin()->first  + 1);

      // --- 4. Write one frame per window ---
      podio::ROOTWriter writer(m_outputFile.value());
      long long totalContribs = 0;

      // Returns hits collection, contributions collection, and parallel source_id vector.
      auto buildWindowColl =
          [&](int wIdx,
              const std::map<int, std::map<uint64_t, std::vector<ContribData>>>& windows,
              const std::unordered_map<uint64_t, int>& sourceMap)
          -> std::tuple<edm4hep::SimCalorimeterHitCollection,
                        edm4hep::CaloHitContributionCollection,
                        std::vector<int>> {
        edm4hep::SimCalorimeterHitCollection   hits;
        edm4hep::CaloHitContributionCollection contribs;
        std::vector<int> sourceIDs;

        auto it = windows.find(wIdx);
        if (it == windows.end()) return {std::move(hits), std::move(contribs), std::move(sourceIDs)};

        for (const auto& [cellID, contribVec] : it->second) {
          auto hit = hits.create();
          hit.setCellID(cellID);
          hit.setPosition({contribVec[0].hitX, contribVec[0].hitY, contribVec[0].hitZ});

          const int hitSourceID = sourceMap.count(cellID) ? sourceMap.at(cellID) : 0;
          sourceIDs.push_back(hitSourceID);

          float totalE = 0.0f;
          for (const auto& cd : contribVec) {
            auto contrib = contribs.create();
            contrib.setEnergy(cd.energy);
            contrib.setTime(cd.time);
            contrib.setPDG(cd.source_id);
            contrib.setStepPosition({cd.stepX, cd.stepY, cd.stepZ});
            hit.addToContributions(contrib);
            totalE += cd.energy;
            ++totalContribs;
          }
          hit.setEnergy(totalE);
        }
        return {std::move(hits), std::move(contribs), std::move(sourceIDs)};
      };

      for (int wIdx = 0; wIdx < numWindows; ++wIdx) {
        podio::Frame windowFrame;

        auto [siTargetColl, siTargetContribs, windowSourceIDsSiTarget] =
            buildWindowColl(wIdx, windowsSiTarget, siTargetSourceMap);
        auto [siPixelColl, siPixelContribs, windowSourceIDsSiPixel] =
            buildWindowColl(wIdx, windowsSiPixel, siPixelSourceMap);

        windowFrame.put(std::move(siTargetColl),    m_outTargetName.value());
        windowFrame.put(std::move(siTargetContribs), m_outTargetName.value() + "_Contributions");
        windowFrame.put(std::move(siPixelColl),     m_outPixelName.value());
        windowFrame.put(std::move(siPixelContribs),  m_outPixelName.value()  + "_Contributions");

        windowFrame.putParameter("SiTargetSourceIDs", windowSourceIDsSiTarget);
        windowFrame.putParameter("SiPixelSourceIDs",  windowSourceIDsSiPixel);

        writer.writeFrame(windowFrame, "events");
      }
      writer.finish();

      // --- 5. Summary ---
      const long long totalInputContribs = static_cast<long long>(siTargetList.size())
                                         + static_cast<long long>(siPixelList.size());
      const double avgPerWindow = (numWindows > 0)
          ? static_cast<double>(totalInputContribs) / numWindows
          : 0.0;

      info() << "[EventWindowSplitter] Total windows written: " << numWindows << endmsg;
      info() << "[EventWindowSplitter] Time range: " << t_start << " ns to " << t_end_global << " ns" << endmsg;
      info() << "[EventWindowSplitter] Avg contributions per window: " << avgPerWindow << endmsg;

      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[EventWindowSplitter] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventWindowSplitter] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override {
    try {
      return Gaudi::Algorithm::finalize();
    } catch (const std::exception& e) {
      error() << "[EventWindowSplitter] Exception in finalize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventWindowSplitter] Unknown exception in finalize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

private:
  Gaudi::Property<std::string> m_inputFile{
      this, "InputFile", "shuffled.root", "Input super-event file"};
  Gaudi::Property<std::string> m_inTargetName{
      this, "InputCollectionSiTarget", "SiTargetHitsMerged", "Input SiTarget collection"};
  Gaudi::Property<std::string> m_inPixelName{
      this, "InputCollectionSiPixel", "SiPixelHitsMerged", "Input SiPixel collection"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "timewindows.edm4hep.root", "Output ROOT file for time-windowed events"};
  Gaudi::Property<std::string> m_outTargetName{
      this, "OutputCollectionSiTarget", "SiTargetHitsWindowed", "Output SiTarget collection name"};
  Gaudi::Property<std::string> m_outPixelName{
      this, "OutputCollectionSiPixel", "SiPixelHitsWindowed", "Output SiPixel collection name"};
  Gaudi::Property<double> m_windowSize{
      this, "WindowSize", 25.0, "Time window size (ns)"};
};

DECLARE_COMPONENT(EventWindowSplitter)
