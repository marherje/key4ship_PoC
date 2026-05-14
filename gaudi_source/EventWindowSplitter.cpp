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
      ok &= validateCollection(frame, m_inPadName.value(),    "EventWindowSplitter", msgStream());
      ok &= validateCollection(frame, m_inSciFiName.value(),  "EventWindowSplitter", msgStream());
      ok &= validateCollection(frame, m_inScintName.value(),  "EventWindowSplitter", msgStream());
      if (!ok) return StatusCode::FAILURE;

      const auto& inTarget = frame.get<edm4hep::SimCalorimeterHitCollection>(m_inTargetName.value());
      const auto& inPad    = frame.get<edm4hep::SimCalorimeterHitCollection>(m_inPadName.value());
      const auto& inSciFi  = frame.get<edm4hep::SimCalorimeterHitCollection>(m_inSciFiName.value());
      const auto& inScint  = frame.get<edm4hep::SimCalorimeterHitCollection>(m_inScintName.value());

      validateParameter(frame, "SiTargetSourceIDs", inTarget.size(),
                        "EventWindowSplitter", msgStream());
      validateParameter(frame, "SiPadSourceIDs",    inPad.size(),
                        "EventWindowSplitter", msgStream());
      validateParameter(frame, "MTCSciFiSourceIDs", inSciFi.size(),
                        "EventWindowSplitter", msgStream());
      validateParameter(frame, "MTCScintSourceIDs", inScint.size(),
                        "EventWindowSplitter", msgStream());

      // Read source_id parameters written by EventShuffler
      std::vector<int> sourceIDsSiTarget =
          frame.getParameter<std::vector<int>>("SiTargetSourceIDs").value_or(std::vector<int>{});
      std::vector<int> sourceIDsSiPad =
          frame.getParameter<std::vector<int>>("SiPadSourceIDs").value_or(std::vector<int>{});
      std::vector<int> sourceIDsMTCSciFi =
          frame.getParameter<std::vector<int>>("MTCSciFiSourceIDs").value_or(std::vector<int>{});
      std::vector<int> sourceIDsMTCScint =
          frame.getParameter<std::vector<int>>("MTCScintSourceIDs").value_or(std::vector<int>{});

      auto makeSourceMap = [](const edm4hep::SimCalorimeterHitCollection& coll,
                              const std::vector<int>& ids) {
        std::unordered_map<uint64_t, int> m;
        size_t idx = 0;
        for (const auto& hit : coll)
          m[hit.getCellID()] = (idx < ids.size()) ? ids[idx++] : 0;
        return m;
      };

      // Build cellID -> source_id maps for fast lookup during window construction
      auto siTargetSourceMap = makeSourceMap(inTarget, sourceIDsSiTarget);
      auto siPadSourceMap    = makeSourceMap(inPad,    sourceIDsSiPad);
      auto mtcSciFiSourceMap = makeSourceMap(inSciFi,  sourceIDsMTCSciFi);
      auto mtcScintSourceMap = makeSourceMap(inScint,  sourceIDsMTCScint);

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
      const auto siPadList    = buildList(inPad);
      const auto mtcSciFiList = buildList(inSciFi);
      const auto mtcScintList = buildList(inScint);

      if (siTargetList.empty() && siPadList.empty() &&
          mtcSciFiList.empty() && mtcScintList.empty()) {
        info() << "[EventWindowSplitter] No contributions found; nothing to write." << endmsg;
        return StatusCode::SUCCESS;
      }

      // --- 2. Determine global t_start ---
      float t_start = std::numeric_limits<float>::max();
      for (const auto* lst : {&siTargetList, &siPadList, &mtcSciFiList, &mtcScintList})
        if (!lst->empty()) t_start = std::min(t_start, lst->front().data.time);

      float t_end_global = std::numeric_limits<float>::lowest();
      for (const auto* lst : {&siTargetList, &siPadList, &mtcSciFiList, &mtcScintList})
        if (!lst->empty()) t_end_global = std::max(t_end_global, lst->back().data.time);

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

      const auto windowsSiTarget  = distribute(siTargetList);
      const auto windowsSiPad     = distribute(siPadList);
      const auto windowsMTCSciFi  = distribute(mtcSciFiList);
      const auto windowsMTCScint  = distribute(mtcScintList);

      int numWindows = 0;
      for (const auto* w : {&windowsSiTarget, &windowsSiPad, &windowsMTCSciFi, &windowsMTCScint})
        if (!w->empty()) numWindows = std::max(numWindows, w->rbegin()->first + 1);

      // --- 4. Write one frame per window ---
      podio::ROOTWriter writer(m_outputFile.value());
      long long totalContribs = 0;

      // Returns hits collection, contributions collection, and parallel source_id vector.
      auto buildWindowColl =
          [&](int wIdx,
              float t_window_start,
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
            contrib.setTime(cd.time - t_window_start);
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
        const float t_window_start =
            t_start + static_cast<float>(wIdx) * static_cast<float>(m_windowSize);
        podio::Frame windowFrame;

        auto [siTargetColl,  siTargetContribs,  wSrcTarget]  =
            buildWindowColl(wIdx, t_window_start, windowsSiTarget,  siTargetSourceMap);
        auto [siPadColl,    siPadContribs,    wSrcPad]    =
            buildWindowColl(wIdx, t_window_start, windowsSiPad,     siPadSourceMap);
        auto [mtcSciFiColl, mtcSciFiContribs, wSrcSciFi]  =
            buildWindowColl(wIdx, t_window_start, windowsMTCSciFi,  mtcSciFiSourceMap);
        auto [mtcScintColl, mtcScintContribs, wSrcScint]  =
            buildWindowColl(wIdx, t_window_start, windowsMTCScint,  mtcScintSourceMap);

        windowFrame.put(std::move(siTargetColl),    m_outTargetName.value());
        windowFrame.put(std::move(siTargetContribs), m_outTargetName.value() + "_Contributions");
        windowFrame.put(std::move(siPadColl),       m_outPadName.value());
        windowFrame.put(std::move(siPadContribs),   m_outPadName.value()  + "_Contributions");
        windowFrame.put(std::move(mtcSciFiColl),    m_outSciFiName.value());
        windowFrame.put(std::move(mtcSciFiContribs), m_outSciFiName.value() + "_Contributions");
        windowFrame.put(std::move(mtcScintColl),    m_outScintName.value());
        windowFrame.put(std::move(mtcScintContribs), m_outScintName.value() + "_Contributions");

        windowFrame.putParameter("SiTargetSourceIDs", wSrcTarget);
        windowFrame.putParameter("SiPadSourceIDs",    wSrcPad);
        windowFrame.putParameter("MTCSciFiSourceIDs", wSrcSciFi);
        windowFrame.putParameter("MTCScintSourceIDs", wSrcScint);
        windowFrame.putParameter("t_window_start",    t_window_start);

        writer.writeFrame(windowFrame, "events");
      }
      writer.finish();

      // --- 5. Summary ---
      const long long totalInputContribs = static_cast<long long>(siTargetList.size())
                                         + static_cast<long long>(siPadList.size())
                                         + static_cast<long long>(mtcSciFiList.size())
                                         + static_cast<long long>(mtcScintList.size());
      const double avgPerWindow = (numWindows > 0)
          ? static_cast<double>(totalInputContribs) / numWindows : 0.0;

      info() << "[EventWindowSplitter] Total windows written: " << numWindows << endmsg;
      info() << "[EventWindowSplitter] Time range: " << t_start << " to " << t_end_global << " ns" << endmsg;
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
      this, "InputCollectionSiTarget", "SiTargetHitsMerged",   "Input SiTarget collection"};
  Gaudi::Property<std::string> m_inPadName{
      this, "InputCollectionSiPad",    "SiPadHitsMerged",      "Input SiPad collection"};
  Gaudi::Property<std::string> m_inSciFiName{
      this, "InputCollectionMTCSciFi", "MTCSciFiHitsMerged",   "Input MTC SciFi collection"};
  Gaudi::Property<std::string> m_inScintName{
      this, "InputCollectionMTCScint", "MTCScintHitsMerged",   "Input MTC Scint collection"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "timewindows.edm4hep.root", "Output ROOT file for time-windowed events"};
  Gaudi::Property<std::string> m_outTargetName{
      this, "OutputCollectionSiTarget", "SiTargetHitsWindowed", "Output SiTarget collection name"};
  Gaudi::Property<std::string> m_outPadName{
      this, "OutputCollectionSiPad",    "SiPadHitsWindowed",    "Output SiPad collection name"};
  Gaudi::Property<std::string> m_outSciFiName{
      this, "OutputCollectionMTCSciFi", "MTCSciFiHitsWindowed", "Output MTC SciFi collection name"};
  Gaudi::Property<std::string> m_outScintName{
      this, "OutputCollectionMTCScint", "MTCScintHitsWindowed", "Output MTC Scint collection name"};
  Gaudi::Property<double> m_windowSize{
      this, "WindowSize", 25.0, "Time window size (ns)"};
};

DECLARE_COMPONENT(EventWindowSplitter)
