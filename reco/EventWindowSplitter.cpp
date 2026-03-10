#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "podio/ROOTWriter.h"
#include "podio/Frame.h"

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <vector>

struct ContribData {
  float    energy;
  float    time;
  int      pdg;
  // Exact Geant4 step position inside the sensitive volume
  float    stepX, stepY, stepZ;
  // Centre of the parent SimCalorimeterHit sensitive volume
  float    hitX, hitY, hitZ;
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
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;
    m_inTargetHandle = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inTargetName.value(), Gaudi::DataHandle::Reader, this);
    m_inPixelHandle  = std::make_unique<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inPixelName.value(),  Gaudi::DataHandle::Reader, this);
    return sc;
  }

  // execute() runs exactly once on the single super-event produced by EventMerger.
  StatusCode execute(const EventContext&) const override {
    const auto* inTarget = m_inTargetHandle->get();
    const auto* inPixel  = m_inPixelHandle->get();

    // --- 1. Build flat, sorted lists of timed contributions ---
    auto buildList = [](const edm4hep::SimCalorimeterHitCollection& coll)
        -> std::vector<TimedContrib> {
      std::vector<TimedContrib> list;
      for (const auto& hit : coll) {
        const auto& hp = hit.getPosition();
        for (const auto& contrib : hit.getContributions()) {
          const auto& sp = contrib.getStepPosition();
          list.push_back({hit.getCellID(),
                          {contrib.getEnergy(), contrib.getTime(), contrib.getPDG(),
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

    const auto siTargetList = buildList(*inTarget);
    const auto siPixelList  = buildList(*inPixel);

    if (siTargetList.empty() && siPixelList.empty()) {
      info() << "[EventWindowSplitter] No contributions found; nothing to write." << endmsg;
      return StatusCode::SUCCESS;
    }

    // --- 2. Determine global t_start ---
    float t_start = std::numeric_limits<float>::max();
    if (!siTargetList.empty()) t_start = std::min(t_start, siTargetList.front().data.time);
    if (!siPixelList.empty())  t_start = std::min(t_start, siPixelList.front().data.time);

    const float t_end_global = [] (const auto& a, const auto& b) -> float {
      float t = std::numeric_limits<float>::lowest();
      if (!a.empty()) t = std::max(t, a.back().data.time);
      if (!b.empty()) t = std::max(t, b.back().data.time);
      return t;
    }(siTargetList, siPixelList);

    // --- 3. Assign each contribution to a window index ---
    // Window k covers [t_start + k*W, t_start + (k+1)*W).
    // Boundary contributions (time == t_start + k*W for k>0) go to window k (next).
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

    // Determine total number of windows
    int numWindows = 0;
    if (!windowsSiTarget.empty()) numWindows = std::max(numWindows, windowsSiTarget.rbegin()->first + 1);
    if (!windowsSiPixel.empty())  numWindows = std::max(numWindows, windowsSiPixel.rbegin()->first  + 1);

    // --- 4. Write one frame per window ---
    podio::ROOTWriter writer(m_outputFile.value());
    long long totalContribs = 0;

    auto buildWindowColl =
        [&](int wIdx,
            const std::map<int, std::map<uint64_t, std::vector<ContribData>>>& windows)
        -> std::pair<edm4hep::SimCalorimeterHitCollection,
                     edm4hep::CaloHitContributionCollection> {
      edm4hep::SimCalorimeterHitCollection  hits;
      edm4hep::CaloHitContributionCollection contribs;

      auto it = windows.find(wIdx);
      if (it == windows.end()) return {std::move(hits), std::move(contribs)};

      for (const auto& [cellID, contribVec] : it->second) {
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
          ++totalContribs;
        }
        hit.setEnergy(totalE);
      }
      return {std::move(hits), std::move(contribs)};
    };

    for (int wIdx = 0; wIdx < numWindows; ++wIdx) {
      podio::Frame frame;

      auto [siTargetColl,   siTargetContribs] = buildWindowColl(wIdx, windowsSiTarget);
      auto [siPixelColl,    siPixelContribs]  = buildWindowColl(wIdx, windowsSiPixel);

      frame.put(std::move(siTargetColl),    m_outTargetName.value());
      frame.put(std::move(siTargetContribs), m_outTargetName.value() + "_Contributions");
      frame.put(std::move(siPixelColl),     m_outPixelName.value());
      frame.put(std::move(siPixelContribs),  m_outPixelName.value()  + "_Contributions");

      writer.writeFrame(frame, "events");
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
  }

  StatusCode finalize() override {
    m_inTargetHandle.reset();
    m_inPixelHandle.reset();
    return Gaudi::Algorithm::finalize();
  }

private:
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

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inTargetHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>> m_inPixelHandle;
};

DECLARE_COMPONENT(EventWindowSplitter)
