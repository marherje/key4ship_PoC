#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "podio/ROOTReader.h"
#include "podio/ROOTWriter.h"
#include "podio/Frame.h"
#include "ValidationUtils.h"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

// ---------------------------------------------------------------------------
// EventOverlay — event-by-event replacement for EventShuffler+EventWindowSplitter.
//
// Overlays event i of every input source into output event i: no inter-event
// delays, no time stream, no window splitting. Contribution times are kept as
// the raw simulation times. The output uses the same collection names
// (*Windowed) and frame parameters (SiTargetSourceIDs, *ContribPDGs,
// t_window_start) as EventWindowSplitter, so job3/job4/job5 run unchanged.
// ---------------------------------------------------------------------------

class EventOverlay : public Gaudi::Algorithm {
public:
  EventOverlay(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    try {
      StatusCode sc = Gaudi::Algorithm::initialize();
      if (sc.isFailure()) return sc;

      const auto& files   = m_inputFiles.value();
      const auto& ids     = m_sourceIDs.value();
      const auto& colsT   = m_collsSiTarget.value();
      const auto& colsP   = m_collsSiPad.value();
      const auto& colsMTC = m_collsMTC.value();

      if (files.empty() ||
          files.size() != ids.size()   ||
          files.size() != colsT.size() ||
          files.size() != colsP.size() ||
          files.size() != colsMTC.size()) {
        error() << "[EventOverlay] InputFiles, SourceIDs, CollectionsSiTarget,"
                << " CollectionsSiPad, CollectionsMTC must all have the same"
                << " (nonzero) size. Got: "
                << files.size() << ", " << ids.size() << ", " << colsT.size()
                << ", " << colsP.size() << ", " << colsMTC.size() << endmsg;
        return StatusCode::FAILURE;
      }

      info() << "[EventOverlay] Configured with " << files.size()
             << " input sources (event-by-event overlay):" << endmsg;
      for (size_t i = 0; i < files.size(); ++i) {
        info() << "[EventOverlay]   source_id=" << ids[i]
               << "  file=" << files[i] << endmsg;
      }
      return sc;
    } catch (const std::exception& e) {
      error() << "[EventOverlay] Exception in initialize(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventOverlay] Unknown exception in initialize()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // execute() runs exactly once (EvtMax=1). All work is done here by reading
  // the input files directly via podio::ROOTReader (bypassing IOSvc).
  StatusCode execute(const EventContext&) const override {
    try {
      const auto& files   = m_inputFiles.value();
      const auto& ids     = m_sourceIDs.value();
      const auto& colsT   = m_collsSiTarget.value();
      const auto& colsP   = m_collsSiPad.value();
      const auto& colsMTC = m_collsMTC.value();
      const int   planeOff = m_mtcPlaneOffset.value();

      // Open all sources and determine the number of overlaid events.
      std::vector<std::unique_ptr<podio::ROOTReader>> readers;
      size_t nEventsCommon = std::numeric_limits<size_t>::max();
      for (size_t i = 0; i < files.size(); ++i) {
        auto r = std::make_unique<podio::ROOTReader>();
        r->openFile(files[i]);
        const size_t n = r->getEntries("events");
        if (n == 0) {
          error() << "[EventOverlay] Source source_id=" << ids[i]
                  << " file=" << files[i] << " has 0 events." << endmsg;
          return StatusCode::FAILURE;
        }
        info() << "[EventOverlay] source_id=" << ids[i]
               << ": events=" << n << "  file=" << files[i] << endmsg;
        nEventsCommon = std::min(nEventsCommon, n);
        readers.push_back(std::move(r));
      }
      if (m_maxEventsPerSource.value() > 0) {
        nEventsCommon = std::min(
            nEventsCommon, static_cast<size_t>(m_maxEventsPerSource.value()));
      }

      // Validate collections on the first frame of every source.
      for (size_t i = 0; i < files.size(); ++i) {
        auto firstFrame = podio::Frame(readers[i]->readEntry("events", 0));
        bool ok = true;
        ok &= validateCollection(firstFrame, colsT[i],   "EventOverlay", msgStream());
        ok &= validateCollection(firstFrame, colsP[i],   "EventOverlay", msgStream());
        ok &= validateCollection(firstFrame, colsMTC[i], "EventOverlay", msgStream());
        if (!ok) {
          error() << "[EventOverlay] Validation failed for file: " << files[i]
                  << ". Aborting." << endmsg;
          return StatusCode::FAILURE;
        }
      }

      podio::ROOTWriter writer(m_outputFile.value());
      long long totalContribs = 0;

      for (size_t ev = 0; ev < nEventsCommon; ++ev) {
        // Per-subsystem accumulation buffers for THIS event only.
        std::map<uint64_t, std::vector<ContribData>> bufSiTarget;
        std::map<uint64_t, std::vector<ContribData>> bufSiPad;
        std::map<uint64_t, std::vector<ContribData>> bufMTCSciFi;
        std::map<uint64_t, std::vector<ContribData>> bufMTCScint;

        for (size_t i = 0; i < files.size(); ++i) {
          auto frame = podio::Frame(readers[i]->readEntry("events", ev));
          const int sourceID = ids[i];

          const auto& hitsTarget =
              frame.get<edm4hep::SimCalorimeterHitCollection>(colsT[i]);
          const auto& hitsPad =
              frame.get<edm4hep::SimCalorimeterHitCollection>(colsP[i]);
          const auto& hitsMTC =
              frame.get<edm4hep::SimCalorimeterHitCollection>(colsMTC[i]);

          auto accum = [&](const edm4hep::SimCalorimeterHitCollection& hits,
                           std::map<uint64_t, std::vector<ContribData>>& buf) {
            for (const auto& hit : hits) {
              auto& vec = buf[hit.getCellID()];
              const auto& hp = hit.getPosition();
              for (const auto& contrib : hit.getContributions()) {
                const auto& sp = contrib.getStepPosition();
                ContribData cd;
                cd.energy = contrib.getEnergy();
                cd.time   = contrib.getTime();  // raw sim time, no offsets
                // ddsim leaves CaloHitContribution.PDG = 0; pull it from the
                // linked MCParticle instead.
                const auto mcp = contrib.getParticle();
                cd.pdg       = mcp.isAvailable() ? mcp.getPDG() : 0;
                cd.source_id = sourceID;
                cd.stepX = sp.x; cd.stepY = sp.y; cd.stepZ = sp.z;
                cd.hitX  = hp.x; cd.hitY  = hp.y; cd.hitZ  = hp.z;
                vec.push_back(cd);
              }
            }
          };

          accum(hitsTarget, bufSiTarget);
          accum(hitsPad,    bufSiPad);

          // Route MTC hits to SciFi or Scint buffer based on 'plane' field in
          // the cell ID (plane 0/1 = SciFi U/V, plane 2 = Scint).
          for (const auto& hit : hitsMTC) {
            const int plane =
                static_cast<int>((hit.getCellID() >> planeOff) & 0x3);
            auto& buf = (plane < 2) ? bufMTCSciFi : bufMTCScint;
            auto& vec = buf[hit.getCellID()];
            const auto& hp = hit.getPosition();
            for (const auto& contrib : hit.getContributions()) {
              const auto& sp = contrib.getStepPosition();
              ContribData cd;
              cd.energy = contrib.getEnergy();
              cd.time   = contrib.getTime();
              const auto mcp = contrib.getParticle();
              cd.pdg       = mcp.isAvailable() ? mcp.getPDG() : 0;
              cd.source_id = sourceID;
              cd.stepX = sp.x; cd.stepY = sp.y; cd.stepZ = sp.z;
              cd.hitX  = hp.x; cd.hitY  = hp.y; cd.hitZ  = hp.z;
              vec.push_back(cd);
            }
          }
        }

        // Build output collections (same conventions as EventShuffler /
        // EventWindowSplitter: source_id encoded in the contribution PDG
        // field, real PDGs in a parallel frame parameter).
        auto buildColl =
            [&totalContribs](const std::map<uint64_t, std::vector<ContribData>>& buf)
            -> std::tuple<edm4hep::SimCalorimeterHitCollection,
                          edm4hep::CaloHitContributionCollection,
                          std::vector<int>,
                          std::vector<int>> {
          edm4hep::SimCalorimeterHitCollection   hits;
          edm4hep::CaloHitContributionCollection contribs;
          std::vector<int> sourceIDs;
          std::vector<int> pdgs;
          for (const auto& [cellID, vec] : buf) {
            auto outHit = hits.create();
            outHit.setCellID(cellID);
            outHit.setPosition({vec[0].hitX, vec[0].hitY, vec[0].hitZ});

            // Hit-level source_id: common value if uniform, else 0 (mixed).
            int hitSourceID = vec[0].source_id;
            for (const auto& cd : vec) {
              if (cd.source_id != hitSourceID) { hitSourceID = 0; break; }
            }
            sourceIDs.push_back(hitSourceID);

            float totalE = 0.0f;
            for (const auto& cd : vec) {
              auto contrib = contribs.create();
              contrib.setEnergy(cd.energy);
              contrib.setTime(cd.time);
              contrib.setPDG(cd.source_id);
              contrib.setStepPosition({cd.stepX, cd.stepY, cd.stepZ});
              outHit.addToContributions(contrib);
              pdgs.push_back(cd.pdg);
              totalE += cd.energy;
              ++totalContribs;
            }
            outHit.setEnergy(totalE);
          }
          return {std::move(hits), std::move(contribs),
                  std::move(sourceIDs), std::move(pdgs)};
        };

        auto [outSiTarget, outSiTargetContribs, srcSiTarget, pdgSiTarget] =
            buildColl(bufSiTarget);
        auto [outSiPad, outSiPadContribs, srcSiPad, pdgSiPad] =
            buildColl(bufSiPad);
        auto [outMTCSciFi, outMTCSciFiContribs, srcMTCSciFi, pdgMTCSciFi] =
            buildColl(bufMTCSciFi);
        auto [outMTCScint, outMTCScintContribs, srcMTCScint, pdgMTCScint] =
            buildColl(bufMTCScint);

        podio::Frame outFrame;
        outFrame.put(std::move(outSiTarget),         m_outTargetName.value());
        outFrame.put(std::move(outSiTargetContribs), m_outTargetName.value() + "_Contributions");
        outFrame.put(std::move(outSiPad),            m_outPadName.value());
        outFrame.put(std::move(outSiPadContribs),    m_outPadName.value() + "_Contributions");
        outFrame.put(std::move(outMTCSciFi),         m_outSciFiName.value());
        outFrame.put(std::move(outMTCSciFiContribs), m_outSciFiName.value() + "_Contributions");
        outFrame.put(std::move(outMTCScint),         m_outScintName.value());
        outFrame.put(std::move(outMTCScintContribs), m_outScintName.value() + "_Contributions");

        outFrame.putParameter("SiTargetSourceIDs", srcSiTarget);
        outFrame.putParameter("SiPadSourceIDs",    srcSiPad);
        outFrame.putParameter("MTCSciFiSourceIDs", srcMTCSciFi);
        outFrame.putParameter("MTCScintSourceIDs", srcMTCScint);
        outFrame.putParameter("SiTargetContribPDGs", pdgSiTarget);
        outFrame.putParameter("SiPadContribPDGs",    pdgSiPad);
        outFrame.putParameter("MTCSciFiContribPDGs", pdgMTCSciFi);
        outFrame.putParameter("MTCScintContribPDGs", pdgMTCScint);
        outFrame.putParameter("t_window_start", 0.0f);

        writer.writeFrame(outFrame, "events");
      }
      writer.finish();

      info() << "[EventOverlay] Events written: " << nEventsCommon
             << "  total contributions: " << totalContribs
             << "  output: " << m_outputFile.value() << endmsg;
      return StatusCode::SUCCESS;
    } catch (const std::exception& e) {
      error() << "[EventOverlay] Exception in execute(): " << e.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (...) {
      error() << "[EventOverlay] Unknown exception in execute()." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  StatusCode finalize() override { return Gaudi::Algorithm::finalize(); }

private:
  struct ContribData {
    float energy;
    float time;        // raw simulation time in ns
    int   source_id;
    int   pdg;         // real Geant4 PDG code from the input ddsim file
    float stepX, stepY, stepZ;
    float hitX, hitY, hitZ;
  };

  Gaudi::Property<std::vector<std::string>> m_inputFiles{
      this, "InputFiles", {}, "List of input edm4hep ROOT files (one per source)"};
  Gaudi::Property<std::vector<int>> m_sourceIDs{
      this, "SourceIDs", {}, "Source ID for each input file (same order as InputFiles)"};
  Gaudi::Property<std::vector<std::string>> m_collsSiTarget{
      this, "CollectionsSiTarget", {}, "SiTarget collection name for each input file"};
  Gaudi::Property<std::vector<std::string>> m_collsSiPad{
      this, "CollectionsSiPad", {}, "SiPad collection name for each input file"};
  Gaudi::Property<std::vector<std::string>> m_collsMTC{
      this, "CollectionsMTC", {}, "MTCDetHits collection name for each input file"};
  Gaudi::Property<int> m_mtcPlaneOffset{
      this, "MTCPlaneOffset", 22,
      "Bit offset of 'plane' field in MTCDetHits cell ID (system:8,station:2,layer:8,slice:4 = 22)"};
  Gaudi::Property<std::string> m_outputFile{
      this, "OutputFile", "events.edm4hep.root", "Output ROOT file (one frame per overlaid event)"};
  Gaudi::Property<std::string> m_outTargetName{
      this, "OutputCollectionSiTarget", "SiTargetHitsWindowed", "Output SiTarget collection name"};
  Gaudi::Property<std::string> m_outPadName{
      this, "OutputCollectionSiPad", "SiPadHitsWindowed", "Output SiPad collection name"};
  Gaudi::Property<std::string> m_outSciFiName{
      this, "OutputCollectionMTCSciFi", "MTCSciFiHitsWindowed", "Output MTC SciFi collection name"};
  Gaudi::Property<std::string> m_outScintName{
      this, "OutputCollectionMTCScint", "MTCScintHitsWindowed", "Output MTC Scint collection name"};
  Gaudi::Property<int> m_maxEventsPerSource{
      this, "MaxEventsPerSource", 0,
      "Maximum number of events to overlay (0 = all common events)"};
};

DECLARE_COMPONENT(EventOverlay)
