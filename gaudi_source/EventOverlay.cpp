#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/MCParticleCollection.h"
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
//
// It also copies every source's MCParticles into a single output MCParticles
// collection (WriteMCParticles, default true) and records where each source
// starts in it as the frame parameter MCParticleSourceOffsets. That is what
// makes *ContribParticleIDs usable downstream: the index stored there is local
// to its own source, so the particle a contribution came from is
// MCParticles[ContribParticleIDs[k] + MCParticleSourceOffsets[source_index]].
// Parent/daughter relations are preserved, remapped through the same offsets.
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

      // MCParticle collection name per source; empty property => "MCParticles"
      // everywhere, which is what ddsim writes.
      std::vector<std::string> colsMC = m_collsMCParticles.value();
      if (colsMC.empty()) colsMC.assign(files.size(), "MCParticles");
      if (m_writeMCParticles.value() && colsMC.size() != files.size()) {
        error() << "[EventOverlay] CollectionsMCParticles has " << colsMC.size()
                << " entries but there are " << files.size() << " sources."
                << endmsg;
        return StatusCode::FAILURE;
      }

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
        // Merged MCParticles of this event and the index at which each source
        // starts in them (parallel to InputFiles/SourceIDs).
        std::vector<MCData> bufMC;
        std::vector<int>    mcOffsets;
        std::vector<int>    mcSourceIDs;

        for (size_t i = 0; i < files.size(); ++i) {
          auto frame = podio::Frame(readers[i]->readEntry("events", ev));
          const int sourceID = ids[i];

          if (m_writeMCParticles.value()) {
            const int base = static_cast<int>(bufMC.size());
            mcOffsets.push_back(base);
            mcSourceIDs.push_back(sourceID);
            // Frame::get<T>() throws when the collection is absent, so go
            // through the CollectionBase overload, which returns nullptr.
            const auto* mcColl = dynamic_cast<const edm4hep::MCParticleCollection*>(
                frame.get(colsMC[i]));
            if (!mcColl) {
              if (ev == 0)
                warning() << "[EventOverlay] Source source_id=" << sourceID
                          << " has no collection '" << colsMC[i]
                          << "'; its contributions cannot be resolved to an "
                             "MCParticle." << endmsg;
            } else {
              for (const auto& mcp : *mcColl) {
                MCData m;
                m.pdg       = mcp.getPDG();
                m.genStatus = mcp.getGeneratorStatus();
                m.simStatus = mcp.getSimulatorStatus();
                m.charge    = mcp.getCharge();
                m.time      = mcp.getTime();
                m.mass      = mcp.getMass();
                const auto& v  = mcp.getVertex();
                const auto& e  = mcp.getEndpoint();
                const auto& pm = mcp.getMomentum();
                const auto& pe = mcp.getMomentumAtEndpoint();
                m.vx = v.x;  m.vy = v.y;  m.vz = v.z;
                m.ex = e.x;  m.ey = e.y;  m.ez = e.z;
                m.px = pm.x; m.py = pm.y; m.pz = pm.z;
                m.epx = pe.x; m.epy = pe.y; m.epz = pe.z;
                m.helicity = mcp.getHelicity();
                // Relations are stored by index; the source's own indices are
                // shifted by `base` so they address the merged collection.
                for (const auto& par : mcp.getParents())
                  m.parents.push_back(base + par.getObjectID().index);
                for (const auto& dau : mcp.getDaughters())
                  m.daughters.push_back(base + dau.getObjectID().index);
                bufMC.push_back(std::move(m));
              }
            }
          }

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
                cd.pdg         = mcp.isAvailable() ? mcp.getPDG() : 0;
                cd.particle_id = mcp.isAvailable() ? mcp.getObjectID().index : -1;
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
              cd.pdg         = mcp.isAvailable() ? mcp.getPDG() : 0;
              cd.particle_id = mcp.isAvailable() ? mcp.getObjectID().index : -1;
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
                          std::vector<int>,
                          std::vector<int>> {
          edm4hep::SimCalorimeterHitCollection   hits;
          edm4hep::CaloHitContributionCollection contribs;
          std::vector<int> sourceIDs;
          std::vector<int> pdgs;
          std::vector<int> particleIDs;
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
              particleIDs.push_back(cd.particle_id);
              totalE += cd.energy;
              ++totalContribs;
            }
            outHit.setEnergy(totalE);
          }
          return {std::move(hits), std::move(contribs),
                  std::move(sourceIDs), std::move(pdgs),
                  std::move(particleIDs)};
        };

        auto [outSiTarget, outSiTargetContribs, srcSiTarget, pdgSiTarget,
              partSiTarget] = buildColl(bufSiTarget);
        auto [outSiPad, outSiPadContribs, srcSiPad, pdgSiPad,
              partSiPad] = buildColl(bufSiPad);
        auto [outMTCSciFi, outMTCSciFiContribs, srcMTCSciFi, pdgMTCSciFi,
              partMTCSciFi] = buildColl(bufMTCSciFi);
        auto [outMTCScint, outMTCScintContribs, srcMTCScint, pdgMTCScint,
              partMTCScint] = buildColl(bufMTCScint);

        // Merged MCParticles: create every particle first, then wire the
        // relations, so an index can point forward as well as backward.
        edm4hep::MCParticleCollection outMC;
        if (m_writeMCParticles.value()) {
          std::vector<edm4hep::MutableMCParticle> created;
          created.reserve(bufMC.size());
          for (const auto& m : bufMC) {
            auto p = outMC.create();
            p.setPDG(m.pdg);
            p.setGeneratorStatus(m.genStatus);
            p.setSimulatorStatus(m.simStatus);
            p.setCharge(m.charge);
            p.setTime(m.time);
            p.setMass(m.mass);
            p.setVertex({m.vx, m.vy, m.vz});
            p.setEndpoint({m.ex, m.ey, m.ez});
            p.setMomentum({m.px, m.py, m.pz});
            p.setMomentumAtEndpoint({m.epx, m.epy, m.epz});
            p.setHelicity(m.helicity);
            created.push_back(p);
          }
          const int nMC = static_cast<int>(created.size());
          for (int k = 0; k < nMC; ++k) {
            for (int ip : bufMC[k].parents)
              if (ip >= 0 && ip < nMC) created[k].addToParents(created[ip]);
            for (int id : bufMC[k].daughters)
              if (id >= 0 && id < nMC) created[k].addToDaughters(created[id]);
          }
        }

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
        // Same layout as the ContribPDGs vectors above (one entry per
        // contribution, in the same order), but identifying the PARTICLE
        // rather than its species — see ContribData::particle_id.
        outFrame.putParameter("SiTargetContribParticleIDs", partSiTarget);
        outFrame.putParameter("SiPadContribParticleIDs",    partSiPad);
        outFrame.putParameter("MTCSciFiContribParticleIDs", partMTCSciFi);
        outFrame.putParameter("MTCScintContribParticleIDs", partMTCScint);
        if (m_writeMCParticles.value()) {
          outFrame.put(std::move(outMC), m_outMCName.value());
          // Index at which each source's particles start in MCParticles, in
          // the order of InputFiles/SourceIDs. Together with
          // *ContribParticleIDs (which is source-local) this identifies the
          // particle behind every contribution.
          outFrame.putParameter("MCParticleSourceOffsets", mcOffsets);
          // Source id of each offset above, so a contribution (whose source id
          // is stored in CaloHitContribution::PDG) can pick its own offset.
          outFrame.putParameter("MCParticleSourceIDs", mcSourceIDs);
        }
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
    // Index of the MCParticle this contribution came from, within its own
    // source file's MCParticle collection. Unlike pdg it identifies the
    // PARTICLE, so two same-species particles of one shower stay distinct —
    // which is what a track-purity study needs and what pdg cannot give.
    // Unique only together with source_id, since each source counts from 0.
    int   particle_id;
    float stepX, stepY, stepZ;
    float hitX, hitY, hitZ;
  };

  // One source MCParticle, flattened so it outlives the podio::Frame it was
  // read from. Parent/daughter indices are already offset into the merged
  // output collection.
  struct MCData {
    int    pdg = 0;
    int    genStatus = 0;
    int    simStatus = 0;
    float  charge = 0.f;
    float  time = 0.f;
    double mass = 0.0;
    double vx = 0, vy = 0, vz = 0;
    double ex = 0, ey = 0, ez = 0;
    double px = 0, py = 0, pz = 0;
    double epx = 0, epy = 0, epz = 0;
    int    helicity = 0;
    std::vector<int> parents;    // indices into the merged collection
    std::vector<int> daughters;  // indices into the merged collection
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
  Gaudi::Property<bool> m_writeMCParticles{
      this, "WriteMCParticles", true,
      "Copy every source's MCParticles into the output frame and write the "
      "MCParticleSourceOffsets parameter"};
  Gaudi::Property<std::vector<std::string>> m_collsMCParticles{
      this, "CollectionsMCParticles", {},
      "MCParticle collection name for each input file (empty = 'MCParticles' "
      "for all sources)"};
  Gaudi::Property<std::string> m_outMCName{
      this, "OutputCollectionMCParticles", "MCParticles",
      "Output MCParticle collection name"};
};

DECLARE_COMPONENT(EventOverlay)
