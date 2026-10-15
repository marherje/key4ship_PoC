#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "DDRec/CellIDPositionConverter.h"

#include "SNDCellDecoder.h"
#include "SNDSamplingWeights.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoNode.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDCaloHitBuilder — first stage of SND-PFA.
//
// Turns the overlay's SimCalorimeterHits into calibrated CalorimeterHits that
// every later stage consumes:
//
//   energy       deposited energy x the plane's sampling weight x a global
//                scale per detector, i.e. an estimate of what the shower
//                carried through this plane rather than what the thin active
//                layer absorbed (see SNDSamplingWeights.h);
//   energyError  the same deposit expressed in MIPs, which is the quantity the
//                clusterer and the muon corridor cut on. EDM4hep has no field
//                for it, and carrying it here keeps every later stage from
//                re-deriving the plane's MIP scale;
//   position     the cell centre from the segmentation, NEVER the sim hit's
//                position: SND_SciFiAction overwrites the y of every MTC hit
//                with an energy-weighted step average for its fibre
//                attenuation model;
//   type         det*10 + view, the convention in docs/snd_pfa_plan.md.
//
// Hits below ThresholdMIP are dropped. The survivors keep the order of the
// input collection and get one CaloHitSimCaloHitLink each, which is what lets
// a consumer reach the *ContribParticleIDs frame parameter -- that parameter
// is indexed by SIM hit, so the link, not the output index, is the way back.
// ---------------------------------------------------------------------------

class SNDCaloHitBuilder : public Gaudi::Algorithm {
public:
  SNDCaloHitBuilder(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    const auto& colls    = m_inputCollections.value();
    const auto& readouts = m_readoutDescriptors.value();
    const auto& dets     = m_detectorIDs.value();
    if (colls.empty() || colls.size() != readouts.size() ||
        colls.size() != dets.size()) {
      error() << "[SNDCaloHitBuilder] InputCollections, ReadoutDescriptors and "
                 "DetectorIDs must have the same (nonzero) size. Got "
              << colls.size() << ", " << readouts.size() << ", " << dets.size()
              << endmsg;
      return StatusCode::FAILURE;
    }

    const auto& outColls = m_outputCollections.value();
    const auto& linkColls = m_outputLinkCollections.value();
    if (outColls.size() != colls.size() || linkColls.size() != colls.size()) {
      error() << "[SNDCaloHitBuilder] OutputCollections and "
                 "OutputLinkCollections must match InputCollections in size."
              << endmsg;
      return StatusCode::FAILURE;
    }

    std::map<std::string, std::string> readoutByColl;
    std::map<std::string, snd::Det>    detByColl;
    for (std::size_t i = 0; i < colls.size(); ++i) {
      readoutByColl[colls[i]] = readouts[i];
      detByColl[colls[i]]     = static_cast<snd::Det>(dets[i]);
      m_handles.push_back(
          std::make_unique<DataHandle<edm4hep::SimCalorimeterHitCollection>>(
              colls[i], Gaudi::DataHandle::Reader, this));
      m_outHandles.push_back(
          std::make_unique<DataHandle<edm4hep::CalorimeterHitCollection>>(
              outColls[i], Gaudi::DataHandle::Writer, this));
      m_linkHandles.push_back(
          std::make_unique<DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection>>(
              linkColls[i], Gaudi::DataHandle::Writer, this));
    }
    m_decoder = std::make_unique<snd::SNDCellDecoder>(
        readoutByColl, detByColl, m_mtcLayersPerStation.value());

    info() << "[SNDCaloHitBuilder] " << colls.size() << " input collections, "
           << "MTC layers per station = " << m_mtcLayersPerStation.value()
           << endmsg;

    if (m_weightMode.value() == "sampling" || m_needGeometry) {
      StatusCode wsc = buildSamplingWeights();
      if (wsc.isFailure()) return wsc;
    }
    if (m_weightMode.value() != "sampling") {
      info() << "[SNDCaloHitBuilder] WeightMode='" << m_weightMode.value()
             << "': every plane weighted 1." << endmsg;
    }

    m_posConv = std::make_unique<dd4hep::rec::CellIDPositionConverter>(
        dd4hep::Detector::getInstance());

    for (const auto& entry : m_mipOverride.value()) {
      // "det:view:MeV", e.g. "2:2:3.5" to put the MTC scintillator MIP back on
      // the value the tracking report measured.
      int det = -1, view = -1;
      double mev = 0.0;
      if (std::sscanf(entry.c_str(), "%d:%d:%lf", &det, &view, &mev) == 3 &&
          mev > 0.0) {
        m_mipOverrides[{det, view}] = mev;
        info() << "[SNDCaloHitBuilder] MIP override det=" << det
               << " view=" << view << " -> " << mev << " MeV" << endmsg;
      } else {
        error() << "[SNDCaloHitBuilder] MIPEnergyOverride entry '" << entry
                << "' is not 'det:view:MeV'." << endmsg;
        return StatusCode::FAILURE;
      }
    }
    return sc;
  }

  // Walks the DD4hep geometry for the sensitive planes of the three
  // detectors and asks SNDSamplingWeights what lies in front of each.
  StatusCode buildSamplingWeights() {
    auto& detector = dd4hep::Detector::getInstance();
    if (detector.state() != dd4hep::Detector::READY) {
      if (m_compactFile.value().empty()) {
        error() << "[SNDCaloHitBuilder] No geometry is loaded and CompactFile "
                   "is empty; cannot compute sampling weights." << endmsg;
        return StatusCode::FAILURE;
      }
      detector.fromXML(m_compactFile.value());
      info() << "[SNDCaloHitBuilder] DD4hep geometry loaded from "
             << m_compactFile.value() << endmsg;
    }

    std::vector<snd::SensitivePlane> planes = snd::collectSensitivePlanes();
    if (planes.empty()) {
      error() << "[SNDCaloHitBuilder] No sensitive planes found in the "
                 "geometry." << endmsg;
      return StatusCode::FAILURE;
    }

    // Where the passive budget of each detector's first plane starts: the
    // upstream face of its own envelope, so the absorber in front of layer 0
    // is not lost. Taken as the first plane's z minus one layer pitch, which
    // the geometry gives through the spacing of the first two layers.
    std::map<int, double> startZ;
    for (int det = 0; det <= 2; ++det) {
      std::vector<double> zs;
      for (const auto& p : planes)
        if (p.det == det) zs.push_back(p.z);
      if (zs.size() < 2) continue;
      std::sort(zs.begin(), zs.end());
      double pitch = zs[1] - zs[0];
      for (std::size_t i = 1; i + 1 < zs.size(); ++i)
        pitch = std::max(pitch, zs[i + 1] - zs[i]);
      startZ[det] = zs.front() - pitch;
    }

    std::set<std::string> unknown;
    m_weights = snd::computePlaneWeights(detector, planes, startZ, &unknown);
    if (!unknown.empty()) {
      std::string list;
      for (const auto& u : unknown) list += " " + u;
      warning() << "[SNDCaloHitBuilder] materials with no dE/dx entry (using "
                << snd::kDefaultDEdxMin << " MeV cm2/g):" << list << endmsg;
    }
    m_planes = std::move(planes);
    info() << "[SNDCaloHitBuilder] sampling weights for " << m_weights.size()
           << " planes." << endmsg;
    return StatusCode::SUCCESS;
  }


  StatusCode execute(const EventContext&) const override {
    const auto& colls = m_inputCollections.value();
    for (std::size_t i = 0; i < m_handles.size(); ++i) {
      auto outColl = std::make_unique<edm4hep::CalorimeterHitCollection>();
      auto linkColl =
          std::make_unique<edm4hep::CaloHitSimCaloHitLinkCollection>();

      const auto* hits = m_handles[i]->get();
      if (!hits) {
        m_outHandles[i]->put(std::move(outColl));
        m_linkHandles[i]->put(std::move(linkColl));
        continue;
      }
      for (const auto& hit : *hits) {
        const auto a = m_decoder->decode(colls[i], hit.getCellID());
        if (!a.valid) {
          ++m_nUndecoded;
          continue;
        }
        const Key key{static_cast<int>(a.det), static_cast<int>(a.view),
                      a.layerOrdinal};
        ++m_counts[key];
        ++m_nHits;

        const double eMip = mipEnergyMeV(a);
        const double eDepMeV = hit.getEnergy() * 1000.0;  // EDM4hep GeV -> MeV
        const double inMips = (eMip > 0.0) ? eDepMeV / eMip : 0.0;
        if (inMips < m_thresholdMip.value()) {
          ++m_nBelowThreshold;
          continue;
        }

        auto outHit = outColl->create();
        outHit.setCellID(hit.getCellID());
        outHit.setEnergy(static_cast<float>(
            hit.getEnergy() * weightOf(a) * globalScale(a)));
        outHit.setEnergyError(static_cast<float>(inMips));
        outHit.setType(static_cast<int>(a.det) * 10 +
                       static_cast<int>(a.view));
        outHit.setTime(earliestTime(hit));
        outHit.setPosition(cellPosition(hit, a));

        auto link = linkColl->create();
        link.setFrom(outHit);
        link.setTo(hit);
        ++m_nOut;
        if (m_dumpCells.value() && m_nDumped < m_maxDumpedCells.value()) {
          ++m_nDumped;
          info() << "[SNDCaloHitBuilder][cell] coll=" << colls[i]
                 << " cellID=" << hit.getCellID()
                 << " det=" << static_cast<int>(a.det)
                 << " station=" << a.station << " layer=" << a.layer
                 << " plane=" << a.plane
                 << " view=" << static_cast<int>(a.view)
                 << " ordinal=" << a.layerOrdinal
                 << " channel=" << m_decoder->channel(colls[i], hit.getCellID())
                 << endmsg;
        }
      }
      m_outHandles[i]->put(std::move(outColl));
      m_linkHandles[i]->put(std::move(linkColl));
    }
    ++m_nEvents;
    return StatusCode::SUCCESS;
  }

  // ---- per-hit quantities -------------------------------------------------

  // MIP energy of this plane's active layer, i.e. what the sampling scan
  // found inside it. An override wins, so a measured MIP scale can replace
  // the dE/dx estimate without touching the code.
  double mipEnergyMeV(const snd::CellAddress& a) const {
    auto ov = m_mipOverrides.find({static_cast<int>(a.det),
                                   static_cast<int>(a.view)});
    if (ov != m_mipOverrides.end()) return ov->second;
    auto it = m_weights.find({static_cast<int>(a.det),
                              static_cast<int>(a.view), a.layerOrdinal});
    return (it != m_weights.end()) ? it->second.activeMeV : 0.0;
  }

  double weightOf(const snd::CellAddress& a) const {
    if (m_weightMode.value() != "sampling") return 1.0;
    auto it = m_weights.find({static_cast<int>(a.det),
                              static_cast<int>(a.view), a.layerOrdinal});
    return (it != m_weights.end()) ? it->second.weight : 1.0;
  }

  double globalScale(const snd::CellAddress& a) const {
    return (a.det == snd::Det::MTC) ? m_globalScaleMTC.value()
                                    : m_globalScaleECAL.value();
  }

  // The earliest contribution: a shower's later deposits in the same cell say
  // nothing about when the particle crossed it.
  static float earliestTime(const edm4hep::SimCalorimeterHit& hit) {
    float t = std::numeric_limits<float>::max();
    for (const auto& c : hit.getContributions()) t = std::min(t, c.getTime());
    return (t == std::numeric_limits<float>::max()) ? 0.f : t;
  }

  edm4hep::Vector3f cellPosition(const edm4hep::SimCalorimeterHit& hit,
                                 const snd::CellAddress& a) const {
    if (m_posConv) {
      try {
        const auto p = m_posConv->position(hit.getCellID());
        // CellIDPositionConverter answers in dd4hep units (cm in this build).
        return edm4hep::Vector3f{static_cast<float>(p.x() / dd4hep::mm),
                                 static_cast<float>(p.y() / dd4hep::mm),
                                 static_cast<float>(p.z() / dd4hep::mm)};
      } catch (const std::exception&) {
        // The VolumeManager does not know every cellID this geometry can
        // produce; rather than lose the hit, fall back on the position
        // Geant4 recorded. It is trustworthy everywhere except in the y of
        // the MTC, which SND_SciFiAction overwrites with an energy-weighted
        // step average -- and which the stereo strips do not measure anyway,
        // so it is zeroed to keep it from being mistaken for a measurement.
        ++m_nPosFallback;
      }
    }
    const auto& p = hit.getPosition();
    const bool mtcStrip = (a.det == snd::Det::MTC && a.view != snd::View::TwoD);
    return edm4hep::Vector3f{p.x, mtcStrip ? 0.f : p.y, p.z};
  }

  StatusCode finalize() override {
    info() << "[SNDCaloHitBuilder] events=" << m_nEvents
           << " hits=" << m_nHits << " written=" << m_nOut
           << " below threshold=" << m_nBelowThreshold
           << " undecoded=" << m_nUndecoded << endmsg;
    if (m_nPosFallback > 0) {
      warning() << "[SNDCaloHitBuilder] " << m_nPosFallback
                << " hits used the simulated position because the "
                   "VolumeManager could not resolve their cellID" << endmsg;
    }
    if (m_nEvents > 0) {
      info() << "[SNDCaloHitBuilder] per event: in="
             << static_cast<double>(m_nHits) / m_nEvents
             << " out=" << static_cast<double>(m_nOut) / m_nEvents << endmsg;
    }
    // One line per (detector, view): how many layers were hit and how often.
    std::map<std::pair<int, int>, std::pair<long, int>> perView;
    for (const auto& [key, n] : m_counts) {
      auto& e = perView[{key.det, key.view}];
      e.first += n;
      e.second += 1;
    }
    for (const auto& [dv, e] : perView) {
      info() << "[SNDCaloHitBuilder]   det=" << dv.first
             << " view=" << dv.second << ": hits=" << e.first
             << " layers hit=" << e.second << endmsg;
    }

    if (!m_weights.empty()) {
      info() << "[SNDCaloHitBuilder] sampling weights (det, view, ordinal, "
                "z [mm], passive [MeV], active [MeV], w):" << endmsg;
      for (const auto& p : m_planes) {
        auto it = m_weights.find({p.det, p.view, p.layerOrdinal});
        if (it == m_weights.end()) continue;
        const auto& w = it->second;
        info() << "[SNDCaloHitBuilder][w] det=" << p.det
               << " view=" << p.view << " ordinal=" << p.layerOrdinal
               << " z=" << p.z << " mat=" << p.material
               << " thick=" << p.thickness
               << " passive=" << w.passiveMeV
               << " active=" << w.activeMeV
               << " layerPassive=" << w.layerPassiveMeV
               << " layerActive=" << w.layerActiveMeV
               << " w=" << w.weight << endmsg;
      }
    }
    return Gaudi::Algorithm::finalize();
  }

private:
  struct Key {
    int det, view, ordinal;
    bool operator<(const Key& o) const {
      return std::tie(det, view, ordinal) < std::tie(o.det, o.view, o.ordinal);
    }
  };

  Gaudi::Property<std::vector<std::string>> m_inputCollections{
      this, "InputCollections", {},
      "SimCalorimeterHit collections to read (SiTarget, SiPad, MTC SciFi, "
      "MTC scint)"};
  Gaudi::Property<std::vector<std::string>> m_readoutDescriptors{
      this, "ReadoutDescriptors", {},
      "cellID descriptor of each input collection's readout, in the same "
      "order (parse_geometry.SNDGeometry.bitfields)"};
  Gaudi::Property<std::vector<int>> m_detectorIDs{
      this, "DetectorIDs", {},
      "Detector of each input collection: 0=SiTarget, 1=SiPad, 2=MTC"};
  Gaudi::Property<int> m_mtcLayersPerStation{
      this, "MTCLayersPerStation", 14,
      "MTC layers per station (DD4hep constant MTC_NLayers), used to make the "
      "MTC layer ordinal contiguous across stations"};
  Gaudi::Property<bool> m_dumpCells{
      this, "DumpDecodedCells", false,
      "Log the decoded address of the first MaxDumpedCells hits, for the "
      "step-0.4 cross-check against an independent Python decode. "
      "Diagnostic only."};
  Gaudi::Property<std::vector<std::string>> m_outputCollections{
      this, "OutputCollections", {},
      "CalorimeterHit collection to write for each input, same order"};
  Gaudi::Property<std::vector<std::string>> m_outputLinkCollections{
      this, "OutputLinkCollections", {},
      "CaloHitSimCaloHitLink collection to write for each input, same order. "
      "It is the only way back to the *ContribParticleIDs frame parameter, "
      "which is indexed by sim hit."};
  Gaudi::Property<double> m_thresholdMip{
      this, "ThresholdMIP", 0.5,
      "Drop hits below this many MIPs, as BasicDigitizer does"};
  Gaudi::Property<double> m_globalScaleECAL{
      this, "GlobalScaleECAL", 1.0,
      "Global energy scale for SiTarget+SiPad, fitted in phase 1.4"};
  Gaudi::Property<double> m_globalScaleMTC{
      this, "GlobalScaleMTC", 1.0,
      "Global energy scale for the MTC, fitted in phase 1.4"};
  Gaudi::Property<std::vector<std::string>> m_mipOverride{
      this, "MIPEnergyOverride", {},
      "Entries 'det:view:MeV' replacing the scanned MIP energy of a plane "
      "(e.g. '2:2:3.5' for the MTC scintillator). Empty = use the scan."};
  Gaudi::Property<std::string> m_weightMode{
      this, "WeightMode", "sampling",
      "'sampling' computes a per-plane weight from the material the geometry "
      "actually has in front of each sensitive layer; 'none' weights every "
      "plane 1 (raw deposited energy)."};
  Gaudi::Property<std::string> m_compactFile{
      this, "CompactFile", "",
      "DD4hep compact used to compute the sampling weights. Only read when no "
      "geometry is loaded yet, so it can be left empty in a job that already "
      "runs ACTSGeoSvc."};
  Gaudi::Property<int> m_maxDumpedCells{
      this, "MaxDumpedCells", 1000,
      "Upper bound on the lines DumpDecodedCells produces. It has to be large "
      "enough to reach every collection: the SiTarget contributes few hits to "
      "an event whose neutrino interacted downstream of it, so a small budget "
      "is spent entirely on the MTC."};

  std::vector<std::unique_ptr<DataHandle<edm4hep::SimCalorimeterHitCollection>>>
      m_handles;
  std::vector<std::unique_ptr<DataHandle<edm4hep::CalorimeterHitCollection>>>
      m_outHandles;
  std::vector<
      std::unique_ptr<DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection>>>
      m_linkHandles;
  std::unique_ptr<dd4hep::rec::CellIDPositionConverter> m_posConv;
  std::map<std::pair<int, int>, double> m_mipOverrides;
  // The position converter needs the geometry even when weights are off.
  static constexpr bool m_needGeometry = true;
  std::unique_ptr<snd::SNDCellDecoder> m_decoder;

  std::vector<snd::SensitivePlane> m_planes;
  // (det, view, layerOrdinal) -> sampling weight
  std::map<std::tuple<int, int, int>, snd::PlaneWeight> m_weights;

  mutable std::map<Key, long> m_counts;
  mutable long m_nHits = 0;
  mutable long m_nOut = 0;
  mutable long m_nBelowThreshold = 0;
  mutable long m_nPosFallback = 0;
  mutable long m_nUndecoded = 0;
  mutable long m_nEvents = 0;
  mutable int  m_nDumped = 0;
};

DECLARE_COMPONENT(SNDCaloHitBuilder)
