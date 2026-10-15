#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/DataHandle.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/TrackCollection.h"

#include "SNDSamplingWeights.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDMuonFinder — finds the muon of a neutrino interaction and measures it.
//
// This is CyberPFA's TrackMatchingAlg without the helices: a track is walked
// layer by layer and the calorimeter clusters lying along it are collected
// into a CORRIDOR. A muon is the one object in this detector that leaves a
// thin, minimum-ionising trail through metres of iron while everything else
// showers, so the corridor's occupancy is what identifies it.
//
// MOMENTUM COMES FROM THE CURVATURE, with the range only as a lower bound
// (decision of 2026-09-27, replacing the opposite choice made on 2026-09-22).
// The reason the first choice was wrong is worth keeping: the whole MTC holds
// 2352.75 mm of iron, whose stopping power reaches 2.69 GeV, while the muons
// of a SHiP neutrino interaction run from 2 to over 100 GeV. NOT ONE of them
// stops inside. What the iron along a corridor then measures is not a range
// but how many layers the tracking managed to reconstruct -- median 476 mm at
// 2 GeV rising to 1960 mm above 20 GeV -- which is a tracking-quality number
// wearing a momentum's clothes, and it biased p by -55 to -85 %.
//
// So: p = 1/|q/p| from the fit, with a sigma68 of 40-50 % that is honest and
// unbiased, and p_range = 0.0011427 GeV/mm x iron as a LOWER BOUND, used when
// the fit never measured the momentum (a track that crossed no field keeps its
// seed q/p, recognisable by var(q/p)/(q/p)^2 sitting at the seed value). The
// charge always comes from the sign of q/p, which the fit gets right 84-100 %
// of the time.
//
// The corridor is built from the track's AtOther states, which carry the
// fitted position on every surface it crossed, and extrapolated straight
// beyond the last one -- a muon that exits has no more measurements, but the
// clusters it made on the way out are still there to be collected.
// ---------------------------------------------------------------------------

namespace {

// Energy loss of a muon per millimetre of iron, in the flat part of the
// Bethe curve: dE/dx_min 1.451 MeV cm^2/g x 7.874 g/cm^3 / 10 mm/cm.
constexpr double kIronMeVPerMm = 1.1427;

// Momentum a muon needs to cross a given range of iron, from the same
// relation. It is deliberately the MIP approximation: above a few GeV a muon
// loses energy by radiation too, which this underestimates, and the honest
// place to fix that is a fit to the guns rather than a fudge factor here.
constexpr double kRangeGeVPerMm = 0.0011427;

struct Corridor {
  std::vector<edm4hep::Cluster> clusters;
  int    nLayersMTC = 0;      // MTC layers with a corridor cluster
  int    nLayersSpanned = 0;  // MTC layers between the first and the last
  int    lastLayerMTC = -1;   // ordinal of the last MTC layer it reached
  double sumMip = 0.0;
  double maxMip = 0.0;
  double ironMm = 0.0;
  // Ionisation loss along the corridor, over every detector it crossed. The
  // fit measures q/p at the LAST hit, deep in the iron, so this is what has
  // to be added back to recover the momentum at the vertex.
  double eLossGeV = 0.0;
  int    nRejectedShowerLike = 0;
};

}  // namespace

class SNDMuonFinder : public Gaudi::Algorithm {
public:
  SNDMuonFinder(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;
    m_tanStereo = std::tan(m_mtcStereoDeg.value() * M_PI / 180.0);

    for (const auto& c : m_trackCollections.value())
      m_trackHandles.push_back(
          std::make_unique<DataHandle<edm4hep::TrackCollection>>(
              c, Gaudi::DataHandle::Reader, this));
    for (const auto& c : m_clusterCollections.value())
      m_clusterHandles.push_back(
          std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
              c, Gaudi::DataHandle::Reader, this));
    if (m_trackHandles.empty() || m_clusterHandles.empty()) {
      error() << "[SNDMuonFinder] TrackCollections and ClusterCollections "
                 "must both be set." << endmsg;
      return StatusCode::FAILURE;
    }

    m_outMuons = std::make_unique<
        DataHandle<edm4hep::ReconstructedParticleCollection>>(
        m_outputMuons.value(), Gaudi::DataHandle::Writer, this);
    m_outCorridor = std::make_unique<DataHandle<edm4hep::ClusterCollection>>(
        m_outputCorridor.value(), Gaudi::DataHandle::Writer, this);
    m_outPID = std::make_unique<DataHandle<edm4hep::ParticleIDCollection>>(
        m_outputPID.value(), Gaudi::DataHandle::Writer, this);

    if (StatusCode gsc = loadIronProfile(); gsc.isFailure()) return gsc;

    info() << "[SNDMuonFinder] " << m_trackHandles.size()
           << " track collections, " << m_clusterHandles.size()
           << " cluster collections; MinMTCLayers="
           << m_minMTCLayers.value() << " MinFracLayers="
           << m_minFracLayers.value() << " MaxMeanMIP="
           << m_maxMeanMip.value() << endmsg;
    return sc;
  }

  // Iron thickness of every MTC layer, derived from the same material scan
  // the sampling weights use, so a geometry variant changes it by itself.
  StatusCode loadIronProfile() {
    auto& detector = dd4hep::Detector::getInstance();
    if (detector.state() != dd4hep::Detector::READY) {
      if (m_compactFile.value().empty()) {
        error() << "[SNDMuonFinder] no geometry loaded and CompactFile empty; "
                   "the iron range cannot be computed." << endmsg;
        return StatusCode::FAILURE;
      }
      detector.fromXML(m_compactFile.value());
    }
    const auto planes = snd::collectSensitivePlanes();
    std::map<int, double> startZ;
    std::map<int, std::vector<double>> zs;
    for (const auto& p : planes) zs[p.det].push_back(p.z);
    for (auto& [det, v] : zs) {
      std::sort(v.begin(), v.end());
      if (v.size() < 2) continue;
      double pitch = v[1] - v[0];
      for (std::size_t i = 1; i + 1 < v.size(); ++i)
        pitch = std::max(pitch, v[i + 1] - v[i]);
      startZ[det] = v.front() - pitch;
    }
    const auto weights = snd::computePlaneWeights(detector, planes, startZ);

    // Two profiles from the same table. The iron of each MTC layer, in mm,
    // for the lower bound; and the MIP energy loss of every layer of every
    // detector, for the energy-loss correction. Taking both from the weight
    // table keeps "the material of this layer" one definition, not two.
    for (const auto& [key, w] : weights) {
      const auto& [det, view, ordinal] = key;
      if (view != 0 && !(det == 1 && view == 2)) continue;  // once per layer
      m_eLossMeVByLayer[{det, ordinal}] = w.layerPassiveMeV + w.layerActiveMeV;
      if (det == 2) m_ironMmByLayer[ordinal] = w.layerPassiveMeV / kIronMeVPerMm;
    }
    double total = 0.0;
    for (const auto& [ordinal, mm] : m_ironMmByLayer) total += mm;
    m_lastMTCLayer = m_lastMTCLayerProp.value() >= 0
        ? m_lastMTCLayerProp.value()
        : (m_ironMmByLayer.empty() ? 0 : m_ironMmByLayer.rbegin()->first);
    info() << "[SNDMuonFinder] iron profile: " << m_ironMmByLayer.size()
           << " MTC layers, " << total << " mm of iron in total" << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override {
    auto muons = std::make_unique<edm4hep::ReconstructedParticleCollection>();
    auto corridorOut = std::make_unique<edm4hep::ClusterCollection>();
    auto pids = std::make_unique<edm4hep::ParticleIDCollection>();

    // All layer clusters, with the layer and view their shape parameters
    // carry, and a flag so no cluster is claimed by two corridors.
    struct ClusterRef {
      edm4hep::Cluster cluster;
      int det = -1, view = -1, layer = -1;
      double mips = 0.0;
      bool used = false;
    };
    std::vector<ClusterRef> clusters;
    for (const auto& handle : m_clusterHandles) {
      const auto* coll = handle->get();
      if (!coll) continue;
      for (const auto& cl : *coll) {
        const auto sp = cl.getShapeParameters();
        if (sp.size() < 5) continue;
        ClusterRef r;
        r.cluster = cl;
        r.det = cl.getType() / 10;
        r.view = cl.getType() % 10;
        r.layer = static_cast<int>(sp[0]);
        r.mips = sp[4];
        clusters.push_back(r);
      }
    }

    // Tracks, strongest first: a corridor built around a long track has first
    // claim on the clusters it contains.
    std::vector<edm4hep::Track> tracks;
    for (const auto& handle : m_trackHandles) {
      const auto* coll = handle->get();
      if (!coll) continue;
      for (const auto& t : *coll) tracks.push_back(t);
    }
    std::sort(tracks.begin(), tracks.end(),
              [](const edm4hep::Track& a, const edm4hep::Track& b) {
                return a.trackerHits_size() > b.trackerHits_size();
              });

    for (const auto& track : tracks) {
      const Corridor corridor = buildCorridor(track, clusters);
      const bool isMuon = identify(track, corridor);
      if (!isMuon) continue;

      // --- momentum -----------------------------------------------------
      double qOverP = 0.0, varQoP = 0.0, phi = 0.0, tanLambda = 0.0;
      bool haveState = false;
      for (const auto& st : track.getTrackStates()) {
        if (st.location != edm4hep::TrackState::AtLastHit) continue;
        qOverP = st.omega;      // q/p [1/GeV], see ACTSSequentialTracker
        varQoP = st.getCovMatrix(edm4hep::TrackParams::omega,
                                 edm4hep::TrackParams::omega);
        phi = st.phi;
        tanLambda = st.tanLambda;
        haveState = true;
        break;
      }

      // Lower bound from the iron the corridor crossed: a muon that got this
      // far cannot have had less momentum than this.
      const double pRangeBound = kRangeGeVPerMm * corridor.ironMm
                                 * m_rangeScale.value();
      const double relVar = (qOverP != 0.0) ? varQoP / (qOverP * qOverP)
                                            : 1e9;
      const double pSpec = (qOverP != 0.0) ? 1.0 / std::abs(qOverP) : 0.0;
      const bool specUsable = haveState && relVar < m_maxRelQopVar.value() &&
                              pSpec > 0.5 && pSpec < 200.0;
      // The fit measures the momentum; the bound only stands in when it did
      // not, and then it is raised to the bound rather than pretending to a
      // measurement. When both exist the bound still wins if the fit fell
      // below it, since the muon demonstrably crossed that much iron.
      // The fit reports the momentum where it last measured it, after the
      // muon has paid its ionisation loss crossing the silicon and the iron
      // upstream of that point -- some 2-4 GeV, which is why q/p compared to
      // the truth at the vertex came out 15-30 % low in step 0.3. Adding the
      // corridor's own material back is the correction, and it needs no new
      // constant: the same per-layer MIP budget that gives the sampling
      // weights gives this.
      const double pSpecAtVertex = (pSpec > 0.0)
          ? pSpec + m_eLossScale.value() * corridor.eLossGeV : 0.0;
      double p = specUsable ? pSpecAtVertex : pRangeBound;
      if (specUsable && p < pRangeBound) p = pRangeBound;
      const bool useSpec = specUsable && p == pSpecAtVertex;

      const double charge = (qOverP < 0.0) ? -1.0 : (qOverP > 0.0 ? 1.0 : 0.0);
      const double theta = M_PI / 2.0 - std::atan(tanLambda);
      const double px = p * std::sin(theta) * std::cos(phi);
      const double py = p * std::sin(theta) * std::sin(phi);
      const double pz = p * std::cos(theta);

      auto mu = muons->create();
      mu.setPDG(charge < 0 ? 13 : -13);
      mu.setCharge(static_cast<float>(charge));
      // ACTS global is (beam, ddY, ddX): the momentum is written in DD4hep
      // axes, so the beam component goes last, as everything else does.
      mu.setMomentum({static_cast<float>(py), static_cast<float>(pz),
                      static_cast<float>(px)});
      mu.setEnergy(static_cast<float>(
          std::sqrt(p * p + 0.10566 * 0.10566)));
      mu.addToTracks(track);

      auto pid = pids->create();
      pid.setPDG(mu.getPDG());
      pid.setLikelihood(static_cast<float>(
          corridor.nLayersSpanned > 0
              ? double(corridor.nLayersMTC) / corridor.nLayersSpanned : 0.0));
      pid.addToParameters(static_cast<float>(corridor.nLayersMTC));
      pid.addToParameters(static_cast<float>(corridor.nLayersSpanned));
      pid.addToParameters(static_cast<float>(pid.getLikelihood()));
      pid.addToParameters(static_cast<float>(
          corridor.clusters.empty() ? 0.0
                                    : corridor.sumMip / corridor.clusters.size()));
      pid.addToParameters(static_cast<float>(corridor.maxMip));
      pid.addToParameters(static_cast<float>(corridor.ironMm));
      pid.addToParameters(static_cast<float>(
          track.getNdf() > 0 ? track.getChi2() / track.getNdf() : 0.0));
      // p_source: 0 = curvature, 1 = the iron lower bound stood in for it.
      pid.addToParameters(static_cast<float>(useSpec ? 0.0 : 1.0));
      pid.addToParameters(static_cast<float>(relVar));
      pid.addToParameters(static_cast<float>(pRangeBound));
      pid.addToParameters(static_cast<float>(pSpec));
      pid.addToParameters(static_cast<float>(corridor.eLossGeV));
      pid.setParticle(mu);

      for (const auto& cl : corridor.clusters) {
        mu.addToClusters(cl);
        corridorOut->push_back(cl.clone(false));
      }
      ++m_nMuons;
    }

    ++m_nEvents;
    m_outMuons->put(std::move(muons));
    m_outCorridor->put(std::move(corridorOut));
    m_outPID->put(std::move(pids));
    return StatusCode::SUCCESS;
  }

  StatusCode finalize() override {
    info() << "[SNDMuonFinder] " << m_nMuons << " muon candidates in "
           << m_nEvents << " events";
    if (m_nEvents > 0)
      info() << " (" << static_cast<double>(m_nMuons) / m_nEvents
             << " per event)";
    info() << endmsg;
    return Gaudi::Algorithm::finalize();
  }

private:
  struct ClusterRefFwd;   // see execute(); the corridor works on its vector

  template <typename RefVec>
  Corridor buildCorridor(const edm4hep::Track& track, RefVec& clusters) const {
    Corridor c;

    // The trajectory as the fit left it: one point per crossed surface.
    struct Pt { double beamZ, ddX, ddY; };
    std::vector<Pt> path;
    for (const auto& st : track.getTrackStates()) {
      if (st.location != edm4hep::TrackState::AtOther) continue;
      const auto& rp = st.referencePoint;
      path.push_back({rp.z, rp.x, rp.y});
    }
    if (path.size() < 2) return c;
    std::sort(path.begin(), path.end(),
              [](const Pt& a, const Pt& b) { return a.beamZ < b.beamZ; });

    // Position at an arbitrary beam z: interpolated between the two nearest
    // states, and extrapolated along the end slope outside them, which is how
    // a muon that left the detector is followed to its last cluster.
    auto at = [&path](double z, double* x, double* y) {
      if (z <= path.front().beamZ) {
        const auto& a = path.front();
        const auto& b = path[std::min<std::size_t>(1, path.size() - 1)];
        const double dz = b.beamZ - a.beamZ;
        const double t = (dz != 0.0) ? (z - a.beamZ) / dz : 0.0;
        *x = a.ddX + t * (b.ddX - a.ddX);
        *y = a.ddY + t * (b.ddY - a.ddY);
        return;
      }
      if (z >= path.back().beamZ) {
        const auto& b = path.back();
        const auto& a = path[path.size() >= 2 ? path.size() - 2 : 0];
        const double dz = b.beamZ - a.beamZ;
        const double t = (dz != 0.0) ? (z - b.beamZ) / dz : 0.0;
        *x = b.ddX + t * (b.ddX - a.ddX);
        *y = b.ddY + t * (b.ddY - a.ddY);
        return;
      }
      for (std::size_t i = 1; i < path.size(); ++i) {
        if (z <= path[i].beamZ) {
          const auto& a = path[i - 1];
          const auto& b = path[i];
          const double dz = b.beamZ - a.beamZ;
          const double t = (dz != 0.0) ? (z - a.beamZ) / dz : 0.0;
          *x = a.ddX + t * (b.ddX - a.ddX);
          *y = a.ddY + t * (b.ddY - a.ddY);
          return;
        }
      }
    };

    const double zLo = path.front().beamZ - m_extrapolateMm.value();
    const double zHi = path.back().beamZ + m_extrapolateMm.value();

    std::map<int, bool> mtcLayerHit;
    for (auto& r : clusters) {
      if (r.used) continue;
      const auto sp = r.cluster.getShapeParameters();
      const double zPlane = (sp.size() > 3) ? sp[3] : r.cluster.getPosition().z;
      if (zPlane < zLo || zPlane > zHi) continue;

      double ex = 0.0, ey = 0.0;
      at(zPlane, &ex, &ey);
      const auto& p = r.cluster.getPosition();

      // Compare only the coordinate this view measures. A stereo SciFi strip
      // knows where the track crossed it along one tilted direction and
      // nothing about the other; asking for a 2D distance there would reject
      // perfectly good clusters. Its cluster x is the fibre coordinate at
      // y = 0, s = x -+ y tan(a) for U / V (checked on truth muons on
      // 2026-10-14), so the track is compared at that same s: x alone is off
      // by |y| tan(a), 9-17 mm at |y| = 100-200 mm.
      double dist = 0.0;
      if (r.det == 2 && r.view != 2)
        dist = std::abs(p.x - (ex + (r.view == 0 ? -1.0 : 1.0) * ey *
                                         m_tanStereo));
      else if (r.view == 0) dist = std::abs(p.x - ex);
      else if (r.view == 1) dist = std::abs(p.y - ey);
      else dist = std::hypot(p.x - ex, p.y - ey);

      const double halfWidth = corridorHalfWidth(r.det, r.view);
      if (dist > halfWidth) continue;

      if (r.mips > m_corridorMaxMip.value()) {
        ++c.nRejectedShowerLike;
        continue;
      }

      r.used = true;
      c.clusters.push_back(r.cluster);
      c.sumMip += r.mips;
      c.maxMip = std::max(c.maxMip, r.mips);
      if (r.det == 2) mtcLayerHit[r.layer] = true;
    }

    // Energy loss along the corridor: every layer it touched, once.
    std::map<std::pair<int, int>, bool> layersTouched;
    for (const auto& cl : c.clusters) {
      const auto sp = cl.getShapeParameters();
      if (sp.empty()) continue;
      layersTouched[{cl.getType() / 10, static_cast<int>(sp[0])}] = true;
    }
    for (const auto& [key, _] : layersTouched) {
      auto it = m_eLossMeVByLayer.find(key);
      if (it != m_eLossMeVByLayer.end()) c.eLossGeV += it->second / 1000.0;
    }

    if (!mtcLayerHit.empty()) {
      c.nLayersMTC = static_cast<int>(mtcLayerHit.size());
      const int first = mtcLayerHit.begin()->first;
      const int last = mtcLayerHit.rbegin()->first;
      c.nLayersSpanned = last - first + 1;
      c.lastLayerMTC = last;
      for (int l = first; l <= last; ++l) {
        auto it = m_ironMmByLayer.find(l);
        if (it != m_ironMmByLayer.end()) c.ironMm += it->second;
      }
    }
    return c;
  }

  double corridorHalfWidth(int det, int view) const {
    if (det == 0) return m_halfWidthSiTarget.value();
    if (det == 1) return m_halfWidthSiPad.value();
    return (view == 2) ? m_halfWidthMTCScint.value()
                       : m_halfWidthMTCSciFi.value();
  }

  bool identify(const edm4hep::Track& track, const Corridor& c) const {
    if (!(track.getType() & 4)) return false;          // must reach the MTC
    // A muon from a vertex deep in the MTC leaves through the back after a
    // few layers; one that reaches the last layer has shown all the MTC
    // could show of it, so it needs fewer layers than one that stopped.
    const bool exits = c.lastLayerMTC >= m_lastMTCLayer - m_exitTolLayers.value();
    const int minLayers = exits ? std::min(m_minMTCLayers.value(),
                                           m_minMTCLayersAtExit.value())
                                : m_minMTCLayers.value();
    if (c.nLayersMTC < minLayers) return false;
    const double frac = (c.nLayersSpanned > 0)
        ? double(c.nLayersMTC) / c.nLayersSpanned : 0.0;
    if (frac < m_minFracLayers.value()) return false;
    const double meanMip = c.clusters.empty()
        ? 0.0 : c.sumMip / c.clusters.size();
    if (meanMip > m_maxMeanMip.value()) return false;
    return true;
  }

  Gaudi::Property<std::vector<std::string>> m_trackCollections{
      this, "TrackCollections", {},
      "Track collections to search, most complete first (GlobalTracks, then "
      "the locals that were never spliced)"};
  Gaudi::Property<std::vector<std::string>> m_clusterCollections{
      this, "ClusterCollections", {}, "Layer cluster collections"};
  Gaudi::Property<std::string> m_outputMuons{
      this, "OutputMuons", "MuonCandidates", "ReconstructedParticles to write"};
  Gaudi::Property<std::string> m_outputCorridor{
      this, "OutputCorridorClusters", "MuonCorridorClusters",
      "Clusters the corridor claimed, so phase 4 can exclude them from E_had"};
  Gaudi::Property<std::string> m_outputPID{
      this, "OutputPID", "MuonPID", "ParticleID with the corridor quantities"};
  Gaudi::Property<std::string> m_compactFile{
      this, "CompactFile", "",
      "Compact used for the iron profile; only read if no geometry is loaded"};

  Gaudi::Property<double> m_halfWidthSiTarget{
      this, "CorridorHalfWidthSiTarget", 2.0, "[mm]"};
  Gaudi::Property<double> m_halfWidthSiPad{
      this, "CorridorHalfWidthSiPad", 10.0, "[mm]"};
  Gaudi::Property<double> m_halfWidthMTCSciFi{
      this, "CorridorHalfWidthMTCSciFi", 15.0, "[mm]"};
  Gaudi::Property<double> m_halfWidthMTCScint{
      this, "CorridorHalfWidthMTCScint", 40.0, "[mm]"};
  Gaudi::Property<double> m_mtcStereoDeg{
      this, "MTCStereoAngleDeg", 5.0,
      "SciFi stereo angle: U and V clusters measure x -+ y tan(angle)"};
  Gaudi::Property<double> m_corridorMaxMip{
      this, "CorridorMaxMIP", 3.0,
      "A corridor cluster above this is shower, not a MIP trail"};
  Gaudi::Property<double> m_extrapolateMm{
      this, "ExtrapolateMm", 400.0,
      "How far past the fitted states the straight extrapolation collects "
      "clusters [mm]"};
  Gaudi::Property<int> m_minMTCLayers{
      this, "MinMTCLayers", 8, "MTC layers the corridor must have"};
  Gaudi::Property<int> m_minMTCLayersAtExit{
      this, "MinMTCLayersAtExit", 999,
      "MTC layers needed instead when the corridor reaches the last MTC "
      "layer (the muon leaves through the back); >= MinMTCLayers turns it off"};
  Gaudi::Property<int> m_exitTolLayers{
      this, "ExitToleranceLayers", 1,
      "How many layers short of the last MTC layer still counts as leaving "
      "through the back"};
  Gaudi::Property<double> m_minFracLayers{
      this, "MinFracLayers", 0.6,
      "Fraction of the spanned MTC layers that must be hit"};
  Gaudi::Property<double> m_maxMeanMip{
      this, "MaxMeanMIP", 2.0, "Mean MIPs per corridor cluster"};
  Gaudi::Property<double> m_maxRelQopVar{
      this, "MaxRelQopVar", 0.30,
      "var(q/p)/(q/p)^2 above which the fit never measured the momentum. A "
      "field-free track keeps the seed value and lands at exactly 0.36."};
  Gaudi::Property<double> m_rangeScale{
      this, "RangeScale", 1.0,
      "Scale on the iron lower bound, for the phase-3 tuning"};
  Gaudi::Property<int> m_lastMTCLayerProp{
      this, "LastMTCLayer", -1,
      "Ordinal of the last MTC layer; -1 takes it from the iron profile"};

  std::vector<std::unique_ptr<DataHandle<edm4hep::TrackCollection>>>
      m_trackHandles;
  std::vector<std::unique_ptr<DataHandle<edm4hep::ClusterCollection>>>
      m_clusterHandles;
  std::unique_ptr<DataHandle<edm4hep::ReconstructedParticleCollection>>
      m_outMuons;
  std::unique_ptr<DataHandle<edm4hep::ClusterCollection>> m_outCorridor;
  std::unique_ptr<DataHandle<edm4hep::ParticleIDCollection>> m_outPID;

  Gaudi::Property<double> m_eLossScale{
      this, "ELossScale", 1.0,
      "Scale on the ionisation loss added back to the fitted momentum. 1 is "
      "the MIP estimate from the material; above a few GeV a muon also "
      "radiates, which this does not model."};

  std::map<int, double> m_ironMmByLayer;
  // (detector, layer ordinal) -> MIP energy loss of that layer [MeV]
  std::map<std::pair<int, int>, double> m_eLossMeVByLayer;
  int m_lastMTCLayer = 0;
  double m_tanStereo = 0.0;

  mutable long m_nMuons = 0;
  mutable long m_nEvents = 0;
};

DECLARE_COMPONENT(SNDMuonFinder)
