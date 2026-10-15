// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/DataHandle.h"

// edm4hep input
#include "edm4hep/SimCalorimeterHitCollection.h"

// edm4hep output
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/MutableTrackerHit3D.h"
#include "edm4hep/CovMatrix3f.h"
#include "edm4hep/Vector3d.h"

// DD4hep segmentation
#include "DDSegmentation/BitFieldCoder.h"

// ACTS
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

// SND geometry service
#include "ISNDGeoSvc.h"

// Standard
#include <atomic>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace {
// One decoded SciFi hit. plane is restricted to 0 (U) or 1 (V); scintillator
// (plane=2) is filtered out before this struct is filled.
struct DecodedHit {
  std::size_t index;    // position in input collection
  int         station;
  int         layer;
  int         plane;    // 0 = U, 1 = V
  double      s_raw;    // strip stereo coordinate [mm] (pos.x, already mm)
  double      z;        // beam coord [mm]
  float       time;
  double      edep;
  uint64_t    cellID;
};

// Build groups of hit indices according to the preselection strategy.
//   group.size() == 1 → emit as 1D measurement on the U or V stereo surface.
//   group.size() == 2 → emit as 2D measurement on the combined (flat) surface.
// To add a new strategy (clustering, ML, …): add a `method == "<Name>"` branch
// below. The converter's I/O contract does not change.
std::vector<std::vector<std::size_t>>
groupHits(const std::vector<DecodedHit>& decoded,
          const std::string&             method,
          double                         maxDz) {
  std::vector<std::vector<std::size_t>> groups;
  if (method == "None") {
    groups.reserve(decoded.size());
    for (std::size_t i = 0; i < decoded.size(); ++i)
      groups.push_back({i});
    return groups;
  }
  if (method == "AllPairs") {
    std::map<std::tuple<int, int, int>, std::vector<std::size_t>> bucket;
    for (std::size_t i = 0; i < decoded.size(); ++i) {
      const auto& d = decoded[i];
      bucket[{d.station, d.layer, d.plane}].push_back(i);
    }
    std::set<std::pair<int, int>> stationLayers;
    for (const auto& kv : bucket)
      stationLayers.insert({std::get<0>(kv.first), std::get<1>(kv.first)});
    for (auto [station, layer] : stationLayers) {
      auto itU = bucket.find({station, layer, 0});
      auto itV = bucket.find({station, layer, 1});
      const std::vector<std::size_t>* Us =
          (itU != bucket.end()) ? &itU->second : nullptr;
      const std::vector<std::size_t>* Vs =
          (itV != bucket.end()) ? &itV->second : nullptr;
      if (Us && Vs) {
        std::vector<bool> usedU(Us->size(), false);
        std::vector<bool> usedV(Vs->size(), false);
        for (std::size_t iu = 0; iu < Us->size(); ++iu) {
          for (std::size_t iv = 0; iv < Vs->size(); ++iv) {
            const double dz =
                std::abs(decoded[(*Us)[iu]].z - decoded[(*Vs)[iv]].z);
            if (dz <= maxDz) {
              groups.push_back({(*Us)[iu], (*Vs)[iv]});
              usedU[iu] = true;
              usedV[iv] = true;
            }
          }
        }
        for (std::size_t iu = 0; iu < Us->size(); ++iu)
          if (!usedU[iu]) groups.push_back({(*Us)[iu]});
        for (std::size_t iv = 0; iv < Vs->size(); ++iv)
          if (!usedV[iv]) groups.push_back({(*Vs)[iv]});
      } else if (Us) {
        for (auto i : *Us) groups.push_back({i});
      } else if (Vs) {
        for (auto i : *Vs) groups.push_back({i});
      }
    }
    return groups;
  }
  // Unknown method → fail-safe singleton emission (= "None").
  for (std::size_t i = 0; i < decoded.size(); ++i)
    groups.push_back({i});
  return groups;
}
}  // namespace

// ---------------------------------------------------------------------------
// MTCSciFiMeasConverter
//
// Converts MTCSciFiHitsWindowed (SimCalorimeterHitCollection of raw windowed
// hits) into TrackerHit3D measurements for ACTS tracking.
//
// Unit note: SND_SciFiAction stores position.x as seg_p.x()*10 (ROOT cm→mm),
// so all three position components arrive in mm.  position.y (truth avg-Y for
// light propagation) is NOT stored here — it must never be used in tracking.
//
// Surface matching uses ACTSGeoSvc::surfaceByAddress(detID=2, station, layer,
// plane) decoded from the hit cellID — no Z-proximity arithmetic.
//
// Output Quality field encodes the kind of measurement:
//   0 = U plane 1D strip  → ACTS surface plane=0, 1D loc0
//   1 = V plane 1D strip  → ACTS surface plane=1, 1D loc0
//   3 = combined U+V pair → ACTS surface plane=3, 2D (loc0,loc1) = (DD4hep X,Y)
//
// Strategy is selected by the PairMethod property:
//   "None"     — emit U and V hits as separate 1D measurements (legacy)
//   "AllPairs" — group hits by (station, layer), emit every U×V combination
//                as a 2D measurement on the combined surface; fall back to 1D
//                for layers with only one plane firing.
// Future methods (clustering, ML preselection) plug in by adding a strategy
// case to groupHits() — the converter's I/O contract stays unchanged.
// ---------------------------------------------------------------------------

class MTCSciFiMeasConverter : public Gaudi::Algorithm {
public:
  MTCSciFiMeasConverter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  // ---- Properties ----
  Gaudi::Property<std::string> m_inputCollection{
      this, "InputCollection", "MTCSciFiHitsWindowed",
      "Input SimCalorimeterHitCollection name (windowed raw hits)"};

  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "MTCSciFiMeasurements",
      "Output TrackerHit3DCollection name"};

  Gaudi::Property<std::string> m_bitField{
      this, "BitField",
      "system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9",
      "DD4hep BitField string for MTC cellID decoding"};

  Gaudi::Property<double> m_stripPitch{
      this, "StripPitch", 1.0,
      "SciFi strip pitch [mm] — determines measurement variance"};

  Gaudi::Property<double> m_stereoAngleDeg{
      this, "StereoAngleDeg", 5.0,
      "Fiber stereo angle magnitude [degrees] — used to compute eBoundLoc0 = cos(α)·x_stereo"};

  Gaudi::Property<std::string> m_pairMethod{
      this, "PairMethod", "AllPairs",
      "U/V preselection strategy: 'None' emits separate 1D U,V hits; "
      "'AllPairs' emits every U×V combination per (station,layer) as a 2D "
      "measurement on the combined surface. Future: 'Clustering', 'ML', ..."};

  Gaudi::Property<double> m_pairMaxDz{
      this, "PairMaxDz", 25.0,
      "Maximum |z_U − z_V| [mm] within a (station,layer) to allow pairing. "
      "Layers with separation > this fall back to 1D output."};

  // A U/V pair is a GHOST when its two fibres were crossed by different
  // particles: the reconstructed crossing then sits where no particle ever
  // was. With N particles in a layer, AllPairs yields N^2 crossings of which
  // N are real — the MTC counterpart of the SiTarget's StripX x StripY ghost
  // problem that ACTSSequentialTracker's GhostResolveWindow addresses using
  // SiPad points. Here the confirmer is the MTC's own scintillator plane
  // (plane=2), which the reconstruction otherwise ignores entirely.
  //
  // Measured on 500 GENIE events: a particle lights a median of 2 scint cells,
  // and in 53% of multi-particle layers two particles light COMPLETELY
  // disjoint cell sets (median Jaccard overlap 0.00), so the scint does
  // separate particles spatially — which is what makes this veto possible.
  Gaudi::Property<std::string> m_scintCollection{
      this, "ScintCollection", "",
      "MTC scintillator hit collection used to confirm U/V crossings. Empty "
      "(default) disables the ghost veto entirely, leaving every pipeline that "
      "does not set it byte-for-byte unchanged."};
  Gaudi::Property<double> m_scintConfirmWindow{
      this, "ScintConfirmWindow", 0.0,
      "Half-width [mm] of the window around a scintillator cell centre within "
      "which a U/V crossing counts as confirmed. Cells are "
      "MTC_ScintCellSize=50 mm wide and the crossing resolution is a few mm, "
      "so ~30 mm covers the cell plus slack. A pair whose crossing lands in no "
      "lit cell of its own layer is dropped as a ghost. 0.0 disables. "
      "Fail-open by construction: a layer with NO scint hits confirms "
      "everything, so scintillator inefficiency can never delete real hits."};

  Gaudi::Property<double> m_acceptanceTolerance{
      this, "AcceptanceTolerance", 2.0,
      "A U/V crossing is kept only if it lies inside the combined surface's "
      "bounds (the SciFi plane, +-300 mm) plus this margin [mm]. AllPairs "
      "crosses every U with every V of a layer; two fibres far apart in x "
      "intersect at |y| = |dx| / (2 tan 5 deg), i.e. a 300 mm separation "
      "'lands' at |y| = 1.7 m, well outside any physical plane. Those points "
      "are pure ghosts, cost seeds and can be drawn outside the detector. "
      "Negative disables the check."};

  // ---- DataHandles ----
  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>
      m_inputHandle;

  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>
      m_scintHandle;

  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_outputHandle;

  // ---- Services ----
  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  // ---- Internal state ----
  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::atomic<long long> m_eventCount{0};
};

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode MTCSciFiMeasConverter::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    m_inputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
        m_inputCollection.value(), Gaudi::DataHandle::Reader, this);

    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    const bool ghostVeto =
        !m_scintCollection.value().empty() && m_scintConfirmWindow.value() > 0.0;
    if (ghostVeto) {
      m_scintHandle = std::make_unique<
          k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection>>(
          m_scintCollection.value(), Gaudi::DataHandle::Reader, this);
      info() << "[MTCSciFiMeasConverter] U/V ghost veto ON: confirming "
             << "crossings against '" << m_scintCollection.value()
             << "' within " << m_scintConfirmWindow.value() << " mm." << endmsg;
    } else if (!m_scintCollection.value().empty() ||
               m_scintConfirmWindow.value() > 0.0) {
      // Half-configured is a silent no-op otherwise, and the symptom (ghosts
      // still there) looks identical to the veto not working.
      warning() << "[MTCSciFiMeasConverter] U/V ghost veto INACTIVE: it needs "
                << "BOTH ScintCollection (now '" << m_scintCollection.value()
                << "') and ScintConfirmWindow > 0 (now "
                << m_scintConfirmWindow.value() << ")." << endmsg;
    }

    m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
        m_bitField.value());

    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[MTCSciFiMeasConverter] Failed to retrieve ACTSGeoSvc."
              << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "[MTCSciFiMeasConverter] Initialized. GeoSvc has "
           << m_geoSvc->allSurfaces().size() << " surfaces." << endmsg;

    return sc;
  } catch (const std::exception& e) {
    error() << "[MTCSciFiMeasConverter] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[MTCSciFiMeasConverter] Unknown exception in initialize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode MTCSciFiMeasConverter::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    auto* output = m_outputHandle->createAndPut();

    const auto* hits = m_inputHandle->get();
    if (!hits || hits->size() == 0) {
      debug() << "[MTCSciFiMeasConverter] evt=" << evtNum
              << " no MTC SciFi hits." << endmsg;
      return StatusCode::SUCCESS;
    }

    const double pitch     = m_stripPitch.value();     // mm
    const double alpha_rad = m_stereoAngleDeg.value() * M_PI / 180.0;
    const double cos_a     = std::cos(alpha_rad);
    const double tan_a     = std::tan(alpha_rad);

    // 1D loc0 variance: eBoundLoc0 = cos(α)·x_stereo → σ² scales with cos²α.
    const double var1D = cos_a * cos_a * (pitch * pitch) / 12.0;       // mm²
    // 2D combined variances on the flat (non-stereo) surface.
    // Canonical convention (see ACTSProtoTracker comments around toGlobalX/Y):
    //   s_U =  cos(α)·X − sin(α)·Y   (plane=0, stereo +α)
    //   s_V =  cos(α)·X + sin(α)·Y   (plane=1, stereo −α)
    // ⇒ X_g = (s_U + s_V)/(2 cos α);  Y_g = (s_V − s_U)/(2 sin α).
    //   σ²(X_g) = pitch²/(24 cos²α),  σ²(Y_g) = pitch²/(24 sin²α).
    //   Off-diagonal = (σ²(s_U) − σ²(s_V)) / (4 sin α cos α) = 0 since σ_U = σ_V.
    const double var2D_x = (pitch * pitch) / 24.0;
    const double var2D_y = (pitch * pitch) / (24.0 * tan_a * tan_a);

    // ---- Decode hits, skip scintillator -----------------------------------
    std::vector<DecodedHit> decoded;
    decoded.reserve(hits->size());
    std::size_t nScint = 0;
    for (std::size_t i = 0; i < hits->size(); ++i) {
      const auto& hit = (*hits)[i];
      const uint64_t cid = hit.getCellID();
      const int plane   = static_cast<int>(m_decoder->get(cid, "plane"));
      if (plane != 0 && plane != 1) {
        ++nScint;
        continue;
      }
      DecodedHit d;
      d.index   = i;
      d.station = static_cast<int>(m_decoder->get(cid, "station"));
      d.layer   = static_cast<int>(m_decoder->get(cid, "layer"));
      d.plane   = plane;
      d.s_raw   = hit.getPosition().x;
      d.z       = hit.getPosition().z;
      d.time    = (hit.contributions_size() > 0)
                  ? static_cast<float>(hit.getContributions(0).getTime())
                  : 0.0f;
      d.edep    = hit.getEnergy();
      d.cellID  = cid;
      decoded.push_back(d);
    }

    // ---- Scintillator confirmation map ------------------------------------
    // (station, layer) -> transverse centres of the lit scint cells. Built once
    // per event; empty when the veto is off, which short-circuits every check
    // below and leaves the output identical to before.
    std::map<std::pair<int, int>, std::vector<std::pair<double, double>>> scintCells;
    const double scintWin =
        m_scintHandle ? m_scintConfirmWindow.value() : 0.0;
    if (m_scintHandle) {
      const auto* scintHits = m_scintHandle->get();
      if (scintHits) {
        for (const auto& sh : *scintHits) {
          const uint64_t scid = sh.getCellID();
          const int st = static_cast<int>(m_decoder->get(scid, "station"));
          const int ly = static_cast<int>(m_decoder->get(scid, "layer"));
          const auto& sp = sh.getPosition();
          scintCells[{st, ly}].emplace_back(sp.x, sp.y);
        }
      }
    }

    // ---- Apply preselection strategy --------------------------------------
    const auto groups =
        groupHits(decoded, m_pairMethod.value(), m_pairMaxDz.value());

    // ---- Emit output ------------------------------------------------------
    std::size_t n1D = 0, n2D = 0, nSkipped = 0, nGhostVetoed = 0, nOutside = 0;
    for (const auto& group : groups) {
      if (group.size() == 1) {
        const auto& d = decoded[group[0]];
        const Acts::Surface* surf =
            m_geoSvc->surfaceByAddress(2, d.station, d.layer, d.plane);
        if (!surf) {
          warning() << "[MTCSciFiMeasConverter] evt=" << evtNum
                    << " 1D hit station=" << d.station
                    << " layer=" << d.layer
                    << " plane=" << d.plane
                    << " has no matching ACTS surface — skipping." << endmsg;
          ++nSkipped;
          continue;
        }
        const double localCoord = cos_a * d.s_raw;
        edm4hep::CovMatrix3f cov{};
        cov[0] = static_cast<float>(var1D);
        auto mhit = output->create();
        mhit.setCellID(d.cellID);
        mhit.setType(2);
        mhit.setQuality(d.plane);   // 0 = U 1D, 1 = V 1D
        mhit.setTime(d.time);
        mhit.setEDep(d.edep);
        mhit.setEDepError(0.0f);
        mhit.setPosition(edm4hep::Vector3d{localCoord, 0.0, d.z});
        mhit.setCovMatrix(cov);
        ++n1D;
      } else if (group.size() == 2) {
        const auto& d0 = decoded[group[0]];
        const auto& d1 = decoded[group[1]];
        const auto& dU = (d0.plane == 0) ? d0 : d1;
        const auto& dV = (d0.plane == 1) ? d0 : d1;
        const Acts::Surface* surf =
            m_geoSvc->surfaceByAddress(2, dU.station, dU.layer, 3);
        if (!surf) {
          warning() << "[MTCSciFiMeasConverter] evt=" << evtNum
                    << " 2D pair station=" << dU.station
                    << " layer=" << dU.layer
                    << " has no combined ACTS surface (plane=3) — skipping."
                    << endmsg;
          ++nSkipped;
          continue;
        }
        // s_raw is what CartesianStripXStereo::position() returns: the strip
        // centre at y = 0, i.e. x − y·tanα for U (+α) and x + y·tanα for V (−α).
        // The stereo coordinate proper is s = cosα·s_raw, so inverting in
        // s_raw directly: X = (s_U + s_V)/2,  Y = (s_V − s_U)/(2 tan α).
        // (Dividing by 2cosα / 2sinα, as before, scaled both by 1/cosα —
        // a 0.4% bias, 1 mm at the acceptance edge, larger than the pitch.)
        const double X_g = 0.5 * (dU.s_raw + dV.s_raw);
        const double Y_g = (dV.s_raw - dU.s_raw) / (2.0 * tan_a);
        const double zMid = 0.5 * (dU.z + dV.z);

        // Acceptance: the crossing must sit on the plane it claims to be on.
        // Surface frame is loc0 = -X, loc1 = +Y (see SeqCalibrator's header).
        if (m_acceptanceTolerance.value() >= 0.0 &&
            !surf->bounds().inside(
                Acts::Vector2(-X_g, Y_g),
                Acts::BoundaryTolerance::AbsoluteEuclidean(
                    m_acceptanceTolerance.value()))) {
          ++nOutside;
          continue;
        }

        // Ghost veto: the crossing must land on a lit scintillator cell of its
        // OWN layer. Fail-open — a layer whose scintillator recorded nothing
        // confirms everything, so an inefficient or absent scint plane can
        // never remove a real hit, it only stops vetoing.
        if (scintWin > 0.0) {
          auto itCells = scintCells.find({dU.station, dU.layer});
          if (itCells != scintCells.end() && !itCells->second.empty()) {
            bool confirmed = false;
            for (const auto& [cx, cy] : itCells->second) {
              if (std::abs(X_g - cx) <= scintWin &&
                  std::abs(Y_g - cy) <= scintWin) { confirmed = true; break; }
            }
            if (!confirmed) { ++nGhostVetoed; continue; }
          }
        }
        edm4hep::CovMatrix3f cov{};
        cov[0] = static_cast<float>(var2D_x);   // xx
        cov[3] = static_cast<float>(var2D_y);   // yy
        auto mhit = output->create();
        mhit.setCellID(dU.cellID);              // station,layer carried via cellID
        mhit.setType(2);
        mhit.setQuality(3);                     // 3 = combined U+V (2D)
        mhit.setTime(0.5f * (dU.time + dV.time));
        mhit.setEDep(dU.edep + dV.edep);
        mhit.setEDepError(0.0f);
        mhit.setPosition(edm4hep::Vector3d{X_g, Y_g, zMid});
        mhit.setCovMatrix(cov);
        ++n2D;
      }
    }

    debug() << "[MTCSciFiMeasConverter] evt=" << evtNum
            << " input=" << hits->size()
            << " method=" << m_pairMethod.value()
            << " emitted 1D=" << n1D
            << " 2D=" << n2D
            << " scintSkipped=" << nScint
            << " ghostVetoed=" << nGhostVetoed
            << " outsideAcceptance=" << nOutside
            << " noSurface=" << nSkipped << endmsg;

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[MTCSciFiMeasConverter] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[MTCSciFiMeasConverter] Unknown exception in execute()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode MTCSciFiMeasConverter::finalize() {
  try {
    m_inputHandle.reset();
    m_outputHandle.reset();
    info() << "[MTCSciFiMeasConverter] Done. Events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[MTCSciFiMeasConverter] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[MTCSciFiMeasConverter] Unknown exception in finalize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(MTCSciFiMeasConverter)
