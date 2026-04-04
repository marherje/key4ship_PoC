// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/DataHandle.h"

// edm4hep input
#include "edm4hep/TrackerHit3DCollection.h"

// edm4hep output
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MutableTrack.h"

// ACTS
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

// ACTS propagator
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"

// ACTS track fitting
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"

// ACTS track containers
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"

// ACTS magnetic field
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

// ACTS calibration context
#include "Acts/Utilities/CalibrationContext.hpp"

// ACTS measurement helpers
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"

// SND geometry service
#include "ISNDGeoSvc.h"

// Standard
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
// SNDMeasurement
// ---------------------------------------------------------------------------

struct SNDMeasurement {
  const Acts::Surface* surface     = nullptr;
  double               localCoord  = 0.0;
  double               localCoord2 = 0.0;
  double               variance    = 0.0;
  double               variance2   = 0.0;
  bool                 is2D        = false;
  int                  detectorID  = -1;
  int                  plane       = -1;
  float                time        = 0.0f;
  float                eDep        = 0.0f;

  SNDMeasurement(const Acts::Surface* sf,
                 double lc, double lc2,
                 double var, double var2,
                 bool twod, int detID, int pl,
                 float t, float e)
      : surface(sf), localCoord(lc), localCoord2(lc2),
        variance(var), variance2(var2), is2D(twod),
        detectorID(detID), plane(pl), time(t), eDep(e) {}
};

// ---------------------------------------------------------------------------
// SNDSourceLink
// ---------------------------------------------------------------------------

struct SNDSourceLink {
  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
  std::size_t index = 0;
  void setGeometryId(Acts::GeometryIdentifier gid) { m_geometryId = gid; }
private:
  Acts::GeometryIdentifier m_geometryId;
};

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

using SNDStepper         = Acts::EigenStepper<>;
using SNDDirectNavigator = Acts::DirectNavigator;
using SNDPropagator      = Acts::Propagator<SNDStepper, SNDDirectNavigator>;

using SNDTrackContainer = Acts::TrackContainer<
    Acts::VectorTrackContainer,
    Acts::VectorMultiTrajectory,
    std::shared_ptr>;

using SNDKalmanFitter = Acts::KalmanFitter<SNDPropagator,
                                           Acts::VectorMultiTrajectory>;

// ---------------------------------------------------------------------------
// ACTSProtoTracker
// ---------------------------------------------------------------------------

class ACTSProtoTracker : public Gaudi::Algorithm {
public:
  ACTSProtoTracker(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  Gaudi::Property<std::string> m_inputSiTarget{
      this, "InputSiTarget", "SiTargetMeasurements",
      "SiTarget TrackerHit3DCollection from SiTargetMeasConverter"};
  Gaudi::Property<std::string> m_inputSiPad{
      this, "InputSiPad", "SiPadMeasurements",
      "SiPad TrackerHit3DCollection from SiPadMeasConverter"};
  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "ACTSTracks",
      "Output edm4hep::TrackCollection name"};
  Gaudi::Property<double> m_bFieldX{this, "BFieldX", 0.0, "BField X [T]"};
  Gaudi::Property<double> m_bFieldY{this, "BFieldY", 0.0, "BField Y [T]"};
  Gaudi::Property<double> m_bFieldZ{this, "BFieldZ", 0.0, "BField Z [T]"};

  // ---- Seed configuration (DD4hep convention: Z=beam) ----
  // SeedPositions: flat list of (x, y, z) triplets in mm.
  //   x = transverse X [mm]
  //   y = transverse Y [mm]
  //   z = beam position [mm] — used only to select starting surface
  // Default: one seed at the beam center.
  Gaudi::Property<std::vector<double>> m_seedPositions{
      this, "SeedPositions", {0.0, 0.0, 0.0},
      "Seed positions in DD4hep convention (Z=beam) as flat (x,y,z) triplets [mm]."};

  // SeedDirections: flat list of (dx, dy, dz) unit vectors.
  //   dz = beam direction component (dominant for forward tracks)
  // The algorithm internally swaps to ACTS convention before fitting.
  // Default: one seed along +Z (beam direction).
  Gaudi::Property<std::vector<double>> m_seedDirections{
      this, "SeedDirections", {0.0, 0.0, 1.0},
      "Seed directions in DD4hep convention (Z=beam) as flat (dx,dy,dz) triplets."};

  // Seed momentum magnitude [GeV], same for all seeds.
  Gaudi::Property<double> m_seedMomentum{
      this, "SeedMomentum", 10.0,
      "Seed momentum magnitude [GeV]."};

  // Enable automatic Hough Transform seeding (ignores SeedPositions/SeedDirections).
  Gaudi::Property<bool> m_autoSeed{
      this, "AutoSeed", false,
      "If true, use Hough Transform automatic seeding instead of manual seeds."};

  // Maximum number of seeds from auto-seeding.
  Gaudi::Property<int> m_maxSeeds{
      this, "MaxSeeds", 5,
      "Maximum number of seeds returned by auto-seeding."};

  // Hough Transform histogram bin size [mm]
  Gaudi::Property<double> m_houghBinSize{
      this, "HoughBinSize", 5.0,
      "Bin size for Hough Transform 2D histogram [mm]. "
      "Should be ~1 SiPad pixel size."};

  // Hough Transform detector half-size [mm] — sets histogram range
  Gaudi::Property<double> m_houghHalfSize{
      this, "HoughHalfSize", 200.0,
      "Half-size of transverse detector for Hough histogram [mm]."};

  // Minimum number of votes for a Hough peak to be considered a seed
  Gaudi::Property<int> m_houghMinVotes{
      this, "HoughMinVotes", 3,
      "Minimum number of hits voting for a Hough peak to create a seed."};

  // Compatibility radius: hit is compatible with seed if within this [mm]
  Gaudi::Property<double> m_seedCompatRadius{
      this, "SeedCompatRadius", 10.0,
      "Radius [mm] within which a hit is considered compatible with a seed."};

  // Strip pitch used for seed position refinement [mm]
  Gaudi::Property<double> m_stripPitch{
      this, "SeedStripPitch", 0.0755,
      "Strip pitch [mm] used for most-frequent-strip seed refinement."};

  Gaudi::Property<double> m_maxChi2PerMeas{
    this, "MaxChi2PerMeas", 500.0,
    "Maximum chi2/nMeas threshold for track acceptance. "
    "Tracks above this threshold are rejected as false seeds."};
      
  // ---- Hough-based auto-seeding -------------------------------------------
  struct SeedCandidate {
    double x;          // transverse X [mm] (DD4hep convention)
    double y;          // transverse Y [mm]
    double z_start;    // beam Z of first compatible layer [mm]
    int    nVotes;     // number of hits supporting this seed
  };

  std::vector<SeedCandidate> findSeeds(
      const std::vector<SNDMeasurement>& measurements,
      const Acts::GeometryContext& gctx) const;

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siTargetHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siPadHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackCollection>>
      m_outputHandle;

  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  mutable std::atomic<long long> m_eventCount{0};
  mutable Acts::MagneticFieldContext m_mctx;
  mutable Acts::CalibrationContext   m_cctx;

  // Dynamically computed in initialize(): number of SiTarget surfaces.
  // SiPad surfaces start at this index in allSurfaces().
  mutable std::size_t m_nSiTargetSurfaces{0};

};

// ---------------------------------------------------------------------------
// findSeeds() — Hough Transform automatic seeder
// ---------------------------------------------------------------------------

std::vector<ACTSProtoTracker::SeedCandidate> ACTSProtoTracker::findSeeds(
    const std::vector<SNDMeasurement>& measurements,
    const Acts::GeometryContext& gctx) const
{
  std::vector<SeedCandidate> seeds;
  if (measurements.empty()) return seeds;

  // =========================================================================
  // STEP 1: Build list of unambiguous 2D points
  // =========================================================================
  // Source A: SiPad hits (directly 2D)
  // Source B: SiTarget strip crossings (StripX × StripY per station)

  struct Point2D {
    double x;      // transverse X [mm]
    double y;      // transverse Y [mm]
    double z;      // beam Z [mm] (ACTS X coordinate of surface)
    int    weight; // 2 for crossings (two hits used), 1 for SiPad
  };
  std::vector<Point2D> points2D;

  // --- Source A: SiPad 2D hits ---
  for (const auto& m : measurements) {
    if (!m.is2D) continue;
    double beamZ = m.surface->center(gctx).x();  // ACTS X = DD4hep beam Z
    points2D.push_back({m.localCoord, m.localCoord2, beamZ, 1});
  }

  // --- Source B: SiTarget strip crossings ---
  // Group SiTarget hits by surface beam-Z position (rounded to 1mm).
  // StripX (plane=0) and StripY (plane=1) are at slightly different Z
  // but belong to the same station. Group by station (every 2 consecutive
  // surfaces in SiTarget share a station).
  //
  // Strategy: group by surface geometry ID prefix — surfaces in the same
  // station differ only in plane (0 vs 1). Since surfaces are sorted by
  // beam Z, consecutive surfaces at similar Z form a station pair.
  // Use beam Z tolerance of 10mm to pair StripX and StripY planes.

  const double stationTolerance = 10.0;  // mm — max Z gap within a station

  // Collect SiTarget measurements grouped by approximate beam Z
  // (round to nearest 10mm to group station pairs)
  std::map<int, std::vector<const SNDMeasurement*>> stationGroups;
  for (const auto& m : measurements) {
    if (m.is2D) continue;  // skip SiPad
    double beamZ = m.surface->center(gctx).x();
    // Group key: round to nearest 10mm
    int key = static_cast<int>(std::round(beamZ / stationTolerance));
    stationGroups[key].push_back(&m);
  }

  // For each station group, form all StripX × StripY crossings
  for (const auto& [key, stationMeas] : stationGroups) {
    // Separate into plane=0 (StripX) and plane=1 (StripY)
    std::vector<const SNDMeasurement*> stripX, stripY;
    for (const auto* m : stationMeas) {
      if (m->plane == 0) stripX.push_back(m);
      else if (m->plane == 1) stripY.push_back(m);
    }

    if (stripX.empty() || stripY.empty()) continue;

    // Form all (StripX, StripY) crossing pairs
    // localCoord of StripX = X position [mm]
    // localCoord of StripY = Y position [mm]
    for (const auto* mx : stripX) {
      for (const auto* my : stripY) {
        double crossX = mx->localCoord;  // X from StripX
        double crossY = my->localCoord;  // Y from StripY
        // Use Z midpoint of the two planes
        double zX = mx->surface->center(gctx).x();
        double zY = my->surface->center(gctx).x();
        double crossZ = 0.5 * (zX + zY);
        points2D.push_back({crossX, crossY, crossZ, 2});
      }
    }
  }

  if (points2D.empty()) {
    SeedCandidate sc;
    sc.x = 0.0; sc.y = 0.0; sc.nVotes = 0;
    sc.z_start = measurements.empty() ? -370.0
               : measurements.front().surface->center(gctx).x();
    seeds.push_back(sc);
    return seeds;
  }

  // =========================================================================
  // STEP 2: Hough Transform on 2D points only
  // =========================================================================
  // Each point2D votes for one (x,y) bin — no strip ambiguity!

  const double halfSize = m_houghHalfSize.value();
  const double binSize  = m_houghBinSize.value();
  const int    nBins    = static_cast<int>(2.0 * halfSize / binSize) + 1;

  std::vector<std::vector<int>> histo(nBins, std::vector<int>(nBins, 0));

  auto toBin = [&](double coord) -> int {
    int bin = static_cast<int>((coord + halfSize) / binSize);
    return std::max(0, std::min(nBins - 1, bin));
  };
  auto fromBin = [&](int bin) -> double {
    return -halfSize + (bin + 0.5) * binSize;
  };

  for (const auto& p : points2D) {
    int ix = toBin(p.x);
    int iy = toBin(p.y);
    histo[ix][iy] += p.weight;
  }

  // =========================================================================
  // STEP 3: Find peaks with non-maximum suppression
  // =========================================================================

  const int minVotes    = m_houghMinVotes.value();
  const int suppressRad = static_cast<int>(
      std::ceil(15.0 / binSize));  // suppress within 15mm

  struct Peak { int ix; int iy; int votes; };
  std::vector<Peak> peaks;

  for (int ix = 0; ix < nBins; ++ix) {
    for (int iy = 0; iy < nBins; ++iy) {
      int v = histo[ix][iy];
      if (v < minVotes) continue;
      bool isMax = true;
      for (int dx = -1; dx <= 1 && isMax; ++dx) {
        for (int dy = -1; dy <= 1 && isMax; ++dy) {
          if (dx == 0 && dy == 0) continue;
          int nx = ix + dx, ny = iy + dy;
          if (nx < 0 || nx >= nBins || ny < 0 || ny >= nBins) continue;
          if (histo[nx][ny] > v) isMax = false;
        }
      }
      if (isMax) peaks.push_back({ix, iy, v});
    }
  }

  std::sort(peaks.begin(), peaks.end(),
            [](const Peak& a, const Peak& b) { return a.votes > b.votes; });

  std::vector<bool> suppressed(peaks.size(), false);
  for (std::size_t i = 0; i < peaks.size(); ++i) {
    if (suppressed[i]) continue;
    for (std::size_t j = i + 1; j < peaks.size(); ++j) {
      if (suppressed[j]) continue;
      int dx = std::abs(peaks[i].ix - peaks[j].ix);
      int dy = std::abs(peaks[i].iy - peaks[j].iy);
      if (dx <= suppressRad && dy <= suppressRad)
        suppressed[j] = true;
    }
  }

  // =========================================================================
  // STEP 4: Refine seed position using most-frequent strip
  // =========================================================================
  // The centroid of compatible 2D points is biased by ghost crossings.
  // Instead, find the most frequently occurring strip value (rounded to
  // strip pitch) in X and Y separately. The real muon strip dominates
  // in frequency — ghost crossings scatter to different strip positions.
  // This gives strip-level precision (~75μm) for the seed position.

  const double compatR  = m_seedCompatRadius.value();
  const double compatR2 = compatR * compatR;
  const int    maxS     = m_maxSeeds.value();

  // Strip pitch for frequency binning [mm]
  const double stripPitch = m_stripPitch.value();

  for (std::size_t pi = 0; pi < peaks.size() && (int)seeds.size() < maxS; ++pi) {
    if (suppressed[pi]) continue;

    const double peakX = fromBin(peaks[pi].ix);
    const double peakY = fromBin(peaks[pi].iy);

    // Collect compatible 2D points and find most frequent X and Y strips
    // Use std::map with rounded strip index as key
    std::map<int, int> xStripFreq;  // strip_index → total weight
    std::map<int, int> yStripFreq;
    double firstZ = std::numeric_limits<double>::max();
    int    nPts   = 0;

    for (const auto& p : points2D) {
      double dx = p.x - peakX;
      double dy = p.y - peakY;
      if (dx*dx + dy*dy < compatR2) {
        // Round to nearest strip index
        int xStrip = static_cast<int>(std::round(p.x / stripPitch));
        int yStrip = static_cast<int>(std::round(p.y / stripPitch));
        xStripFreq[xStrip] += p.weight;
        yStripFreq[yStrip] += p.weight;
        nPts += p.weight;
        if (p.z < firstZ) firstZ = p.z;
      }
    }

    if (nPts == 0) continue;

    // Find most frequent X strip
    double refinedX = peakX;
    int maxXFreq = 0;
    for (const auto& [strip, freq] : xStripFreq) {
      if (freq > maxXFreq) {
        maxXFreq = freq;
        refinedX = strip * stripPitch;
      }
    }

    // Find most frequent Y strip
    double refinedY = peakY;
    int maxYFreq = 0;
    for (const auto& [strip, freq] : yStripFreq) {
      if (freq > maxYFreq) {
        maxYFreq = freq;
        refinedY = strip * stripPitch;
      }
    }

    SeedCandidate sc;
    sc.x       = refinedX;
    sc.y       = refinedY;
    sc.z_start = firstZ;
    sc.nVotes  = peaks[pi].votes;
    seeds.push_back(sc);
  }

  return seeds;
}

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    m_siTargetHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiTarget.value(), Gaudi::DataHandle::Reader, this);
    m_siPadHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiPad.value(), Gaudi::DataHandle::Reader, this);
    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[ACTSProtoTracker] Failed to retrieve ACTSGeoSvc." << endmsg;
      return StatusCode::FAILURE;
    }

    // Compute the number of SiTarget surfaces dynamically.
    // Find the largest gap in X between consecutive surfaces —
    // that gap separates SiTarget from SiPad.
    const auto& allSurf = m_geoSvc->allSurfaces();
    const auto& gctx    = m_geoSvc->geometryContext();

    if (allSurf.size() < 2) {
      error() << "[ACTSProtoTracker] Too few surfaces in geometry." << endmsg;
      return StatusCode::FAILURE;
    }

    double maxGap      = 0.0;
    std::size_t gapIdx = 0;
    for (std::size_t i = 1; i < allSurf.size(); ++i) {
      // Surfaces are positioned along X in ACTS (beam coord swapped X<->Z)
      double gap = allSurf[i]->center(gctx).x() -
                   allSurf[i-1]->center(gctx).x();
      if (gap > maxGap) {
        maxGap  = gap;
        gapIdx  = i;
      }
    }
    m_nSiTargetSurfaces = gapIdx;

    info() << "[ACTSProtoTracker] Initialized. GeoSvc has "
           << allSurf.size() << " surfaces. "
           << "SiTarget surfaces: " << m_nSiTargetSurfaces
           << " SiPad surfaces: " << (allSurf.size() - m_nSiTargetSurfaces)
           << " (gap=" << maxGap << " mm)" << endmsg;

    // Validate seed properties
    const std::size_t nPosVals = m_seedPositions.value().size();
    const std::size_t nDirVals = m_seedDirections.value().size();
    if (nPosVals % 3 != 0 || nDirVals % 3 != 0) {
      error() << "[ACTSProtoTracker] SeedPositions and SeedDirections must "
              << "have sizes that are multiples of 3." << endmsg;
      return StatusCode::FAILURE;
    }
    if (nPosVals / 3 != nDirVals / 3) {
      error() << "[ACTSProtoTracker] SeedPositions and SeedDirections must "
              << "have the same number of (x,y,z) triplets." << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "[ACTSProtoTracker] Configured with "
           << nPosVals / 3 << " seed(s)." << endmsg;

    return sc;
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in initialize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    auto* output = m_outputHandle->createAndPut();

    const auto& allSurfaces = m_geoSvc->allSurfaces();
    const auto& gctx        = m_geoSvc->geometryContext();

    // =========================================================================
    // STEP 1: Collect measurements
    // =========================================================================
    std::vector<SNDMeasurement> measurements;

    // ---- SiTarget (1D strips) -----------------------------------------------
    const auto* stHits = m_siTargetHandle->get();
    if (stHits) {
      for (std::size_t i = 0; i < stHits->size(); ++i) {
        const auto& hit   = (*stHits)[i];
        const int   plane = hit.getQuality();
        const auto& pos   = hit.getPosition();

        const Acts::Surface* bestSurface = nullptr;
        double bestDist = 1e9;
        for (std::size_t j = 0; j < m_nSiTargetSurfaces; ++j) {
          double dz = std::abs(allSurfaces[j]->center(gctx).x() - pos.z);
          if (dz < bestDist) {
            bestDist    = dz;
            bestSurface = allSurfaces[j];
          }
        }

        if (!bestSurface) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " SiTarget hit " << i << " has no matching surface."
                    << endmsg;
          continue;
        }

        const auto& cov   = hit.getCovMatrix();
        // Z=beam: plane=0 measures X, plane=1 measures Y
        double localCoord = (plane == 0) ? pos.x : pos.y;
        double var        = (plane == 0) ? cov[0] : cov[3];
        // cov[0]=xx for plane=0 (StripX), cov[3]=yy for plane=1 (StripY)

        measurements.emplace_back(bestSurface, localCoord, 0.0,
                                   var, 0.0, false, 0, plane,
                                   hit.getTime(), hit.getEDep());
      }
    }

    // ---- SiPad (2D pixels) --------------------------------------------------
    const auto* spHits = m_siPadHandle->get();
    if (spHits) {
      for (std::size_t i = 0; i < spHits->size(); ++i) {
        const auto& hit = (*spHits)[i];
        const auto& pos = hit.getPosition();

        const Acts::Surface* bestSurface = nullptr;
        double bestDist = 1e9;
        for (std::size_t j = m_nSiTargetSurfaces; j < allSurfaces.size(); ++j) {
          double dz = std::abs(allSurfaces[j]->center(gctx).x() - pos.z);
          if (dz < bestDist) {
            bestDist    = dz;
            bestSurface = allSurfaces[j];
          }
        }

        if (!bestSurface) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " SiPad hit " << i << " has no matching surface."
                    << endmsg;
          continue;
        }

        const auto& cov = hit.getCovMatrix();
        // Z=beam: SiPad measures X and Y (transverse coordinates)
        measurements.emplace_back(bestSurface, pos.x, pos.y,
                                   cov[0], cov[3], true, 1, -1,
                                   hit.getTime(), hit.getEDep());
      }
    }

    if (measurements.empty()) {
      debug() << "[ACTSProtoTracker] evt=" << evtNum
              << " no measurements." << endmsg;
      return StatusCode::SUCCESS;
    }

    std::sort(measurements.begin(), measurements.end(),
              [&](const SNDMeasurement& a, const SNDMeasurement& b) {
                return a.surface->center(gctx).x() <
                       b.surface->center(gctx).x();
              });

    debug() << "[ACTSProtoTracker] evt=" << evtNum
            << " SiTarget=" << (stHits ? stHits->size() : 0)
            << " SiPad=" << (spHits ? spHits->size() : 0)
            << " total measurements=" << measurements.size() << endmsg;

    // =========================================================================
    // STEP 4: Build shared KF components (once per event)
    // =========================================================================

    auto bField = std::make_shared<Acts::ConstantBField>(
        Acts::Vector3(m_bFieldX.value() * Acts::UnitConstants::T,
                      m_bFieldY.value() * Acts::UnitConstants::T,
                      m_bFieldZ.value() * Acts::UnitConstants::T));

    // ---- SurfaceAccessor: maps SourceLink → Surface -------------------------
    struct SNDSurfaceAccessor {
      const std::vector<SNDMeasurement>* meas = nullptr;
      const Acts::Surface* operator()(const Acts::SourceLink& sl) const {
        return (*meas)[sl.get<SNDSourceLink>().index].surface;
      }
    };

    // ---- Calibrator ---------------------------------------------------------
    struct SNDCalibrator {
      const std::vector<SNDMeasurement>* meas = nullptr;
      void operator()(const Acts::GeometryContext& /*gctx*/,
                      const Acts::CalibrationContext& /*cctx*/,
                      const Acts::SourceLink& sl,
                      Acts::VectorMultiTrajectory::TrackStateProxy ts) const {
        const auto& ssl = sl.get<SNDSourceLink>();
        const auto& m   = (*meas)[ssl.index];
        if (m.is2D) {
          constexpr std::array<Acts::BoundIndices, 2> indices = {
              Acts::eBoundLoc0, Acts::eBoundLoc1};
          ts.allocateCalibrated(2);
          ts.template calibrated<2>() =
              Acts::ActsVector<2>(m.localCoord, m.localCoord2);
          ts.template calibratedCovariance<2>() =
              Acts::ActsSquareMatrix<2>{{m.variance, 0.0},
                                        {0.0, m.variance2}};
          ts.setProjectorSubspaceIndices(indices);
        } else {
          // 1D SiTarget strip measurement.
          // plane=0 (StripX): DD4hep X → global Z → eBoundLoc0
          // plane=1 (StripY): DD4hep Y → global Y → eBoundLoc1
          ts.allocateCalibrated(1);
          ts.template calibrated<1>() =
              Acts::ActsVector<1>(m.localCoord);
          ts.template calibratedCovariance<1>() =
              Acts::ActsSquareMatrix<1>{{m.variance}};
          if (m.plane == 1) {
            constexpr std::array<Acts::BoundIndices, 1> indices = {
                Acts::eBoundLoc1};
            ts.setProjectorSubspaceIndices(indices);
          } else {
            constexpr std::array<Acts::BoundIndices, 1> indices = {
                Acts::eBoundLoc0};
            ts.setProjectorSubspaceIndices(indices);
          }
        }
      }
    };

    using KFExtensions = Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory>;
    Acts::GainMatrixUpdater gainMatrixUpdater;
    SNDSurfaceAccessor surfaceAccessor;
    surfaceAccessor.meas = &measurements;
    SNDCalibrator calibrator;
    calibrator.meas = &measurements;

    KFExtensions extensions;
    extensions.updater
        .connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
            &gainMatrixUpdater);
    extensions.calibrator
        .template connect<&SNDCalibrator::operator()>(&calibrator);
    extensions.surfaceAccessor
        .template connect<&SNDSurfaceAccessor::operator()>(&surfaceAccessor);

    // Shared seed covariance (loose — same for all seeds)
    Acts::BoundSquareMatrix seedCov = Acts::BoundSquareMatrix::Zero();
    seedCov(Acts::eBoundLoc0,   Acts::eBoundLoc0)   = 1e6;
    seedCov(Acts::eBoundLoc1,   Acts::eBoundLoc1)   = 1e6;
    seedCov(Acts::eBoundPhi,    Acts::eBoundPhi)     = 1.0;
    seedCov(Acts::eBoundTheta,  Acts::eBoundTheta)   = 1.0;
    seedCov(Acts::eBoundQOverP, Acts::eBoundQOverP)  = 10.0;
    seedCov(Acts::eBoundTime,   Acts::eBoundTime)    = 1e9;

    const double seedQoverP =
        1.0 / (m_seedMomentum.value() * Acts::UnitConstants::GeV);

    // ---- Build seed list (auto or manual) -----------------------------------
    // Each entry is (dd_x, dd_y, dd_z) in DD4hep convention.
    std::vector<std::array<double, 3>> seedPositions3D;
    std::vector<std::array<double, 3>> seedDirections3D;

    if (m_autoSeed.value()) {
      // Hough Transform seeding using all measurements
      auto autoSeeds = findSeeds(measurements, gctx);

      for (const auto& sc : autoSeeds) {
        seedPositions3D.push_back({sc.x, sc.y, sc.z_start});
        seedDirections3D.push_back({0.0, 0.0, 1.0});
      }

      if (evtNum < 3) {
        debug() << "[ACTSProtoTracker] evt=" << evtNum
                << " Hough seeding: " << autoSeeds.size()
                << " seed(s) found" << endmsg;
        for (std::size_t si = 0; si < autoSeeds.size(); ++si) {
          debug() << "[ACTSProtoTracker]   seed " << si
                  << ": x=" << autoSeeds[si].x
                  << " y=" << autoSeeds[si].y
                  << " z_start=" << autoSeeds[si].z_start
                  << " votes=" << autoSeeds[si].nVotes << endmsg;
        }
      }

      if (autoSeeds.empty()) {
        // Fallback: single seed at first surface
        double zFirst = allSurfaces.empty() ? -370.0
                      : allSurfaces.front()->center(gctx).x();
        seedPositions3D.push_back({0.0, 0.0, zFirst});
        seedDirections3D.push_back({0.0, 0.0, 1.0});
        if (evtNum < 3) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " Hough seeding found no seeds, using fallback." << endmsg;
        }
      }
    } else {
      // Manual seeding from Gaudi properties
      const auto& seedPosCfg = m_seedPositions.value();
      const auto& seedDirCfg = m_seedDirections.value();
      const std::size_t nManual = seedPosCfg.size() / 3;
      for (std::size_t i = 0; i < nManual; ++i) {
        seedPositions3D.push_back({seedPosCfg[3*i+0],
                                   seedPosCfg[3*i+1],
                                   seedPosCfg[3*i+2]});
        seedDirections3D.push_back({seedDirCfg[3*i+0],
                                    seedDirCfg[3*i+1],
                                    seedDirCfg[3*i+2]});
      }
    }

    const std::size_t nSeeds = seedPositions3D.size();

    std::size_t nTracks = 0;

    // Fingerprints store sets of hit indices (not surface geoIDs) for
    // accurate duplicate detection across seeds with different hit selection.
    std::vector<std::set<std::size_t>> acceptedFingerprints;
    const double kDuplicateOverlapFraction = 0.7;  // reject if >70% shared hits

    // =========================================================================
    // STEPS 5-7: Loop over seeds
    // =========================================================================
    for (std::size_t iSeed = 0; iSeed < nSeeds; ++iSeed) {

      // ---- Build fresh propagator for this seed ----------------------------
      // (propagator is moved into KalmanFitter, so rebuild each iteration)
      SNDStepper        stepper_i(bField);
      SNDDirectNavigator navigator_i;
      SNDPropagator propagator_i(std::move(stepper_i),
                                  std::move(navigator_i));

      // ---- Seed position (DD4hep → ACTS coordinate swap) -------------------
      // DD4hep convention: x=transverse X, y=transverse Y, z=beam.
      const double dd_x = seedPositions3D[iSeed][0];  // transverse X [mm]
      const double dd_y = seedPositions3D[iSeed][1];  // transverse Y [mm]
      // dd_z (beam) is used to find the starting surface
      const double dd_z = seedPositions3D[iSeed][2];  // beam Z [mm]

      // Always use the first surface as the KF seed reference surface.
      // The Hough z_start is used only to inform hit selection, not the KF seed.
      // Starting from the first surface is safe — the KF handles empty layers
      // (holes) gracefully with loose covariance.
      const Acts::Surface* sfSeed = allSurfaces.front();

      Acts::Vector4 seedPos4;
      // Beam coordinate goes into ePos0 (ACTS X, after X<->Z swap in geometry).
      // Surface rotation rot90Y maps: local X = global Z, local Y = global Y.
      // Therefore: eBoundLoc0 = global Z, eBoundLoc1 = global Y.
      // BoundTrackParameters maps: ePos1 → eBoundLoc1, ePos2 → eBoundLoc0.
      // DD4hep convention: dd_x = transverse X (→ global Z → eBoundLoc0 → ePos2)
      //                    dd_y = transverse Y (→ global Y → eBoundLoc1 → ePos1)
      seedPos4[Acts::ePos0] = sfSeed->center(gctx).x();  // beam coord (ACTS X)
      seedPos4[Acts::ePos1] = dd_y;  // DD4hep Y → global Y → eBoundLoc1
      seedPos4[Acts::ePos2] = dd_x;  // DD4hep X → global Z → eBoundLoc0
      seedPos4[Acts::eTime] = 0.0;

      // ---- Seed direction (DD4hep → ACTS coordinate swap) ------------------
      // DD4hep convention: (dx, dy, dz), Z=beam.
      // ACTS swap: ACTS_x = dz (beam), ACTS_y = dx, ACTS_z = dy
      const double dd_dx = seedDirections3D[iSeed][0];
      const double dd_dy = seedDirections3D[iSeed][1];
      const double dd_dz = seedDirections3D[iSeed][2];

      Acts::Vector3 seedDir(dd_dz, dd_dx, dd_dy);
      if (seedDir.norm() < 1e-6) seedDir = Acts::Vector3(1.0, 0.001, 0.001);
      seedDir = seedDir.normalized();

      // ---- Seed-dependent hit selection ------------------------------------
      // Select the best hit per surface relative to THIS seed's transverse
      // position (dd_x, dd_y in DD4hep convention).
      // For 1D measurements: minimize |localCoord - seed_ref|
      //   plane=0 (StripX): seed_ref = dd_x
      //   plane=1 (StripY): seed_ref = dd_y
      // For 2D measurements (SiPad): minimize 2D distance to (dd_x, dd_y)
      const double seed_ref_x = dd_x;  // transverse X reference [mm]
      const double seed_ref_y = dd_y;  // transverse Y reference [mm]

      std::unordered_map<Acts::GeometryIdentifier, std::size_t> bestHitPerSurface;
      for (std::size_t i = 0; i < measurements.size(); ++i) {
        const auto& m   = measurements[i];
        auto        gid = m.surface->geometryId();

        double dist = 0.0;
        if (m.is2D) {
          double dx = m.localCoord  - seed_ref_x;
          double dy = m.localCoord2 - seed_ref_y;
          dist = std::sqrt(dx * dx + dy * dy);
        } else {
          if (m.plane == 0)
            dist = std::abs(m.localCoord - seed_ref_x);
          else
            dist = std::abs(m.localCoord - seed_ref_y);
        }

        auto it = bestHitPerSurface.find(gid);
        if (it == bestHitPerSurface.end()) {
          bestHitPerSurface[gid] = i;
        } else {
          const auto& mExist = measurements[it->second];
          double existDist = 0.0;
          if (mExist.is2D) {
            double dx = mExist.localCoord  - seed_ref_x;
            double dy = mExist.localCoord2 - seed_ref_y;
            existDist = std::sqrt(dx * dx + dy * dy);
          } else {
            if (mExist.plane == 0)
              existDist = std::abs(mExist.localCoord - seed_ref_x);
            else
              existDist = std::abs(mExist.localCoord - seed_ref_y);
          }
          if (dist < existDist) it->second = i;
        }
      }

      // Build source link list for this seed (one per surface, sorted by X)
      std::vector<Acts::SourceLink> allSourceLinks;
      allSourceLinks.reserve(bestHitPerSurface.size());
      for (const auto& [gid, idx] : bestHitPerSurface) {
        SNDSourceLink ssl;
        ssl.index = idx;
        ssl.setGeometryId(gid);
        allSourceLinks.emplace_back(ssl);
      }
      std::sort(allSourceLinks.begin(), allSourceLinks.end(),
                [&](const Acts::SourceLink& a, const Acts::SourceLink& b) {
                  const auto& ma = measurements[a.get<SNDSourceLink>().index];
                  const auto& mb = measurements[b.get<SNDSourceLink>().index];
                  return ma.surface->center(gctx).x() <
                         mb.surface->center(gctx).x();
                });

      // ---- Create seed parameters -------------------------------------------
      auto seedParamsResult = Acts::BoundTrackParameters::create(
          gctx, sfSeed->getSharedPtr(), seedPos4, seedDir, seedQoverP,
          seedCov, Acts::ParticleHypothesis::muon());

      if (!seedParamsResult.ok()) {
        warning() << "[ACTSProtoTracker] evt=" << evtNum
                  << " seed=" << iSeed
                  << " seed creation failed: "
                  << seedParamsResult.error() << endmsg;
        continue;
      }
      const auto& seedParams = *seedParamsResult;

      // ---- KF options -------------------------------------------------------
      Acts::PropagatorPlainOptions pOptions(gctx, m_mctx);
      pOptions.direction = Acts::Direction::Forward();
      pOptions.stepping.maxStepSize = 10.0;
      pOptions.maxSteps = 10000;

      Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(
          gctx, m_mctx, std::cref(m_cctx),
          extensions,
          pOptions,
          nullptr,
          false,
          false);

      // ---- Run KF -----------------------------------------------------------
      auto trackBackend = std::make_shared<Acts::VectorTrackContainer>();
      auto trajBackend  = std::make_shared<Acts::VectorMultiTrajectory>();
      SNDTrackContainer kfTracks(trackBackend, trajBackend);

      SNDKalmanFitter kf(std::move(propagator_i),
                         Acts::getDefaultLogger("KF", Acts::Logging::WARNING));

      // Always use the full surface sequence.
      // The KF will mark surfaces before the first hit as holes — this is fine.
      auto kfResult = kf.fit(
          allSourceLinks.begin(), allSourceLinks.end(),
          seedParams,
          kfOptions,
          allSurfaces,        // full surface sequence, always from first surface
          kfTracks);

      if (!kfResult.ok()) {
        warning() << "[ACTSProtoTracker] evt=" << evtNum
                  << " seed=" << iSeed
                  << " KF failed: " << kfResult.error() << endmsg;
        continue;
      }

      const auto& kfTrack = *kfResult;
      const std::size_t nMeas = kfTrack.nMeasurements();

      debug() << "[ACTSProtoTracker] evt=" << evtNum
              << " seed=" << iSeed
              << " KF nMeas=" << nMeas
              << " nHoles=" << kfTrack.nHoles()
              << " chi2=" << kfTrack.chi2() << endmsg;

      // ---- Write output if track is good ------------------------------------
      if (nMeas >= 3) {
        // ---- Duplicate rejection -------------------------------------------
        // Build fingerprint: set of measurement (hit) indices used by this fit.
        // Two tracks are duplicates if they share too many of the same hits,
        // not just the same surfaces (which would always be 100% overlap since
        // both seeds visit all surfaces).
        std::set<std::size_t> fingerprint;
        for (const auto& [gid, idx] : bestHitPerSurface) {
          fingerprint.insert(idx);
        }

        // Check overlap with already accepted tracks.
        bool isDuplicate = false;
        for (const auto& accepted : acceptedFingerprints) {
          std::size_t nShared = 0;
          for (const auto& idx : fingerprint) {
            if (accepted.count(idx)) ++nShared;
          }
          const double smaller = static_cast<double>(
              std::min(fingerprint.size(), accepted.size()));
          const double overlapFraction = (smaller > 0) ? nShared / smaller : 0.0;

          if (overlapFraction > kDuplicateOverlapFraction) {
            isDuplicate = true;
            debug() << "[ACTSProtoTracker] evt=" << evtNum
                    << " seed=" << iSeed
                    << " rejected as duplicate (hit overlap="
                    << overlapFraction << ")" << endmsg;
            break;
          }
        }

        if (isDuplicate) continue;  // skip to next seed

        const double chi2PerMeas = kfTrack.chi2() / std::max(1.0, (double)nMeas);
        if (chi2PerMeas > m_maxChi2PerMeas.value()) {
          debug() << "[ACTSProtoTracker] evt=" << evtNum
          << " seed=" << iSeed
          << " rejected: chi2/nMeas=" << chi2PerMeas << endmsg;
          continue;
        }

        // Accept this track
        acceptedFingerprints.push_back(fingerprint);

        // ---- Write output --------------------------------------------------
        auto track = output->create();
        track.setType(1);
        track.setChi2(static_cast<float>(kfTrack.chi2()));
        track.setNdf(static_cast<int>(kfTrack.nDoF()));

        // Store seed transverse position as the first TrackState (location=AtIP).
        // D0 = seed transverse X [mm] (DD4hep convention)
        // Z0 = seed transverse Y [mm] (DD4hep convention)
        // This allows the event display to reconstruct per-track hit assignment.
        {
          edm4hep::TrackState seedState{};
          seedState.location = edm4hep::TrackState::AtIP;
          seedState.D0       = static_cast<float>(dd_x);
          seedState.Z0       = static_cast<float>(dd_y);
          track.addToTrackStates(seedState);
        }

        try {
          auto tipIdx = kfTrack.tipIndex();
          auto& mutableTraj = kfTracks.trackStateContainer();
          while (true) {
            auto ts = mutableTraj.getTrackState(tipIdx);
            if (ts.hasCalibrated()) {
              edm4hep::TrackState edm4ts;
              edm4ts.location  = edm4hep::TrackState::AtOther;
              edm4ts.D0        = 0.0f;
              edm4ts.Z0        = 0.0f;
              if (ts.hasSmoothed()) {
                edm4ts.phi       = static_cast<float>(
                    ts.smoothed()[Acts::eBoundPhi]);
                edm4ts.tanLambda = static_cast<float>(
                    std::tan(M_PI / 2.0 - ts.smoothed()[Acts::eBoundTheta]));
                edm4ts.omega     = static_cast<float>(
                    ts.smoothed()[Acts::eBoundQOverP]);
              } else if (ts.hasFiltered()) {
                edm4ts.phi       = static_cast<float>(
                    ts.filtered()[Acts::eBoundPhi]);
                edm4ts.tanLambda = 0.0f;
                edm4ts.omega     = 0.0f;
              } else {
                edm4ts.phi       = 0.0f;
                edm4ts.tanLambda = 0.0f;
                edm4ts.omega     = 0.0f;
              }
              track.addToTrackStates(edm4ts);
            }
            if (!ts.hasPrevious()) break;
            tipIdx = ts.previous();
          }
        } catch (const std::exception& e) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " seed=" << iSeed
                    << " trackStates iteration failed: " << e.what() << endmsg;
        }
        ++nTracks;

      }

    }  // end loop over seeds

    info() << "[ACTSProtoTracker] evt=" << evtNum
           << " measurements=" << measurements.size()
           << " KF tracks=" << nTracks << endmsg;

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in execute()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::finalize() {
  try {
    m_siTargetHandle.reset();
    m_siPadHandle.reset();
    m_outputHandle.reset();
    info() << "[ACTSProtoTracker] Done. Total events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in finalize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(ACTSProtoTracker)