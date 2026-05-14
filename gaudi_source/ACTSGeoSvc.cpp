#include "ACTSGeoSvc.h"
#include "DD4hep/Detector.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <string>

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

ACTSGeoSvc::ACTSGeoSvc(const std::string& name, ISvcLocator* svcLoc)
    : extends(name, svcLoc) {}

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode ACTSGeoSvc::initialize() {
  StatusCode sc = Service::initialize();
  if (sc.isFailure()) return sc;

  if (m_compactFile.value().empty()) {
    error() << "[ACTSGeoSvc] CompactFile property must be set." << endmsg;
    return StatusCode::FAILURE;
  }

  // Load DD4hep geometry
  auto& desc = dd4hep::Detector::getInstance();
  desc.fromXML(m_compactFile.value());
  info() << "[ACTSGeoSvc] DD4hep geometry loaded from: "
         << m_compactFile.value() << endmsg;

  // Struct holding extracted plane info.
  // TGeo positions are in cm; we convert to mm (* 10) immediately.
  struct PlaneInfo {
    double z;           // global Z center [mm]  // beam axis is now Z
    double halfX;       // half-size in X [mm]   // transverse
    double halfY;       // half-size in Y [mm]   // transverse
    double thickness;   // half-thickness in Z [mm]
    int    plane;       // 0=U/StripX, 1=V/StripY, -1=pixel (SiPad)
    int    detID;       // 0=SiTarget, 1=SiPad, 2=MTC
    int    station;     // station index; -1 for non-MTC
    int    layerInDet;  // layer index within station (MTC) or detector
  };

  // Walk TGeo tree and collect sensitive volumes for the requested detector.
  auto extractPlanes = [&](const std::string& detName)
      -> std::vector<PlaneInfo> {
    std::vector<PlaneInfo> planes;
    TGeoManager* mgr = gGeoManager;
    if (!mgr) return planes;

    std::function<void(TGeoNode*, const std::string&)> walk =
        [&](TGeoNode* node, const std::string& path) {
      TGeoVolume* vol = node->GetVolume();
      if (!vol) return;

      const std::string nodeName = node->GetName();
      const std::string volName  = vol->GetName();
      const std::string curPath  = path + "/" + nodeName;

      bool inDetector  = (curPath.find(detName) != std::string::npos);
      bool isSensitive = false;
      if (detName == "SiTarget") {
        isSensitive = (volName.find("SiTarget_layer_") != std::string::npos) &&
                (volName.find("_slice_2") != std::string::npos ||
                 volName.find("_slice_4") != std::string::npos);
      } else if (detName == "SiPad") {
        isSensitive = (volName.find("SiPad_layer_") != std::string::npos) &&
                (volName.find("_slice_4") != std::string::npos);
      } else if (detName == "MTC") {
        isSensitive = (volName.find("_slice_2") != std::string::npos ||
                       volName.find("_slice_4") != std::string::npos);
      }

            bool isSliceContainer = (volName.find("_col") == std::string::npos);
      if (inDetector && isSensitive && isSliceContainer) {
        if (mgr->cd(curPath.c_str())) {
          TGeoMatrix* mat  = mgr->GetCurrentMatrix();
          const double* tr = mat->GetTranslation();
          TGeoBBox* box    = dynamic_cast<TGeoBBox*>(vol->GetShape());
          if (box) {
            PlaneInfo pi;
            pi.z         = tr[2]        * 10.0;
            pi.halfX     = box->GetDX() * 10.0;
            pi.halfY     = box->GetDY() * 10.0;
            pi.thickness = box->GetDZ() * 10.0;
            pi.plane      = -1;
            pi.detID      = -1;
            pi.station    = -1;
            pi.layerInDet = -1;

            // Helper: parse integer after tag, stopping at next '_'
            auto parseIntAfter = [&](const std::string& s,
                                     const std::string& tag) -> int {
              auto pos = s.find(tag);
              if (pos == std::string::npos) return -1;
              pos += tag.size();
              auto end = s.find('_', pos);
              return std::stoi(s.substr(pos, end == std::string::npos
                                              ? std::string::npos
                                              : end - pos));
            };

            if (detName == "SiTarget") {
              pi.detID      = 0;
              pi.station    = -1;
              pi.layerInDet = parseIntAfter(volName, "SiTarget_layer_");
              if (volName.find("slice_2") != std::string::npos)      pi.plane = 0;
              else if (volName.find("slice_4") != std::string::npos) pi.plane = 1;
            } else if (detName == "SiPad") {
              pi.detID      = 1;
              pi.station    = -1;
              pi.layerInDet = parseIntAfter(volName, "SiPad_layer_");
              pi.plane      = -1;
            } else if (detName == "MTC") {
              pi.detID   = 2;
              pi.plane   = (volName.find("_slice_2") != std::string::npos) ? 0 : 1;
              pi.station = volName.find("MTC40") != std::string::npos ? 0
                         : volName.find("MTC50") != std::string::npos ? 1 : 2;
              auto pos_l = volName.find("_layer_");
              auto pos_s = volName.find("_slice_", pos_l);
              if (pos_l != std::string::npos && pos_s != std::string::npos)
                pi.layerInDet = std::stoi(volName.substr(pos_l + 7,
                                                          pos_s - (pos_l + 7)));
            }

            planes.push_back(pi);
          }
        }
      }

      for (int i = 0; i < node->GetNdaughters(); ++i)
        walk(node->GetDaughter(i), curPath);
    };

    walk(mgr->GetTopNode(), "");

    std::sort(planes.begin(), planes.end(),
              [](const PlaneInfo& a, const PlaneInfo& b) {
                return a.z < b.z;
              });
    return planes;
  };

  // Detector names in Z order (SiTarget most upstream, MTC most downstream)
  const std::vector<std::string> detNames = {"SiTarget", "SiPad", "MTC"};

  // Collect ALL planes from ALL detectors into a single flat list.
  // This avoids CuboidVolumeBuilder multi-volume offset bugs.
  std::vector<PlaneInfo> allPlanes;
  double globalHalfX = 0.0;
  double globalHalfY = 0.0;

  for (const std::string& detName : detNames) {
    auto planes = extractPlanes(detName);

    if (planes.empty()) {
      error() << "[ACTSGeoSvc] No sensitive planes found for detector: "
              << detName << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "[ACTSGeoSvc] Detector " << detName
           << ": found " << planes.size() << " sensitive planes." << endmsg;

    for (const auto& pi : planes) {
      info() << "[ACTSGeoSvc]   plane z=" << pi.z
             << " halfX=" << pi.halfX << " halfY=" << pi.halfY
             << " thickness=" << pi.thickness
             << " detID=" << pi.detID
             << " station=" << pi.station
             << " layer=" << pi.layerInDet
             << " plane=" << pi.plane << endmsg;
      // Track largest detector half-sizes for volume bounds
      globalHalfX = std::max(globalHalfX, pi.halfX);
      globalHalfY = std::max(globalHalfY, pi.halfY);
      allPlanes.push_back(pi);
    }
  }

  // Sort all planes by X (should already be sorted since SiTarget < SiPad in X)
  std::sort(allPlanes.begin(), allPlanes.end(),
            [](const PlaneInfo& a, const PlaneInfo& b) {
              return a.z < b.z;
            });

  // =========================================================================
  // Build TrackingGeometry using PlaneLayer directly.
  // This avoids CuboidVolumeBuilder's position offset bug.
  // Each plane becomes a PlaneSurface + PlaneLayer with absolute position.
  // All layers are placed in a single CuboidVolumeBounds TrackingVolume.
  // =========================================================================

  // ---- Step 2: Create one PlaneLayer per surface --------------------------
  // Pass nullptr for surfaceArray: the PlaneLayer IS the sensitive surface.
  // The Navigator will resolve it directly via the layer surface representation.
  std::vector<Acts::LayerPtr> allLayers;
  allLayers.reserve(allPlanes.size());

  // Rotation 90° around Y: maps local Z (surface normal) to global X.
  // This avoids the theta=0 singularity in ACTS bound coordinates
  // when the beam (and track direction) is along global Z.
  // After this rotation, the surface normal points in X, and
  // seedDir=(1,0,0) gives theta=π/2 — far from the coordinate singularity.
  Acts::RotationMatrix3 rot90Y =
      Eigen::AngleAxisd(M_PI / 2.0,
                        Acts::Vector3::UnitY()).toRotationMatrix();

  for (std::size_t i = 0; i < allPlanes.size(); ++i) {
    const auto& pi = allPlanes[i];

    Acts::Transform3 transform = Acts::Transform3::Identity();
    // For SiTarget/SiPad: plain rot90Y maps surface normal to beam (ACTS X).
    // For MTC SciFi: R_X(∓α)·rot90Y tilts eBoundLoc0 to stereo direction
    // x∓y·tan(α), matching CartesianStripXStereo's strip_centre_at_y0.
    Acts::RotationMatrix3 surfaceRot = rot90Y;
    if (pi.detID == 2) {
      const double alpha_rad = m_mtcStereoAngle.value() * M_PI / 180.0;
      surfaceRot = Eigen::AngleAxisd(pi.plane == 0 ? +alpha_rad : -alpha_rad,
                       Acts::Vector3::UnitX()).toRotationMatrix() * rot90Y;
    }
    transform.rotate(surfaceRot);
    transform.pretranslate(Acts::Vector3(pi.z, 0.0, 0.0));

    // Bounds in the local XY plane (transverse to beam)
    auto bounds = std::make_shared<Acts::RectangleBounds>(
        pi.halfX, pi.halfY);

    const double layerThickness = 2.0 * pi.thickness + 0.01;

    auto layer = Acts::PlaneLayer::create(
        transform,
        bounds,
        std::unique_ptr<Acts::SurfaceArray>{},
        layerThickness,
        nullptr,
        Acts::active);

    allLayers.push_back(layer);
  }

  // ---- Step 3: Create navigation layers at boundaries --------------------
  // ACTS requires NavigationLayer objects at the volume boundaries for
  // correct extrapolation. Create one at xMin and one at xMax.

  const double zMin = allPlanes.front().z - allPlanes.front().thickness - 5.0;
  const double zMax = allPlanes.back().z  + allPlanes.back().thickness  + 5.0;

  // Navigation layer at zMin
  {
    Acts::Transform3 t = Acts::Transform3::Identity();
    t.rotate(rot90Y);
    t.pretranslate(Acts::Vector3(zMin, 0.0, 0.0));
    auto navBounds = std::make_shared<Acts::RectangleBounds>(
        globalHalfX + 5.0, globalHalfY + 5.0);
    auto navLayer = Acts::PlaneLayer::create(t, navBounds, nullptr, 0.0,
                                             nullptr,
                                             Acts::navigation);
    allLayers.insert(allLayers.begin(), navLayer);
  }

  // Navigation layer at zMax
  {
    Acts::Transform3 t = Acts::Transform3::Identity();
    t.rotate(rot90Y);
    t.pretranslate(Acts::Vector3(zMax, 0.0, 0.0));
    auto navBounds = std::make_shared<Acts::RectangleBounds>(
        globalHalfX + 5.0, globalHalfY + 5.0);
    auto navLayer = Acts::PlaneLayer::create(t, navBounds, nullptr, 0.0,
                                             nullptr,
                                             Acts::navigation);
    allLayers.push_back(navLayer);
  }

  // ---- Step 4: Build LayerArray -------------------------------------------
  // The LayerArrayCreator bins the layers along X (beam axis).

  Acts::LayerArrayCreator::Config lacCfg{};
  Acts::LayerArrayCreator lac(
      lacCfg,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::WARNING));

  auto layerArray = lac.layerArray(
      m_gctx,
      allLayers,
      zMin - 1.0,
      zMax + 1.0,
      Acts::arbitrary,          // use arbitrary binning (respects actual positions)
      Acts::AxisDirection::AxisX);

  // ---- Step 5: Build TrackingVolume with CuboidVolumeBounds ---------------
  // The volume must fully contain all layers.
  // halfX must be > abs(xMin) and > abs(xMax) since volume is centered at 0.
  // Beware of the rotation X<->Z: beam axis is now X, so we use zMin/zMax for volume halfX.

  const double volumeHalfZ = globalHalfX + 15.0;
  const double volumeHalfY = globalHalfY + 15.0;
  const double volumeHalfX = std::max(std::abs(zMin), std::abs(zMax)) + 2.0;

  auto volumeBounds = std::make_shared<Acts::CuboidVolumeBounds>(
      volumeHalfX, volumeHalfY, volumeHalfZ);

  // Volume is centered at the origin
  Acts::Transform3 volumeTransform = Acts::Transform3::Identity();

  auto trackingVolume = std::make_shared<Acts::TrackingVolume>(
    volumeTransform,
    volumeBounds,
    nullptr,                    // no volume material
    std::move(layerArray),      // staticLayerArray
    nullptr,                    // no contained volume array
    Acts::MutableTrackingVolumeVector{},  // no dense volumes
    "SNDVolume");

  // ---- Step 6: Build TrackingGeometry -------------------------------------
  Acts::GeometryIdentifierHook hook{};
  m_trackingGeometry = std::make_shared<Acts::TrackingGeometry>(
    trackingVolume,   // shared_ptr<TrackingVolume>
    nullptr,          // no material decorator
    hook,
    Acts::getDummyLogger());

  if (!m_trackingGeometry) {
    error() << "[ACTSGeoSvc] Failed to build ACTS TrackingGeometry."
            << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "[ACTSGeoSvc] TrackingGeometry built with PlaneLayer."
         << endmsg;

  // Collect sensitive surfaces from confined layers.
  // With nullptr SurfaceArray, the layer's surfaceRepresentation() IS
  // the sensitive surface that the Navigator will visit.
  const Acts::TrackingVolume* world =
      m_trackingGeometry->highestTrackingVolume();
  if (world && world->confinedLayers()) {
    for (const auto& layer : world->confinedLayers()->arrayObjects()) {
      if (layer && layer->layerType() == Acts::active) {
        m_allSurfaces.push_back(&layer->surfaceRepresentation());
      }
    }
  }
  std::sort(m_allSurfaces.begin(), m_allSurfaces.end(),
            [&](const Acts::Surface* a, const Acts::Surface* b) {
              return a->center(m_gctx).x() < b->center(m_gctx).x();
            });

  // allPlanes sorted by z matches m_allSurfaces sorted by ACTS-X — build map.
  if (allPlanes.size() == m_allSurfaces.size()) {
    for (std::size_t i = 0; i < allPlanes.size(); ++i) {
      const auto& pi = allPlanes[i];
      m_surfaceByAddressMap[{pi.detID, pi.station, pi.layerInDet, pi.plane}]
          = m_allSurfaces[i];
    }
  } else {
    warning() << "[ACTSGeoSvc] allPlanes.size()=" << allPlanes.size()
              << " != m_allSurfaces.size()=" << m_allSurfaces.size()
              << " — surfaceByAddress map not populated." << endmsg;
  }

  info() << "[ACTSGeoSvc] TrackingGeometry built. Total surfaces: "
         << m_allSurfaces.size() << endmsg;

  return StatusCode::SUCCESS;
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode ACTSGeoSvc::finalize() {
  return Service::finalize();
}

// ---------------------------------------------------------------------------
// Interface method implementations
// ---------------------------------------------------------------------------

const Acts::TrackingGeometry& ACTSGeoSvc::trackingGeometry() const {
  return *m_trackingGeometry;
}

const std::vector<const Acts::Surface*>& ACTSGeoSvc::allSurfaces() const {
  return m_allSurfaces;
}

const Acts::Surface* ACTSGeoSvc::surfaceByAddress(
    int detID, int station, int layer, int plane) const {
  auto it = m_surfaceByAddressMap.find({detID, station, layer, plane});
  return (it != m_surfaceByAddressMap.end()) ? it->second : nullptr;
}

const Acts::GeometryContext& ACTSGeoSvc::geometryContext() const {
  return m_gctx;
}

DECLARE_COMPONENT(ACTSGeoSvc)