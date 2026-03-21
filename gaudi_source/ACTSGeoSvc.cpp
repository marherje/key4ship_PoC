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
#include "Acts/Surfaces/PlaneSurface.hpp"
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
    double x;          // global X center [mm]
    double halfY;      // half-size in Y [mm]
    double halfZ;      // half-size in Z [mm]
    double thickness;  // half-thickness in X [mm]
    int    plane;      // 0=StripY, 1=StripZ, -1=pixel (SiPad)
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
      }

      if (inDetector && isSensitive) {
        if (mgr->cd(curPath.c_str())) {
          TGeoMatrix* mat  = mgr->GetCurrentMatrix();
          const double* tr = mat->GetTranslation();
          TGeoBBox* box    = dynamic_cast<TGeoBBox*>(vol->GetShape());
          if (box) {
            PlaneInfo pi;
            // TGeo uses cm; convert to mm (* 10)
            pi.x         = tr[0]        * 10.0;
            pi.halfY     = box->GetDY() * 10.0;
            pi.halfZ     = box->GetDZ() * 10.0;
            pi.thickness = box->GetDX() * 10.0;

            if (detName == "SiTarget") {
              if (volName.find("slice_2") != std::string::npos)
                pi.plane = 0;
              else if (volName.find("slice_4") != std::string::npos)
                pi.plane = 1;
            } else if (detName == "SiPad") {
              if (volName.find("slice_4") != std::string::npos)
                pi.plane = -1;
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
                return a.x < b.x;
              });
    return planes;
  };

  // Detector names in X order (SiTarget more negative, SiPad less negative)
  const std::vector<std::string> detNames = {"SiTarget", "SiPad"};

  Acts::RotationMatrix3 rotIdentity = Acts::RotationMatrix3::Identity();
  Acts::RotationMatrix3 rotStripZ;
  rotStripZ = Eigen::AngleAxisd(M_PI / 2.0, Acts::Vector3::UnitX());

  // Collect ALL planes from ALL detectors into a single flat list.
  // This avoids CuboidVolumeBuilder multi-volume offset bugs.
  std::vector<PlaneInfo> allPlanes;
  double globalHalfY = 0.0;
  double globalHalfZ = 0.0;

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
      info() << "[ACTSGeoSvc]   plane x=" << pi.x
             << " halfY=" << pi.halfY << " halfZ=" << pi.halfZ
             << " thickness=" << pi.thickness
             << " plane=" << pi.plane << endmsg;
      // Track largest detector half-sizes for volume bounds
      globalHalfY = std::max(globalHalfY, pi.halfY);
      globalHalfZ = std::max(globalHalfZ, pi.halfZ);
      allPlanes.push_back(pi);
    }
  }

  // Sort all planes by X (should already be sorted since SiTarget < SiPad in X)
  std::sort(allPlanes.begin(), allPlanes.end(),
            [](const PlaneInfo& a, const PlaneInfo& b) {
              return a.x < b.x;
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

  for (std::size_t i = 0; i < allPlanes.size(); ++i) {
    const auto& pi = allPlanes[i];
    Acts::RotationMatrix3 rot = (pi.plane == 1) ? rotStripZ : rotIdentity;

    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.rotate(rot);
    transform.pretranslate(Acts::Vector3(pi.x, 0.0, 0.0));

    auto bounds = std::make_shared<Acts::RectangleBounds>(
        pi.halfY, pi.halfZ);

    const double layerThickness = 2.0 * pi.thickness + 0.01;

    // nullptr SurfaceArray: the layer surface representation itself
    // is the sensitive surface that the CKF will associate hits with.
    auto layer = Acts::PlaneLayer::create(
        transform,
        bounds,
        std::unique_ptr<Acts::SurfaceArray>{}, 
        layerThickness,
        nullptr,           // no approach descriptor
        Acts::active);     // mark as active (sensitive)

    allLayers.push_back(layer);
  }

  // ---- Step 3: Create navigation layers at boundaries --------------------
  // ACTS requires NavigationLayer objects at the volume boundaries for
  // correct extrapolation. Create one at xMin and one at xMax.

  const double xMin = allPlanes.front().x - allPlanes.front().thickness - 5.0;
  const double xMax = allPlanes.back().x  + allPlanes.back().thickness  + 5.0;

  // Navigation layer at xMin
  {
    Acts::Transform3 t = Acts::Transform3::Identity();
    t.pretranslate(Acts::Vector3(xMin, 0.0, 0.0));
    auto navBounds = std::make_shared<Acts::RectangleBounds>(
        globalHalfY + 5.0, globalHalfZ + 5.0);
    auto navLayer = Acts::PlaneLayer::create(t, navBounds, nullptr, 0.0,
                                             nullptr,
                                             Acts::navigation);
    allLayers.insert(allLayers.begin(), navLayer);
  }

  // Navigation layer at xMax
  {
    Acts::Transform3 t = Acts::Transform3::Identity();
    t.pretranslate(Acts::Vector3(xMax, 0.0, 0.0));
    auto navBounds = std::make_shared<Acts::RectangleBounds>(
        globalHalfY + 5.0, globalHalfZ + 5.0);
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
      xMin - 1.0,
      xMax + 1.0,
      Acts::arbitrary,          // use arbitrary binning (respects actual positions)
      Acts::AxisDirection::AxisX);

  // ---- Step 5: Build TrackingVolume with CuboidVolumeBounds ---------------
  // The volume must fully contain all layers.
  // halfX must be > abs(xMin) and > abs(xMax) since volume is centered at 0.

  const double volumeHalfX = std::max(std::abs(xMin), std::abs(xMax)) + 2.0;
  const double volumeHalfY = globalHalfY + 15.0;
  const double volumeHalfZ = globalHalfZ + 15.0;

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

  // Populate m_allSurfaces: all rectangle surfaces sorted by X
  const Acts::TrackingVolume* world = m_trackingGeometry->highestTrackingVolume();
  if (world && world->confinedLayers()) {
    for (const auto& layer : world->confinedLayers()->arrayObjects()) {
      if (layer && layer->layerType() == Acts::active) {
        m_allSurfaces.push_back(&layer->surfaceRepresentation());
      }
    }
  }
  // Sort by X
  std::sort(m_allSurfaces.begin(), m_allSurfaces.end(),
          [&](const Acts::Surface* a, const Acts::Surface* b) {
            return a->center(m_gctx).x() < b->center(m_gctx).x();
          });
  std::sort(m_allSurfaces.begin(), m_allSurfaces.end(),
            [&](const Acts::Surface* a, const Acts::Surface* b) {
              return a->center(m_gctx).x() < b->center(m_gctx).x();
            });

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

const Acts::GeometryContext& ACTSGeoSvc::geometryContext() const {
  return m_gctx;
}

DECLARE_COMPONENT(ACTSGeoSvc)