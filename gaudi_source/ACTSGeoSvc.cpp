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
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <string>

// ---------------------------------------------------------------------------
// SNDDetectorElement
// A minimal DetectorElementBase so the CKF actor classifies each plane as
// sensitive (CKF checks associatedDetectorElement() != nullptr at line 442
// of CombinatorialKalmanFilter.hpp; without this, every surface is "passive").
// ---------------------------------------------------------------------------

class SNDDetectorElement : public Acts::DetectorElementBase {
public:
  SNDDetectorElement(std::shared_ptr<const Acts::PlanarBounds> bounds,
                     Acts::Transform3 transform, double thickness)
      : m_transform(std::move(transform)), m_thickness(thickness) {
    // PlaneSurface(bounds, detElement) sets m_associatedDetectorElement = this
    m_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(bounds, *this);
  }

  const Acts::Transform3& transform(const Acts::GeometryContext&) const override {
    return m_transform;
  }
  const Acts::Surface& surface() const override { return *m_surface; }
  Acts::Surface&       surface() override       { return *m_surface; }
  double thickness() const override { return m_thickness; }

private:
  Acts::Transform3 m_transform;
  double m_thickness;
  std::shared_ptr<Acts::PlaneSurface> m_surface;
};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

ACTSGeoSvc::ACTSGeoSvc(const std::string& name, ISvcLocator* svcLoc)
    : extends(name, svcLoc) {}

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

namespace {
// Query DD4hep for radiation length and nuclear interaction length of a
// named material, then pack into an Acts::MaterialSlab of the given thickness.
// DD4hep radLength()/intLength() return cm; ACTS wants mm.
Acts::MaterialSlab makeSlab(const std::string& matName, double thicknessMm) {
    auto& desc = dd4hep::Detector::getInstance();
    dd4hep::Material mat = desc.material(matName);
    float X0  = static_cast<float>(mat.radLength() * 10.0);   // cm → mm
    float L0  = static_cast<float>(mat.intLength() * 10.0);   // cm → mm
    float rho = static_cast<float>(mat.density());             // g/cm³
    float A   = static_cast<float>(mat.A());
    float Z   = static_cast<float>(mat.Z());
    return Acts::MaterialSlab(
        Acts::Material::fromMassDensity(X0, L0, A, Z, rho),
        static_cast<float>(thicknessMm));
}
} // namespace

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

  // ---- Insert combined MTC surfaces (plane=3) at z_mid of each (station, layer)
  // These are flat (non-stereo) planes used for 2D paired-strip measurements.
  // Without them, paired U+V hits would have to land on a stereo surface where
  // the covariance is non-diagonal. The 1D unpaired hits still target the U/V
  // surfaces.
  {
    std::map<std::pair<int,int>, std::pair<const PlaneInfo*, const PlaneInfo*>> mtcUV;
    for (const auto& pi : allPlanes) {
      if (pi.detID != 2) continue;
      if (pi.plane != 0 && pi.plane != 1) continue;
      auto& slot = mtcUV[{pi.station, pi.layerInDet}];
      if (pi.plane == 0) slot.first  = &pi;
      else               slot.second = &pi;
    }
    std::vector<PlaneInfo> combined;
    combined.reserve(mtcUV.size());
    for (const auto& [key, uv] : mtcUV) {
      if (!uv.first || !uv.second) continue;
      PlaneInfo pc = *uv.first;
      pc.z         = 0.5 * (uv.first->z + uv.second->z);
      pc.plane     = 3;                          // combined / paired
      // Virtual measurement surface: zero physical extent so the layer array
      // does not flag a bounds overlap with the closely-spaced U/V layers.
      // Material lives on U (iron+scint) and V (scint); the combined surface
      // contributes nothing to MS/energyLoss.
      pc.thickness = 0.0;
      combined.push_back(pc);
    }
    for (const auto& pc : combined) allPlanes.push_back(pc);
    info() << "[ACTSGeoSvc] Added " << combined.size()
           << " combined MTC surfaces (plane=3) at U/V midpoints." << endmsg;
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
  // Each active layer gets a PlaneSurface module + SurfaceArray so the ACTS
  // Navigator can find it via resolveSensitive=true. Without a non-null
  // SurfaceArray, Layer::resolve() always returns false and the CKF visits
  // zero surfaces.
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
    if (pi.detID == 2 && (pi.plane == 0 || pi.plane == 1)) {
      // MTC SciFi U (plane=0) / V (plane=1): apply ±α stereo tilt.
      // Combined MTC surfaces (plane=3) keep rot90Y only (flat).
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

    // Detector element: sets associatedDetectorElement() on the surface so the
    // CKF actor treats it as sensitive (not passive). Must outlive the surface.
    auto detElem = std::make_shared<SNDDetectorElement>(
        bounds, transform, pi.thickness);
    m_detectorElements.push_back(detElem);

    auto surfaceArray = std::make_unique<Acts::SurfaceArray>(
        detElem->surface().getSharedPtr());

    auto layer = Acts::PlaneLayer::create(
        transform,
        bounds,
        std::move(surfaceArray),
        layerThickness,
        nullptr,
        Acts::active);

    // Attach surface material from DD4hep — never hard-code X0/rho.
    // Layout confirmed from SND_compact.xml layer slice order.
    {
        std::shared_ptr<const Acts::ISurfaceMaterial> surfMat;
        if (pi.detID == 0) {
            // SiTarget: TungstenDens1910(3.5mm) upstream of plane=0 only
            if (pi.plane == 0) {
                surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                    Acts::MaterialSlab::combineLayers(
                        makeSlab("TungstenDens1910", 3.5),
                        makeSlab("Silicon", 0.3)));
            } else {
                surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                    makeSlab("Silicon", 0.3));
            }
        } else if (pi.detID == 1) {
            // SiPad: TungstenDens1910(3.5mm) + Silicon(0.65mm); skip thin Cu/CF
            surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                Acts::MaterialSlab::combineLayers(
                    makeSlab("TungstenDens1910", 3.5),
                    makeSlab("Silicon", 0.65)));
        } else if (pi.detID == 2) {
            if (pi.plane == 0) {
                // MTC U: 50mm outer + 3mm inner iron upstream, then scifi
                surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                    Acts::MaterialSlab::combineLayers(
                        makeSlab("Iron", 53.0),
                        makeSlab("Scintillator", 1.35)));
            } else if (pi.plane == 1) {
                // MTC V: only 1mm air between U and V — just attach scifi
                surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                    makeSlab("Scintillator", 1.35));
            } else if (pi.plane == 2) {
                // MTC Scint: 3mm inner iron upstream, then scintillator
                surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
                    Acts::MaterialSlab::combineLayers(
                        makeSlab("Iron", 3.0),
                        makeSlab("Scintillator", 15.0)));
            }
        }
        if (surfMat) {
            detElem->surface().assignSurfaceMaterial(surfMat);
        }
    }

    allLayers.push_back(layer);
  }

  // Material sanity dump — eyeball one surface: Si X0~93.7mm, Fe X0~17.6mm, W X0~3.5mm
  for (const auto& de : m_detectorElements) {
    const auto* mat = de->surface().surfaceMaterial();
    if (!mat) continue;
    auto slab = mat->materialSlab(Acts::Vector2{0, 0});
    info() << "[ACTSGeoSvc] MatDump: X0=" << slab.material().X0()
           << " mm  L0=" << slab.material().L0()
           << " mm  t=" << slab.thickness() << " mm" << endmsg;
    break;
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
  // We use the module surfaces stored in each layer's SurfaceArray — these
  // are the same surfaces the Navigator visits via resolveSensitive=true.
  // Using surfaceRepresentation() instead would give the wrong geometry IDs
  // (the layer envelope, not the sensitive module surface).
  const Acts::TrackingVolume* world =
      m_trackingGeometry->highestTrackingVolume();
  if (world && world->confinedLayers()) {
    for (const auto& layer : world->confinedLayers()->arrayObjects()) {
      if (layer && layer->layerType() == Acts::active) {
        const Acts::SurfaceArray* sa = layer->surfaceArray();
        if (sa) {
          for (const Acts::Surface* sf : sa->surfaces()) {
            if (sf) m_allSurfaces.push_back(sf);
          }
        } else {
          m_allSurfaces.push_back(&layer->surfaceRepresentation());
        }
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