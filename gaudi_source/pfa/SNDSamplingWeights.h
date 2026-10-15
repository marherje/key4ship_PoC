#pragma once

#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SNDSamplingWeights — how much energy a sampling plane did NOT see.
//
// SiTarget, SiPad and the MTC are sampling calorimeters: a hit records the
// energy left in a thin active layer, while most of the shower is absorbed in
// the tungsten or iron in front of it. Turning a hit into an energy therefore
// needs the ratio
//
//     w = (dE_passive + dE_active) / dE_active
//
// with dE the energy a minimum-ionising particle would lose crossing each, so
// that E_hit = E_deposited * w estimates what the shower actually carried.
//
// The weight is computed PER LAYER, not per sensitive plane: dE_passive is
// all the absorber of the layer and dE_active is the sum over ITS planes, so
// every plane of a layer gets the same w and the layer's energy comes out as
// w * sum(deposits). Attributing the absorber to the first plane behind it
// instead -- the obvious reading, and what this file did until 2026-09-23 --
// leaves the second view weighing ~1 and makes it worthless for energy: on
// 50 GENIE events the MTC V plane held MORE hits than U (2307 vs 2251) and
// carried 0.45 % of the energy against U's 55 %. Two things go wrong with
// that. The orthogonal views sample the same particle crossing, so averaging
// them is a sqrt(2) on the sampling fluctuation, thrown away; and if the
// heavy plane's hit falls below threshold while its partner's does not, the
// layer loses essentially all of its energy. Sharing the absorber across the
// layer keeps the scale (SiTarget: 33.9 * (E_X + E_Y) = 67.8 * E_X when both
// planes see the same) and fixes both.
//
// The weights are COMPUTED FROM THE GEOMETRY, never hard-coded: the MOBO
// search varies the tungsten thickness of both silicon detectors, so a literal
// would silently decalibrate every variant. dd4hep::rec::MaterialManager walks
// the real material between consecutive sensitive planes along the beam and
// this header multiplies each slab by the minimum-ionising dE/dx of its
// material.
//
// Limits worth knowing:
//   * the walk is along (x, y) = (0, 0), the beam axis, so it assumes the
//     layers are uniform in the transverse plane -- true for this detector,
//     where the only transverse structure is the MTC's magnetised iron core
//     and its steel filler, which the axis crosses at its centre;
//   * dE/dx_min is a flat value per material, i.e. the weight is a MIP-scale
//     ratio, not a shower-dependent one. That is the usual sampling-fraction
//     approximation and the reason phase 1.4 still fits one global scale per
//     detector on top of it (e/h and the EM scale).
//
// HOW THE SCAN IS DONE. The obvious tool, dd4hep::rec::MaterialManager's
// materialsBetween(), was tried first and dropped: on this geometry it
// returned only the last 3 mm of the 53 mm of iron between two MTC layers
// (measured 2026-09-23, scan of z = 81..134 mm), while covering other spans
// correctly. Rather than work around a tool whose answer cannot be trusted
// slab by slab, the walk here is TGeo's own -- InitTrack, FindNextBoundary,
// Step -- which is what the material managers use underneath and which
// reproduces the compact exactly: 53.000 mm of Iron for that same span,
// 3.5 mm W + 0.6 mm Si for a SiTarget layer, 10 mm W + 0.65 mm Si for a
// SiPad one.
//
// UNITS. TGeo works in centimetres and this header speaks millimetres, so
// positions are divided by 10 going in and step lengths multiplied by 10
// coming out. TGeoMaterial::GetDensity() is g/cm^3 already.
// ---------------------------------------------------------------------------

namespace snd {

// Minimum ionising dE/dx [MeV cm^2 / g] per material name, as DD4hep spells
// it in SND_compact.xml / materials.xml. PDG values for the elements; the
// compounds are close enough to 1.7-1.9 that the exact choice moves a weight
// by well under a percent, since they are all thin.
inline const std::map<std::string, double>& dEdxMinTable() {
  static const std::map<std::string, double> kTable = {
      {"TungstenDens1910", 1.145}, {"Tungsten", 1.145},
      {"Iron", 1.451},             {"Steel", 1.451},
      {"StainlessSteel", 1.451},   {"Silicon", 1.664},
      {"Copper", 1.403},           {"Air", 1.815},
      {"Scintillator", 1.956},     {"Polystyrene", 1.936},
      {"PVT", 1.956},              {"Carbon", 1.745},
      {"CarbonFiber", 1.745},      {"CarbonFiber_25percent", 1.745},
      {"Kapton", 1.700},           {"PCB", 1.700},
      {"G10", 1.700},              {"Epoxy", 1.900},
      {"Glue", 1.900},             {"Vacuum", 0.0},
  };
  return kTable;
}

inline constexpr double kDefaultDEdxMin = 1.70;  // MeV cm^2/g

// A sensitive plane, as the geometry describes it.
struct SensitivePlane {
  int    det = -1;           // 0=SiTarget, 1=SiPad, 2=MTC
  int    view = -1;          // 0=X/U, 1=Y/V, 2=2D
  int    layerOrdinal = -1;  // contiguous along the beam within the detector
  double z = 0.0;            // centre along the beam [mm]
  double thickness = 0.0;    // full thickness [mm]
  std::string material;      // active material name
};

// What the material scan found for one plane.
struct PlaneWeight {
  double passiveMeV = 0.0;  // MIP energy loss upstream, since the previous plane
  double activeMeV = 0.0;   // MIP energy loss inside the plane itself
  // The same two sums over every plane of this plane's LAYER; the weight is
  // built from these, so all planes of a layer share it.
  double layerPassiveMeV = 0.0;
  double layerActiveMeV = 0.0;
  double weight = 1.0;      // (layerPassive + layerActive) / layerActive
  bool   valid = false;
};

inline double dEdxMin(const std::string& material,
                      std::set<std::string>* unknown = nullptr) {
  const auto& table = dEdxMinTable();
  auto it = table.find(material);
  if (it != table.end()) return it->second;
  if (unknown) unknown->insert(material);
  return kDefaultDEdxMin;
}

// MIP energy loss of one slab: dE/dx [MeV cm^2/g] * rho [g/cm^3] * L [mm]/10.
inline double slabMeV(const std::string& material, double density,
                      double thicknessMm,
                      std::set<std::string>* unknown = nullptr) {
  return dEdxMin(material, unknown) * density * (thicknessMm / 10.0);
}

// The sensitive volumes, identified by DD4hep's own sensitivity flag rather
// than by name, with det/view read off the slice they belong to.
inline std::vector<SensitivePlane> collectSensitivePlanes() {
  std::vector<SensitivePlane> planes;
  std::set<std::tuple<int, int, long>> seenPlane;
  TGeoManager* mgr = gGeoManager;
  if (!mgr) return planes;

  std::function<void(TGeoNode*, const std::string&)> walk =
      [&](TGeoNode* node, const std::string& path) {
    TGeoVolume* vol = node->GetVolume();
    if (!vol) return;
    const std::string volName = vol->GetName();
    const std::string curPath = path + "/" + node->GetName();

    // Where the sensitivity actually sits differs per detector: the MTC
    // flags the slice itself, while the SiTarget flags each
    // <slice>_colN_rowM tile and the SiPad each <slice>_wafer_pads. Taking
    // only the flagged volumes and de-duplicating them by (det, view, z)
    // afterwards handles all three without a per-detector rule -- 1920
    // SiTarget volumes collapse to 240 planes, 792 SiPad ones to 22.
    if (dd4hep::Volume(vol).isSensitive()) {
      int det = -1, view = -1;
      if (volName.find("SiTarget_layer_") != std::string::npos) {
        det = 0;
        if (volName.find("_slice_2") != std::string::npos) view = 0;
        else if (volName.find("_slice_4") != std::string::npos) view = 1;
      } else if (volName.find("SiPad_layer_") != std::string::npos) {
        det = 1;
        view = 2;
      } else if (curPath.find("MTC") != std::string::npos) {
        det = 2;
        if (volName.find("_slice_2") != std::string::npos) view = 0;
        else if (volName.find("_slice_4") != std::string::npos) view = 1;
        else if (volName.find("_slice_6") != std::string::npos) view = 2;
      }

      if (det >= 0 && view >= 0 && mgr->cd(curPath.c_str())) {
        const double* tr = mgr->GetCurrentMatrix()->GetTranslation();
        if (auto* box = dynamic_cast<TGeoBBox*>(vol->GetShape())) {
          SensitivePlane p;
          p.det = det;
          p.view = view;
          p.z = tr[2] * 10.0;                 // TGeo cm -> mm
          p.thickness = box->GetDZ() * 20.0;  // half-size cm -> full mm
          if (vol->GetMaterial()) p.material = vol->GetMaterial()->GetName();
          // One entry per physical plane, not per tile.
          const long zKey = std::lround(p.z * 1000.0);
          if (seenPlane.emplace(det, view, zKey).second)
            planes.push_back(p);
        }
      }
    }
    for (int i = 0; i < node->GetNdaughters(); ++i)
      walk(node->GetDaughter(i), curPath);
  };

  if (mgr->GetTopNode()) walk(mgr->GetTopNode(), "");

  // Sort along the beam and number the layers per detector. The ordinal has
  // to match SNDCellDecoder's, which counts physical layers, so planes of
  // the same layer (SiTarget X/Y, MTC U/V/scint) share one.
  std::sort(planes.begin(), planes.end(),
            [](const SensitivePlane& a, const SensitivePlane& b) {
              return (a.det != b.det) ? (a.det < b.det) : (a.z < b.z);
            });
  const std::map<int, int> planesPerLayer{{0, 2}, {1, 1}, {2, 3}};
  std::map<int, int> seen;
  for (auto& p : planes) {
    const int n = planesPerLayer.at(p.det);
    p.layerOrdinal = seen[p.det] / n;
    ++seen[p.det];
  }
  return planes;
}

// Millimetres of one material along the beam axis between two z positions,
// which is what a muon's range is measured in: the MTC's iron.
inline double scanMaterialMm(TGeoManager* mgr, double z0mm, double z1mm,
                             const std::string& materialName) {
  if (!mgr || z1mm <= z0mm) return 0.0;
  double total = 0.0;
  constexpr double kNudge = 1e-7;
  double z = z0mm / 10.0 + kNudge;
  const double zEnd = z1mm / 10.0;
  mgr->InitTrack(0.0, 0.0, z, 0.0, 0.0, 1.0);
  int guard = 0;
  while (z < zEnd - 1e-9 && ++guard < 100000) {
    TGeoNode* node = mgr->GetCurrentNode();
    const TGeoMaterial* mat =
        (node && node->GetVolume()) ? node->GetVolume()->GetMaterial() : nullptr;
    mgr->FindNextBoundary();
    double step = mgr->GetStep();
    if (step > zEnd - z) step = zEnd - z;
    if (step <= kNudge) {
      z += kNudge;
      mgr->SetCurrentPoint(0.0, 0.0, z);
      mgr->FindNode();
      continue;
    }
    if (mat && materialName == mat->GetName()) total += step * 10.0;
    mgr->Step();
    z += step;
    if (mgr->IsOutside()) break;
  }
  return total;
}

// MIP energy loss along the beam axis between two z positions [mm], walking
// the real geometry boundary by boundary.
inline double scanMeV(TGeoManager* mgr, double z0mm, double z1mm,
                      std::set<std::string>* unknown = nullptr) {
  if (!mgr || z1mm <= z0mm) return 0.0;
  double total = 0.0;
  const double zEnd = z1mm / 10.0;
  // Nudge [cm]: both ends of the interval are volume boundaries by
  // construction (they are the faces of a sensitive plane), and TGeo answers
  // a zero step when a track starts exactly on one. Treating that as the end
  // of the scan -- which this did until 2026-09-23 -- silently truncated some
  // layers: the SiTarget came out with two different weights, 33.91 on 63
  // layers and 47.06 on 54, where the geometry is periodic and they must all
  // be equal. Starting just inside, and stepping over a null boundary instead
  // of stopping at it, makes every layer of a detector agree.
  constexpr double kNudge = 1e-7;   // 1 nm
  double z = z0mm / 10.0 + kNudge;
  mgr->InitTrack(0.0, 0.0, z, 0.0, 0.0, 1.0);
  int guard = 0;
  while (z < zEnd - 1e-9 && ++guard < 100000) {
    TGeoNode* node = mgr->GetCurrentNode();
    const TGeoMaterial* mat =
        (node && node->GetVolume()) ? node->GetVolume()->GetMaterial() : nullptr;
    mgr->FindNextBoundary();
    double step = mgr->GetStep();
    if (step > zEnd - z) step = zEnd - z;
    if (step <= kNudge) {
      // Sitting on a boundary: cross it and carry on.
      z += kNudge;
      mgr->SetCurrentPoint(0.0, 0.0, z);
      mgr->FindNode();
      continue;
    }
    if (mat)
      total += slabMeV(mat->GetName(), mat->GetDensity(), step * 10.0, unknown);
    mgr->Step();
    z += step;
    if (mgr->IsOutside()) break;
  }
  return total;
}

// Computes the weight of every plane in `planes` (which must be sorted along
// the beam within each detector). The passive budget of a plane is everything
// between the downstream face of the previous sensitive plane and its own
// upstream face; for the first plane of a detector it is measured from
// `upstreamStartZ` for that detector, so the absorber in front of layer 0
// counts too.
inline std::map<std::tuple<int, int, int>, PlaneWeight> computePlaneWeights(
    dd4hep::Detector& detector,
    const std::vector<SensitivePlane>& planes,
    const std::map<int, double>& upstreamStartZ,
    std::set<std::string>* unknownMaterials = nullptr) {
  std::map<std::tuple<int, int, int>, PlaneWeight> out;
  if (planes.empty()) return out;

  TGeoManager* mgr = detector.world().volume().ptr()
      ? gGeoManager : gGeoManager;
  if (!mgr) return out;

  std::map<int, double> prevEnd;  // per detector: end of the previous plane
  for (const auto& p : planes) {
    PlaneWeight w;
    const double halfT = 0.5 * p.thickness;
    const double zFront = p.z - halfT;
    const double zBack  = p.z + halfT;

    double zStart = zFront;
    auto itPrev = prevEnd.find(p.det);
    if (itPrev != prevEnd.end()) {
      zStart = itPrev->second;
    } else {
      auto itStart = upstreamStartZ.find(p.det);
      if (itStart != upstreamStartZ.end()) zStart = itStart->second;
    }

    if (zFront - zStart > 1e-3)
      w.passiveMeV = scanMeV(mgr, zStart, zFront, unknownMaterials);

    // The active slab is taken from the geometry too rather than from the
    // plane's nominal material alone, so a plane built of several sub-slabs
    // is still accounted correctly.
    w.activeMeV = scanMeV(mgr, zFront, zBack, unknownMaterials);
    if (w.activeMeV <= 0.0) {
      // Fall back on the nominal material so a plane never gets weight 0.
      w.activeMeV = slabMeV(p.material, 1.0, p.thickness, unknownMaterials);
    }

    // Per-plane numbers for now; the per-layer sharing follows below.
    w.weight = (w.activeMeV > 0.0)
        ? (w.passiveMeV + w.activeMeV) / w.activeMeV
        : 1.0;
    w.valid = (w.activeMeV > 0.0);
    out[{p.det, p.view, p.layerOrdinal}] = w;
    prevEnd[p.det] = zBack;
  }

  // Share each layer's absorber over its own sensitive planes: every plane of
  // a layer ends up with the same weight, so summing w * deposit over the
  // planes reconstructs the layer's energy once, using all of its views.
  std::map<std::pair<int, int>, std::pair<double, double>> perLayer;  // passive, active
  for (const auto& [key, w] : out) {
    const auto& [det, view, ordinal] = key;
    (void)view;
    auto& acc = perLayer[{det, ordinal}];
    acc.first += w.passiveMeV;
    acc.second += w.activeMeV;
  }
  for (auto& [key, w] : out) {
    const auto& [det, view, ordinal] = key;
    (void)view;
    const auto& acc = perLayer.at({det, ordinal});
    w.layerPassiveMeV = acc.first;
    w.layerActiveMeV = acc.second;
    w.weight = (acc.second > 0.0) ? (acc.first + acc.second) / acc.second : 1.0;
    w.valid = (acc.second > 0.0);
  }
  return out;
}

}  // namespace snd
