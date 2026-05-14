//====================================================================
//  DD4hep detector plugin for MTCDetector
//  Three stations (MTC40/50/60) built from a single shared layer
//  template. A parent Assembly anchors the "MTC" DetElement in the
//  world; each station envelope lives inside it at its absolute Z.
//  Layer structure per station:
//    outer iron | inner iron | SciFi U (+alpha) | gap | SciFi V (-alpha)
//    | inner iron | scintillator
//
//  If the detector XML has a field_y attribute (e.g., 1.0*tesla),
//  a uniform By field is added to the global DD4hep field for each
//  outer iron absorber slab (slice 0, 50 mm thick per layer).
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/FieldTypes.h"
#include "DD4hep/Printout.h"
#include "XML/Utilities.h"
#include <cmath>
#include <string>
#include <vector>

using namespace dd4hep;

// ---------------------------------------------------------------------------
// Field object: returns By only inside the registered outer iron slabs
// ---------------------------------------------------------------------------
namespace {
  struct MTCIronField : public CartesianField::Object {
    struct Slab { double zlo, zhi, xhalf, yhalf; };
    double            m_by = 0.0;
    std::vector<Slab> m_slabs;

    void fieldComponents(const double* pos, double* field) override {
      double z  = pos[2];
      double ax = std::abs(pos[0]);
      double ay = std::abs(pos[1]);
      for (const auto& s : m_slabs) {
        if (z >= s.zlo && z <= s.zhi && ax <= s.xhalf && ay <= s.yhalf) {
          field[1] += m_by;
          return;
        }
      }
    }
  };
}

// ---------------------------------------------------------------------------
// Detector factory
// ---------------------------------------------------------------------------
dd4hep::Ref_t create_MTCDetector(dd4hep::Detector& description,
                                  xml_h element,
                                  dd4hep::SensitiveDetector sens) {

  xml_det_t   x_det  = element;
  std::string name   = x_det.nameStr();
  int         det_id = x_det.id();
  DetElement  sdet(name, det_id);

  sens.setType("calorimeter");
  xml::setDetectorTypeFlag(element, sdet);

  Material mat_air = description.air();

  // -- Read optional magnetic field strength (0 = no field) -----------------
  bool   hasField = x_det.hasAttr(_Unicode(field_y));
  double bFieldY  = hasField ? x_det.attr<double>(_Unicode(field_y)) : 0.0;

  // Outer iron slabs for the field map (populated during build)
  std::vector<MTCIronField::Slab> ironSlabs;

  // -- Parent assembly gives sdet a valid placement -------------------------
  // Assembly auto-sizes to contain all station envelopes placed inside it.
  Assembly     parent_asm(name);
  Volume       motherVol = description.pickMotherVolume(sdet);
  PlacedVolume parent_pv = motherVol.placeVolume(parent_asm, Position(0.0, 0.0, 0.0));
  parent_pv.addPhysVolID("system", det_id);
  sdet.setPlacement(parent_pv);

  // -- Pre-compute total Z extent from shared layer template ----------------
  double total_z = 0.0;
  for (xml_coll_t lColl(x_det, _Unicode(layer)); lColl; ++lColl) {
    xml_comp_t x_layer = lColl;
    int repeat = x_layer.hasAttr(_Unicode(repeat))
                   ? x_layer.attr<int>(_Unicode(repeat)) : 1;
    double layer_thick = 0.0;
    for (xml_coll_t sColl(x_layer, _Unicode(slice)); sColl; ++sColl)
      layer_thick += xml_comp_t(sColl).attr<double>(_Unicode(thickness));
    total_z += repeat * layer_thick;
  }

  // -- Loop over <station> children -----------------------------------------
  for (xml_coll_t stColl(x_det, _Unicode(station)); stColl; ++stColl) {
    xml_comp_t  x_st    = stColl;
    std::string st_name = x_st.attr<std::string>(_Unicode(name));
    int         st_id   = x_st.attr<int>(_Unicode(station_id));
    double      env_w   = x_st.attr<double>(_Unicode(env_width));
    double      env_h   = x_st.attr<double>(_Unicode(env_height));
    double      z0      = x_st.attr<double>(_Unicode(z_position));
    double      ang_deg = x_st.attr<double>(_Unicode(fiber_angle));

    // Para (parallelepiped) parameters for SciFi U/V slices.
    // alpha_para is the TGeoPara shear angle (degrees): encodes the stereo tilt
    // without rotating the volume, so the slice stays within the station envelope.
    double cos_a      = std::cos(ang_deg * M_PI / 180.0);
    double sin_a      = std::sin(ang_deg * M_PI / 180.0);
    double alpha_para = -std::atan(sin_a) * 180.0 / M_PI;  // degrees, U plane
    double scifi_dx   = env_w / 2.0 - (env_h / 2.0) * sin_a;             // half-width of SciFi slice

    dd4hep::printout(dd4hep::INFO, "MTCDetector",
        "%s: env=(%.0f x %.0f mm) z0=%.1f mm fiber_angle=%.1f deg",
        st_name.c_str(), env_w, env_h, z0, ang_deg);

    // -- Build station envelope and place inside parent assembly -------------
    const double safety = 0.1;  // mm
    Box    env_shape(env_w / 2.0 + safety, env_h / 2.0 + safety, total_z / 2.0 + safety);
    Volume env_vol(st_name + "_envelope", env_shape, mat_air);
    env_vol.setVisAttributes(description, "InvisibleNoDaughters");

    // Station envelope placed at its absolute Z inside the assembly (assembly is at origin).
    PlacedVolume env_pv = parent_asm.placeVolume(env_vol, Position(0.0, 0.0, z0));
    env_pv.addPhysVolID("station", st_id);

    DetElement st_de(sdet, st_name, st_id);
    st_de.setPlacement(env_pv);

    // -- Layer and slice loop ------------------------------------------------
    int    layer_idx = 0;
    double cur_z     = -total_z / 2.0;

    for (xml_coll_t lColl(x_det, _Unicode(layer)); lColl; ++lColl) {
      xml_comp_t x_layer = lColl;
      int repeat = x_layer.hasAttr(_Unicode(repeat))
                     ? x_layer.attr<int>(_Unicode(repeat)) : 1;

      double layer_thick = 0.0;
      for (xml_coll_t sColl(x_layer, _Unicode(slice)); sColl; ++sColl)
        layer_thick += xml_comp_t(sColl).attr<double>(_Unicode(thickness));

      for (int rep = 0; rep < repeat; rep++, layer_idx++) {

        std::string  layer_name = st_name + "_layer_" + std::to_string(layer_idx);
        Box          layer_box(env_w / 2.0, env_h / 2.0, layer_thick / 2.0);
        Volume       layer_vol(layer_name, layer_box, mat_air);
        layer_vol.setVisAttributes(description, "InvisibleWithDaughters");

        double       layer_center_z = cur_z + layer_thick / 2.0;
        PlacedVolume layer_pv = env_vol.placeVolume(layer_vol,
                                  Position(0.0, 0.0, layer_center_z));
        layer_pv.addPhysVolID("layer", layer_idx);

        DetElement layer_de(st_de, layer_name, layer_idx);
        layer_de.setPlacement(layer_pv);

        // -- Slices ----------------------------------------------------------
        double local_z        = -layer_thick / 2.0;
        int    slice_in_layer = 0;

        for (xml_coll_t sColl(x_layer, _Unicode(slice)); sColl; ++sColl, ++slice_in_layer) {
          xml_comp_t x_slice = sColl;

          double      thick    = x_slice.attr<double>(_Unicode(thickness));
          std::string mat_name = x_slice.attr<std::string>(_Unicode(material));
          bool        is_sens  = x_slice.hasAttr(_Unicode(sensitive))
                                   ? x_slice.attr<bool>(_Unicode(sensitive)) : false;

          int plane_id = -1;
          if (is_sens && x_slice.hasAttr(_Unicode(plane)))
            plane_id = x_slice.attr<int>(_Unicode(plane));

          // Slice 0 = outer iron absorber: register for magnetic field
          if (slice_in_layer == 0 && bFieldY != 0.0) {
            double slab_lo = z0 + cur_z;           // global Z start of outer iron
            double slab_hi = slab_lo + thick;       // global Z end
            ironSlabs.push_back({slab_lo, slab_hi, env_w / 2.0 + 1.0, env_h / 2.0 + 1.0});
          }

          std::string slice_name   = layer_name + "_slice_" + std::to_string(slice_in_layer);
          Material    mat          = description.material(mat_name);
          double      sl_center_z  = local_z + thick / 2.0;

          // SciFi planes (plane 0 = U, plane 1 = V): use Para to encode the stereo
          // angle as a shear, avoiding envelope overflow from a pure Z-rotation.
          // All other slices remain rectangular.
          Solid sl_solid;
          if (plane_id == 0)
            sl_solid = Solid(new TGeoPara((slice_name + "_shape").c_str(),
                                          scifi_dx, env_h / 2.0, thick / 2.0,
                                          alpha_para, 0.0, 0.0));
          else if (plane_id == 1)
            sl_solid = Solid(new TGeoPara((slice_name + "_shape").c_str(),
                                          scifi_dx, env_h / 2.0, thick / 2.0,
                                         -alpha_para, 0.0, 0.0));
          else
            sl_solid = Box(env_w / 2.0, env_h / 2.0, thick / 2.0);

          Volume sl_vol(slice_name, sl_solid, mat);

          if (x_slice.hasAttr(_Unicode(vis)))
            sl_vol.setVisAttributes(description, x_slice.attr<std::string>(_Unicode(vis)));

          if (is_sens)
            sl_vol.setSensitiveDetector(sens);

          PlacedVolume sl_pv = layer_vol.placeVolume(sl_vol,
                                 Transform3D(Position(0.0, 0.0, sl_center_z)));
          sl_pv.addPhysVolID("slice", slice_in_layer);
          if (plane_id >= 0)
            sl_pv.addPhysVolID("plane", plane_id);

          DetElement slice_de(layer_de, slice_name,
                              layer_idx * 10 + slice_in_layer);
          slice_de.setPlacement(sl_pv);

          local_z += thick;
        }

        cur_z += layer_thick;
      }
    }
  }

  // -- Add magnetic field in outer iron slabs to the global DD4hep field ----
  if (bFieldY != 0.0 && !ironSlabs.empty()) {
    auto* obj        = new MTCIronField();
    obj->field_type  = CartesianField::MAGNETIC;
    obj->m_by        = bFieldY;
    obj->m_slabs     = std::move(ironSlabs);

    CartesianField cf;
    cf.assign(obj, name + "_BField", "magnetic");
    description.field().add(cf);

    dd4hep::printout(dd4hep::INFO, "MTCDetector",
        "%s: By=%.2f T magnetic field registered in %zu outer iron slabs",
        name.c_str(), bFieldY / dd4hep::tesla, obj->m_slabs.size());
  }

  return sdet;
}

DECLARE_DETELEMENT(MTCDetector, create_MTCDetector)
