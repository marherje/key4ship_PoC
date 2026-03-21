#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "DD4hep/Printout.h"
#include <cmath>
#include <string>

using namespace dd4hep;

dd4hep::Ref_t create_SiTarget(dd4hep::Detector& description,
                               xml_h element,
                               dd4hep::SensitiveDetector sens) {

  xml_det_t   x_det  = element;
  std::string name   = x_det.nameStr();
  int         det_id = x_det.id();
  DetElement  sdet(name, det_id);

  sens.setType("calorimeter");
  xml::setDetectorTypeFlag(element, sdet);   // needed by ACTS convertDD4hepDetector

  // -- Global envelope parameters -------------------------------------------
  xml_comp_t x_par = x_det.child(_Unicode(parameters));
  double env_w = x_par.attr<double>(_Unicode(env_width));
  double env_h = x_par.attr<double>(_Unicode(env_height));
  double x0    = x_par.attr<double>(_Unicode(x_position));
  Material mat_air = description.air();

  // -- First pass: compute total Y extent ------------------------------------
  double total_x = 0.0;
  for (xml_coll_t lColl(x_det, _Unicode(layer)); lColl; ++lColl) {
    xml_comp_t x_layer = lColl;
    int    repeat  = x_layer.hasAttr(_Unicode(repeat))
                       ? x_layer.attr<int>(_Unicode(repeat)) : 1;
    double spacing = x_layer.hasAttr(_Unicode(spacing))
                       ? x_layer.attr<double>(_Unicode(spacing)) : 0.0;
    double layer_thick = 0.0;
    for (xml_coll_t sColl(x_layer, _Unicode(slice)); sColl; ++sColl)
      layer_thick += xml_comp_t(sColl).attr<double>(_Unicode(thickness));
    total_x += repeat * (layer_thick + spacing);
  }

  // -- Global envelope -------------------------------------------------------
  double env_cx = x0;
  Box    env_shape( total_x / 2.0, env_w / 2.0, env_h / 2.0);
  Volume env_vol(name + "_envelope", env_shape, mat_air);
  env_vol.setVisAttributes(description, "InvisibleNoDaughters");

  Volume       motherVol = description.pickMotherVolume(sdet);
  PlacedVolume env_pv   = motherVol.placeVolume(env_vol, Position(env_cx, 0.0, 0.0));
  env_pv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_pv);

  // -- Global counters -------------------------------------------------------
  int    layer_idx = 0;
  double cur_x     = -total_x / 2.0;

  // -- Loop over <layer> -----------------------------------------------------
  for (xml_coll_t lColl(x_det, _Unicode(layer)); lColl; ++lColl) {
    xml_comp_t x_layer = lColl;

    int    repeat  = x_layer.hasAttr(_Unicode(repeat))
                       ? x_layer.attr<int>(_Unicode(repeat)) : 1;
    double spacing = x_layer.hasAttr(_Unicode(spacing))
                       ? x_layer.attr<double>(_Unicode(spacing)) : 0.0;

    double layer_thick = 0.0;
    for (xml_coll_t sColl(x_layer, _Unicode(slice)); sColl; ++sColl)
      layer_thick += xml_comp_t(sColl).attr<double>(_Unicode(thickness));

    for (int rep = 0; rep < repeat; rep++, layer_idx++) {

      std::string  layer_name = name + "_layer_" + std::to_string(layer_idx);
      Box          layer_box(layer_thick / 2.0, env_w / 2.0, env_h / 2.0);
      Volume       layer_vol(layer_name, layer_box, mat_air);
      layer_vol.setVisAttributes(description, "InvisibleWithDaughters");

      double       layer_center_x = cur_x + layer_thick / 2.0;
      PlacedVolume layer_pv = env_vol.placeVolume(layer_vol,
                                Position(layer_center_x, 0.0, 0.0));
      layer_pv.addPhysVolID("layer", layer_idx);

      DetElement layer_de(sdet, layer_name, layer_idx);
      layer_de.setPlacement(layer_pv);

      // -- Loop over <slice> -----------------------------------------------
      double local_x        = -layer_thick / 2.0;
      int    slice_in_layer = 0;

      for (xml_coll_t sColl(x_layer, _Unicode(slice)); sColl; ++sColl, ++slice_in_layer) {
        xml_comp_t x_slice = sColl;

        double      thick    = x_slice.attr<double>(_Unicode(thickness));
        std::string mat_name = x_slice.attr<std::string>(_Unicode(material));
        bool        is_sens  = x_slice.hasAttr(_Unicode(sensitive))
                                 ? x_slice.attr<bool>(_Unicode(sensitive))
                                 : false;

        int plane_id = 0;
        if (is_sens && x_slice.hasAttr(_Unicode(plane)))
          plane_id = x_slice.attr<int>(_Unicode(plane));

        std::string  slice_name = layer_name + "_slice_" + std::to_string(slice_in_layer);
        Material     mat        = description.material(mat_name);
        Box          sl_box(thick / 2.0, env_w / 2.0, env_h / 2.0);
        Volume       sl_vol(slice_name, sl_box, mat);

        if (x_slice.hasAttr(_Unicode(vis)))
          sl_vol.setVisAttributes(description, x_slice.attr<std::string>(_Unicode(vis)));

        if (is_sens)
          sl_vol.setSensitiveDetector(sens);

        double       sl_center_x = local_x + thick / 2.0;
        PlacedVolume sl_pv;
        if (is_sens && plane_id == 1) {
          // Rotate Y plane 90 degrees around X so that local Z -> global Y.
          // CartesianStripX will then measure strips along global Y.
          sl_pv = layer_vol.placeVolume(sl_vol,
                    Transform3D(RotationZYX(0.0, 0.0, M_PI / 2.0),
                                Position(sl_center_x,0.0, 0.0)));
        } else {
          sl_pv = layer_vol.placeVolume(sl_vol,
                    Position(sl_center_x,0.0, 0.0));
        }

        sl_pv.addPhysVolID("slice",  slice_in_layer);
       
        if (is_sens) {
          sl_pv.addPhysVolID("plane", plane_id);
          DetElement sl_de(layer_de,
                           "plane_" + std::to_string(plane_id),
                           layer_idx * 10 + plane_id);
          sl_de.setPlacement(sl_pv);
        }

        local_x += thick;
      }

      cur_x += layer_thick + spacing;
    }
  }

  return sdet;
}

DECLARE_DETELEMENT(SiTargetDetector, create_SiTarget)
