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
  double z0    = x_par.attr<double>(_Unicode(z_position));
  Material mat_air = description.air();

  // -- First pass: compute total Y extent ------------------------------------
  double total_z = 0.0;
  for (xml_coll_t lColl(x_det, _Unicode(layer)); lColl; ++lColl) {
    xml_comp_t x_layer = lColl;
    int    repeat  = x_layer.hasAttr(_Unicode(repeat))
                       ? x_layer.attr<int>(_Unicode(repeat)) : 1;
    double spacing = x_layer.hasAttr(_Unicode(spacing))
                       ? x_layer.attr<double>(_Unicode(spacing)) : 0.0;
    double layer_thick = 0.0;
    for (xml_coll_t sColl(x_layer, _Unicode(slice)); sColl; ++sColl)
      layer_thick += xml_comp_t(sColl).attr<double>(_Unicode(thickness));
    total_z += repeat * (layer_thick + spacing);
  }

  // -- Global envelope -------------------------------------------------------
  double env_cz = z0;
  // Z=beam: envelope is thin along Z (beam), wide in X and Y (transverse).
  // Add a small safety margin so that daughter volumes do not touch the
  // envelope boundary exactly (TGeo flags exact-boundary as "outside").
  const double safety = 0.1;  // mm
  Box    env_shape( env_w / 2.0 + safety, env_h / 2.0 + safety, total_z / 2.0 + safety);
  Volume env_vol(name + "_envelope", env_shape, mat_air);
  env_vol.setVisAttributes(description, "InvisibleNoDaughters");

  Volume       motherVol = description.pickMotherVolume(sdet);
  PlacedVolume env_pv   = motherVol.placeVolume(env_vol, Position(0.0, 0.0, env_cz));
  env_pv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_pv);

  // -- Global counters -------------------------------------------------------
  int    layer_idx = 0;
  double cur_z     = -total_z / 2.0;

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
      // Z=beam: layer box is wide in X and Y (transverse), thin along Z (beam).
      Box          layer_box(env_w / 2.0, env_h / 2.0, layer_thick / 2.0);
      Volume       layer_vol(layer_name, layer_box, mat_air);
      layer_vol.setVisAttributes(description, "InvisibleWithDaughters");

      double       layer_center_z = cur_z + layer_thick / 2.0;
      PlacedVolume layer_pv = env_vol.placeVolume(layer_vol,
                                Position(0.0, 0.0, layer_center_z));
      layer_pv.addPhysVolID("layer", layer_idx);

      DetElement layer_de(sdet, layer_name, layer_idx);
      layer_de.setPlacement(layer_pv);

      // -- Loop over <slice> -----------------------------------------------
      double local_z        = -layer_thick / 2.0;
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
        // Z=beam: slice box wide in X and Y (transverse), thin along Z (beam).
        Box          sl_box(env_w / 2.0, env_h / 2.0, thick / 2.0);
        Volume       sl_vol(slice_name, sl_box, mat);

        if (x_slice.hasAttr(_Unicode(vis)))
          sl_vol.setVisAttributes(description, x_slice.attr<std::string>(_Unicode(vis)));

        if (is_sens)
          sl_vol.setSensitiveDetector(sens);

        double       sl_center_z = local_z + thick / 2.0;
        // Z=beam: no rotation needed. Both plane=0 (StripX) and plane=1
        // (StripY) are transverse slices; normal is already along Z.
        PlacedVolume sl_pv = layer_vol.placeVolume(sl_vol,
                    Position(0.0, 0.0, sl_center_z));

        sl_pv.addPhysVolID("slice",  slice_in_layer);
       
        if (is_sens) {
          sl_pv.addPhysVolID("plane", plane_id);
          DetElement sl_de(layer_de,
                           "plane_" + std::to_string(plane_id),
                           layer_idx * 10 + plane_id);
          sl_de.setPlacement(sl_pv);
        }

        local_z += thick;
      }

      cur_z += layer_thick + spacing;
    }
  }

  return sdet;
}

DECLARE_DETELEMENT(SiTargetDetector, create_SiTarget)
