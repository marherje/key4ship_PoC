//====================================================================
//  DD4hep detector plugin for SiTarget (Silicon Target)
//  Based on FairShip SiliconTarget geometry.
//
//  Layout per station:
//    W absorber
//    X-strip Si plane  (2 rows × 4 columns of sensors, no rotation)
//    Y-strip Si plane  (2 rows × 4 columns of sensors, rotated 90° around z)
//
//  Both planes share one readout "SiTargetHits" (CartesianGridXY).
//  The plane physVolID (0=X, 1=Y) distinguishes them in the cellID.
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include <cmath>
#include <string>

using namespace dd4hep;

dd4hep::Ref_t create_SiTarget(dd4hep::Detector& description,
                               xml_h element,
                               dd4hep::SensitiveDetector sens) {

  xml_det_t   x_det = element;
  std::string name  = x_det.nameStr();
  int         id    = x_det.id();
  DetElement  sdet(name, id);

  // Read parameters from XML
  xml_comp_t x_par       = x_det.child(_Unicode(parameters));
  int    n_layers         = x_par.attr<int>(_Unicode(n_layers));
  double absorber_thick   = x_par.attr<double>(_Unicode(target_thickness));
  double layer_pitch      = x_par.attr<double>(_Unicode(target_spacing));
  double module_offset    = x_par.attr<double>(_Unicode(module_offset));
  double sensor_w         = x_par.attr<double>(_Unicode(sensor_width));
  double sensor_l         = x_par.attr<double>(_Unicode(sensor_length));
  double sensor_thick      = x_par.attr<double>(_Unicode(sensor_thickness));
  double env_w            = x_par.attr<double>(_Unicode(env_width));
  double env_h            = x_par.attr<double>(_Unicode(env_height));
  double z0               = x_par.attr<double>(_Unicode(z_position));
  std::string mat_W_name  = x_par.attr<std::string>(_Unicode(absorber_material));
  std::string mat_Si_name = x_par.attr<std::string>(_Unicode(sensor_material));

  // planeSpacing = targetSpacing - targetThickness - 2*moduleOffset
  const double plane_gap       = layer_pitch - absorber_thick - 2.0 * module_offset;
  // Full active area: 4 columns × 2 rows with 1 mm inter-sensor gaps
  const double full_plane_w    = 4.0 * sensor_w + 3.0 * dd4hep::mm;
  const double full_plane_h    = 2.0 * sensor_l + 1.0 * dd4hep::mm;

  // Materials (names come from XML parameters)
  Material mat_W   = description.material(mat_W_name);
  Material mat_Si  = description.material(mat_Si_name);
  Material mat_air = description.air();

  sens.setType("calorimeter");

  // --- Envelope ---
  double total_z  = (n_layers - 1) * layer_pitch
                    + absorber_thick + module_offset
                    + sensor_thick + plane_gap + sensor_thick
                    + 2.0 * dd4hep::mm;
  double env_cz   = z0 + total_z / 2.0 - 1.0 * dd4hep::mm;

  Box    env_shape(env_w / 2.0, env_h / 2.0, total_z / 2.0);
  Volume env_vol(name + "_env", env_shape, mat_air);

  Volume motherVol    = description.pickMotherVolume(sdet);
  PlacedVolume env_pv = motherVol.placeVolume(env_vol, Position(0.0, 0.0, env_cz));
  env_pv.addPhysVolID("system", id);
  sdet.setPlacement(env_pv);

  // --- Reusable volumes ---
  // W absorber spans the full envelope width/height
  Box    w_box(env_w / 2.0, env_h / 2.0, absorber_thick / 2.0);
  Volume vol_W(name + "_W", w_box, mat_W);

  // One sensitive slab per view covering the full active area (standard k4geo approach)
  Box    si_plane_box(full_plane_w / 2.0, full_plane_h / 2.0, sensor_thick / 2.0);
  Volume vol_si_plane(name + "_si", si_plane_box, mat_Si);
  vol_si_plane.setSensitiveDetector(sens);

  // --- Station loop ---
  for (int i = 0; i < n_layers; i++) {
    double gz_W_start  = z0 + i * layer_pitch;
    double gz_W_center = gz_W_start + absorber_thick / 2.0;
    double gz_X        = gz_W_start + absorber_thick + module_offset + sensor_thick / 2.0;
    double gz_Y        = gz_X + sensor_thick / 2.0 + plane_gap + sensor_thick / 2.0;

    double lz_W = gz_W_center - env_cz;
    double lz_X = gz_X        - env_cz;
    double lz_Y = gz_Y        - env_cz;

    // W absorber (not sensitive — no DetElement needed)
    PlacedVolume pv_w = env_vol.placeVolume(vol_W, Position(0.0, 0.0, lz_W));
    pv_w.addPhysVolID("layer", i);

    // One sensitive slab per view — avoids repeated-logical-volume ambiguity
    for (int plane = 0; plane < 2; plane++) {
      double sz = (plane == 0) ? lz_X : lz_Y;

      // Y plane: rotate 90° around z so the slab's local x stays the measurement axis
      Transform3D tf;
      if (plane == 0) {
        tf = Transform3D(RotationZYX(0, 0, 0),           Position(0.0, 0.0, sz));
      } else {
        tf = Transform3D(RotationZYX(M_PI / 2.0, 0, 0), Position(0.0, 0.0, sz));
      }

      PlacedVolume pv = env_vol.placeVolume(vol_si_plane, tf);
      pv.addPhysVolID("system", id);
      pv.addPhysVolID("layer",  i);
      pv.addPhysVolID("plane",  plane);

      DetElement plane_de(sdet, "plane_" + std::to_string(i * 2 + plane), i * 2 + plane);
      plane_de.setPlacement(pv);
    }
  }

  return sdet;
}

DECLARE_DETELEMENT(SiTargetDetector, create_SiTarget)
