#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include <string>

using namespace dd4hep;

// Passive emulsion target: port of FairShip SND/EmulsionTarget/Target.cxx
// (Design 4 default: 5 walls x (36 W plates + 37 emulsion films), passive).
// No sensitive volumes, no readout, no segmentation, no digitization.
dd4hep::Ref_t create_EmulsionTarget(dd4hep::Detector& description,
                                    xml_h element,
                                    dd4hep::SensitiveDetector /*sens*/) {

  xml_det_t   x_det  = element;
  std::string name   = x_det.nameStr();
  int         det_id = x_det.id();
  DetElement  sdet(name, det_id);

  // -- Envelope, brick and film parameters ------------------------------------
  xml_comp_t x_par    = x_det.child(_Unicode(parameters));
  double env_x    = x_par.attr<double>(_Unicode(env_width));
  double env_y    = x_par.attr<double>(_Unicode(env_height));
  double env_z    = x_par.attr<double>(_Unicode(env_length));
  double z0       = x_par.attr<double>(_Unicode(z_position));
  double brick_x  = x_par.attr<double>(_Unicode(brick_width));
  double brick_y  = x_par.attr<double>(_Unicode(brick_height));
  double brick_z  = x_par.attr<double>(_Unicode(brick_length));
  int    n_walls  = x_par.attr<int>   (_Unicode(n_walls));
  double tt_gap   = x_par.attr<double>(_Unicode(tt_gap));
  int    n_plates = x_par.attr<int>   (_Unicode(n_plates));
  double w_thick  = x_par.attr<double>(_Unicode(tungsten_thickness));
  double em_thick = x_par.attr<double>(_Unicode(emulsion_thickness));
  double pb_thick = x_par.attr<double>(_Unicode(base_thickness));

  double film_thick = 2.0 * em_thick + pb_thick;   // FairShip EmPlateWidth
  double cell_thick = film_thick + w_thick;         // FairShip AllPlateWidth

  Material mat_air = description.air();
  Material mat_w    = description.material("TungstenDens1910");
  Material mat_em   = description.material("NuclearEmulsion");
  Material mat_pb    = description.material("PlasticBase");

  // -- Envelope ----------------------------------------------------------------
  Box    env_shape(env_x / 2.0, env_y / 2.0, env_z / 2.0);
  Volume env_vol(name + "_envelope", env_shape, mat_air);
  env_vol.setVisAttributes(description, "InvisibleNoDaughters");

  Volume       motherVol = description.pickMotherVolume(sdet);
  PlacedVolume env_pv    = motherVol.placeVolume(env_vol, Position(0.0, 0.0, z0));
  env_pv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_pv);

  // -- Shared brick volume: 36 tungsten plates interleaved with 37 emulsion
  //    films (each film = emulsion + plastic base + emulsion), replicating
  //    the FairShip Target::ConstructGeometry layout. -----------------------
  Volume brick_vol(name + "_brick",
                   Box(brick_x / 2.0, brick_y / 2.0, brick_z / 2.0), mat_air);
  brick_vol.setVisAttributes(description, "InvisibleWithDaughters");

  Volume abs_vol(name + "_absorber",
                 Box(brick_x / 2.0, brick_y / 2.0, w_thick / 2.0), mat_w);
  abs_vol.setVisAttributes(description, "TungstenVis");

  Volume em_vol(name + "_emulsion",
                Box(brick_x / 2.0, brick_y / 2.0, em_thick / 2.0), mat_em);
  em_vol.setVisAttributes(description, "SiVis");

  Volume base_vol(name + "_plasticbase",
                  Box(brick_x / 2.0, brick_y / 2.0, pb_thick / 2.0), mat_pb);
  base_vol.setVisAttributes(description, "BaseVis");

  for (int n = 0; n < n_plates; ++n)
    brick_vol.placeVolume(abs_vol,
        Position(0, 0, -brick_z / 2.0 + film_thick + w_thick / 2.0 + n * cell_thick));

  for (int n = 0; n < n_plates + 1; ++n) {
    brick_vol.placeVolume(em_vol,             // bottom emulsion
        Position(0, 0, -brick_z / 2.0 + em_thick / 2.0 + n * cell_thick));
    brick_vol.placeVolume(base_vol,           // plastic base
        Position(0, 0, -brick_z / 2.0 + em_thick + pb_thick / 2.0 + n * cell_thick));
    brick_vol.placeVolume(em_vol,             // top emulsion
        Position(0, 0, -brick_z / 2.0 + em_thick + pb_thick + em_thick / 2.0
                          + n * cell_thick));
  }

  // -- Place n_walls bricks along Z, separated by target-tracker gaps -------
  double zc = -env_z / 2.0 + tt_gap;
  for (int w = 0; w < n_walls; ++w) {
    env_vol.placeVolume(brick_vol, Position(0, 0, zc + brick_z / 2.0));
    zc += brick_z + tt_gap;
  }

  dd4hep::printout(dd4hep::INFO, "EmulsionTarget",
      "%s: %d walls, %d W plates/brick, passive, z-centre %.1f mm",
      name.c_str(), n_walls, n_plates, z0 / dd4hep::mm);

  return sdet;
}

DECLARE_DETELEMENT(EmulsionTargetDetector, create_EmulsionTarget)
