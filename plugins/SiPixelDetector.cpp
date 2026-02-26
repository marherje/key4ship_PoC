//====================================================================
//  DD4hep detector plugin for SiPixelDetector
//  Based on CaloPrototype_v02 pattern, updated for DD4hep v01-35
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"

using namespace dd4hep;
using namespace dd4hep::detail;

dd4hep::Ref_t create_SiPixelDetector(dd4hep::Detector& description,
                             xml_h element,
                             dd4hep::SensitiveDetector sens) {

  xml_det_t   x_det    = element;
  std::string det_name = x_det.nameStr();
  DetElement  sdet(det_name, x_det.id());

  Layering   layering(x_det);
  xml_dim_t  dim = x_det.dimensions();

  // Create envelope volume and place it into the world
  Volume envelope = xml::createPlacedEnvelope(description, element, sdet);
  xml::setDetectorTypeFlag(element, sdet);

  if (description.buildType() == BUILD_ENVELOPE) return sdet;

  Material air = description.air();
  sens.setType("calorimeter");

  double cal_hx = dim.x() / 2.0;
  double cal_hy = dim.y() / 2.0;
  double cal_hz = dim.z() / 2.0;

  int    layer_num   = 0;
  int    layerType   = 0;
  double layer_pos_z = -cal_hz;

  for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
    xml_comp_t x_layer        = c;
    int        repeat         = x_layer.repeat();
    const Layer* lay          = layering.layer(layer_num);
    double       layer_thickness = lay->thickness();
    std::string  layer_type_name = _toString(layerType, "layerType%d");

    for (int j = 0; j < repeat; j++) {
      std::string layer_name = _toString(layer_num, "layer%d");
      DetElement  layer(layer_name, layer_num);

      Volume layer_vol(layer_type_name, Box(cal_hx, cal_hy, layer_thickness / 2.0), air);

      double slice_pos_z  = -(layer_thickness / 2.0);
      int    slice_number = 0;

      for (xml_coll_t k(x_layer, _U(slice)); k; ++k) {
        xml_comp_t  x_slice        = k;
        std::string slice_name     = _toString(slice_number, "slice%d");
        double      slice_thickness = x_slice.thickness();
        Material    slice_mat      = description.material(x_slice.materialStr());
        DetElement  slice(layer, slice_name, slice_number);

        slice_pos_z += slice_thickness / 2.0;

        Volume slice_vol(slice_name, Box(cal_hx, cal_hy, slice_thickness / 2.0), slice_mat);

        if (x_slice.isSensitive()) {
          slice_vol.setSensitiveDetector(sens);
        }

        slice_vol.setAttributes(description,
                                x_slice.regionStr(),
                                x_slice.limitsStr(),
                                x_slice.visStr());

        PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol,
                                                       Position(0.0, 0.0, slice_pos_z));
        slice_phv.addPhysVolID("slice", slice_number);
        slice.setPlacement(slice_phv);

        slice_pos_z += slice_thickness / 2.0;
        ++slice_number;
      }

      layer_vol.setAttributes(description,
                              x_layer.regionStr(),
                              x_layer.limitsStr(),
                              x_layer.visStr());

      layer_pos_z += layer_thickness / 2.0;

      PlacedVolume layer_phv = envelope.placeVolume(layer_vol,
                                                    Position(0.0, 0.0, layer_pos_z));
      layer_phv.addPhysVolID("layer", layer_num);
      layer.setPlacement(layer_phv);

      layer_pos_z += layer_thickness / 2.0;
      ++layer_num;
    }
    ++layerType;
  }

  return sdet;
}

DECLARE_DETELEMENT(SiPixelDetector, create_SiPixelDetector)
