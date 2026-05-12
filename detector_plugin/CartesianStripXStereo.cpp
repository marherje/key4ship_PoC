#include <DDSegmentation/CartesianStripX.h>
#include <DDSegmentation/SegmentationParameter.h>
#include <DD4hep/Factories.h>
#include <DD4hep/detail/SegmentationsInterna.h>
#include <cmath>

// How DECLARE_SEGMENTATION(Name, func) works:
//   DD4HEP_OPEN_PLUGIN(dd4hep, Name) creates a stub "struct Name {}" inside
//   dd4hep::{anonymous}, making 'Name' the plugin lookup key.
//   'func' must be a factory that creates the REAL implementation object.
//   Therefore the real class must have a DIFFERENT name from the plugin key.

namespace dd4hep {

/// CartesianStripX variant that computes the strip from the stereo-corrected
/// coordinate (local_x - local_y * tan(stereo_angle_deg)).
/// position() returns the strip centre at y=0 (volume reference plane).
class CartesianStripXStereoImpl : public DDSegmentation::CartesianStripX {
public:
  CartesianStripXStereoImpl(const std::string& cellEncoding = "")
      : DDSegmentation::CartesianStripX(cellEncoding) { init(); }
  CartesianStripXStereoImpl(const DDSegmentation::BitFieldCoder* decoder)
      : DDSegmentation::CartesianStripX(decoder) { init(); }
  virtual ~CartesianStripXStereoImpl() = default;

  virtual DDSegmentation::CellID
  cellID(const DDSegmentation::Vector3D& localPosition,
         const DDSegmentation::Vector3D& /*globalPosition*/,
         const DDSegmentation::VolumeID& volumeID) const override {
    const double tan_a    = std::tan(_stereoAngle * M_PI / 180.0);
    const double x_stereo = localPosition.X + localPosition.Y * tan_a;
    DDSegmentation::CellID cID = volumeID;
    _decoder->set(cID, _xId, positionToBin(x_stereo, _stripSizeX, _offsetX));
    return cID;
  }

  virtual DDSegmentation::Vector3D
  position(const DDSegmentation::CellID& cID) const override {
    return { binToPosition(_decoder->get(cID, _xId), _stripSizeX, _offsetX),
             0., 0. };
  }

private:
  void init() {
    _type        = "CartesianStripXStereo";
    _description = "CartesianStripX with stereo tilt: strip from (x - y*tan(angle_deg))";
    registerParameter("stereo_angle",
                      "Stereo tilt in degrees (+U plane, -V plane)",
                      _stereoAngle, 0.0,
                      DDSegmentation::SegmentationParameter::NoUnit);
  }
  double _stereoAngle = 0.0;
};

} // namespace dd4hep

// Plugin name "CartesianStripXStereo" is what the XML type= attribute resolves to.
// create_segmentation<CartesianStripXStereoImpl> creates the real implementation.
// (Unqualified create_segmentation reaches dd4hep::{anonymous} from inside namespace dd4hep.)
DECLARE_SEGMENTATION(CartesianStripXStereo,
    create_segmentation<dd4hep::CartesianStripXStereoImpl>)
