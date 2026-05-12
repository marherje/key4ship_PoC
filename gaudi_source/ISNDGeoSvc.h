#pragma once
#include "k4ActsTracking/IActsGeoSvc.h"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <vector>

class ISNDGeoSvc : virtual public IActsGeoSvc {
public:
  DeclareInterfaceID(ISNDGeoSvc, 1, 0);

  // Returns all rectangle surfaces sorted by X, for all detectors combined.
  virtual const std::vector<const Acts::Surface*>& allSurfaces() const = 0;

  // Returns the ACTS geometry context (needed by algorithms in execute()).
  virtual const Acts::GeometryContext& geometryContext() const = 0;

  // Returns the ACTS surface for the given detector address, or nullptr.
  // detID: 0=SiTarget, 1=SiPad, 2=MTC
  // station: station index (-1 for detectors with no stations)
  // layer:   layer index within station (or within detector for SiTarget/SiPad)
  // plane:   0/1 for strips, -1 for pixels
  virtual const Acts::Surface* surfaceByAddress(
      int detID, int station, int layer, int plane) const = 0;

  virtual ~ISNDGeoSvc() {}
};
