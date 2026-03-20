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

  virtual ~ISNDGeoSvc() {}
};
