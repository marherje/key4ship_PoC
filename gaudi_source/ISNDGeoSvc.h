#pragma once
#include "k4ActsTracking/IActsGeoSvc.h"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <vector>

class ISNDGeoSvc : virtual public IActsGeoSvc {
public:
  DeclareInterfaceID(ISNDGeoSvc, 1, 2);

  // Returns all rectangle surfaces sorted by X, for all detectors combined.
  virtual const std::vector<const Acts::Surface*>& allSurfaces() const = 0;

  // Returns the surfaces of ONE detector, sorted by X (beam) like allSurfaces().
  // detID: 0=SiTarget, 1=SiPad, 2=MTC; anything else yields an empty vector.
  //
  // This is the authoritative per-detector partition: it is recorded while the
  // TGeo tree is walked, where the detector identity is known from the volume
  // name. Deriving it afterwards from the largest z-gaps in allSurfaces() is
  // wrong in general — the MTC inter-station gaps can exceed the inter-detector
  // ones — so any algorithm needing "just this detector's surfaces" must ask
  // here rather than re-deriving it.
  virtual const std::vector<const Acts::Surface*>& surfacesOf(int detID) const = 0;

  // Returns the ACTS geometry context (needed by algorithms in execute()).
  virtual const Acts::GeometryContext& geometryContext() const = 0;

  // Returns the ACTS surface for the given detector address, or nullptr.
  // detID: 0=SiTarget, 1=SiPad, 2=MTC
  // station: station index (-1 for detectors with no stations)
  // layer:   layer index within station (or within detector for SiTarget/SiPad)
  // plane:   0/1 for strips, -1 for pixels
  virtual const Acts::Surface* surfaceByAddress(
      int detID, int station, int layer, int plane) const = 0;

  // Returns the ordinal of the layer this surface belongs to, counted within
  // its own detector along the beam (0 = most upstream), or -1 if unknown.
  //
  // Not the same as the raw `layer` of surfaceByAddress(): that one restarts at
  // 0 in every MTC station, so it cannot answer "how many layers lie between
  // these two surfaces". The ordinal is contiguous across a whole detector, so
  // subtracting two of them is a layer count. All surfaces of one physical
  // layer share an ordinal — the MTC's U (plane=0), V (plane=1) and combined
  // (plane=3) surfaces map to a single layer, which is what makes a hit count
  // per layer independent of whether the U/V pair was merged.
  virtual int layerOrdinalOf(const Acts::Surface* surface) const = 0;

  // Number of distinct layers in a detector, i.e. one past the largest ordinal
  // layerOrdinalOf() can return for it. 0 for an unknown detID.
  virtual int nLayersOf(int detID) const = 0;

  virtual ~ISNDGeoSvc() {}
};
