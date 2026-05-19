# ACTS Track Reconstruction Integration

## Overview

ACTS (A Common Tracking Software) is integrated via `ACTSGeoSvc` (a Gaudi Service) and `ACTSProtoTracker` (a Gaudi Algorithm), both compiled into `libSND_reco.so`.

ACTS version: 44.3.0 (from key4hep 2026-02-01 CVMFS stack).

---

## ACTSGeoSvc

**Files:** `gaudi_source/ACTSGeoSvc.h`, `gaudi_source/ACTSGeoSvc.cpp`  
**Interface:** `gaudi_source/ISNDGeoSvc.h` (extends `IActsGeoSvc` from k4ActsTracking)

### Purpose
Loads DD4hep geometry, walks the TGeo tree, and builds ACTS `TrackingGeometry` consisting of `PlaneSurface` objects for each sensitive detector plane.

### Configuration
```python
from Configurables import ACTSGeoSvc
geo = ACTSGeoSvc("ACTSGeoSvc")
geo.CompactFile    = "simulation/geometry/SND_compact.xml"  # required
geo.MTCStereoAngle = 5.0  # degrees, must match CartesianStripXStereo XML value
ApplicationMgr(..., ExtSvc=[geo])  # MUST be in ExtSvc, not TopAlg
```

### Interface methods
```cpp
// In algorithm initialize():
ServiceHandle<ISNDGeoSvc> m_geoSvc{this, "ACTSGeoSvc", "ACTSGeoSvc"};

const Acts::TrackingGeometry& geo  = m_geoSvc->trackingGeometry();
const Acts::GeometryContext&  gctx = m_geoSvc->geometryContext();
const std::vector<const Acts::Surface*>& surfs = m_geoSvc->allSurfaces();
const Acts::Surface* surf = m_geoSvc->surfaceByAddress(detID, station, layer, plane);
```

`surfaceByAddress(detID, station, layer, plane)` is the primary lookup used by measurement converters.

### Geometry construction (ACTSGeoSvc.cpp)
1. Loads DD4hep geometry via `dd4hep::Detector::getInstance().fromXML(...)`
2. Walks `gGeoManager` TGeo tree, finds sensitive volumes by detector name
3. Extracts global Z position and half-sizes from TGeo matrices (TGeo in **cm** → multiply ×10 for ACTS mm)
4. Constructs `Acts::PlaneSurface` with `Acts::RectangleBounds` for each plane
5. Wraps surfaces in `Acts::PlaneLayer` → `Acts::TrackingVolume` → `Acts::TrackingGeometry`

### Unit notes
- TGeo positions are in **cm**; all values multiplied by 10 before passing to ACTS
- ACTS uses **mm** and **GeV** throughout
- `Acts::UnitConstants::mm` = 1.0, `Acts::UnitConstants::GeV` = 1.0

---

## Measurement Converters

### SiTargetMeasConverter (`gaudi_source/SiTargetMeasConverter.cpp`)
Converts `SiTargetHitsWindowed` (EDM4HEP `TrackerHit3D`) to ACTS measurements on SiTarget surfaces.

### SiPadMeasConverter (`gaudi_source/SiPadMeasConverter.cpp`)
Converts `SiPadHitsWindowed` to measurements on SiPad surfaces.

### MTCSciFiMeasConverter (`gaudi_source/MTCSciFiMeasConverter.cpp`)
Converts `MTCSciFiHitsWindowed` (plane 0 = U, plane 1 = V) to 1D `TrackerHit3D` measurements on MTC SciFi surfaces. Scintillator hits (plane 2) are skipped. Strip coordinate is stored in `position.x` [mm]; variance = pitch²/12.

**Status:** Fully implemented and active in `muon_pipeline`. Not wired in other pipelines (see [Pipeline differences](#pipeline-differences--mtc-and-b-field)).

### SNDMeasurement struct
```cpp
struct SNDMeasurement {
  const Acts::Surface* surface;
  double localCoord;    // 1D strip coordinate (SiTarget/SciFi) or X (SiPad)
  double localCoord2;   // Y coordinate (SiPad only, is2D=true)
  double variance;      // measurement uncertainty²
  double variance2;     // Y uncertainty² (SiPad only)
  bool   is2D;          // true for SiPad (pad detector)
  int    detectorID;    // 0=SiTarget, 1=SiPad, 2=MTC
  int    plane;         // strip plane index
  float  time;          // hit time [ns]
};
```

---

## ACTSProtoTracker

**File:** `gaudi_source/ACTSProtoTracker.cpp`

### Algorithm flow
1. `initialize()` — retrieves `ACTSGeoSvc`; dynamically determines surface group sizes (SiTarget / SiPad / MTC) by finding the two largest Z-gaps in the surface list
2. `execute()` — per event:
   a. Collect measurements from SiTarget + SiPad + MTC converters (MTC active when `InputMTC` collection is present)
   b. `findSeeds()` — 2D Hough transform over all measurements (SiPad 2D hits, SiTarget X×Y crossings, MTC U×V stereo crossings)
   c. For each seed: run ACTS Kalman fitter → `Acts::TrackContainer`
   d. Write output tracks to `ACTSTracks` EDM4HEP collection
3. `finalize()` — log statistics

### Hough seeding
2D Hough transform over (x, z) or (y, z) plane. Parameters:
- `HoughBinSize` — accumulator bin width
- `HoughHalfSize` — accumulator range
- `HoughMinVotes` — min hits per bin to form a seed
- `AutoSeed` — if true, automatically constructs seeds from Hough peaks
- `MaxSeeds` — cap on number of seeds per event

### Kalman fitter setup

The B-field provider is chosen per event based on the `IronFieldRanges` property:

- **`IronFieldRanges` set** (muon_pipeline): uses `IronSlabBField`, a custom `Acts::MagneticFieldProvider` that returns a configurable `By` inside each registered rectangular slab (the MTC outer iron absorbers) and zero everywhere else.
- **`IronFieldRanges` empty** (other pipelines): falls back to `Acts::ConstantBField(BFieldX, BFieldY, BFieldZ)`, all defaulting to 0.

```cpp
// muon_pipeline: IronSlabBField with geometry-driven slab ranges
// other pipelines: ConstantBField(0, 0, 0)
auto stepper    = Acts::EigenStepper<>(bField);
auto navigator  = Acts::DirectNavigator{};
auto propagator = Acts::Propagator(std::move(stepper), std::move(navigator));
// KalmanFitter constructed per seed (propagator is moved in)
```

### Track output
Tracks are stored in an `edm4hep::TrackCollection` named `ACTSTracks`.  
Seeds (Hough candidates before fitting) are also stored in `ShipHits.root` via job5.

---

## Magnetic Field

The MTC contains iron absorber slabs that carry a magnetic field (for muon momentum measurement). `ACTSProtoTracker` models this with a custom `IronSlabBField` provider that returns `By` inside each iron slab volume and zero elsewhere.

The slab ranges and field strength are passed via `IronFieldRanges` — a flat list of `[xlo, xhi, ylo, yhi, zlo, zhi, by]` septets in ACTS coordinates (mm, T). In `muon_pipeline/job4_tracking.py` these are computed automatically from `parse_geometry.SNDGeometry` (which reads `SND_compact.xml`), so the field map stays in sync with the geometry.

The SiTarget and SiPad regions have no magnetic field; `IronSlabBField` returns zero there automatically.

Other pipelines (no `IronFieldRanges`) fall back to `Acts::ConstantBField(0, 0, 0)`.

---

## Pipeline differences — MTC and B-field

`muon_pipeline` is the complete reference pipeline. It is the only one that wires MTC and the iron B-field:

| Feature | `muon_pipeline` | other pipelines |
|---|---|---|
| `MTCSciFiMeasConverter` in `TopAlg` | Yes | No |
| `InputMTC` wired to `ACTSProtoTracker` | Yes (`MTCSciFiMeasurements`) | No (default collection absent) |
| MTC U×V crossings in Hough seeding | Yes | No |
| MTC 1D hits fed into Kalman fitter | Yes | No |
| `IronFieldRanges` (MTC iron B-field) | Yes, from `parse_geometry` | No |
| Effective B-field | `IronSlabBField` (non-zero By in iron) | `ConstantBField(0,0,0)` |

To bring another pipeline up to full feature parity, follow `gaudi_jobs/muon_pipeline/job4_tracking.py` as the reference.
