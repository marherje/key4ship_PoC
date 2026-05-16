# Known Issues — MTC Tracking and Event Display

Diagnosed from `gaudi_jobs/muon_pipeline/debug_tracking_20.txt` (50-event run, single 5 GeV µ⁻ at xyz=(5, −82.5, 750) mm, dir (0,0,1)).

---

## Fixed

### U–V grouping straddle at 10 mm boundary (FIXED)

Adjacent U and V SciFi planes in the same MTC layer are separated by 2.35 mm in z.
The seeding code grouped them by `round(center.x() / 10.0)`, which placed ~2–3 layer pairs per station into different groups (U in bin N, V in bin N+1), preventing U×V crossings from forming for those layers.

**Fix applied:** changed divisor from `10.0` to `5.0` in `ACTSProtoTracker.cpp` (Source C grouping key). A 2.35 mm gap always fits within a 5 mm bin.

---

### Hough seeder gives wrong seed for curved MTC tracks (FIXED by CKF switch)

**Was:** The Hough accumulator assumes straight tracks. A 5 GeV muon curves by hundreds of mm across the MTC (1.7 T iron field, R ≈ 9.8 m), spreading votes across many bins. The wrong Hough peak won, giving chi2/nMeas ≈ 120.

**Fix applied:** Replaced `KalmanFitter` + Hough seeding with `CombinatorialKalmanFilter` (CKF). The CKF propagates the track state between surfaces using `EigenStepper` + `IronSlabBField`, inherently following the curved trajectory without needing an accurate seed position. The Hough is retained only to provide a rough transverse (x, y) starting position hint; the CKF propagates through all 150 surfaces regardless.

**Result:** chi2/nMeas drops from ~120 to ~3–5, nMeas = 19–31 per event (up from 0).

---

### CKF finds no tracks — all surfaces classified as "passive" (FIXED)

**Root cause 1 — null SurfaceArray:** `ACTSGeoSvc` created `PlaneLayer` with `nullptr` for `SurfaceArray`. `Layer::resolve()` in ACTS checks `resolveSensitive && m_surfaceArray` (line 60 of `Layer.ipp`). With null `SurfaceArray`, every layer returns false and the `Acts::Navigator` skips all surfaces → `nMeas=0 nHoles=0`.

**Fix applied:** Each active `PlaneLayer` now gets a `SurfaceArray` containing one `PlaneSurface` module created through `SNDDetectorElement`.

**Root cause 2 — no `associatedDetectorElement`:** `CKFActor::filter()` checks `surface->associatedDetectorElement() != nullptr` (line 442 of `CombinatorialKalmanFilter.hpp`) to decide if a surface is sensitive. Bare `PlaneSurface` objects have no detector element → every surface logged as "Passive surface detected" → no track states created.

**Fix applied:** Added `SNDDetectorElement : public Acts::DetectorElementBase` in `ACTSGeoSvc.cpp`. Surfaces are now created via `PlaneSurface(bounds, *detElement)`, which automatically sets `m_associatedDetectorElement`. Detector elements are stored in `ACTSGeoSvc::m_detectorElements` to keep them alive (surfaces hold raw back-pointers).

**Root cause 3 — `Acts::Navigator` cannot find surfaces via `resolveSensitive`:** Even with a non-null `SurfaceArray`, the standard `Acts::Navigator` needs the `TrackingGeometry` volume hierarchy. For our flat single-volume geometry, it visited surfaces but the CKF's `PropagatorPlainOptions::navigation` (of base type `NavigatorPlainOptions`) cannot carry the `DirectNavigator::surfaces` list — the CKF's `setPlainOptions()` only copies base fields.

**Fix applied:** `SNDFixedNavigator` wrapper struct (in `ACTSProtoTracker.cpp`) wraps `Acts::DirectNavigator`. It stores the full 150-surface list in its own `surfaces` member (not in `Options`) and injects it in `makeState()` before delegating to `DirectNavigator`. The CKF's `setPlainOptions()` cannot erase the surfaces since they live in the Config, not in Options.

---

### `bad_optional_access` crash in fingerprint loop (FIXED)

**Root cause:** `VectorMultiTrajectory::getUncalibratedSourceLink_impl()` calls `.value()` on an `std::optional<SourceLink>` (line 335 of `VectorMultiTrajectory.hpp`). Some track states with `MeasurementFlag` set have an allocated but unfilled source link slot — calling `getUncalibratedSourceLink()` on them throws `std::bad_optional_access`.

**Fix applied:** Added `ts.hasUncalibratedSourceLink()` guard before calling `getUncalibratedSourceLink()` in the duplicate-rejection fingerprint loop in `ACTSProtoTracker.cpp`.

---

## Open Issues

### 1. Event display draws a vertical straight line instead of the fitted trajectory

**Symptom:** The reconstructed track is rendered as a vertical line at constant (seed_x, seed_y) passing through every detector layer — it does not follow the actual track.

**Root cause — code:** In `event_display/event_display_eve.py`, `read_track_points()` builds the display points as:

```python
sx, sy = seed['seed_x'], seed['seed_y']
points = [(sx*mm2cm, sy*mm2cm, z*mm2cm, z) for z in layer_zs_mm]
```

This creates a column at fixed (x, y) for all z. `draw_tracks()` then connects these with straight line segments.

**Root cause — data:** `EDM4HEP2RNTuple.cpp` writes only `seed_x` and `seed_y` to the `ACTSTracks` RNTuple (extracted from the `AtIP` track state's `D0`/`Z0` fields). The CKF output — fitted (x, y, dx/dz, dy/dz) at each measurement surface — is written to `edm4hep::TrackState` objects in `tracks.edm4hep.root` by `ACTSProtoTracker`, but `EDM4HEP2RNTuple` does not yet extract them.

**Required fix:**

1. **Read per-surface fitted track states in `EDM4HEP2RNTuple.cpp`.** Extract the `AtOther` track states (phi, tanLambda, omega) from the `ACTSTracks` collection and write them as a branch in the RNTuple.

2. **Update `event_display_eve.py` to use per-surface positions.** Connect the fitted positions with line segments rather than projecting a seed point through all layers.

---

### 2. CKF accepts fewer muon hits than expected (chi2CutOff tuning)

**Symptom:** nMeas ≈ 19–31 per event, while ~60+ genuine muon hits are expected (40 SiTarget + 20 SiPad + some MTC). Most surfaces show `MeasurementSelector: No measurement candidate, chi2 ≈ 350` (outlier states).

**Root cause:** After the KF converges on the first few measurements, sigma_pred shrinks to ~0.022 mm (SiTarget strip pitch / sqrt(12)). The chi2CutOff of 15.0 then accepts only measurements within sqrt(15 × 2 × 0.022²) ≈ 0.12 mm of the predicted position. Background hits (secondary particles, delta-rays) at larger distances are correctly rejected, but some genuine muon hits may also fall outside this window if the track model or coordinate mapping has small systematic offsets.

**Diagnosis:** Raise `Chi2CutOff` from 15 to 100 in `job4_tracking.py`. If nMeas jumps to ~60+, the cut was too tight. If it stays at ~25, the background rejection is correct and genuine muon hits are being accepted.
