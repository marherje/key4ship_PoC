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

## Fixed (continued)

### Event display draws a vertical straight line instead of the fitted trajectory (FIXED)

**Was:** `read_track_points()` built display points as `(seed_x, seed_y, z)` for all layer z values — a vertical line. `EDM4HEP2RNTuple` wrote only the `AtIP` seed state.

**Fix applied:**

1. **`ACTSProtoTracker.cpp`** — Each `AtOther` TrackState now has `referencePoint` filled:
   - `.x = smoothed[eBoundLoc0]` (DD4hep transverse X in mm)
   - `.y = smoothed[eBoundLoc1]` (DD4hep transverse Y in mm)
   - `.z = surface.center(gctx).x()` (beam Z in mm)
   Falls back to filtered parameters if no smoothed state.

2. **`EDM4HEP2RNTuple.cpp`** — Added `ACTSTrackStates` RNTuple (one row per surface per track): `window_id, track_id, state_idx, x, y, z` (all in mm, upstream→downstream order).

3. **`event_display_eve.py`** — `read_track_points()` now calls `read_track_states()` to read `ACTSTrackStates` and connects the actual fitted positions. Falls back to seed projection only if no states are found (e.g. old files).

---

### 2. CKF accepts fewer muon hits than expected (chi2CutOff tuning)

**Symptom:** nMeas ≈ 19–31 per event, while ~60+ genuine muon hits are expected (40 SiTarget + 20 SiPad + some MTC). Most surfaces show `MeasurementSelector: No measurement candidate, chi2 ≈ 350` (outlier states).

**Root cause:** After the KF converges on the first few measurements, sigma_pred shrinks to ~0.022 mm (SiTarget strip pitch / sqrt(12)). The chi2CutOff of 15.0 then accepts only measurements within sqrt(15 × 2 × 0.022²) ≈ 0.12 mm of the predicted position. Background hits (secondary particles, delta-rays) at larger distances are correctly rejected, but some genuine muon hits may also fall outside this window if the track model or coordinate mapping has small systematic offsets.

**Diagnosis:** Raise `Chi2CutOff` from 15 to 100 in `job4_tracking.py`. If nMeas jumps to ~60+, the cut was too tight. If it stays at ~25, the background rejection is correct and genuine muon hits are being accepted.

---

### MTC SciFi measurement value: cos(α) scale and stereo encoding (FIXED)

The ACTS surface for each MTC SciFi plane has rotation `R_X(±α)·R_Y(π/2)` (`ACTSGeoSvc.cpp`), making:

```
eBoundLoc0_U = cos(α)·dd_x − sin(α)·dd_y   (plane 0)
eBoundLoc0_V = cos(α)·dd_x + sin(α)·dd_y   (plane 1)
```

The 1D strip read-out from `CartesianStripXStereo` gives `x_stereo = dd_x ∓ dd_y·tan(α)`, which equals `eBoundLoc0 / cos(α)`. The stereo contribution (`±sin(α)·dd_y`) is already encoded in `x_stereo`; no separate recovery of `dd_y` is needed. The original issue (plain `CartesianStripX` stored only `dd_x`, missing the stereo term entirely) was resolved when the segmentation was changed to `CartesianStripXStereo`.

The remaining cos(α) scale factor (`x_stereo ≠ eBoundLoc0`) was corrected in `MTCSciFiMeasConverter.cpp`:

```cpp
const double localCoord = cos_a * pos.x;                         // eBoundLoc0 = cos(α)·x_stereo
const double variance   = cos_a * cos_a * (pitch * pitch) / 12.0; // σ²(eBoundLoc0)
```

For α = 5° the correction is 0.4%, but the explicit formula is correct for any stereo angle.

---

## Fixed (continued)

### MTC SciFi B-field in ACTS CKF (FIXED)

The `ACTSProtoTracker` `BFieldX/Y/Z` properties default to zero. The fix was to add `IronSlabBField` (a custom `Acts::MagneticFieldProvider`) and populate it from `SND_compact.xml` constants via `parse_geometry.py` in `job4_tracking.py`.

`IronFieldRanges` is set to 45 slabs (3 stations × 15 layers), each covering the 50 mm outer iron absorber at `By = 1.7 T`. Confirmed active in `debug_tracking_31.txt` (line 30: `ACTSProtoTracker.IronFieldRanges [245.0, 295.0, … 3583.8, -301.0, 301.0, -301.0, 301.0, 1.7]`).

`MTCSciFiMeasConverter` and `proto.InputMTC` are correctly wired in `muon_pipeline/job4_tracking.py`.

---

## Open issues

### PropagatorError:2 — propagator hits step-count limit in ~6% of events

**Symptom:** `Propagation reached the configured maximum number of steps with the initial parameters` (PropagatorError:2). Observed in `debug_tracking_31.txt`: events 3 and 43 (all 3 seeds fail → 0 tracks); event 24 seed 1 fails while seeds 0 and 2 succeed.

**Root cause:** `pOptions.maxSteps` was hardcoded to 10000 with `maxStepSize = 10 mm`. The CKF visits 150 surfaces across ~4 m of detector. For off-axis Hough seeds (event 24 seed 1: loc0 = −93 mm, loc1 = −200 mm, 1 mm inside the ±201 mm slab boundary; events 3/43: all three seeds near the correct muon position) the stepper exhausts its budget before completing the track. The exact trigger per event is not fully understood without step-level diagnostics.

**Fix applied:** `MaxPropSteps` and `MaxStepSize` are now Gaudi properties (defaults: 100 000 steps, 100 mm). At 100 mm steps in 1.7 T / 10 GeV iron the positional error per step is ~0.06 mm — well within the 1 mm SciFi pitch. The old 10 mm cap was unnecessarily tight and contributed to step-count exhaustion.

**To verify:** Re-run `muon_pipeline` after rebuild. Events 3, 24, 43 should no longer show PropagatorError:2.

---

### Event display: stereo U/V assignment was swapped (FIXED — but see "Latent risks" below)

> **Update 2026-05-17:** the diagnostic plots in
> `gaudi_jobs/muon_pipeline/plots_curvature/` show that with the **current**
> geometry the reco-y reconstructed by `(loc0_V − loc0_U)/(2·tan 5°)` has the
> opposite sign of truth-y again. See the latent-risk note below for why this
> "fix" is convention-dependent and how to test it.

**Symptom:** The reconstructed y coordinate from MTC SciFi stereo pairing had wrong sign in the event display.

**Root cause:** The stereo pairing code in `read_track_states()` assigned `tilt > 0` → U plane and `tilt < 0` → V plane. But the ACTSGeoSvc rotation `R_X(+α)·R_Y(π/2)` for plane 0 gives `stereoTiltZ = +sin(α)` and yields the V-type measurement (`eBoundLoc0 = cos·dd_x + sin·dd_y`), while `tilt < 0` (plane 1, `R_X(−α)`) gives U-type (`cos·dd_x − sin·dd_y`). The assignment was therefore backwards, flipping the sign of `gy = (loc0_V − loc0_U) / (2·tan 5°)`.

Verified: with `tilt > 0`, loc0 = −5.623 mm ≈ cos5°·1.5 + sin5°·(−82) = −5.66 mm → V plane. With `tilt < 0`, loc0 = +8.548 mm ≈ cos5°·1.5 − sin5°·(−82) = +8.64 mm → U plane.

**Fix applied** (`event_display_eve.py`):
```python
# Before (wrong):
loc0_U = loc0  if tilt  > 0 else loc02
loc0_V = loc02 if tilt  > 0 else loc0

# After (correct):
loc0_V = loc0  if tilt  > 0 else loc02   # positive tilt → V plane
loc0_U = loc02 if tilt  > 0 else loc0    # negative tilt → U plane
```

---

## Latent risks (silent failure modes)

### MTC SciFi stereo `loc0` sign convention — touches 4 unrelated files

The mapping from `eBoundLoc0` on a stereo SciFi plane to the physical
DD4hep (dd_x, dd_y) depends on **four conventions that must all agree** but
live in four unrelated files. Flipping any one of them silently flips the
reconstructed Y on every MTC track. No build error, no runtime warning — only
the diagnostic plots will reveal it.

**The four conventions:**

| # | Where | What sets the convention |
|--:|-------|--------------------------|
| 1 | `simulation/geometry/MTCDetector.xml` | sign of `stereo_angle="…"` on each `CartesianStripXStereo` `key_value` |
| 2 | `detector_plugin/segmentations/CartesianStripXStereo.cpp` | which perpendicular direction is "+strip-index" for a given `stereo_angle` |
| 3 | `gaudi_source/ACTSGeoSvc.cpp:273` | sign of `α` in `R_X(±α) · R_Y(π/2)` per `pi.plane` |
| 4 | `gaudi_source/MTCSciFiMeasConverter.cpp:181` | `cos α · pos.x` — assumes `pos.x` is already the strip centre at y=0 in mm, in the "natural" direction for this α |

If (1)/(2)/(3)/(4) collectively give
`loc0_plane0 = cos α · dd_x + sin α · dd_y` (the convention `docs/known_issues.md`
recorded in May 2025), the display formula
`gy = (loc0_V − loc0_U) / (2 tan α)` with `V = (tilt > 0)` returns `+dd_y`.

If any one of them is flipped, the same formula returns `−dd_y` and the bug
re-appears on the muon plots without any other symptom — X still tracks
correctly because the **sum** `loc0_V + loc0_U` is invariant to the sign of the
`sin α · dd_y` term.

**Empirical test (90 seconds, run after any change to the four files above):**

```bash
cd gaudi_jobs/muon_pipeline
python diagnose_mtc_curvature.py --max-windows 5
# look at plots_curvature/mtc_curvature_w*.pdf
# if reco y on MTC stations is the *mirror image* of truth y (same magnitude,
# opposite sign), one of the four conventions has flipped.
```

**Where to fix when the test fails:**

The single-file remedy is `event_display_eve.py:464-465` and
`gaudi_jobs/muon_pipeline/diagnose_mtc_curvature.py:115-118` (swap `U`/`V`
labels, or equivalently flip the subtraction order in `gy`). The root-cause
remedy is to identify which of the four files changed and revert that change.

Note also that `gaudi_source/ACTSProtoTracker.cpp:1474-1479` (the
`toGlobalX/toGlobalY` lambdas that write the `x`/`y` RNTuple columns)
assumes its own convention, documented in the comment block at lines 1467-1473:
`loc0 = cos α · dd_x − sin α · dd_y`. The diagnostic does **not** read those
columns (it stereo-pairs `loc0` itself), so the diagnostic sign and the
RNTuple `y` sign can disagree. Whatever swap is applied to fix the display
must also be reflected as a sign flip on the `loc0 · stereoTiltZ` term in
`toGlobalY`, otherwise downstream consumers of the `y` column will be off.

### No surface material on ACTS tracking geometry — CKF MS/EL flags are no-ops (FIXED — Path B, 2026-05-17)

**Was:** `ACTSGeoSvc.cpp` constructed each `PlaneLayer` with no
`Acts::HomogeneousSurfaceMaterial` attached. The CKF MS/EL flags at
`ACTSProtoTracker.cpp:1340-1341` were no-ops; 50 mm MTC iron per layer was
treated as vacuum, biasing q/p and producing 3–10× over-curvature.

**Fix applied (Path B from [`acts_material_migration.md`](acts_material_migration.md)):**
A `makeSlab` helper in `ACTSGeoSvc.cpp` (anonymous namespace) queries DD4hep's
already-loaded material database by name and constructs an `Acts::MaterialSlab`.
After each `PlaneLayer::create(...)` call, `detElem->surface().assignSurfaceMaterial()`
attaches a `HomogeneousSurfaceMaterial` keyed on `pi.detID` and `pi.plane`:

| Subdetector | Plane | Material assigned |
|---|---|---|
| SiTarget (detID=0) | plane=0 (X) | TungstenDens1910(3.5mm) + Silicon(0.3mm) |
| SiTarget (detID=0) | plane=1 (Y) | Silicon(0.3mm) |
| SiPad (detID=1) | plane=−1 | TungstenDens1910(3.5mm) + Silicon(0.65mm) |
| MTC SciFi (detID=2) | plane=0 (U) | Iron(53mm) + Scintillator(1.35mm) |
| MTC SciFi (detID=2) | plane=1 (V) | Scintillator(1.35mm) |
| MTC Scint (detID=2) | plane=2 | Iron(3mm) + Scintillator(15mm) |

Material properties (X₀, L₀, ρ, A, Z) come from `dd4hep::Detector::material(name)`;
DD4hep `radLength()`/`intLength()` return cm, multiplied by 10 for ACTS mm.
`MaterialSlab::combineLayers` is used to fold absorber + sensor into one effective slab.

**Verification:** Static log at GeoSvc init: `MatDump: X0=3.83mm L0=111mm t=3.8mm`
(first SiTarget surface; W X₀≈3.5mm dominates, total t=3.5+0.3=3.8mm — matches PDG).

**Result (50-event muon pipeline, post-fix):**
- Median |R_truth/R_reco| = **0.94** (was 3–10×); reco within 6% of truth
- Median reco |dx| through MTC: **1–6 mm mean, 10–20 mm max** (was 1600–3400 mm)
- Sign agreement: 22/26 = 85% (4 outliers; seeding ghosts, not material-related)

**Follow-up (Path A):** Full `Acts::convertDD4hepDetector` migration remains the
long-term goal (single source of truth, passive iron layers as genuine ACTS layers).
See [`acts_material_migration.md`](acts_material_migration.md).
