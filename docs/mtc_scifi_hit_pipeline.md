# MTC SciFi Hit Pipeline: Detector → Gaudi

## Geometry — `simulation/geometry/MTCDetector.xml`

Three stations (MTC40/50/60, `station_id` 0/1/2) each with 15 identical layers:

| slice | material | `plane` field | sensitive |
|-------|----------|---------------|-----------|
| Iron absorber (×2) | — | — | no |
| SciFi U | Silicon | 0 | yes |
| Air gap | — | — | no |
| SciFi V | Silicon | 1 | yes |
| Scintillator | Scintillator | 2 | yes |

All sensitive slices share one readout `MTCDetHits`:

```xml
<readout name="MTCDetHits">
  <segmentation type="MultiSegmentation" key="plane">
    <segmentation type="CartesianStripXStereo" key_value="0"
                  strip_size_x="MTC_SciFiChannelSize"
                  offset_x="-(MTC_max_env_width/2.0)+0.001*mm"
                  stereo_angle="5.0"/>    <!-- U: strip from x - y*tan(+5°) -->
    <segmentation type="CartesianStripXStereo" key_value="1"
                  strip_size_x="MTC_SciFiChannelSize"
                  offset_x="-(MTC_max_env_width/2.0)+0.001*mm"
                  stereo_angle="-5.0"/>   <!-- V: strip from x - y*tan(-5°) -->
    <segmentation type="CartesianGridXY" key_value="2" .../>   <!-- scint pads -->
  </segmentation>
  <id>system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9</id>
</readout>
```

The `plane` field at bit offset 22 is the routing key used by all downstream code.

### CartesianStripXStereo segmentation — `detector_plugin/CartesianStripXStereo.cpp`

Custom DD4hep segmentation plugin replacing `CartesianStripX` for the MTC SciFi planes.

**Strip assignment:** `strip = floor((local_x - local_y * tan(stereo_angle_deg * π/180)) / pitch)`

- U plane (stereo_angle=+5°): `x_stereo = x - y·tan(5°)` — strip number shifts right as y increases
- V plane (stereo_angle=−5°): `x_stereo = x + y·tan(5°)` — strip number shifts left as y increases

**`position(cellID)`** returns the strip centre at y=0 (reference plane), i.e. `x_0 = strip_centre_at_y=0`. This is the correct stereo coordinate for tracking.

**Registration:** The `DECLARE_SEGMENTATION` macro creates a stub struct with the plugin name, so the real implementation class must have a different name (`CartesianStripXStereoImpl` vs plugin name `CartesianStripXStereo`). The `stereo_angle` XML attribute must be a **literal number** (e.g. `"5.0"`) — DD4hep does not evaluate constant expressions for custom segmentation parameters (it uses raw string parsing via `stringstream`).

**Verified:** For a muon at y≈−42.8mm, U and V strips in the same layer differ by ~7 strips. Expected: `2 × 4.28cm × tan(5°) / 0.1cm ≈ 7.5` — consistent within 1-strip quantization. ✓

---

## DDG4 sensitive action — `simulation/ddg4/SND_SciFiAction.cpp`

`SND_SciFiAction` replaces the default DDG4 calorimeter action for the MTC detector.

**What it does, per G4 event:**
1. `process()` — on each step, gets the cellID via `m_segmentation.cellID(step)`. With `CartesianStripXStereo`, the cellID now encodes the stereo-corrected strip number. Sets `hit->position.x = seg_p.x()` (strip centre at y=0, in ROOT cm units). Accumulates energy-weighted step-Y: `sumEY[cellID] += edep * y_mid`
2. `end()` — overwrites `hit.position.y = sumEY[cellID] / sumE[cellID]` (energy-weighted average Y of all steps in that strip)

**Why position.y is repurposed:** `CartesianStripXStereo::position()` returns y=0 (no Y information in a strip). The EDM4HEP mapper does not transfer `MonteCarloContrib::y` to `CaloHitContribution.stepPosition` in key4hep-2026-02-01. Repurposing `position.y` is the only field reliably transferred that is otherwise unused for this detector type. `SciFiDigitizer` reads it as the transverse hit coordinate for light-propagation distance computation.

**Unit note (important for downstream):** `hit.position.x` is set from `m_segmentation.position(cid).x()`, which returns the raw ROOT-cm value from `CartesianStripXStereoImpl::position()` without the cm→mm conversion that DDG4 would normally apply. In contrast, `hit.position.y` and `hit.position.z` are in Geant4/DDG4 mm. Therefore:
- `position.x` is in **ROOT cm** (e.g. 0.45 for a strip at 4.5mm transverse)
- `position.y` is in **Geant4 mm** (e.g. −42.8 for a hit at y=−42.8mm)
- `position.z` is in **Geant4 mm**

`SciFiDigitizer` is unaffected because `x_hit` cancels in the `dx = x_hit − x_sipm` term. Any future tracking code reading `position.x` must multiply by 10 to convert to mm.

**Registration:** activated via `SIM.action.mapActions["MTC"] = "SND_SciFiAction"` in the ddsim steering file.

**Build note:** `libSND_SciFiAction_plugin.so` links `DD4hep::DDG4` (which needs Geant4 transitive libs). It requires `INSTALL_RPATH_USE_LINK_PATH TRUE` in `CMakeLists.txt` so the installed `.so` embeds RPATH to `geant4/11.4.0-wr5h4o/lib64`.

---

## Collection name chain

```
ddsim  (Geant4 + DD4hep + SND_SciFiAction + CartesianStripXStereo)
  MTCDetHits              SimCalorimeterHit   output_*.edm4hep.root
    cellID   = stereo strip number (x - y*tan(±5°) / pitch), correct for tracking
    position.x = strip centre at y=0  [ROOT cm — NOT mm, see unit note above]
    position.y = energy-weighted avg step-Y  [Geant4 mm] ← written by SND_SciFiAction
    position.z = step midpoint Z  [Geant4 mm]

job1  EventShuffler
  routes hits by plane field (bits 22–23 of cellID):
    plane 0,1 → MTCSciFiHitsMerged   SimCalorimeterHit   shuffled.edm4hep.root
    plane 2   → MTCScintHitsMerged   SimCalorimeterHit

job2  EventWindowSplitter  (25 ns windows)
  MTCSciFiHitsWindowed    SimCalorimeterHit   timewindows.edm4hep.root
  MTCScintHitsWindowed    SimCalorimeterHit

job3  SciFiDigitizer        (SciFi only — see physics below)
  MTCSciFiHitsDigi        SimCalorimeterHit   digitized.edm4hep.root
      .energy = QDC counts
      .position copied from MTCSciFiHitsWindowed (same unit caveat applies)
job3  GeV2MIPConversion + BasicDigitizer  (scintillator)
  MTCScintHitsMIP → MTCScintHitsDigi

job4  ACTSProtoTracker  (SiTarget + SiPad only; MTC not yet used for tracking)
  iosvc "keep *" forwards MTCSciFiHitsWindowed + MTCScintHitsWindowed
  → tracks.edm4hep.root

job5  EDM4HEP2RNTuple
  MTCSciFi (RNTuple, reads MTCSciFiHitsWindowed)   ShipHits.root
  MTCScint (RNTuple, reads MTCScintHitsWindowed)
```

**Note:** job5 writes the **windowed** (pre-digitization) collections to the RNTuple. The digitized outputs (`MTCSciFiHitsDigi`, `MTCScintHitsDigi`) are currently not persisted to `ShipHits.root`.

---

## SciFiDigitizer physics — `gaudi_source/SciFiDigitizer.cpp`

Processes `MTCSciFiHitsWindowed` (plane 0 and 1 only; plane 2 is skipped). Each hit represents one strip cell.

### Step 1 — propagation distance

```
hit.position.x  → strip centre x [ROOT cm — see unit note; x_sipm formula absorbs this]
hit.position.y  → energy-weighted avg step-Y [mm] ← from SND_SciFiAction
station_id      → from bits 8–9 of cellID
h_half          → EnvHeightHalf[station_id]  (MTC40=200, MTC50=250, MTC60=300 mm)
sign            → +1 for U plane (plane=0), -1 for V plane (plane=1)
```

The fiber is tilted by `FiberAngleDeg` (5°) from the Y axis. The SiPM is at the `+y` end (`SiPMSide = +1`):

```
x_sipm = x_hit + sign * y_sipm * tan(angle)   # SiPM x for THIS fiber [same unit as x_hit]
y_sipm = SiPMSide * h_half                    # y of SiPM end [mm]
dx = x_hit - x_sipm = -sign * h_half * tan(angle)   # x_hit cancels — independent of strip
dy = y_hit - y_sipm
d  = sqrt(dx^2 + dy^2)                        # 2D distance hit → SiPM
d_cm = d / 10
fiber_len_cm = (2 * h_half / cos(angle)) / 10
```

`x_hit` cancels in `dx`, so the propagation distance depends only on `y_hit` and the fiber geometry — the unit inconsistency in `position.x` has no effect here.

### Step 2 — photon yield with attenuation

```
n_photons = E_GeV * PhotonsPerGeV
n_direct   = n_photons * exp(-(d_cm - AttenuationOffset) / AttenuationLength)
n_reflect  = n_photons * MirrorReflectivity
           * exp(-((fiber_len_cm - d_cm) - AttenuationOffset) / AttenuationLength)
n_total = n_direct + n_reflect
```

Parameters (from `job3_digitize.py`): `AttenuationLength=300 cm`, `AttenuationOffset=20 cm`, `PhotonsPerGeV=1.6e5`, `MirrorReflectivity=0.9`.

### Step 3 — SiPM saturation + Poisson

```python
n_detected = Poisson(n_total)
n_fired = SiPMNpixels * (1 - exp(-n_detected / SiPMNpixels))  # saturation
```

Threshold cut: if `n_fired < PhotonThreshold (3.5)`, hit is dropped.

### Step 4 — QDC conversion

```python
A = Gauss(QDC_A, QDC_sigmaA)   # 0.172 ± 0.006
B = Gauss(QDC_B, QDC_sigmaB)   # -1.31 ± 0.33
signal = A * n_fired + B        # QDC counts
```

Output: `hit.setEnergy(signal)` where energy now means QDC counts, not GeV.

---

## RNTuple columns (`MTCSciFi` in `ShipHits.root`)

From `job5_rntuple.py` BitField `"system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9"`:

| column | content |
|--------|---------|
| `E` | energy deposit (GeV, pre-digitization windowed hit) |
| `x` | strip centre at y=0 [ROOT cm, ~0.1–1 range]; multiply ×10 for mm |
| `y` | energy-weighted avg step-Y [Geant4 mm] |
| `z` | step midpoint Z [Geant4 mm] |
| `bf_system` | detector system ID (3 for MTC) |
| `bf_station` | station index: 0=MTC40, 1=MTC50, 2=MTC60 |
| `bf_layer` | layer index within station (0–14) |
| `bf_slice` | slice index within layer |
| `bf_plane` | 0=SciFi U, 1=SciFi V |
| `bf_strip` | stereo strip channel: encodes `x - y*tan(±5°)` quantised to 1mm pitch |
| `bf_x`, `bf_y` | pad coordinates (unused for SciFi; used for scintillator) |
| `source_id` | muon source index (from EventShuffler) |
| `t` | hit time (ns) |
| `t_window_start`, `window_id` | time window bookkeeping |


## Problem

The values after the digitization are ~10^3 lower than the values in the fairroot simulation. Why if the gev to number of photon conversion is correct?
