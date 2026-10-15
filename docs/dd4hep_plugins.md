# DD4hep Detector Plugins

## Overview

Three sub-detectors, each a C++ DD4hep plugin (`dd4hep_add_plugin` in CMakeLists.txt):

| Plugin target | Source | Readout | Segmentation |
|--------------|--------|---------|--------------|
| `SiTargetDetector_plugin` | `detector_plugin/SiTargetDetector.cpp` | `SiTargetHits` | MultiSegmentation(CartesianStripX/Y keyed on `plane`) |
| `SiPadDetector_plugin` | `detector_plugin/SiPadDetector.cpp` | `SiPadHits` | MultiSegmentation(36 × TiledLayerGridXY keyed on `wafer`) |
| `MTCDetector_plugin` | `detector_plugin/MTCDetector.cpp` + `CartesianStripXStereo.cpp` | `MTCDetHits` | MultiSegmentation(CartesianStripXStereo U/V + CartesianGridXY keyed on `plane`) |

---

## Compact XML Structure

Top-level: `simulation/geometry/SND_compact.xml`  
Includes sub-detectors via `<include ref="..."/>`.

Detector ordering (Z along beam axis, particle gun enters at Z ≈ −1000 mm):
1. `SiTarget.xml` — upstream silicon strip tracker (**120 layers**, X/Y strips), z = −1690…−370 mm
2. `SiPadDetector.xml` — silicon pad ECAL (middle, **22 layers**), z = −370…0 mm
3. `MTCDetector.xml` — downstream muon/tau catcher (**3 stations × 14 layers**), z = +1…+3361 mm

Layer counts are the values in `SND_compact.xml`, which are also their ceilings
(`SiTarget_NLayers_max`, `SiPad_NLayers_max`); the MOBO geometry variants vary the
count below them. **Do not take layer counts from anywhere else**: the GENIE sample
`output_genie_allflavours_0_2000.edm4hep.root` was simulated against an earlier
166-layer SiTarget and its hits address layers this geometry does not have.

### Readout declaration pattern
```xml
<readout name="MTCDetHits">
  <segmentation type="MultiSegmentation" key="plane">
    <segmentation type="CartesianStripXStereo" key_value="0"
                  strip_size_x="MTC_SciFiChannelSize"
                  offset_x="-(MTC_max_env_width/2.0)+0.001*mm"
                  stereo_angle="5.0"/>
    <segmentation type="CartesianStripXStereo" key_value="1"
                  strip_size_x="MTC_SciFiChannelSize"
                  offset_x="-(MTC_max_env_width/2.0)+0.001*mm"
                  stereo_angle="-5.0"/>
    <segmentation type="CartesianGridXY" key_value="2"
                  grid_size_x="MTC_ScintCellSize" grid_size_y="MTC_ScintCellSize"/>
  </segmentation>
  <id>system:8,station:2,layer:8,slice:4,plane:2,strip:14,x:9,y:9</id>
</readout>
```

**Critical:** Custom segmentation attributes (`stereo_angle`, etc.) must be **literal numbers** — DD4hep does NOT evaluate constant expressions for custom plugins (uses raw `stringstream` parsing).

---

## SiTarget Detector (`SiTargetDetector.cpp`)

- **120 layers**; each is 3.5 mm of tungsten followed by two sensitive silicon
  planes of 0.3 mm, layer pitch 11.0 mm
- `plane=0` → CartesianStripX (X-measuring strips), `plane=1` → CartesianStripY
- Strip pitch 75.5 µm; sensors tiled 4 × 2 per plane, so the sensitive volumes
  are the `<slice>_colN_rowM` children rather than the slice itself
- **Its own cellID layout**, not the MTC's:
  `system:8,layer:8,slice:4,plane:1,column:2,row:2,strip:14` — `plane` is ONE bit
  at offset 20, not two bits at 22. Code that decodes both detectors needs both
  descriptors (see `gaudi_source/pfa/SNDCellDecoder.h`).

---

## SiPad Detector (`SiPadDetector.cpp`)

- **22 layers**, each 10 mm of tungsten + 0.65 mm of silicon, layer pitch 16.55 mm
- Readout is a **MultiSegmentation over 36 `TiledLayerGridXY`** keyed on `wafer`
  (3 × 3 ASUs × 2 × 2 wafers), pad pitch 5.53 mm, 96 × 96 pads per plane, with the
  physical dead rims between wafers reproduced
- cellID: `system:8,layer:8,slice:4,x:9,y:9,wafer:6` — no `plane` field at all
- The sensitive volumes are the `<slice>_wafer_pads` children, not the slice
- Used as ECAL; also feeds ACTS tracking

---

## MTC Detector (`MTCDetector.cpp`)

Three stations: MTC40, MTC50, MTC60, all 600 × 600 mm envelopes, back to back
(`MTC_inter_station_gap = 0`). They differ in the width of the magnetised iron
core: 400 / 500 / 600 mm, flanked by non-magnetic steel filler. The detector
element carries `field_y="MTC_BFieldY"` (−1.9 T), applied by a custom
`CartesianField` **only inside the slice-0 iron core**, not across the envelope.

Each station: **14 layers**, pitch 80.0 mm. Each layer, in beam order:

| Slice | Material | Thickness | `plane` | Sensitive |
|-------|----------|-----------|---------|-----------|
| 0 outer iron (magnetised) | Iron | 50.0 mm | — | no |
| 1 inner iron | Iron | 3.0 mm | — | no |
| 2 SciFi U (+5° stereo) | **Scintillator** | 1.35 mm | 0 | yes |
| 3 gap | Air | 1.3 mm | — | no |
| 4 SciFi V (−5° stereo) | **Scintillator** | 1.35 mm | 1 | yes |
| 5 inner iron | Iron | 3.0 mm | — | no |
| 6 scintillator tile | Scintillator | 20.0 mm | 2 | yes |
| 7 gap | Air | 1e-6 mm | — | no |

So one station is 14 × 56 mm of iron and the three together hold **2352.75 mm**,
whose stopping power reaches only ~2.7 GeV — the MTC cannot measure a muon's
momentum by range for the spectrum it sees.

The SciFi planes are homogeneous scintillator slabs, not individual fibres, and
their material is **Scintillator** (polyvinyltoluene), not silicon as an earlier
version of this document said.

All sensitive slices share one readout (`MTCDetHits`) routed by `plane` field.

### CellID layout
`system:8, station:2, layer:8, slice:4, plane:2, strip:14, x:9, y:9`

Bit offsets:
- system: 0–7
- station: 8–9
- layer: 10–17
- slice: 18–21
- **plane: 22–23** ← routing key for MultiSegmentation, for `EventOverlay`'s
  split of the MTC into SciFi (plane 0/1) and scintillator (plane 2), and for
  every consumer downstream

---

## CartesianStripXStereo Segmentation (`detector_plugin/CartesianStripXStereo.cpp`)

Custom DD4hep segmentation for stereo SciFi strips.

**Strip assignment:**
```
strip = floor((local_x - local_y * tan(stereo_angle_deg * π/180)) / pitch)
```
- U plane (stereo_angle = +5°): strip shifts right as y increases
- V plane (stereo_angle = −5°): strip shifts left as y increases

**`position(cellID)`** returns strip centre at y=0 (x₀ = strip_centre_at_y=0, y=0, z=0). This gives the reference stereo coordinate; the actual y of a hit is stored separately in `hit.position.y` by `SND_SciFiAction`.

**Implementation class naming:** Plugin name (`CartesianStripXStereo`) must differ from impl class name (`CartesianStripXStereoImpl`) because `DECLARE_SEGMENTATION` creates a stub struct named after the plugin string.

**Unit:** `position()` returns ROOT cm (not mm). See unit convention table in the `key4hep-snd` skill.

---

## DDG4 Sensitive Action (`simulation/ddg4/SND_SciFiAction.cpp`)

Replaces the default DDG4 calorimeter action for the MTC detector.

Per Geant4 step:
1. `process()` — computes stereo-corrected strip cellID via `m_segmentation.cellID(step)`. Sets `hit.position.x = seg_p.x()` (strip centre at y=0, ROOT cm). Accumulates energy-weighted step-Y: `sumEY[cellID] += edep * y_mid`.
2. `end()` — overwrites `hit.position.y = sumEY / sumE` (energy-weighted avg Y in Geant4 mm).

Activated in ddsim steering:
```python
SIM.action.mapActions["MTC"] = "SND_SciFiAction"
```
String `"MTC"` must match the detector name in `SND_compact.xml`.

**Build requirement:** Links `DD4hep::DDG4` (transitive Geant4 deps) → must have `INSTALL_RPATH_USE_LINK_PATH TRUE`.

For full details on position encoding and unit inconsistencies, see `docs/mtc_scifi_hit_pipeline.md`.

---

## Keeping this document true

Every number above is in `simulation/geometry/SND_compact.xml` and its includes,
and `simulation/geometry/parse_geometry.py` reads them:

```python
from parse_geometry import SNDGeometry
geo = SNDGeometry()
geo.constants["SiTarget_NLayers"], geo.mtc_n_layers, geo.bitfields["MTCDetHits"]
```

Prefer that over copying a constant into code or into a slide. This file was
last checked against the compact on 2026-09-28; before that it claimed 20
SiTarget layers, 20 SiPad layers, 15 MTC layers per station and silicon SciFi
planes, none of which had been true for months.

---

## Adding a New Detector Plugin

1. Create `detector_plugin/MyDetector.cpp`
2. Add to CMakeLists.txt:
   ```cmake
   dd4hep_add_plugin(MyDetector_plugin SOURCES detector_plugin/MyDetector.cpp)
   target_link_libraries(MyDetector_plugin DD4hep::DDCore)
   ```
3. Add `<include ref="MyDetector.xml"/>` to `SND_compact.xml`
4. `/build` to rebuild
5. Update `init_key4ship.sh` if a new env variable is needed
