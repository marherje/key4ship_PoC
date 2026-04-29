# SND Event Display

Interactive 3D event display for SND@LHC based on ROOT TEve.
Detector geometry and hit-display parameters are driven entirely by
`detector_config.json` — no positions are hardcoded in the Python script.

## Prerequisites

- key4hep environment sourced
- Simulation pipeline completed: `ShipHits.root` must exist in the working
  directory (produced by `job5_rntuple.py`)
- DD4hep geometry built: `../simulation/geometry/SND_compact.xml` must be present

---

## Step 1 — Extract sensitive-layer Z positions

Run `inspect_SND.py` from the `event_display/` directory.
It loads the DD4hep geometry and walks the ROOT TGeo tree to find every
sensitive silicon slice, printing its global Z position (in cm):

```bash
cd key4ship_PoC/event_display
python inspect_SND.py
```

Example output:

```
SiTarget_X   SiTarget_vol_0_layer_0_slice_2   x= -37.04
SiTarget_Y   SiTarget_vol_0_layer_0_slice_4   x= -36.46
...

SITARGET_X_Z = [
    -37.04, -35.94, -34.84, ...
]
SITARGET_Y_Z = [
    -36.46, -35.36, -34.26, ...
]
SIPAD_Z = [
    -13.95, -12.50, -11.05, ...
]
```

> The positions at the bottom of the output are formatted as Python lists
> ready to be pasted directly into the config file.

---

## Step 2 — Update `detector_config.json`

Open `detector_config.json` and replace the `"layers_z_cm"` arrays with the
values printed by `inspect_SND.py`:

| JSON key               | Value from inspect output |
|------------------------|---------------------------|
| `SiTarget_StripX` → `"layers_z_cm"` | `SITARGET_X_Z` list |
| `SiTarget_StripY` → `"layers_z_cm"` | `SITARGET_Y_Z` list |
| `SiPad`           → `"layers_z_cm"` | `SIPAD_Z` list      |

Also update the `"geometry"` section:
- `SiTarget planes` → union of `SITARGET_X_Z` and `SITARGET_Y_Z` (both planes per layer)
- `SiPad planes`    → `SIPAD_Z`

The `"color"` and `"voxel"` fields in the config control the appearance of
each detector; edit them freely without re-running `inspect_SND.py`.

---

## Step 3 — Run the event display

```bash
python event_display_eve.py \
    --hits  ShipHits.root \
    --config detector_config.json \
    --window 0
```

| Option | Description |
|--------|-------------|
| `--hits` | Path to `ShipHits.root` (RNTuple output of `job5_rntuple.py`) |
| `--config` | Path to the JSON config (default: `detector_config.json`) |
| `--window` | Time-window / event index to display (default: `0`) |

Close the TEve window to exit.

---

## Config file structure (`detector_config.json`)

```json
{
  "detectors": [
    {
      "name":        "SiTarget_StripX",
      "ntuple":      "SiTargetMeas",
      "filter":      {"plane": 0},
      "color":       [1.0, 0.4, 0.7],
      "voxel":       {"x": 0.003775, "y": 4.9, "z": 0.015},
      "layers_z_cm": [ ... ]
    },
    ...
  ],
  "geometry": [
    {
      "name":         "SiTarget planes",
      "color":        [0.2, 0.4, 0.9],
      "transparency": 90,
      "voxel":        {"x": 20.0, "y": 20.0, "z": 0.015},
      "layers_z_cm":  [ ... ]
    },
    ...
  ]
}
```

- **`detectors`** — each entry reads one RNTuple from `ShipHits.root` and draws
  one box per hit; `filter` applies an equality cut on any integer field
  (e.g. `{"plane": 0}` keeps only StripX hits).
- **`geometry`** — transparent outline planes drawn in the global scene as a
  detector reference frame.
- **`voxel`** — half-sizes of the box drawn per hit, in cm.
- **`color`** — RGB triplet `[r, g, b]` in the range `[0, 1]`.
