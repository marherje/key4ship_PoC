# Gaudi Reconstruction Pipeline

## Overview

Five sequential `k4run` jobs transform ddsim output to a ROOT RNTuple.

```
output_*.edm4hep.root
  → job1: EventShuffler      → shuffled.edm4hep.root
  → job2: EventWindowSplitter → timewindows.edm4hep.root
  → job3: SciFiDigitizer + GeV2MIPConversion + BasicDigitizer → digitized.edm4hep.root
  → job4: SiTargetMeasConverter + SiPadMeasConverter + [MTCSciFiMeasConverter] + ACTSProtoTracker → tracks.edm4hep.root
  → job5: EDM4HEP2RNTuple     → ShipHits.root
```

**Reference pipeline** (MTC + B-field fully enabled):
```bash
cd gaudi_jobs/muon_pipeline
source muon_pipeline.sh
```

Other pipelines (`2_mu_pipeline`, `3_mu_pipeline`, `2_mu_ang_pipeline`, etc.) do not wire MTC or the iron B-field and use `ConstantBField(0,0,0)`. See `docs/acts_integration.md` for a full feature comparison.

---

## Job 1 — EventShuffler

**File:** `gaudi_source/EventShuffler.cpp`  
**Config:** `gaudi_jobs/*/job1_shuffler.py`

Merges N simulation files into one super-event. Assigns source IDs and time offsets per file to simulate pile-up.

### Architecture note
All work happens in `finalize()`. The algorithm reads files directly via `podio::ROOTReader`, bypassing Gaudi's IOSvc entirely. `execute()` is a no-op.

```python
ApplicationMgr(TopAlg=[shuffler], ExtSvc=[], EvtSel="NONE", EvtMax=1)
# No IOSvc — would conflict with direct podio I/O
```

### Key properties
| Property | Description |
|----------|-------------|
| `InputFiles` | List of edm4hep ROOT files (one per source) |
| `SourceIDs` | Integer ID for each file (same order) |
| `Delays` | Inter-event time delay in ns for each source |
| `CollectionsSiTarget/SiPad/MTC` | Collection name for each input file |
| `OutputFile` | Output file path |
| `MaxEventsPerSource` | 0 = no limit |

### MTC hit routing
The `plane` bitfield (bits 22–23 of cellID, layout `system:8,station:2,layer:8,slice:4,plane:2,...`) routes:
- plane 0, 1 → `MTCSciFiHitsMerged`
- plane 2 → `MTCScintHitsMerged`

### Source ID encoding
`edm4hep::CaloHitContribution` has no source field. Source ID is stored in the **PDG field** (`contrib.setPDG(source_id)`). PDG is unused elsewhere in this pipeline.

---

## Job 2 — EventWindowSplitter

**File:** `gaudi_source/EventWindowSplitter.cpp`  
**Config:** `gaudi_jobs/*/job2_splitter.py`

Splits the merged super-event into 25 ns time windows. Each window becomes one EDM4HEP frame in the output file. Uses Gaudi IOSvc (normal I/O, not bypass mode).

### Key properties
| Property | Description |
|----------|-------------|
| `WindowSize` | Time window width in ns (default 25) |
| Input/output collections | Configured via IOSvc `keep` rules |

---

## Job 3 — Digitization

**Files:** `gaudi_source/SciFiDigitizer.cpp`, `GeV2MIPConversion.cpp`, `BasicDigitizer.cpp`  
**Config:** `gaudi_jobs/*/job3_digitize.py`

Two parallel paths:

**SciFi (MTC plane 0 and 1):** `SciFiDigitizer`
- Converts GeV deposits to QDC counts using fiber attenuation model
- Computes propagation distance from hit Y (energy-weighted avg) to SiPM
- Output energy field = QDC counts (not GeV)
- See `docs/mtc_scifi_hit_pipeline.md` for physics details

**Scintillator (MTC plane 2) + SiTarget + SiPad:** `GeV2MIPConversion` → `BasicDigitizer`
- `GeV2MIPConversion`: multiplies energy by MIP/GeV factor
- `BasicDigitizer`: applies MIP threshold cut, drops hits below threshold

### SciFiDigitizer key properties
| Property | Default | Description |
|----------|---------|-------------|
| `AttenuationLength` | 300 cm | Fiber attenuation length |
| `AttenuationOffset` | 20 cm | Near-end dead zone |
| `PhotonsPerGeV` | 1.6e5 | Photon yield |
| `MirrorReflectivity` | 0.9 | Far-end mirror |
| `SiPMNpixels` | 1600 | SiPM pixel count (for saturation) |
| `PhotonThreshold` | 3.5 | Min fired pixels to keep hit |
| `QDC_A/B/sigmaA/sigmaB` | — | QDC linearization parameters |

---

## Job 4 — Tracking

**Files:** `gaudi_source/SiTargetMeasConverter.cpp`, `SiPadMeasConverter.cpp`, `MTCSciFiMeasConverter.cpp`, `ACTSProtoTracker.cpp`  
**Config:** `gaudi_jobs/*/job4_tracking.py`

**Step 1 — Measurement converters:** Convert EDM4HEP hit collections to ACTS `SourceLink`/measurement objects keyed to ACTS surfaces from `ACTSGeoSvc`.

**Step 2 — ACTSProtoTracker:**
1. Hough transform seeding (2D, using SiTarget X×Y crossings, SiPad 2D hits, and MTC U×V stereo crossings)
2. Kalman filter track fit using ACTS `KalmanFitter` with `IronSlabBField` (MTC iron) or `ConstantBField(0,0,0)`
3. Outputs `ACTSTracks` collection to EDM4HEP

### Pipeline feature matrix

| Feature | `muon_pipeline` | other pipelines |
|---|---|---|
| `MTCSciFiMeasConverter` active | **Yes** | No |
| MTC hits in Hough seeding and KF | **Yes** | No |
| Iron slab B-field (`IronFieldRanges`) | **Yes** | No — `ConstantBField(0,0,0)` |
| Geometry params from `parse_geometry` | **Yes** | No — hardcoded |

`muon_pipeline` is the reference; use `gaudi_jobs/muon_pipeline/job4_tracking.py` as the template for new configurations.

### ACTSProtoTracker key properties
| Property | Default | Description |
|----------|---------|-------------|
| `AutoSeed` | false | Use Hough seeding automatically |
| `MaxSeeds` | 5 | Maximum number of seeds to fit |
| `HoughBinSize` | 5.0 | Hough accumulator bin size [mm] |
| `HoughHalfSize` | 200.0 | Hough accumulator range [mm] |
| `HoughMinVotes` | 3 | Minimum votes to form a seed |
| `InputMTC` | `"MTCSciFiMeasurements"` | MTC measurement collection (must be produced by `MTCSciFiMeasConverter`) |
| `MTCStereoAngle` | 5.0 | MTC SciFi stereo angle [degrees] |
| `IronFieldRanges` | `[]` | Per-slab field map: `[xlo,xhi, ylo,yhi, zlo,zhi, by] × N` in ACTS coords [mm, T] |
| `IsolationWindow` | 0.0 | 2D distance [mm] for crossing isolation filter (0 = disabled) |
| `IsolationMaxNeighbors` | 2 | Max neighbors within window to be considered track-like |
| `HoughMaxMultiplicity` | 1e9 | Max crossings/station ratio to accept a Hough peak (1e9 = disabled) |
| `MaxChi2PerMeas` | 500.0 | Track acceptance threshold on chi²/nMeas |

### ACTSGeoSvc in job config
```python
from Configurables import ACTSGeoSvc
geo = ACTSGeoSvc("ACTSGeoSvc")
geo.CompactFile = "simulation/geometry/SND_compact.xml"
ApplicationMgr(..., ExtSvc=[geo])  # Service in ExtSvc, not TopAlg
```

---

## Job 5 — EDM4HEP2RNTuple

**File:** `gaudi_source/EDM4HEP2RNTuple.cpp`  
**Config:** `gaudi_jobs/*/job5_rntuple.py`

Converts EDM4HEP collections to a ROOT RNTuple (`ShipHits.root`) for analysis.

### Written collections
- `SiTarget` — from `SiTargetHitsWindowed`
- `SiPad` — from `SiPadHitsWindowed`
- `MTCSciFi` — from `MTCSciFiHitsWindowed` (**pre-digitization windowed hits**, not `MTCSciFiHitsDigi`)
- `MTCScint` — from `MTCScintHitsWindowed`
- `Tracks` — from `ACTSTracks`

**Note:** digitized MTC hits (`*Digi`) are NOT written to the RNTuple. See `docs/mtc_scifi_hit_pipeline.md` for full RNTuple column definitions.

---

## Adding a New Gaudi Algorithm

1. Create `gaudi_source/MyAlgorithm.cpp` with `DECLARE_COMPONENT(MyAlgorithm)` at end
2. Add `.cpp` to `gaudi_add_module(SND_reco SOURCES ...)` in `CMakeLists.txt`
3. `/build` to rebuild
4. Import in Python: `from Configurables import MyAlgorithm`
