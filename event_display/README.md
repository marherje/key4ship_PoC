# SND Event Display

Interactive 3D event display for SND@LHC based on ROOT TEve.

Nothing about the detector is hardcoded in the Python script: the layer
z-positions are read from the DD4hep compact given with `--geometry`
(`extract_z_from_geometry()` walks `gGeoManager` and fills `layers_z_cm`
in place), and everything else — which RNTuple feeds each detector, hit box
sizes, colours — comes from `detector_config.json`. That is what makes the
same display work unchanged on the baseline and on any optimization variant.

## Prerequisites

```bash
source init_key4ship.sh          # from the repo root; needed for the DD4hep plugins
```

- A `ShipHits.root` produced by `job5_rntuple.py` (any pipeline or variant).
- The compact XML that **that same file was simulated with**. Pairing hits with
  a different geometry silently snaps hits to the wrong planes.
- An X11 display — TEve is interactive. On lxplus, connect with `ssh -Y`.

---

## Basic commands

`launch.sh` is the short form; it must be run from this directory.

```bash
cd event_display
./launch.sh <ShipHits.root> <compact.xml> [window]
```

The equivalent direct call, which is what you need for the extra options:

```bash
python event_display_eve.py \
    --hits     <ShipHits.root> \
    --geometry <compact.xml> \
    --config   detector_config.json \
    --window   0 \
    --color-by detector
```

| Option | Default | Description |
|--------|---------|-------------|
| `--hits` | *(required)* | `ShipHits.root`, the RNTuple output of `job5_rntuple.py` |
| `--geometry` | `../simulation/geometry/SND_compact.xml` | DD4hep compact used to extract the layer z-positions. **Must match the geometry the hits were simulated with.** |
| `--config` | `detector_config.json` | Detector/hit display config (resolved against the cwd) |
| `--window` | `0` | Event index to display (see "Choosing an event" below) |
| `--color-by` | `detector` | `detector` = one colour per subdetector; `mc` = one colour per MC-origin particle |

Close the TEve window to exit.

> `launch.sh` has no defaults: both paths are required, and it exits with a
> usage message if either is missing or does not exist. This is deliberate —
> a default geometry would silently pair a `ShipHits.root` with a compact it
> was not simulated with. `window` is the only optional argument; when omitted,
> `event_display_eve.py`'s own default (`0`) applies.

---

## Baseline-geometry pipelines

These were all simulated against the committed baseline compact,
`../simulation/geometry/SND_compact.xml`:

```bash
cd event_display

# two muons                    -- events 0-9, all with tracks
./launch.sh ../gaudi_jobs/2_mu_pipeline/ShipHits.root  ../simulation/geometry/SND_compact.xml 1

# mu+ 10 GeV + pi+ 5 GeV       -- events 0-99
./launch.sh ../gaudi_jobs/mu_pi_pipeline/ShipHits.root ../simulation/geometry/SND_compact.xml 1

# e + pi                       -- events 0-99, only ~20 have tracks
./launch.sh ../gaudi_jobs/e_pi_pipeline/ShipHits.root  ../simulation/geometry/SND_compact.xml 3
```

The other `gaudi_jobs/*_pipeline/` directories hold the Gaudi job scripts but
have no `ShipHits.root` committed — run their `*_pipeline.sh` first.

---

## Custom geometries

### Optimization variants (`controller.py`)

Each variant keeps its hits and its compact in two different trees, so both
paths must be given:

```bash
cd event_display
./launch.sh ../Optimization/Analysis/variants/var00007/ShipHits.root \
            ../Optimization/Simulation/geometry/variants/var00007.xml 3
```

Variants with output committed (see `../Optimization/results.csv` for the full
parameter table):

| Variant | Geometry | Events |
|---------|----------|---------|
| `var00001` (thin_W) | SiPad W 5 mm, 32 layers | 0–4 |
| `var00002` (short_sipad) | SiPad 200 mm / 12 layers, SiTarget 136 | 0–4 |
| `var00004`, `var00006` | SiPad W 8 mm | 0–2 |
| `var00007` (coarse10) | SiPad 500 mm / 10 layers, SiTarget 109 | 0–19 |
| `var00008` (fine35) | SiPad 600 mm / 35 layers, SiTarget 100 | 0–19 |

`var00003` and `var00005` are the controller's deliberate failure cases and
have no output.

### MOBO trials

A MOBO trial directory is self-contained — `geometry.xml` and `ShipHits.root`
sit side by side, so the two arguments share a prefix:

```bash
cd event_display
T=../Optimization/mobo/runs/smoke_local/trials/t00003
./launch.sh $T/ShipHits.root $T/geometry.xml 0
```

Runs with output on disk. `Optimization/mobo/.gitignore` excludes `runs/`, so
none of this is in the repo — it lives in the run directory where it was
produced:

| Run | Trials | Events | Notes |
|-----|--------|--------|-------|
| `snd_proxy` | `t00000`–`t00039` | 0–99 | first full optimization; `t00000` = baseline, `t00001`–`t00016` Sobol, `t00017`–`t00039` qLogNEHVI |
| `smoke_local` | `t00000`–`t00005` | 0–9 | `t00000` = baseline, `t00001`–`t00003` Sobol, `t00004`–`t00005` qLogNEHVI |
| `condor_kill` | `t00000`–`t00005` | 0–9 | `t00000` = baseline, rest Sobol |
| `condor_check` | `t00000` | 0–9 | baseline |

The parameters that produced a trial, plus a printed z-extent summary of every
subdetector, are in that trial's `params.yaml`.

#### Two geometries from the `snd_proxy` run

`snd_proxy` optimizes 7 free parameters against two objectives — maximize
`nhits_sipad`, minimize `cost_proxy`. Comparing an early random-phase trial with
one that ended up on the Pareto front shows what the optimizer actually did to
the detector, and both differences are visible in the display without measuring
anything:

```bash
cd event_display

# t00001 -- Sobol, first point of the random phase
T=../Optimization/mobo/runs/snd_proxy/trials/t00001
./launch.sh $T/ShipHits.root $T/geometry.xml 0

# t00038 -- qLogNEHVI, highest-nhits point of the Pareto front
T=../Optimization/mobo/runs/snd_proxy/trials/t00038
./launch.sh $T/ShipHits.root $T/geometry.xml 0
```

|  | `t00001` (Sobol) | `t00038` (Pareto front) |
|---|---|---|
| Tag | `sobol` | `qlognehvi` |
| SiPad | 16 layers, W 9.75 mm | **47 layers**, W 5.00 mm |
| SiPad z-extent | `[-398.1, 0.0]` mm | `[-500.0, 0.0]` mm |
| SiPad layer pitch | 24.88 mm | **10.64 mm** |
| SiTarget | 47 layers, W 2.94 mm | 40 layers, W 2.47 mm |
| SiTarget z-extent | `[-1700.0, -398.1]` mm | `[-1700.0, -500.0]` mm |
| SiTarget spacing | 14.69 mm | 15.00 mm |
| `nhits_sipad` | 5423 | 15132 |
| `cost_proxy` | 716.2 | 951.3 |

Two things to look for on screen:

1. **The SiPad sampling gets almost three times finer** — 16 planes become 47,
   and the pitch drops from 24.9 mm to 10.6 mm. This is the whole Pareto front:
   `sipad_fill` is the only parameter that varies along it.
2. **The SiTarget/SiPad boundary moves upstream**, from −398 mm to −500 mm.
   `SiTarget_dim_z = 1700 mm − SiPad_dim_z` in the compact, so the two
   subdetectors share a fixed 1700 mm budget: every millimetre the SiPad gains
   is one the SiTarget loses. That is why the optimizer pushed `SiPad_dim_z` to
   its upper bound — it buys SiPad planes *and* removes SiTarget layers, which
   is the dominant term in `cost_proxy`.

`--color-by mc` works on both files, and is the clearer view for seeing how much
more of the muon and pion tracks the finer SiPad samples:

```bash
python event_display_eve.py \
    --hits     ../Optimization/mobo/runs/snd_proxy/trials/t00038/ShipHits.root \
    --geometry ../Optimization/mobo/runs/snd_proxy/trials/t00038/geometry.xml \
    --window   0 \
    --color-by mc
```

Both trials ran the `mu_pi` pipeline over 100 events, so `--window` accepts
0–99. Since each trial carries its own `geometry.xml`, never mix the two: a
`ShipHits.root` from `t00001` paired with `t00038`'s compact would snap hits
onto 47 planes that did not exist when it was simulated.

---

## Choosing an event (a.k.a. window)

`--window` selects one entry of the file, counting from 0.

### Why "window" and "event" mean the same thing here

The chain can turn simulated events into reconstruction input in two ways, and
**both produce exactly the same `*Windowed` collections and frame parameters**
(`t_window_start`, `SiTargetSourceIDs`, `*ContribPDGs`), so job3/job4/job5
cannot tell which one made their input:

**1. Time-stream mode — MC merging + window splitting.**
`job1_shuffler.py` (`EventShuffler`) merges several ddsim files into a single
time stream: event *n* of a source is placed at `n * Delays[source]` ns, so
each source gets its own inter-event spacing, and the result is written as
`*Merged` collections. `job2_splitter.py` (`EventWindowSplitter`) then cuts
that stream into fixed slices of `WindowSize` ns (25 ns by default), writes one
frame per slice as `*Windowed`, and re-references every contribution time to
that frame's `t_window_start`. Here a window is a **readout time window** and
holds however many events happened to fall inside it — the pile-up picture.

**2. Event-by-event mode.** `job1_overlay.py` (`EventOverlay`) does it in one
step instead: output event *i* is built from input event *i*, with no delays
and no time splitting, and `t_window_start` is fixed at 0. Here one window is
exactly one event.

Because everything downstream only ever sees `*Windowed` and `window_id`, the
two words are used interchangeably across the framework — the display flag is
`--window`, the RNTuple field is `window_id`, and this document calls it an
event. They are the same number in mode 2, which is what every committed
pipeline runs today; in mode 1 one window can hold several events, and
`--window` is then a time slice rather than an event.

Both algorithms are compiled (see `CMakeLists.txt`) and the shuffler/splitter
configurations are still committed in several `gaudi_jobs/*_pipeline/`
directories, so going back to time-stream mode means running `job1_shuffler.py`
+ `job2_splitter.py` in place of `job1_overlay.py` — nothing downstream
changes, but a `--window` then stops being an event index.

### Multi-source overlay

Even in event-by-event mode there is superposition when the pipeline reads more
than one file. With `SourceIDs = [1, 2]` or `[1, 2, 3]` (`2_mu`, `2_mu_ang`,
`1_mu_1_e`, `2_mu_1_e`, `3_mu`), `--window 3` shows event 3 of *each* input
file superimposed — still one event per source, matched by index, not a
time-based pile-up. The single-source pipelines (`mu_pi`, `e_pi`, `muon`, and
the optimization chain) are strictly one simulated event per window.

Events with no reconstructed track still draw their hits.

To list the events in a file, and which of them have ACTS tracks:

```bash
python - <<'EOF'
import ROOT
ROOT.gSystem.Load("libROOTNTuple")
ROOT.gInterpreter.Declare(r'''
#include <set>
std::set<int> winset(const char* f, const char* nt) {
  auto r = ROOT::RNTupleReader::Open(nt, f);
  auto v = r->GetView<int>("window_id");
  std::set<int> s;
  for (auto i : r->GetEntryRange()) s.insert(v(i));
  return s;
}''')
f = "../Optimization/Analysis/variants/var00007/ShipHits.root"
print("events:", sorted(ROOT.winset(f, "SiTarget")))
print("with tracks:", sorted(ROOT.winset(f, "ACTSTracks")))
EOF
```

---

## Colouring by MC origin (`--color-by mc`)

By default every hit takes the colour of its subdetector. With
`--color-by mc`, each hit is coloured instead by the **MC particle its energy
deposit is attributed to**, read per hit from the `origin_pdg` field of the
RNTuple. `job5_rntuple.py` fills that field with the PDG of the hit's
highest-energy contribution, taken from the `*ContribPDGs` frame parameters
that `EventOverlay` writes from each contribution's linked `MCParticle`.

In the committed samples this resolves to the **primary**, so the flag
separates the activity of each primary rather than track from shower: in the
`mu_pi` events the muon's hits are blue and the pion's green, each including
whatever it induced downstream. That is what makes it useful for telling apart
two overlaid particles in the calorimeter and the MTC.

**`launch.sh` does not forward this flag** — it only takes the two paths and
the event. To colour by MC origin, call the Python script directly and pass
the same two paths that you would have given to `launch.sh`:

```bash
cd event_display

# instead of:  ./launch.sh <hits> <compact> 3
python event_display_eve.py \
    --hits     ../Optimization/Analysis/variants/var00007/ShipHits.root \
    --geometry ../Optimization/Simulation/geometry/variants/var00007.xml \
    --window   3 \
    --color-by mc
```

Hits are grouped in the TEve browser tree by particle (`<detector> <label>
(<pdg>)`), so each species can be toggled on and off independently, and the
legend of the species present is printed on the terminal:

```
[Hits] MC-origin legend (colour by primary particle):
[Hits]   |PDG|=0     other  rgb=(1.00,1.00,1.00)
[Hits]   |PDG|=13    mu     rgb=(0.20,0.45,1.00)
[Hits]   |PDG|=211   pi     rgb=(0.20,0.85,0.30)
```

`|PDG|=0` is not a particle: it is the fallback written when the winning
contribution has no linked `MCParticle` (or the `*ContribPDGs` parameter is
absent), drawn white like any unlisted species — deliberately not grey, which
is the neutron's colour. A few are normal; a file that is *all* `|PDG|=0` means
the parameters never got written, so the colouring is meaningless and
`--color-by detector` is the honest view.

Default palette — the lookup drops the sign, so a particle and its
antiparticle share a colour:

| PDG | Label | Colour |
|-----|-------|--------|
| 11 | `e` | red |
| 13 | `mu` | blue |
| 211 | `pi` | green |
| 321 | `K` | purple |
| 2212 | `p` | orange |
| 2112 | `n` | grey |
| 22 | `gamma` | yellow |
| *(other)* | `other` | white |

To override it, add a top-level `mc_colors` map to the JSON config, keyed by
the **absolute** PDG value:

```json
"mc_colors": {
  "13":  [0.0, 1.0, 1.0],
  "211": [1.0, 0.0, 1.0]
}
```

All `ShipHits.root` files committed in this repo carry `origin_pdg`, so
`--color-by mc` works on every example above.

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
    }
  ],
  "geometry": [
    {
      "name":         "SiTarget planes",
      "color":        [0.2, 0.4, 0.9],
      "transparency": 90,
      "voxel":        {"x": 20.0, "y": 20.0, "z": 0.015},
      "layers_z_cm":  [ ... ]
    }
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
- **`layers_z_cm`** — placeholder only. It is overwritten at startup from the
  compact passed with `--geometry`; the values committed in the JSON are never
  used. If a group comes out with no layers the display aborts, listing the
  offending entries: it means their `slice_index` no longer points at the
  sensitive slice of that layer in this geometry. Failing there is deliberate —
  otherwise the display would open happily and draw that detector empty, which
  looks exactly like "no hits in this event".
