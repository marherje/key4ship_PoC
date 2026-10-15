# GENIE-driven simulation

Feeds real GENIE neutrino-interaction events into ddsim, instead of the
particle guns used by `launch_PG.sh`/`launch_multiplePG.sh`. This is exactly
the extension point `multiplePG_base.py`'s docstring anticipated: a HepMC
file replaces the gun plugins, nothing downstream changes.

There are two entry points:

* **`launch_genie_gen.sh`** — generates the GENIE events itself (`gevgen_fnal`
  against our own GDML) and simulates them, in one HTCondor job per sample.
  This is the one to use for any new geometry. See "Generating with GENIE".
* **`launch_genie.sh`** — consumes an existing gst sample (the 2026_07_30
  EOS production). Kept for reproducing what was done on that sample; note
  its flux-frame caveat below.

## Generating with GENIE (`launch_genie_gen.sh`)

GENIE (3.4.2, tune `G21_11a_00_000`) comes from the SHiP cvmfs stack
(`/cvmfs/ship.cern.ch/25.09/setUp.sh`), which cannot coexist with key4hep in
one shell: the job script runs `gevgen_fnal` + `gntpc` in a subshell with
that stack, then the HepMC conversion and `ddsim` in a subshell with ours.
Flux and cross-section splines are Antonio Iuliano's
(`/eos/experiment/ship/user/aiuliano/inputs_gevgen_fnal/`). One difference
from his production: that used FairShip's `shipgen/genie_config`, which
swaps GENIE's decayer for `Pythia8Decayer2023`; the cvmfs GENIE is built
against Pythia6 and does not have it (`gevgen_fnal` aborts at start-up), so
the launcher leaves `GXMLPATH` unset and GENIE decays with its default
`PythiaDecayer`. `GENIE_CONFIG_DIR=...` restores the override for a build
that has Pythia8.

```bash
./launch_genie_gen.sh               # one job per line of `samples`
RUN_LOCAL=1 ./launch_genie_gen.sh   # tiny tests only: the spline file is 470 MB
```

Outputs per sample `<label>`: `data/genie/<label>/gntp.<run>.{ghep,gst}.root`
(+ GENIE's `genie-mcjob-<run>.status`, which has the exposure in POT),
`data/genie_<label>.hepmc`, `data/output_genie_<label>.edm4hep.root`, logs in
`log/genie_<label>.{gevgen,gntpc,ddsim}.log`. The job skips generation if
the gst file already exists, so a failed ddsim can be re-run without
regenerating.

### Reference frame: why GENIE sees a shifted GDML

The gsimple flux rays are in the FairShip frame (metres) and `gevgen_fnal`
cannot translate them — `-z` only slides each ray along its own direction,
which is why the 2026_07_30 production has its neutrino cone opening from
z ≈ 0 of our frame (the MTC40 front) instead of from the SHiP target
25.76 m upstream. So the geometry goes to the flux instead: GENIE navigates
`simulation/geometry/gdml/SND_compact_fairship_frame.gdml`, the key4hep
export translated by +2576 cm in z (SiPad front at 2539 cm, as in FairShip)
with the world box grown to ±40 m, and `genie_gst_to_hepmc.py --z-shift-m
-25.76` brings the vertices back before ddsim. After ANY geometry change
regenerate both files:

```bash
source ../../init_key4ship.sh
(cd ../geometry && geoConverter -compact2gdml -input SND_compact.xml -output gdml/SND_compact.gdml)
(cd ../geometry/gdml && ./shift_gdml_frame.py SND_compact.gdml SND_compact_fairship_frame.gdml)
```

`shift_gdml_frame.py` only touches the three world-level placements and the
world box; the detector internals — the material budget GENIE samples
vertices over — are exactly what ddsim simulates.

## Input of `launch_genie.sh` (pre-existing sample)

GENIE production for SND (`gevgen_fnal`, all neutrino flavours) lives on EOS:

```
/eos/experiment/ship/user/aiuliano/GENIE_sims/2026_07_30_key4ship_allflavours/
```

The relevant file is `gntp.0.gst.root`, GENIE's flat "gst" summary TTree (one
row per interaction). It was generated against `SND_compact_key4ship.gdml`,
a GDML export of the same geometry as `simulation/geometry/SND_compact.xml`,
so its vertex position (`vtxx`/`vtxy`/`vtxz`, **metres** — GENIE's own
natural length unit, not ROOT/TGeo's cm, even though it navigated a
GDML/TGeo geometry) is already in our DD4hep reference frame
(`--z-shift-m 0`, the default). Its flux, however, was not shifted into
that frame (see "Reference frame" above): XY distributions and z-dependent
rates of this sample are biased toward the centre and toward MTC40.

## How it works

1. **`genie_gst_to_hepmc.py`** reads a chunk of the gst tree and writes a
   HepMC (v3 ascii) file: one GenEvent per interaction, with the GENIE
   interaction vertex, the incoming neutrino as an untracked beam particle,
   and the outgoing lepton (`fspl`) + post-FSI hadronic system (`nf`
   particles) as status-1 final-state particles — i.e. exactly what Geant4
   needs to propagate.

   ```bash
   ./genie_gst_to_hepmc.py \
       --input /eos/experiment/ship/user/aiuliano/GENIE_sims/2026_07_30_key4ship_allflavours/gntp.0.gst.root \
       --start 0 --nevents 2000 \
       --output data/genie_allflavours_0_2000.hepmc
   # --cc-only / --flavour <PDG> to select a subset
   ```

2. **`genie_base.py`** is `multiplePG_base.py`'s counterpart: `configure()`
   points `SIM.inputFiles` at the HepMC file (ddsim auto-selects
   `Geant4InputAction`/`HEPMC3FileReader` from the `.hepmc` extension) and
   sets `numberOfEvents = -1` to consume the whole file.

3. **`launch_genie.sh`** ties both steps together per configured chunk of
   the gst tree (`start:nevents[:label]`) and submits via
   `generic_condor_multiplePG.sh` (generator-agnostic — it just runs
   `ddsim --steeringFile <script>`), same as `launch_multiplePG.sh`.

   ```bash
   ./launch_genie.sh                 # submit to HTCondor
   RUN_LOCAL=1 ./launch_genie.sh     # run locally
   GENIE_INPUT=/other/gntp.0.gst.root ./launch_genie.sh
   ```

## Output structure

Same as `launch_multiplePG.sh`: `steer/`, `data/genie_<label>.hepmc` +
`data/output_genie_<label>.edm4hep.root`, `log/genie_<label>.ddsim.log`.

## Reconstruction (`gaudi_jobs/genie_pipeline/`)

`job3_digitize.py`/`job4_tracking.py`/`job5_rntuple.py` are unmodified
copies of `e_pi_pipeline`'s — digitization and tracking only consume the
`*Windowed` hit collections and don't care what generated the primaries.
Only `job1_overlay.py` differs, pointing at the GENIE ddsim output as its
single source:

```bash
cd ../../gaudi_jobs/genie_pipeline
INPUT_FILE=../../simulation/run_script/data/output_genie_<label>.edm4hep.root \
    ./genie_pipeline.sh
```

Note: GENIE final states are typically dense hadronic showers (several to a
few dozen final-state particles per interaction), not the clean 1-2 track
topologies the PG samples were designed around — `ACTSProtoTracker`'s CKF
may find few or no accepted tracks per event depending on multiplicity. That
is expected and not a pipeline failure; `ShipHits.root` still records the
full hit-level information for every event.
