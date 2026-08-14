# GENIE-driven simulation

Feeds real GENIE neutrino-interaction events into ddsim, instead of the
particle guns used by `launch_PG.sh`/`launch_multiplePG.sh`. This is exactly
the extension point `multiplePG_base.py`'s docstring anticipated: a HepMC
file replaces the gun plugins, nothing downstream changes.

## Input

GENIE production for SND (`gevgen_fnal`, all neutrino flavours) lives on EOS:

```
/eos/experiment/ship/user/aiuliano/GENIE_sims/2026_07_30_key4ship_allflavours/
```

The relevant file is `gntp.0.gst.root`, GENIE's flat "gst" summary TTree (one
row per interaction). It was generated against `SND_compact_key4ship.gdml`,
a GDML export of the same geometry as `simulation/geometry/SND_compact.xml`,
so its vertex position (`vtxx`/`vtxy`/`vtxz`, cm) is already in our DD4hep
reference frame — no FairShip↔key4ship offset to reapply here (that
conversion only concerned the neutrino flux geometry `gevgen_fnal` itself
was run against; see the EOS directory's own `README.md`).

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
