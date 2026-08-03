# mobo — asynchronous multi-objective Bayesian optimization of the SND geometry

Proposes detector geometries, evaluates them on HTCondor, and learns from the
results. The optimizer is [qLogNEHVI](https://botorch.org/) over Gaussian-process
surrogates; the evaluation is the existing chain of this repository (DD4hep
variant → ddsim → `job1_overlay` → `job4_tracking` → `job5_rntuple` → metrics),
which is exactly what `Optimization/controller.py` does for one variant at a
time.

Everything lives under `Optimization/` — the package, its venv, the run
directories and the generated Condor files. Nothing outside is written to; the
only things read from outside are `init_key4ship.sh` and the shared geometry in
`simulation/geometry/`.

```
                 ask                       submit                 payload.py
  MOBOptimizer ────────► AsyncLoop ────────────────► HTCondor ───────────────► worker
       ▲                    │                                                    │
       │   tell             │  TrialStore (sqlite: ids, status, metrics)         │
       └────────────────────┴────────────────◄───────────────────────────────────┘
                                              DONE / FAILED + metrics.json
```

---

## Setup

```bash
cd Optimization/mobo
bash setup_env.sh                       # ~2 min, one off
```

It creates `.venv` on top of the key4hep python with `--system-site-packages`, so
ROOT, numpy, scipy, pandas, matplotlib **and torch** come from the stack
(key4hep 2026-02-01 already ships torch 2.9); only botorch, gpytorch, hydra and
plotly are downloaded, about 170 MB. If the AFS quota cannot take that, the
script fails with a clear message and you can put the venv elsewhere:

```bash
MOBO_VENV_DIR=/eos/user/${USER:0:1}/$USER/mobo-venv bash setup_env.sh
```

The venv is only used by the driver process. Worker nodes run the payload with
the bare key4hep python, so the venv never has to be visible from them.

Every session afterwards:

```bash
source init_key4ship.sh                 # from the repo root
source Optimization/mobo/.venv/bin/activate
```

> **`init_key4ship.sh`, not just key4hep's `setup.sh`.** The geometry uses this
> project's own DD4hep plugins (the SiPad/SiTarget/MTC builders and the
> `CartesianStripXStereo` segmentation) from `install/`, which only reach dd4hep
> through the `LD_LIBRARY_PATH` that `init_key4ship.sh` prepends. Without it
> ddsim dies with `FAILED to create segmentation ... [Missing factory]`.
> `mobo-run` checks this before submitting anything.

### One gotcha worth knowing

key4hep exports its own site-packages on `PYTHONPATH`, and `PYTHONPATH` is
searched **before** the venv. Any package the stack also ships therefore wins
over the venv's copy, whatever pip says it installed. It only bites for tools
that need a newer version than the stack has — `mypy` is the current example, so
`check.sh` puts the venv first for it.

```bash
bash check.sh            # ruff + mypy + the fast tests
bash check.sh --slow     # ... and the slow ones
```

---

## Running

```bash
# smoke test: 4 trials of 10 events, locally, ~5 minutes
mobo-run experiment=snd_proxy executor=local \
         evaluator.nevents=10 optimizer.n_init=2 \
         experiment.loop.max_trials=4 experiment.loop.max_in_flight=2

# the real thing: 20 trials of 100 events, 4 at a time, on the batch system
mobo-run experiment=snd_proxy executor=htcondor

mobo-status runs/snd_proxy          # trials, Pareto front, hypervolume
mobo-status runs/snd_proxy --all    # every trial, not just the front
mobo-report runs/snd_proxy          # -> runs/snd_proxy/report.html
mobo-resume runs/snd_proxy          # continue after a kill, a crash, a logout
mobo-resume runs/snd_proxy --max-trials 50
```

`mobo-run` is a [Hydra](https://hydra.cc) application: every key under `conf/`
is overridable on the command line (`evaluator.nevents=50`,
`optimizer.batch_size=2`, `executor.flavour=tomorrow`). The *resolved* config is
written to `runs/<name>/config.yaml` before anything starts, and that copy — not
`conf/` — is what `mobo-resume` reads back, so editing a default later cannot
silently change a run in progress.

### What a run directory holds

```
runs/snd_proxy/
├── config.yaml            the resolved config this run is pinned to
├── trials.db              sqlite: the single source of truth
├── report.html            regenerated after every completed trial
└── trials/t00000/         one directory per trial, self-contained
    ├── geometry.xml       the compact that was simulated
    ├── params.yaml        the parameters and constants that produced it
    ├── trial.json         the payload's input
    ├── job.sh, job.sub    what was submitted (Condor only)
    ├── steering.py, ddsim.log, job1.log, job4.log, job5.log
    ├── *.root             ddsim, tracking and RNTuple output
    ├── metrics.json       everything that was measured
    └── DONE | FAILED      the completion sentinel
```

A trial can be re-run, inspected or debugged from its directory alone:

```bash
cd runs/snd_proxy/trials/t00007
bash job.sh                                   # exactly what the worker ran
python3 ../../../../src/mobo/ship/payload.py --workdir $PWD
```

---

## The experiment

`conf/experiment/snd_proxy.yaml`. Seven free parameters, two objectives.

| parameter | range | baseline | note |
|---|---|---|---|
| `SiPad_WThickness` | 5 – 15 mm | 10 | |
| `SiPad_dim_z` | 250 – 500 mm | 370 | also sets `SiTarget_dim_z = 1700 − this` |
| `sipad_fill` | 0.3 – 1.0 | 0.978 | → `SiPad_NLayers` |
| `SiTarget_WThickness` | 2 – 5 mm | 3.5 | |
| `SiTarget_spacing` | 8 – 15 mm | 11 | |
| `sitarget_fill` | 0.5 – 1.0 | 1.0 | → `SiTarget_NLayers` |
| `xy_gap_frac` | 0.1 – 1.0 | 1.0 | → `SiTarget_XY_plane_gap` |

`SiPad_frame_gap` (0.1 mm), `SiTarget_module_offset` (1 mm) and
`SiPad_layer_gap` (`auto`) are held at the baseline under `fixed:`; freeing one
is a two-line config change, not a code change.

**Fill fractions instead of layer counts.** `SiPad_NLayers` is bounded by
`floor(SiPad_dim_z / layer_thickness)`, which itself depends on
`SiPad_WThickness` and `SiPad_dim_z`; `SiTarget_XY_plane_gap` is bounded by a
budget built from the spacing, the absorber and the module offset. A box in
those coordinates is mostly infeasible, and an optimizer that spends its budget
discovering that learns nothing about physics. So the search space asks "what
fraction of the space that is available do we use", and the unit cube becomes
feasible by construction — `tests/ship/test_geometry.py` checks all 128 corners
and a random sample of the interior. `config.write_variant(dry_run=True)` is
still run as a final gate before anything is queued; a geometry that fails it is
recorded `INFEASIBLE` without costing a core-hour.

**Where the baseline already sits at the maximum, the scan runs downward from
it.** `xy_gap_frac` and `sitarget_fill` are both 1.0 at the baseline — the gap
uses its whole budget (4.9 mm of 4.9 mm) and the layers fill the whole envelope
— so their ranges end at 1.0 and explore below, rather than stopping short of
the design they are anchored to. (This is why `xy_gap_frac` goes to 1.0 and not
to the 0.9 first sketched: with 0.9 the baseline would fall outside the cube,
and trial 0 could not reproduce it.)

That makes the top of the box a place the optimizer visits often — acquisition
maxima like boundaries — so it has to be exactly reachable, not merely
reachable in principle. It was not, at first: `config.format_value` renders a
float as `"%g*mm"`, six significant digits, and the template compares
`SiTarget_XY_plane_gap` against its own derived maximum. Six digits round *up*
about half the time, turning `gap == gap_max` into `gap > gap_max`, and roughly
a third of the proposals sitting on that boundary were rejected. The baseline
never showed it, because 4.9 is exact in six digits. `Geometry._length_within`
now asks the renderer what it would write and falls back to a full-precision
verbatim string only when the plain float would overshoot — so the boundary is
reachable and the baseline's XML is unchanged.

### Objectives

```yaml
objectives:
  - { name: nhits_sipad, direction: max, ref_from_baseline: 0.5 }
  - { name: cost_proxy,  direction: min, ref_from_baseline: 2.0 }
```

Trial 0 is always the baseline. Once it completes, the hypervolume reference
point is anchored to it (half its hits, twice its cost) and persisted in the
store; a resume reads it back and never recomputes it, because moving the
reference mid-experiment would silently redefine the metric.

`cost_proxy` is analytic — `A · (silicon area) + B · (tungsten mass)`, summed
over SiTarget and SiPad, with `A = 30 k€/m²` and `B = 0.1 k€/kg` in
`conf/evaluator/snd.yaml`. SiTarget carries two silicon planes per layer, SiPad
one. The baseline comes out at 1585.19 k€ (44.39 m² of silicon, 2535 kg of
tungsten); every factor in it is checked by hand in `tests/ship/test_metrics.py`.

Every run records the *whole* metrics dict — hit counts for all four detectors,
areas, masses, channel counts, CPU and wall time — and the objectives select
from it. Adding an objective later therefore never means re-running anything.

### Adding an objective

1. make the payload record it in `metrics.json` (a new function in
   `ship/metrics.py`, merged in `payload.measure`);
2. name it under `objectives:` with a direction and a reference rule.

That is all. Tracking efficiency and resolution from `ACTSTracks` — already
written to `ShipHits.root` by `job5_rntuple.py` — are the intended next ones.
When they arrive with a per-trial standard error, set `optimizer.fixed_noise:
true` and have the payload also write `<metric>_stderr`; the GPs will use
per-observation noise instead of an inferred homoscedastic term.

---

## Architecture

```
src/mobo/
├── core/        generic: search space, store, GP models, acquisition, loop
├── exec/        how a trial runs: local subprocesses, HTCondor
├── ship/        everything SND-specific, behind the Evaluator interface
├── viz/         Pareto/hypervolume plots and the HTML report
└── cli.py       mobo-run / mobo-resume / mobo-status / mobo-report
```

The optimizer is **detector-agnostic**. It knows a unit cube, a vector of
objectives and an `Evaluator` with three methods:

```python
class Evaluator(ABC):
    def validate(self, trial) -> str | None:          # cheap feasibility gate
    def prepare(self, trial) -> tuple[list[str], Path] # -> command, workdir
    def collect(self, trial) -> Result                 # workdir -> metrics
```

Porting the loop to another detector means writing one of those plus a config;
`core/`, `exec/` and `viz/` are untouched. `tests/test_architecture.py` enforces
it: no import of `ship/` may appear in any of the other layers.

Two conventions worth knowing before reading the code:

* **everything is maximized internally.** An objective declared `min` is negated
  at exactly one place (`ObjectiveSpec.signed`) and nowhere else. Plots and
  reports show physical units and physical directions.
* **completion is a file, not a query.** The payload writes `DONE` or `FAILED`
  as its last act. `condor_q` is only a watchdog (held jobs, vanished jobs),
  because a finished job leaves the queue within seconds and a poll that misses
  that window cannot tell "finished" from "never existed". The local executor
  uses the same protocol, which is what makes "it works locally" evidence about
  the Condor path.

### Failures

| status | meaning | retried? |
|---|---|---|
| `INFEASIBLE` | the geometry would not build; rejected by the dry-run gate | no |
| `FAILED` | the payload ran and reported a failed step (`sim_failed`, `job4_failed`, …) | no |
| `FAILED (lost)` | the job vanished without a sentinel — eviction, node crash | yes, same seed |
| `FAILED (held)` | Condor put the job on hold | yes, same seed |
| `FAILED (timeout)` | no sentinel within `timeout_hours` | yes, same seed |

Retries keep the **same seed** on purpose: if the second attempt fails
identically, the cause is the design and not the cluster. Failed and infeasible
trials are excluded from the GP rather than modelled — inventing a sentinel
value would invent a landscape. Learning the feasible region with a classifier
is a later phase.

### Reproducibility

A global `seed` in the config derives everything: `seed_trial = hash(seed,
trial_id)` feeds the Sobol design, the GP fits, the acquisition optimizer and
**Geant4** (`run_sim.py --seed`, which reaches `SIM.random.seed`). Same seed,
same `nhits`; different seed, different `nhits` — both are asserted in
`tests/ship/test_payload_smoke.py`. The store also records the resolved config,
the git SHA and the torch/botorch versions.

---

## HTCondor

```bash
mobo-run experiment=snd_proxy executor=htcondor
```

Nothing is transferred: the payload, the geometry, the analysis chain and the
trial directory are all on AFS and visible from the worker. The generated
`job.sub` includes `MY.SendCredential = true` so the job can *write* to AFS —
without a forwarded kerberos credential every trial fails on its first write, so
check `klist` before launching a long run and renew with `kinit` if needed.

Useful knobs (`conf/executor/htcondor.yaml`): `flavour` (default `workday`),
`request_cpus`, `request_memory`, `poll_interval`, `timeout_hours`,
`grace_seconds`, `extra_submit` for site-specific lines.

`run_local: true` runs the generated wrapper as a local subprocess instead of
submitting it — the whole Condor path except `condor_submit` itself, which is
how it is tested without a scheduler.

### Manual checklist on the real cluster

The parts a unit test cannot cover. Run these once after any change to
`exec/htcondor.py`:

1. **one job end to end**
   ```bash
   mobo-run experiment=snd_proxy executor=htcondor evaluator.nevents=10 \
            experiment.loop.max_trials=1 run_dir=runs/condor_check
   ```
   → the job appears in `condor_q`, `runs/condor_check/trials/t00000/DONE`
   shows up, and the loop reports `COMPLETED` without intervention.

2. **asynchrony**
   ```bash
   mobo-run experiment=snd_proxy executor=htcondor evaluator.nevents=10 \
            experiment.loop.max_trials=8 experiment.loop.max_in_flight=4 \
            run_dir=runs/condor_async
   ```
   → four jobs run at once (`condor_q` shows four), and each completion is
   followed within one poll by a new submission. In the log, the asks after the
   first four must be conditioned on the ones still running — the timeline plot
   in `report.html` shows four overlapping bars throughout.

3. **kill and resume**
   Kill the driver (Ctrl-C twice) while jobs are in flight, then
   `mobo-resume runs/condor_async`. The jobs still in the queue must be
   re-adopted (`re-adopted t000NN` in the log), not resubmitted, and
   `mobo-status` must show no duplicate or missing trial ids.

All three were run on lxplus/HTCondor on 2026-08-01 and passed:

* one job, 10 events: submitted 13:07:35, `COMPLETED` 13:10:36, picked up from
  its sentinel with no intervention. `nhits_sipad = 420` — identical to the
  local run of the same trial with the same seed, so the Condor path changes
  nothing about the physics;
* four jobs in flight, each freed slot refilled on the next poll;
* driver killed with three jobs running → `mobo-resume` re-adopted all three by
  their original cluster ids (no resubmission, no wasted CPU) and the run
  finished with exactly `max_trials` trials, ids dense and unique.

One accidental finding worth keeping: an earlier attempt left a second driver
alive on the same store. Both kept proposing, and still no trial id was
duplicated or lost — the `BEGIN IMMEDIATE` allocation in `store.py` holds across
processes, not just threads. The run did overshoot `max_trials` (each driver
counts the budget independently), so don't do it on purpose.

---

## Tests

```bash
pytest -q                          # everything runnable here (~2 min)
pytest -q -m "not slow"            # fast subset
pytest -q -m "slow"                # GP quality + the full payload (~30 min)
pytest -q -m "not key4hep"         # outside the stack
pytest tests/exec/test_htcondor_dryrun.py --update-golden   # after changing job.sub
```

The tests that matter most:

| test | claim |
|---|---|
| `core/test_optimizer.py::test_beats_sobol` | qLogNEHVI dominates more hypervolume than pure Sobol at equal budget, on BraninCurrin, in 3 of 4 seeds |
| `core/test_store.py::test_concurrent_creates_*` | trial ids are unique and dense under concurrent writers |
| `core/test_loop.py` | the async loop, resume, retries and `max_in_flight`, over real subprocesses |
| `ship/test_geometry.py` | derived limits are exactly the boundary of what the renderer accepts; the whole cube builds; the baseline round-trips byte-identically |
| `ship/test_metrics.py` | every factor of the cost model, by hand |
| `ship/test_payload_smoke.py` | the real chain, and that the seed reaches Geant4 |
| `test_architecture.py` | `ship/` never leaks into `core/`, `exec/` or `viz/`; the payload imports nothing a worker node lacks |

---

## Changes made outside this package

Deliberately small, all inside `Optimization/` except the last one:

| file | change |
|---|---|
| `Optimization/Simulation/run_scripts/run_sim.py` | `--seed` → `randomSeed` in the steering. Default unchanged (no seed = ddsim's own). |
| `Optimization/Analysis/job4_tracking.py` | `OutputLevel` from `SND_OUTPUT_LEVEL`, default `DEBUG` so manual use is unchanged; the payload sets `INFO`, since DEBUG writes hundreds of MB per trial. |
| `Optimization/Simulation/geometry/config.py` | the envelope-overlap check gained a 1 nm tolerance (`OVERLAP_TOL`). SiTarget is *meant* to end exactly where SiPad begins, and the two z values come from different expressions; for a `SiPad_dim_z` that is not a round number they differ in the last bits, and a strict `>` rejected perfectly valid geometries. Only reachable with arbitrary floats, which is why hand-written variants never hit it. |

`controller.py`, `make_variants.py` and the analysis jobs keep their current
behaviour; `tests/ship/test_geometry.py::test_regression_baseline_xml_is_unchanged`
pins the baseline XML against the old path.

---

## Not yet done (deliberately)

- physics objectives: tracking efficiency and resolution from `ACTSTracks`,
  with per-trial uncertainties → `optimizer.fixed_noise: true`;
- learned feasibility (a GP classifier over failures) and output constraints;
- multi-fidelity (event count as the fidelity);
- optimizing the tracker's own knobs in `job4_tracking.py`.

The hooks for the first are already in place: `metrics.json` is open-ended and
objectives are selected by config.
