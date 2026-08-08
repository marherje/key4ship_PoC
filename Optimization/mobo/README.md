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
| `SiPad_WThickness` | 10 – 30 mm | 10 | |
| `SiPad_dim_z` | 220 – 1000 mm | 220 | also sets `SiTarget_dim_z = SiDetector_total_dim_z − this` (2200 mm) |
| `SiPad_NLayers` | 10 – 20 (int) | 10 | |
| `SiTarget_WThickness` | 2 – 5 mm | 3.5 | |
| `SiTarget_NLayers` | 80 – 120 (int) | 120 | |
| `SiTarget_XY_plane_gap_frac` | 0.0 – 1.0 | 1 | 0 = Y plane against the X plane, 1 = against the next W |

`SiPad_frame_gap` (0.1 mm), `SiTarget_module_offset` (1 mm) and both layer gaps
(`auto`) are held at the baseline under `fixed:`; freeing one is a two-line
config change, not a code change.

**Equidistant layers, so the counts can be searched directly.** Both detectors
carry their closing gap on `auto`, which sizes it so that `NLayers` layers span
`dim_z` exactly. The pitch is therefore `dim_z / NLayers` and the sensitive
planes come out evenly spread over the whole envelope. Two consequences:

* **`NLayers` is the sampling-pitch knob**, and it can be an ordinary integer
  parameter. Its only bound is the rigid content of a layer — the absorber plus
  the sensors and the services — which no longer moves with a second quantity
  being optimized at the same time.
* **There is no `SiTarget_spacing` parameter.** A SiTarget layer *is* the W-to-W
  pitch, so a free spacing was never independent of the count: all it decided on
  its own was how much dead air each layer carried, and since the objectives
  count silicon area and tungsten mass — both proportional to `NLayers`, neither
  to the pitch — it was free in the cost function. A nuisance dimension. It is
  now an output, `SiTarget_dim_z / SiTarget_NLayers`.

Asking for fewer layers than fit does not shorten a detector, it dilutes it: the
gap grows until the layers again span the envelope.

*What still needs a fraction.* The distance between the two Si planes is bounded
by what the pitch leaves over once the rigid content is placed — a budget that
genuinely is derived, and that moves with every other parameter. So it is not
searched as a length at all: `SiTarget_XY_plane_gap_frac` says *where the Y
plane sits* in that air, 0 flush against the X plane and 1 flush against the
tungsten of the next layer. It decides how a layer's air is split between
"between the two Si planes" and "after them"; it never changes the pitch.

The interpolation lives in `SND_compact_template.xml`, next to the budget it is
a fraction of, not in the python. That is why `Geometry.constants_for` is a
pass-through and `baseline_params` is the identity: every search parameter is a
template constant, so nothing outside the template can disagree with it about
what a parameter means. `simulation/geometry/SND_compact.xml` carries the same
parametrization — same constant name, same 0/1 endpoints — over its own pitch,
which there is the `SiTarget_spacing` input rather than `dim_z / NLayers`.

*What this costs.* The box is no longer feasible by construction. Far corners
like "20 SiPad layers of 30 mm of tungsten in a 250 mm envelope" ask for 711 mm
of material in 250 mm and simply do not exist. At the ranges above roughly 80%
of the cube builds; `tests/ship/test_geometry.py` asserts that a majority does,
and — more usefully — that every rejection is a *fitting* rejection, so the
parametrization can never be what mangled a feasible point.
`config.write_variant(dry_run=True)` is the final gate before anything is
queued, and a geometry that fails it is recorded `INFEASIBLE` without costing a
core-hour.

*The two SiTarget bounds trade against each other.* `SiTarget_NLayers_max` is
the count that fits with **no** XY gap, and `SiTarget_XY_plane_gap_max` is the
gap that fits at the count actually asked for. Pinning both at once asks for a
layer thinner than its own contents — which is why the feasibility test probes
them at separate points rather than together.

**The top of the range has to be exactly reachable.** The baseline sits at
`SiTarget_XY_plane_gap_frac = 1` — the Y plane against the next absorber — so
trial 0 is itself on the boundary of the box, and boundaries are where the
optimizer keeps going back, since acquisition maxima like them. Reachable *in
principle* is not enough.

It was not reachable, at first, back when this layer turned the fraction into a
length itself. `config.format_value` renders a float as `"%g*mm"`, six
significant digits, and the template then compares `SiTarget_XY_plane_gap`
against its own derived maximum. Six digits round *up* about half the time,
turning `gap == gap_max` into `gap > gap_max`, and roughly a third of the
proposals on that boundary were rejected.

Two things fix it for good, and both are in the template now. The fraction is
written at full precision (17 significant digits round-trip a float exactly),
and the gap is interpolated as `max − (1 − frac)·(max − min)` rather than the
algebraically identical `min + frac·(max − min)`: only the first form yields
*exactly* `max` at `frac = 1`. There is no longer a length rendered at six
digits anywhere on the path, so the failure mode cannot come back;
`tests/ship/test_geometry.py` asserts the equality holds for every variant.

The other endpoint is not 0 either, for a reason worth knowing: the XY gap is a
`<slice>` of `SiTarget.xml` like any other, and Geant4 aborts on a box of null
extent (`G4Box::G4Box, "Dimensions too small for Solid"`). "The two planes are
touching" is rendered as the same nanometre of air (`AUTO_MIN`) that every other
gap with nothing left to give collapses to.

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

## How qLogNEHVI picks the next geometry

Two objectives in conflict have no single optimum, only a Pareto front, and an
optimizer needs a scalar. That scalar is the **dominated hypervolume**: given a
reference point (here anchored to the baseline at half its hits and twice its
cost), the front dominates a region of the objective plane, and its area is the
metric. Improving the front means enlarging that area.

**1. The surrogate.** `build_model` (`core/models.py`) fits one independent
`SingleTaskGP` per objective on the completed trials. Each GP predicts, at any
untried geometry, not a number but a posterior: a mean and an uncertainty.
Independent on purpose — a simulated hit count and an analytic cost share no
latent structure worth the extra hyperparameters.

**2. The acquisition.** For a candidate *x*, ask: how much would the
hypervolume grow if we simulated it? The answer is a distribution, so take its
expectation over the GP posteriors, estimated by quasi-MC with `mc_samples`
(128) draws. Reading the name backwards:

| piece | meaning |
|---|---|
| **HVI** | hypervolume improvement — the area this point would add |
| **E** | expected, integrated over the GP posterior. This is where exploration comes from: a point with a mediocre mean but a large variance can still have a high *expected* improvement, so uncertainty attracts on its own |
| **N** | *noisy* — the current front is itself uncertain (a hit count from 100 events is an estimate), so observed values are not taken as truth; the posterior at the measured points is integrated over too. That is what `X_baseline` is for |
| **q** | proposes *q* points at once; `batch_size: 1` here |
| **Log** | evaluated in log space. Once a decent front exists, plain EHVI underflows to zero across most of the cube and its gradient flattens to exactly zero in floating point, so the optimizer cannot climb. The log formulation keeps the gradient informative |

**3. Maximizing it.** `optimize_acqf` runs multi-start gradient ascent over the
unit cube: `raw_samples` (512) draws to find promising starts, `num_restarts`
(10) of them refined for up to 200 iterations. This is cheap — it only queries
the GPs, never the simulator. Its maximizer is the next trial.

**4. The asynchronous part.** When the loop asks for a trial there are up to
`max_in_flight - 1` others still running. Their locations are known even though
their outcomes are not, so they are passed as **`X_pending`** and the
acquisition integrates over what they might return, discounting improvement
they may already be delivering. Without it every slot of the batch would get a
near-duplicate of the same point, and `max_in_flight` cores would do one core's
work.

Two lesser knobs: `prune_baseline` drops dominated points from `X_baseline`,
which cannot change the front and only cost time; `sequential` picks a q-batch
greedily rather than jointly, and is irrelevant at q=1.

Everything above is in the maximization convention — the sign of a `min`
objective is flipped in exactly one place (`ObjectiveSpec.signed`), and the
plots put it back.

### Knowing when the front has been found

Strictly, you never do. `pareto_mask` (`core/pareto.py`) is BoTorch's
`is_non_dominated` over the trials that were actually run: it reports which of
*your observations* nothing else dominates. It says nothing about the true
front, which for a black-box simulator over a 7D continuum is not computable —
there is no analytic solution to compare against and no way to evaluate the
whole space.

What is measured instead is **convergence of the dominated hypervolume**.
`hypervolume_trace` recomputes the hypervolume over the first *i* trials in
trial order, giving a curve that is monotone non-decreasing by construction —
`viz/progress.py` calls it "the one plot that says whether the optimization is
working", and notes that a *drop* is impossible unless the reference point
moved, i.e. unless something rewrote history. A plateau means the optimizer has
stopped finding improvements; that is the practical stopping signal, and it is
evidence of diminishing returns, not proof of optimality.

Three caveats worth keeping in mind when reading that curve:

* **A plateau can be local.** In 7D, a few dozen trials is a thin sample. The
  curve flattening says the acquisition found nothing better nearby, not that
  nothing better exists.
* **Compare against the Sobol phase.** The honest test of whether the model
  earns its keep is whether the hypervolume grows faster once qLogNEHVI takes
  over than it did under the quasi-random design. If the two slopes match, the
  GP is not adding anything.
* **The front is as noisy as the metrics.** `pareto_mask` runs on observed
  values, so a point can look non-dominated through Monte-Carlo luck. The
  acquisition accounts for that internally (the *N* in the name); the reported
  front does not.

There is deliberately **no convergence-based stopping rule** in the loop: it
ends on `max_trials` or `time_budget_hours` and nothing else. Stopping early on
a plateau would need a noise model for the plateau itself, which is not worth
it while the objectives are proxies.

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
