#!/usr/bin/env python3
"""What one trial actually runs — on a worker node, alone.

Replays the chain `controller.run_variant()` performs (ddsim -> job1 -> job4 ->
job5 -> metrics) for a single variant -- with, when the trial spec carries a
`genie` block, a generation step in front of it (geoConverter -> GDML frame
shift -> gevgen_fnal -> gst -> HepMC), so that the neutrino interactions are
sampled over THIS variant's material and not a shared one's -- and reports
through the filesystem:
`metrics.json` plus a `DONE` or `FAILED` sentinel written as the very last act.
The executor watches for that sentinel and nothing else, which is what lets the
same script be a local subprocess or an HTCondor job with no changes.

Deliberately standalone: stdlib only, no import of the `mobo` package, no
virtualenv. It runs with the bare key4hep python that `init_key4ship.sh`
provides, because that is all a worker node is guaranteed to have. The two
project modules it does use (`metrics.py` next to it, `parse_geometry.py` from
the shared geometry) are loaded by path, the way `controller.py` loads
`config.py` and for the same reason: key4hep's PYTHONPATH is crowded and a
plain import of a generically-named module is a coin toss.

    python3 payload.py --workdir <trial dir>

The trial dir must contain a `trial.json` written by `ship/evaluator.py`.

Style note: this file keeps %-formatting throughout, to match `controller.py`
and `run_sim.py`, which are its neighbours in every practical sense.
"""

# ruff: noqa: UP031
import argparse
import contextlib
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

DONE = "DONE"
FAILED = "FAILED"
LOG_TAIL_LINES = 30


def load_module(name, path):
    """Import a file by path (see the module docstring for why)."""
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class StepFailed(Exception):
    """One step of the chain failed; carries what the log had to say."""

    def __init__(self, step, message, log=None):
        self.step = step
        self.log = log
        super().__init__(message)


def tail(path, n=LOG_TAIL_LINES):
    try:
        lines = Path(path).read_text(errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-n:])


# Wall seconds per step, accumulated by run_step() and copied into
# metrics.json as time_<step>_s: with GENIE, ddsim and the tracker all in one
# job, "the trial took 4 h" is not an answer to where the time went.
STEP_TIMES: dict[str, float] = {}


def run_step(step, cmd, cwd, env, log_path, tolerate=None):
    """Run one subprocess, both streams into log_path. Raises on failure.

    `tolerate(log_text, returncode)` gets a chance to pardon a non-zero exit —
    ddsim is known to abort during ROOT/cppyy teardown *after* the run finished
    and the output file is complete, and failing a trial for that would throw
    away a perfectly good simulation.
    """
    t0 = time.time()
    with open(log_path, "w") as fh:
        proc = subprocess.run(
            [str(c) for c in cmd],
            cwd=str(cwd),
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
        )
    STEP_TIMES[step] = STEP_TIMES.get(step, 0.0) + (time.time() - t0)
    if proc.returncode == 0:
        return
    if tolerate is not None and tolerate(
        Path(log_path).read_text(errors="replace"), proc.returncode
    ):
        print(
            "WARNING: %s exited with %d but its output is complete; continuing"
            % (step, proc.returncode)
        )
        return
    raise StepFailed(
        step, "%s exited with code %d" % (step, proc.returncode), log=tail(log_path)
    )


def ddsim_teardown_ok(log_text, _returncode):
    """ddsim crashed in teardown after completing the event loop."""
    return "Finished run 0 after" in log_text


def simulate(spec, workdir, env, hepmc=None):
    """ddsim, through the project's own run_sim.py so the steering stays shared.

    With `hepmc` (a per-trial file from generate()) the whole file is consumed
    (--events -1): it holds exactly this geometry's interactions, and asking
    run_sim_genie.py for a fixed count would either truncate it or run into
    EOF.
    """
    sim_file = workdir / spec["sim_file"]
    if spec.get("skip_sim") and sim_file.is_file():
        print(f"reusing the existing simulation {sim_file}")
        return sim_file

    cmd = [
        sys.executable,
        spec["run_sim"],
        "--compact",
        spec["compact"],
        "--var-dir",
        str(workdir),
        "--events",
        "-1" if hepmc is not None else str(spec["nevents"]),
        "--seed",
        str(spec["seed"]),
    ]
    if hepmc is not None:
        cmd += ["--hepmc", str(hepmc)]
    run_step(
        "sim",
        cmd,
        workdir,
        env,
        workdir / "run_sim.log",
        tolerate=lambda text, rc: ddsim_teardown_ok(tail(workdir / "ddsim.log", 200), rc),
    )
    if not sim_file.is_file():
        raise StepFailed(
            "sim", f"ddsim produced no {sim_file.name}", log=tail(workdir / "ddsim.log")
        )
    return sim_file


# Environment variables a clean shell for the SHiP stack keeps. The payload
# runs inside key4hep, and the two stacks cannot share LD_LIBRARY_PATH; so the
# GENIE step gets `env -i` plus exactly this: identity, a system PATH, the
# kerberos cache (condor's forwarded credential, needed to read the flux and
# splines on EOS) and condor's scratch settings.
CLEAN_ENV_KEEP = (
    "HOME",
    "USER",
    "LOGNAME",
    "LANG",
    "TERM",
    "KRB5CCNAME",
    "TMPDIR",
    "_CONDOR_SCRATCH_DIR",
    "_CONDOR_JOB_AD",
    "_CONDOR_MACHINE_AD",
    # run_gevgen.sh's own overrides, when the submitter set them.
    "GENIE_FLUX",
    "GENIE_XSEC",
    "GENIE_TUNE",
    "GENIE_FLAVOURS",
    "SHIP_SETUP",
    "GENIE_CONFIG_DIR",
)
CLEAN_PATH = "/usr/local/bin:/usr/bin:/bin"


def clean_env(env):
    kept = {k: env[k] for k in CLEAN_ENV_KEEP if k in env}
    kept["PATH"] = CLEAN_PATH
    return kept


def count_hepmc_events(path):
    """Events in a HepMC3 ascii file: one `E ` line each."""
    n = 0
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("E "):
                n += 1
    return n


def read_genie_pot(ghep_file, log_file):
    """The POT normalization of a gevgen_fnal sample, or None.

    gevgen_fnal stores it as the weight of the ghep file's gtree; the log has
    the same number ("Normalization for generated sample"). Read the file
    first (with GENIE's dictionaries absent ROOT complains about the branches
    but the tree weight is plain TTree metadata), and fall back to the log.
    """
    try:
        import ROOT

        ROOT.gErrorIgnoreLevel = ROOT.kFatal
        f = ROOT.TFile.Open(str(ghep_file))
        if f and not f.IsZombie():
            tree = f.Get("gtree")
            if tree:
                pot = float(tree.GetWeight())
                f.Close()
                if pot > 0:
                    return pot
    except Exception:  # noqa: BLE001 - the log is the fallback
        pass
    try:
        for line in Path(log_file).read_text(errors="replace").splitlines():
            if "Normalization for generated sample" in line:
                return float(line.split(":")[-1].split()[0])
    except (OSError, ValueError, IndexError):
        pass
    return None


def generate(spec, workdir, env):
    """GENIE on this variant's own geometry -> a HepMC file for ddsim.

    Four steps, each a subprocess with its own log in the trial directory:
      1. geoConverter -compact2gdml     geometry.xml -> geometry.gdml
      2. shift_gdml_frame.py            -> geometry_fairship.gdml
         (the gsimple flux is in the FairShip frame and gevgen_fnal cannot
         translate it, so the geometry goes to the flux; see the script)
      3. run_gevgen.sh (SHiP stack)     -> genie/gntp.<run>.{ghep,gst}.root
      4. genie_gst_to_hepmc.py          -> genie.hepmc, vertices back in the
                                           key4ship frame
    Returns (hepmc path, number of events in it, POT of the sample).
    """
    g = spec["genie"]
    gdml = workdir / "geometry.gdml"
    gdml_shifted = workdir / "geometry_fairship.gdml"
    genie_dir = workdir / "genie"
    hepmc = workdir / "genie.hepmc"
    run = int(spec["trial_id"])
    if spec.get("skip_sim") and hepmc.is_file():
        print(f"reusing the existing GENIE sample {hepmc}")
        return hepmc, count_hepmc_events(hepmc), None

    run_step(
        "gdml",
        [
            "geoConverter",
            "-compact2gdml",
            "-input",
            spec["compact"],
            "-output",
            str(gdml),
        ],
        workdir,
        env,
        workdir / "gdml.log",
    )
    if not gdml.is_file():
        raise StepFailed(
            "gdml", "geoConverter wrote no GDML", log=tail(workdir / "gdml.log")
        )

    run_step(
        "gdml",
        [
            sys.executable,
            g["shift_gdml"],
            str(gdml),
            str(gdml_shifted),
            "--dz-cm",
            str(g["dz_cm"]),
        ],
        workdir,
        env,
        workdir / "gdml_shift.log",
    )

    genie_dir.mkdir(exist_ok=True)
    run_step(
        "genie",
        [
            "env",
            "-i",
            *("%s=%s" % kv for kv in clean_env(env).items()),
            "bash",
            g["run_gevgen"],
            str(gdml_shifted),
            str(genie_dir),
            str(spec["nevents"]),
            str(run),
            str(spec["seed"]),
            str(genie_dir),
        ],
        workdir,
        env,
        workdir / "genie.log",
    )
    gst = genie_dir / ("gntp.%d.gst.root" % run)
    ghep = genie_dir / ("gntp.%d.ghep.root" % run)
    if not gst.is_file():
        raise StepFailed(
            "genie",
            "gevgen_fnal/gntpc produced no %s" % gst.name,
            log=tail(genie_dir / ("gevgen_%d.log" % run)),
        )
    pot = read_genie_pot(ghep, genie_dir / ("gevgen_%d.log" % run))

    run_step(
        "hepmc",
        [
            sys.executable,
            g["convert"],
            "--input",
            str(gst),
            "--output",
            str(hepmc),
            "--z-shift-m",
            str(g["z_shift_m"]),
            "--world-half-size-m",
            str(g["world_half_size_m"]),
        ],
        workdir,
        env,
        workdir / "hepmc.log",
    )
    n = count_hepmc_events(hepmc) if hepmc.is_file() else 0
    if n == 0:
        raise StepFailed(
            "hepmc", "no events in the HepMC file", log=tail(workdir / "hepmc.log")
        )
    print("GENIE: %d interactions in this geometry, %s POT" % (n, pot))
    return hepmc, n, pot


# What each analysis step consumes and nobody after it needs, keyed by the
# job's short name (job1_overlay -> "job1"). Deleted right after the step
# unless keep_intermediate is set: the chain is a pipeline of copies (job1
# rewrites the ddsim output as events.edm4hep.root, job4 rewrites that with
# the tracks added), so without this a trial holds two or three full-size
# copies of its events for most of its life -- ~1.9 GB at 10k events, of
# which the ddsim output alone is 760 MB kept alive through the whole
# tracking step.
CONSUMED_BY = {
    "job1": ("{sim_file}",),
    "job4": ("events.edm4hep.root",),
}


def drop_consumed(step, workdir, spec):
    """Delete the inputs `step` has just used up. Returns the bytes freed."""
    if spec.get("keep_intermediate"):
        return 0
    freed = 0
    for name in CONSUMED_BY.get(step, ()):
        path = Path(workdir) / name.format(sim_file=spec["sim_file"])
        try:
            freed += path.stat().st_size
            path.unlink()
        except OSError:
            continue
    if freed:
        print("%s done; dropped its input (%.0f MB)" % (step, freed / 1048576.0))
    return freed


def analyse(spec, workdir, env, sim_file):
    """job1 -> job4 -> job5, exactly as the controller chains them.

    The jobs use bare filenames, so the working directory has to be the trial's.
    Each job's input is deleted as soon as the job has finished with it (see
    CONSUMED_BY), unless the trial keeps its intermediates.
    """
    env = dict(env)
    env["SND_COMPACT"] = spec["compact"]
    env["INPUT_FILE"] = str(sim_file)
    # job4 logs at DEBUG by default, which is right for one interactive variant
    # and unusable at the scale of an optimization run (hundreds of MB per job).
    env["SND_OUTPUT_LEVEL"] = spec.get("output_level", "INFO")

    analysis_dir = Path(spec["analysis_dir"])
    for job in spec.get("jobs") or ("job1_overlay", "job4_tracking", "job5_rntuple"):
        short = job.split("_")[0]
        run_step(
            short,
            ["k4run", str(analysis_dir / (job + ".py"))],
            workdir,
            env,
            workdir / (short + ".log"),
        )
        drop_consumed(short, workdir, spec)


def measure(spec, workdir, nevents=None):
    """Everything about this trial that can be turned into a number.

    `nevents` is what ddsim actually ran -- for a generated sample that is the
    HepMC's count (interactions outside the world are dropped), not the
    requested one -- and is the denominator of every per-event metric.
    """
    if nevents is None:
        nevents = int(spec["nevents"])
    metrics_mod = load_module("mobo_ship_metrics", HERE / "metrics.py")
    parse_geometry = load_module("mobo_parse_geometry", Path(spec["parse_geometry"]))

    # Analytic metrics come from the constants of the compact that was actually
    # simulated, so they describe the geometry that produced the hits and not
    # what we meant to build.
    constants = parse_geometry.SNDGeometry(spec["compact"]).constants
    cost = metrics_mod.CostModel(**spec.get("cost_model", {}))
    metrics = metrics_mod.analytic_metrics(constants, cost)

    hits_file = workdir / spec.get("hits_file", "ShipHits.root")
    if not hits_file.is_file():
        raise StepFailed("metrics", f"no {hits_file.name}")
    hits = metrics_mod.hit_metrics(str(hits_file))
    if "nhits_sipad" not in hits:
        raise StepFailed(
            "metrics",
            f"could not read the SiPad RNTuple from {hits_file.name}",
            log=tail(workdir / "job5.log"),
        )
    metrics.update(hits)

    # Track metrics are opt-in: only a chain whose job5 writes the track
    # RNTuples has them, and an absent container must not be read as zero.
    if spec.get("track_metrics"):
        tracks = metrics_mod.track_metrics(str(hits_file), int(nevents))
        if "eff_track" not in tracks:
            raise StepFailed(
                "metrics",
                "could not read the GlobalTracks RNTuple from %s" % hits_file.name,
                log=tail(workdir / "job5.log"),
            )
        metrics.update(tracks)
    return metrics


def discard_intermediates(workdir, spec):
    """Delete the ROOT files a finished trial no longer needs.

    A trial writes ~37 MB of them -- the ddsim output, the overlay, the tracks
    and the RNTuple -- and every one of those numbers has already been turned
    into metrics.json by the time this runs. Keeping them costs 7.4 GB over a
    200-trial campaign, which does not fit in a 10 GB AFS quota: measured the
    hard way, when ddsim started failing with `SysError in <TFile::Flush>` on a
    full volume and the loop reported it as `job1_failed` two steps later.

    What survives is what makes the trial reproducible and diagnosable:
    metrics.json, params.yaml, trial.json, geometry.xml and every log. Set
    `keep_intermediate: true` in the evaluator config to keep every ROOT file
    -- this final sweep AND the on-the-fly deletion in analyse() -- e.g. to
    feed one variant to the event display.
    """
    if spec.get("keep_intermediate"):
        return 0
    freed = 0
    doomed = list(Path(workdir).glob("*.root"))
    # A generated trial adds its own: the two GDMLs (2.4 MB each), the GENIE
    # ghep/gst pair and the HepMC. metrics.json already has what they said.
    doomed += list(Path(workdir).glob("*.gdml"))
    doomed += list((Path(workdir) / "genie").glob("*.root"))
    doomed += list(Path(workdir).glob("*.hepmc"))
    for path in sorted(doomed):
        try:
            freed += path.stat().st_size
            path.unlink()
        except OSError:
            continue
    return freed


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workdir", required=True, help="the trial directory")
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    started_wall = time.time()
    started_cpu = os.times()

    # Whatever happens from here on, this process leaves exactly one sentinel.
    try:
        spec = json.loads((workdir / "trial.json").read_text())
        env = os.environ.copy()

        nevents = int(spec["nevents"])
        genie = {}
        hepmc = None
        if spec.get("genie"):
            hepmc, nevents, pot = generate(spec, workdir, env)
            genie["nevents_generated"] = float(nevents)
            if pot:
                # The interaction rate of THIS geometry: what a fixed -n
                # throws away, and the one number that scales with material.
                genie["genie_pot"] = float(pot)
                genie["n_int_per_1e18pot"] = nevents / float(pot) * 1e18
        sim_file = simulate(spec, workdir, env, hepmc=hepmc)
        analyse(spec, workdir, env, sim_file)
        metrics = measure(spec, workdir, nevents=nevents)
        metrics.update(genie)

        finished_cpu = os.times()
        metrics["wall_time_s"] = time.time() - started_wall
        # The work happens in subprocesses, so it is the children's CPU that counts.
        metrics["cpu_time_s"] = (
            finished_cpu.children_user - started_cpu.children_user
        ) + (finished_cpu.children_system - started_cpu.children_system)
        metrics["nevents"] = float(nevents)
        metrics.update({"time_%s_s" % k: v for k, v in STEP_TIMES.items()})

        (workdir / "metrics.json").write_text(
            json.dumps(metrics, indent=2, sort_keys=True)
        )
        # After metrics.json exists, never before: a crash here must not lose
        # the measurement, and the sentinel goes last as always.
        freed = discard_intermediates(workdir, spec)
        if freed:
            print("discarded %.1f MB of intermediate ROOT files" % (freed / 1048576.0))
        (workdir / DONE).write_text("ok\n")
        print(f"metrics: {json.dumps(metrics, sort_keys=True)}")
        return 0

    except StepFailed as exc:
        _fail(workdir, f"{exc.step}_failed: {exc}", exc.log)
        return 1
    except Exception as exc:  # noqa: BLE001 - the sentinel is the error channel
        import traceback

        _fail(workdir, f"{type(exc).__name__}: {exc}", traceback.format_exc())
        return 1


def _fail(workdir, reason, detail=None):
    """One line of cause, then the evidence — the store keeps the first line."""
    text = reason.strip() + "\n"
    if detail:
        text += "\n" + detail.strip() + "\n"
    with contextlib.suppress(OSError):
        (workdir / FAILED).write_text(text)
    print(f"FAILED: {reason}", file=sys.stderr)
    if detail:
        print(detail, file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
