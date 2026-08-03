#!/usr/bin/env python3
"""What one trial actually runs — on a worker node, alone.

Replays the chain `controller.run_variant()` performs (ddsim -> job1 -> job4 ->
job5 -> metrics) for a single variant, and reports through the filesystem:
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


def run_step(step, cmd, cwd, env, log_path, tolerate=None):
    """Run one subprocess, both streams into log_path. Raises on failure.

    `tolerate(log_text, returncode)` gets a chance to pardon a non-zero exit —
    ddsim is known to abort during ROOT/cppyy teardown *after* the run finished
    and the output file is complete, and failing a trial for that would throw
    away a perfectly good simulation.
    """
    with open(log_path, "w") as fh:
        proc = subprocess.run(
            [str(c) for c in cmd],
            cwd=str(cwd),
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
        )
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


def simulate(spec, workdir, env):
    """ddsim, through the project's own run_sim.py so the steering stays shared."""
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
        str(spec["nevents"]),
        "--seed",
        str(spec["seed"]),
    ]
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


def analyse(spec, workdir, env, sim_file):
    """job1 -> job4 -> job5, exactly as the controller chains them.

    The jobs use bare filenames, so the working directory has to be the trial's.
    """
    env = dict(env)
    env["SND_COMPACT"] = spec["compact"]
    env["INPUT_FILE"] = str(sim_file)
    # job4 logs at DEBUG by default, which is right for one interactive variant
    # and unusable at the scale of an optimization run (hundreds of MB per job).
    env["SND_OUTPUT_LEVEL"] = spec.get("output_level", "INFO")

    analysis_dir = Path(spec["analysis_dir"])
    for job in ("job1_overlay", "job4_tracking", "job5_rntuple"):
        short = job.split("_")[0]
        run_step(
            short,
            ["k4run", str(analysis_dir / (job + ".py"))],
            workdir,
            env,
            workdir / (short + ".log"),
        )


def measure(spec, workdir):
    """Everything about this trial that can be turned into a number."""
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
    return metrics


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

        sim_file = simulate(spec, workdir, env)
        analyse(spec, workdir, env, sim_file)
        metrics = measure(spec, workdir)

        finished_cpu = os.times()
        metrics["wall_time_s"] = time.time() - started_wall
        # The work happens in subprocesses, so it is the children's CPU that counts.
        metrics["cpu_time_s"] = (
            finished_cpu.children_user - started_cpu.children_user
        ) + (finished_cpu.children_system - started_cpu.children_system)
        metrics["nevents"] = float(spec["nevents"])

        (workdir / "metrics.json").write_text(
            json.dumps(metrics, indent=2, sort_keys=True)
        )
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
