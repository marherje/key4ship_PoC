#!/usr/bin/env python3
"""Orchestrate the geometry optimization loop: variants -> sim -> analysis -> FOM.

Reads a variants spec (same format as make_variants.py), and for each entry
assigns a run id varNNNNN, renders the compact XML, simulates it with ddsim,
runs the Gaudi analysis chain on it and collects the figure of merit into
results.csv. Local and sequential; a variant that fails does not abort the rest.

Usage (needs `source init_key4ship.sh` for ddsim/k4run):
    python3 controller.py                       # every variant of the spec
    python3 controller.py --dry-run             # only render + validate
    python3 controller.py --only thin_W         # a subset (repeatable)
    python3 controller.py --events 5            # shorter simulations
    python3 controller.py --only thin_W --skip-sim   # reuse an existing sim
"""

import argparse
import csv
import importlib.util
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
GEO_DIR = HERE / "Simulation" / "geometry"
RUN_SCRIPTS = HERE / "Simulation" / "run_scripts"
SIM_OUT = HERE / "Simulation" / "output" / "variants"
ANA_DIR = HERE / "Analysis"
ANA_OUT = ANA_DIR / "variants"
RESULTS = HERE / "results.csv"

SIM_FILE = "output_mu_pi.edm4hep.root"
VAR_RE = re.compile(r"var(\d{5})$")


def _load_module(name, path):
    """Import a file by path.

    `config` and `make_variants` are generic names and key4hep's PYTHONPATH is
    crowded, so a plain `sys.path.insert` + import could pick up something else.
    """
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


geo_config = _load_module("geo_config", GEO_DIR / "config.py")
# make_variants.py does a plain `import config`, which from here would either
# fail or pick up an unrelated config.py from key4hep's PYTHONPATH. Binding the
# name to the module we just loaded makes it resolve to ours, and to the very
# same instance the controller uses.
sys.modules.setdefault("config", geo_config)
make_variants = _load_module("geo_make_variants", GEO_DIR / "make_variants.py")


# ── environment ─────────────────────────────────────────────────────────────

def check_environment():
    """ddsim and k4run must be on PATH; they only are inside key4hep."""
    missing = [t for t in ("ddsim", "k4run") if shutil.which(t) is None]
    if missing:
        sys.exit("No se encuentra %s en el PATH: haz `source init_key4ship.sh` "
                 "antes de lanzar el controller." % ", ".join(missing))


# ── run ids ─────────────────────────────────────────────────────────────────

def next_var_id():
    """var00001, var00002, … — the highest existing id plus one."""
    used = [int(m.group(1)) for d in SIM_OUT.glob("var*") if d.is_dir()
            for m in [VAR_RE.match(d.name)] if m]
    return "var%05d" % ((max(used) + 1) if used else 1)


def existing_var_id(name):
    """Id of an already-simulated variant, matched by name in its params.yaml.

    --skip-sim must not mint new ids: it reuses the simulation of a previous
    run, so it has to find which var dir that run wrote.
    """
    for d in sorted(SIM_OUT.glob("var*")):
        params = d / "params.yaml"
        if not params.is_file():
            continue
        try:
            data = yaml.safe_load(params.read_text()) or {}
        except yaml.YAMLError:
            continue
        if data.get("name") == name:
            return d.name
    return None


def write_params_snapshot(var_dir, var_id, name, constants):
    """Record what this id was, so a run can be reproduced from its output."""
    (var_dir / "params.yaml").write_text(yaml.safe_dump(
        {"var_id": var_id, "name": name, "constants": constants},
        sort_keys=True, default_flow_style=False))


# ── steps ───────────────────────────────────────────────────────────────────

def run_step(cmd, cwd, env, log_path):
    """Run one subprocess, both streams to log_path. True if it succeeded."""
    with open(log_path, "w") as fh:
        proc = subprocess.run([str(c) for c in cmd], cwd=str(cwd), env=env,
                              stdout=fh, stderr=subprocess.STDOUT)
    return proc.returncode == 0


def run_variant(var_id, name, constants, args):
    """Full chain for one variant; returns the results.csv row."""
    row = {"var_id": var_id, "name": name, "constants": constants,
           "nhits": "", "status": "ok"}

    sim_dir = SIM_OUT / var_id
    ana_dir = ANA_OUT / var_id
    if not args.dry_run:
        sim_dir.mkdir(parents=True, exist_ok=True)
        ana_dir.mkdir(parents=True, exist_ok=True)

    # 1. geometry
    compact = GEO_DIR / "variants" / ("%s.xml" % var_id)
    try:
        ok, lines, errors = geo_config.write_variant(constants, compact,
                                                     dry_run=args.dry_run)
    except RuntimeError as exc:
        # Raised before the geometry can even be rendered (a parameter that
        # cannot be sized, an unresolved placeholder…). Same outcome as a failed
        # validation, so report it as such instead of as a generic exception.
        print("  !! ERROR:", exc)
        row["status"] = "geometry_error"
        return row
    print("\n".join(lines))
    if not ok:
        for e in errors:
            print("  !! ERROR:", e)
        row["status"] = "geometry_error"
        return row

    if args.dry_run:
        row["status"] = "dry-run"
        return row

    write_params_snapshot(sim_dir, var_id, name, constants)

    # 2. simulation
    sim_file = sim_dir / SIM_FILE
    if args.skip_sim:
        if not sim_file.is_file():
            print("  !! ERROR: --skip-sim pero no existe %s" % sim_file)
            row["status"] = "missing_sim"
            return row
        print("  (--skip-sim: reutilizando %s)" % sim_file)
    else:
        print("  Simulando (%d eventos)…" % args.events)
        if not run_step([sys.executable, RUN_SCRIPTS / "run_sim.py",
                         "--compact", compact, "--var-dir", sim_dir,
                         "--events", args.events],
                        cwd=sim_dir, env=os.environ.copy(),
                        log_path=sim_dir / "run_sim.log"):
            print("  !! ERROR: la simulacion fallo, ver %s" % (sim_dir / "ddsim.log"))
            row["status"] = "sim_failed"
            return row

    # 3. analysis — bare filenames, so cwd must be the variant dir
    env = os.environ.copy()
    env["SND_COMPACT"] = str(compact)
    env["INPUT_FILE"] = str(sim_file)
    for job in ("job1_overlay", "job4_tracking", "job5_rntuple"):
        print("  %s…" % job)
        if not run_step(["k4run", ANA_DIR / ("%s.py" % job)],
                        cwd=ana_dir, env=env,
                        log_path=ana_dir / ("%s.log" % job.split("_")[0])):
            print("  !! ERROR: %s fallo, ver %s"
                  % (job, ana_dir / ("%s.log" % job.split("_")[0])))
            row["status"] = "%s_failed" % job.split("_")[0]
            return row

    # 4. FOM
    if not run_step([sys.executable, ANA_DIR / "compute_fom.py"],
                    cwd=ana_dir, env=env, log_path=ana_dir / "fom.log"):
        print("  !! ERROR: compute_fom fallo, ver %s" % (ana_dir / "fom.log"))
        row["status"] = "fom_failed"
        return row

    fom = ana_dir / "fom.txt"
    row["nhits"] = int(fom.read_text().strip())
    print("  FOM (SiPad nhits): %d" % row["nhits"])
    return row


# ── results ─────────────────────────────────────────────────────────────────

def append_results(rows):
    """Append to results.csv and print the same table to stdout."""
    if not rows:
        return
    keys = sorted({k for r in rows for k in r["constants"]})
    header = ["var_id", "name"] + keys + ["nhits_SiPad", "status"]

    def as_record(r):
        rec = {"var_id": r["var_id"], "name": r["name"],
               "nhits_SiPad": r["nhits"], "status": r["status"]}
        rec.update({k: r["constants"].get(k, "") for k in keys})
        return rec

    exists = RESULTS.is_file()
    with open(RESULTS, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        if not exists:
            writer.writeheader()
        for r in rows:
            writer.writerow(as_record(r))

    widths = [max(len(str(h)), *(len(str(as_record(r)[h])) for r in rows))
              for h in header]
    print("\nResultados:")
    print("  " + "  ".join(str(h).ljust(w) for h, w in zip(header, widths)))
    for r in rows:
        rec = as_record(r)
        print("  " + "  ".join(str(rec[h]).ljust(w) for h, w in zip(header, widths)))


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spec", default=str(GEO_DIR / "variants.yaml"),
                    help="variants YAML (default: the one next to the template)")
    ap.add_argument("--only", action="append", default=None, metavar="NAME",
                    help="only this variant (repeatable)")
    ap.add_argument("--dry-run", action="store_true",
                    help="render and validate the geometry, run nothing")
    ap.add_argument("--events", type=int, default=100, help="events per variant")
    ap.add_argument("--skip-sim", action="store_true",
                    help="reuse the existing simulation of each variant instead "
                         "of running ddsim. Does NOT mint new ids: the id is "
                         "looked up by name in the params.yaml of each var dir, "
                         "so the variant must have been run before.")
    args = ap.parse_args()

    if not args.dry_run:
        check_environment()

    base_path, variants = make_variants.load_variants(args.spec)
    _, base_constants = geo_config.load_params(base_path)

    if args.only:
        unknown = sorted(set(args.only) - set(variants))
        if unknown:
            print("Variantes desconocidas: %s" % ", ".join(unknown))
            return 1
        variants = {k: v for k, v in variants.items() if k in args.only}

    print("Spec:  %s" % args.spec)
    print("Base:  %s" % base_path)
    print("Modo:  %s" % ("dry-run" if args.dry_run else
                         ("skip-sim" if args.skip_sim else "completo")))

    rows = []
    for name, overrides in variants.items():
        overrides = overrides or {}
        print("\n=== %s ===" % name)

        unknown = sorted(set(overrides) - set(base_constants))
        if unknown:
            print("  !! ERROR: parametro(s) inexistente(s) en la base: %s"
                  % ", ".join(unknown))
            rows.append({"var_id": "-", "name": name, "constants": overrides,
                         "nhits": "", "status": "unknown_param"})
            continue

        constants = dict(base_constants)
        constants.update(overrides)

        if args.skip_sim:
            var_id = existing_var_id(name)
            if var_id is None:
                print("  !! ERROR: --skip-sim pero no hay ninguna ejecucion "
                      "previa de '%s' en %s" % (name, SIM_OUT))
                rows.append({"var_id": "-", "name": name, "constants": overrides,
                             "nhits": "", "status": "missing_sim"})
                continue
        else:
            # Assigned right before the dir is created, so the next scan sees it.
            var_id = next_var_id()
        print("  id: %s" % var_id)

        try:
            rows.append(run_variant(var_id, name, constants, args))
        except Exception as exc:               # one variant must not kill the rest
            print("  !! ERROR:", exc)
            rows.append({"var_id": var_id, "name": name, "constants": overrides,
                         "nhits": "", "status": "exception"})

    if args.dry_run:
        for r in rows:
            print("  %-20s %s" % (r["name"], r["status"]))
    else:
        append_results(rows)

    failed = [r for r in rows if r["status"] not in ("ok", "dry-run")]
    if failed:
        print("\n%d variante(s) fallaron: %s"
              % (len(failed), ", ".join(r["name"] for r in failed)))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
