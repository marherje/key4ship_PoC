#!/usr/bin/env python3
"""Run one local mu_pi ddsim simulation for a geometry variant.

The steering is generated per variant rather than kept as a file, because the
only things that change between variants are the compact XML and the output
path. Everything else reproduces the committed mu_pi sample
(simulation/run_script/steer/QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0.py),
so a variant differs from the baseline sample only by its geometry.

Usage (needs `source init_key4ship.sh` for ddsim):
    python3 run_sim.py --compact <variant.xml> --var-dir <dir> [--events N]
    python3 run_sim.py ... --seed 12345      # reproducible Geant4 run
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OPT = HERE.parent.parent                       # Optimization/
REPO = OPT.parent
RUN_SCRIPT_DIR = REPO / "simulation" / "run_script"
OUTPUT_ROOT = OPT / "Simulation" / "output" / "variants"

# Every path rendered into the steering is absolute: ddsim runs with cwd set to
# the variant directory, so a relative path would resolve against the wrong dir.
STEERING_TEMPLATE = '''import sys
sys.path.insert(0, "{run_script_dir}")

from DDSim.DD4hepSimulation import DD4hepSimulation
from multiplePG_base import configure

SIM = DD4hepSimulation()
configure(SIM, dict(
    physicsList    = "QGSP_BERT",
    numberOfEvents = {nevents},
    particles      = [("mu+", 10.0), ("pi+", 5.0)],
    angleDeg       = 1.0,
    position       = (1, 1, -1000),
    compactFile    = "{compact}",
    outputFile     = "{output}",{seed_line}
))
'''

# multiplePG_base.configure() maps randomSeed -> SIM.random.seed. Omitting the
# key entirely (rather than passing None) keeps the steering byte-identical to
# what it was before seeding existed, so an unseeded run is unchanged.
SEED_TEMPLATE = '\n    randomSeed     = {seed},'


def write_steering(var_dir, compact, nevents, seed=None):
    """Render the steering for one variant -> var_dir/steering.py."""
    var_dir = Path(var_dir).resolve()
    var_dir.mkdir(parents=True, exist_ok=True)
    output = var_dir / "output_mu_pi.edm4hep.root"
    steering = var_dir / "steering.py"
    steering.write_text(STEERING_TEMPLATE.format(
        run_script_dir=RUN_SCRIPT_DIR,
        nevents=int(nevents),
        compact=Path(compact).resolve(),
        output=output,
        seed_line="" if seed is None else SEED_TEMPLATE.format(seed=int(seed))))
    return steering


def run_sim(var_dir, compact, nevents=100, seed=None):
    """Simulate one variant; returns the edm4hep output path.

    Without `seed` Geant4 picks its own, which is fine for a one-off variant but
    not for an optimization run: the same geometry has to give the same answer
    twice, or the optimizer is fitting a model to the random number generator.
    """
    if shutil.which("ddsim") is None:
        raise RuntimeError("ddsim not found — source init_key4ship.sh first")

    var_dir = Path(var_dir).resolve()
    compact = Path(compact).resolve()
    if not compact.is_file():
        raise RuntimeError("compact file not found: %s" % compact)

    steering = write_steering(var_dir, compact, nevents, seed)
    log = var_dir / "ddsim.log"
    with open(log, "w") as fh:
        subprocess.run(["ddsim", "--steeringFile", str(steering)],
                       cwd=str(var_dir), stdout=fh, stderr=subprocess.STDOUT,
                       check=True)
    return var_dir / "output_mu_pi.edm4hep.root"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--compact", required=True, help="variant compact XML")
    ap.add_argument("--var-dir", required=True, help="output dir of this variant")
    ap.add_argument("--events", type=int, default=100, help="events to simulate")
    ap.add_argument("--seed", type=int, default=None,
                    help="Geant4 random seed (default: ddsim's own, i.e. not "
                         "reproducible)")
    args = ap.parse_args()

    try:
        out = run_sim(args.var_dir, args.compact, args.events, args.seed)
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print("  !! ERROR:", exc)
        return 1
    print("Simulado %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
