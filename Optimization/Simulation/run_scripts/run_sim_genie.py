#!/usr/bin/env python3
"""Run one local GENIE ddsim simulation for a geometry variant.

The GENIE counterpart of run_sim.py. Same contract -- render a steering per
variant, run ddsim in the variant directory, return the edm4hep path -- with
one difference that matters for an optimization run:

THE HEPMC FILE IS AN INPUT, NOT SOMETHING THIS SCRIPT MAKES. Two ways to
feed it:

  * shared (the `snd_genie` experiment): one HepMC of GENIE final states,
    replayed through every variant. Fair between trials, but the vertices
    were sampled over ONE geometry, so a variant with less material is
    credited with interactions it would not have had;
  * per variant (the `snd_genie_gen` experiment): GENIE is run on the
    variant's own GDML first (mobo's payload does geoConverter ->
    shift_gdml_frame.py -> run_gevgen.sh -> genie_gst_to_hepmc.py) and
    --hepmc points at the result. Then --events -1 is the right count: the
    file holds exactly the interactions of that geometry, and asking for
    more than it has would make ddsim hit EOF.

Usage (needs `source init_key4ship.sh` for ddsim):
    python3 run_sim_genie.py --compact <variant.xml> --var-dir <dir> [--events N]
    python3 run_sim_genie.py ... --seed 12345      # reproducible Geant4 run
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

# Produced by simulation/run_script/launch_genie.sh's conversion step.
DEFAULT_HEPMC = RUN_SCRIPT_DIR / "data" / "genie_allflavours_0_2000.hepmc"

OUTPUT_NAME = "output_genie.edm4hep.root"

# Every path rendered into the steering is absolute: ddsim runs with cwd set to
# the variant directory, so a relative path would resolve against the wrong dir.
STEERING_TEMPLATE = '''import sys
sys.path.insert(0, "{run_script_dir}")

from DDSim.DD4hepSimulation import DD4hepSimulation
from genie_base import configure

SIM = DD4hepSimulation()
configure(SIM, dict(
    physicsList    = "QGSP_BERT",
    numberOfEvents = {nevents},
    hepmcFile      = "{hepmc}",
    compactFile    = "{compact}",
    outputFile     = "{output}",{seed_line}
))
'''

SEED_TEMPLATE = '\n    randomSeed     = {seed},'


def write_steering(var_dir, compact, nevents, hepmc, seed=None):
    """Render the steering for one variant -> var_dir/steering.py."""
    var_dir = Path(var_dir).resolve()
    var_dir.mkdir(parents=True, exist_ok=True)
    output = var_dir / OUTPUT_NAME
    steering = var_dir / "steering.py"
    steering.write_text(STEERING_TEMPLATE.format(
        run_script_dir=RUN_SCRIPT_DIR,
        nevents=int(nevents),
        hepmc=Path(hepmc).resolve(),
        compact=Path(compact).resolve(),
        output=output,
        seed_line="" if seed is None else SEED_TEMPLATE.format(seed=int(seed))))
    return steering


def run_sim(var_dir, compact, nevents=100, hepmc=DEFAULT_HEPMC, seed=None):
    """Simulate one variant; returns the edm4hep output path."""
    if shutil.which("ddsim") is None:
        raise RuntimeError("ddsim not found -- source init_key4ship.sh first")

    var_dir = Path(var_dir).resolve()
    compact = Path(compact).resolve()
    if not compact.is_file():
        raise RuntimeError("compact file not found: %s" % compact)
    hepmc = Path(hepmc).resolve()
    if not hepmc.is_file():
        raise RuntimeError(
            "HepMC file not found: %s -- generate it with "
            "simulation/run_script/genie_gst_to_hepmc.py" % hepmc)

    steering = write_steering(var_dir, compact, nevents, hepmc, seed)
    log = var_dir / "ddsim.log"
    with open(log, "w") as fh:
        subprocess.run(["ddsim", "--steeringFile", str(steering)],
                       cwd=str(var_dir), stdout=fh, stderr=subprocess.STDOUT,
                       check=True)
    return var_dir / OUTPUT_NAME


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--compact", required=True, help="variant compact XML")
    ap.add_argument("--var-dir", required=True, help="output dir of this variant")
    ap.add_argument("--events", type=int, default=100,
                    help="events to simulate; -1 = the whole HepMC file")
    ap.add_argument("--hepmc", default=str(DEFAULT_HEPMC),
                    help="GENIE HepMC primaries (shared across variants)")
    ap.add_argument("--seed", type=int, default=None,
                    help="Geant4 random seed; omit for a non-reproducible run")
    args = ap.parse_args()

    try:
        out = run_sim(args.var_dir, args.compact, args.events, args.hepmc, args.seed)
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print("  !! ERROR: %s" % exc)
        return 1
    print("simulation written to %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
