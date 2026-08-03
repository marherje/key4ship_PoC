#!/usr/bin/env python3
"""Stand-in for the real payload: BraninCurrin instead of a full simulation.

Same contract as `mobo.ship.payload` — read `trial.json` from the workdir, do
the work, write `metrics.json` and then a `DONE`/`FAILED` sentinel — so the loop
and the executors get exercised end to end without a Geant4 run. Stdlib only,
and no import of the mobo package, exactly like the real one.
"""

import argparse
import json
import math
import random
import sys
import time
from pathlib import Path


def branin(x0, x1):
    """BoTorch's BraninCurrin first objective, on the unit square."""
    a = 15.0 * x0 - 5.0
    b = 15.0 * x1
    t1 = b - 5.1 / (4 * math.pi**2) * a**2 + 5.0 / math.pi * a - 6.0
    t2 = 10.0 * (1.0 - 1.0 / (8.0 * math.pi)) * math.cos(a)
    return t1**2 + t2 + 10.0


def currin(x0, x1):
    """Second objective; x1 = 0 is a removable singularity, pinned as BoTorch does."""
    factor = 1.0 if x1 == 0 else 1.0 - math.exp(-1.0 / (2.0 * x1))
    num = 2300 * x0**3 + 1900 * x0**2 + 2092 * x0 + 60
    den = 100 * x0**3 + 500 * x0**2 + 4 * x0 + 20
    return factor * num / den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    started = time.time()
    try:
        spec = json.loads((workdir / "trial.json").read_text())
        x0, x1 = spec["unit_x"]
        rng = random.Random(spec["seed"])

        if rng.random() < spec.get("failure_rate", 0.0):
            raise RuntimeError("injected failure")
        if spec.get("sleep"):
            time.sleep(float(spec["sleep"]))

        noise = spec.get("noise", 0.0)
        metrics = {
            # negated: the loop maximizes, BraninCurrin is a minimization problem
            "f1": -(branin(x0, x1) + (rng.gauss(0, noise) if noise else 0.0)),
            "f2": -(currin(x0, x1) + (rng.gauss(0, noise) if noise else 0.0)),
            "wall_time_s": time.time() - started,
        }
        (workdir / "metrics.json").write_text(json.dumps(metrics, indent=2))
        (workdir / "DONE").write_text("ok\n")
    except Exception as exc:  # noqa: BLE001 - the sentinel is the error channel
        (workdir / "FAILED").write_text(f"{type(exc).__name__}: {exc}\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
