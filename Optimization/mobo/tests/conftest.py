"""Shared fixtures: a synthetic search space and a synthetic evaluator.

The synthetic problem is BoTorch's BraninCurrin (2 parameters, 2 objectives,
known Pareto front), taken with `negate=True` so that both objectives are
maximized — the convention the optimizer works in. It is the reference the core
is validated against, so that a regression in the optimizer shows up without
running a single simulation.
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import pytest

from mobo.core.evaluator import Evaluator
from mobo.core.search_space import SearchSpace
from mobo.core.types import ObjectiveSpec, ParameterSpec, Result


@pytest.fixture
def space2d() -> SearchSpace:
    return SearchSpace([ParameterSpec("x0", 0.0, 1.0), ParameterSpec("x1", 0.0, 1.0)])


@pytest.fixture
def mixed_space() -> SearchSpace:
    """Floats, an integer and a log dimension, plus a fixed parameter."""
    return SearchSpace(
        [
            ParameterSpec("thickness", 5.0, 15.0),
            ParameterSpec("nlayers", 4, 40, kind="int"),
            ParameterSpec("tol", 1e-4, 1e-1, log=True),
        ],
        fixed={"frame_gap": 0.1, "mode": "auto"},
    )


@pytest.fixture
def objectives() -> list[ObjectiveSpec]:
    """BraninCurrin's own reference point, in the maximization convention."""
    return [
        ObjectiveSpec("f1", "max", -18.0),
        ObjectiveSpec("f2", "max", -6.0),
    ]


class BraninCurrinEvaluator:
    """`metrics = evaluate(params)`, optionally noisy and occasionally failing."""

    def __init__(self, noise: float = 0.0, failure_rate: float = 0.0, seed: int = 0):
        from botorch.test_functions.multi_objective import BraninCurrin

        self.problem = BraninCurrin(negate=True).to(dtype=__import__("torch").double)
        self.noise = noise
        self.failure_rate = failure_rate
        self.rng = random.Random(seed)
        self.n_calls = 0
        self.n_failures = 0

    def __call__(self, unit_x) -> dict[str, float] | None:
        import torch

        self.n_calls += 1
        if self.failure_rate and self.rng.random() < self.failure_rate:
            self.n_failures += 1
            return None
        x = torch.tensor([list(map(float, unit_x))], dtype=torch.double)
        y = self.problem(x).squeeze(0).tolist()
        if self.noise:
            y = [v + self.rng.gauss(0.0, self.noise) for v in y]
        return {"f1": y[0], "f2": y[1]}


@pytest.fixture
def braninc() -> type[BraninCurrinEvaluator]:
    return BraninCurrinEvaluator


@pytest.fixture
def store(tmp_path):
    from mobo.core.store import TrialStore

    s = TrialStore(tmp_path / "trials.db")
    yield s
    s.close()


# ── a full Evaluator over the synthetic problem ──────────────────────────────

SYNTHETIC_PAYLOAD = Path(__file__).resolve().parent / "synthetic_payload.py"


class SyntheticEvaluator(Evaluator):
    """Drives `synthetic_payload.py` — a real subprocess, no physics.

    Lets the loop and both executors be tested for real (workdirs, sentinels,
    retries, resume) in seconds instead of hours.
    """

    def __init__(
        self,
        root: Path,
        failure_rate: float = 0.0,
        noise: float = 0.0,
        sleep: float = 0.0,
        infeasible_below: float | None = None,
    ):
        self.root = Path(root)
        self.failure_rate = failure_rate
        self.noise = noise
        self.sleep = sleep
        self.infeasible_below = infeasible_below

    def workdir(self, trial) -> Path:
        return self.root / "trials" / trial.name

    def validate(self, trial):
        if self.infeasible_below is not None and trial.unit_x[0] < self.infeasible_below:
            return f"x0 = {trial.unit_x[0]:.3f} is outside the feasible region"
        return None

    def prepare(self, trial):
        workdir = self.workdir(trial)
        workdir.mkdir(parents=True, exist_ok=True)
        (workdir / "trial.json").write_text(
            json.dumps(
                {
                    "unit_x": list(trial.unit_x),
                    "seed": trial.seed,
                    "failure_rate": self.failure_rate,
                    "noise": self.noise,
                    "sleep": self.sleep,
                }
            )
        )
        return [
            sys.executable,
            str(SYNTHETIC_PAYLOAD),
            "--workdir",
            str(workdir),
        ], workdir

    def collect(self, trial):
        workdir = self.workdir(trial)
        metrics_path = workdir / "metrics.json"
        if not metrics_path.is_file():
            failed = workdir / "FAILED"
            reason = failed.read_text().strip() if failed.is_file() else "no metrics.json"
            return Result.failed(reason)
        return Result.ok(json.loads(metrics_path.read_text()))

    def baseline_params(self):
        return {"x0": 0.5, "x1": 0.5}


@pytest.fixture
def synthetic_evaluator(tmp_path):
    return SyntheticEvaluator(tmp_path)


def pytest_addoption(parser):
    """`--update-golden` rewrites the recorded job.sub / job.sh fixtures.

    Golden files are there so a change to what we submit to the cluster has to
    be deliberate; regenerating them is a one-flag, reviewable diff.
    """
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="rewrite the golden files instead of comparing against them",
    )
