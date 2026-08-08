"""The real payload, on the real baseline geometry, with 10 events.

Slow (a couple of minutes) and only runnable inside the key4hep environment, but
it is the only test that exercises ddsim, the three Gaudi jobs and the RNTuple
readback together — everything else in `ship/` is tested against the renderer
and against arithmetic.

Two claims:

* the chain runs and produces a `metrics.json` whose analytic half agrees with
  the numbers computed independently in `test_metrics.py`;
* the same seed gives the same physics. Without that the optimizer is fitting a
  GP to Geant4's random number generator.

    pytest tests/ship/test_payload_smoke.py -m "slow and key4hep"
"""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

import pytest

from mobo.core.types import Trial, TrialStatus
from mobo.ship.evaluator import EvaluatorConfig, SNDEvaluator

pytestmark = [pytest.mark.slow, pytest.mark.key4hep]

NEVENTS = 10
SEED = 12345


def run_payload(root: Path, seed: int, nevents: int = NEVENTS, trial_id: int = 0):
    """One complete trial, start to finish. Returns (evaluator, trial, result)."""
    evaluator = SNDEvaluator(
        EvaluatorConfig(nevents=nevents, runs_dir=str(root), output_level="INFO")
    )
    problems = evaluator.check_environment()
    if problems:
        pytest.skip("environment not ready: " + "; ".join(problems))

    trial = Trial(
        trial_id=trial_id,
        params=evaluator.baseline_params(),
        unit_x=[0.0],
        seed=seed,
        tag="baseline",
    )
    assert evaluator.validate(trial) is None, "the baseline must always be feasible"

    cmd, workdir = evaluator.prepare(trial)
    started = time.time()
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=3600)
    trial.workdir = str(workdir)
    result = evaluator.collect(trial)
    elapsed = time.time() - started
    print(f"payload rc={proc.returncode} in {elapsed:.0f} s -> {result.status}")
    return evaluator, trial, result


@pytest.fixture(scope="module")
def baseline_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("payload_smoke")
    return run_payload(root, SEED)


def test_the_chain_completes(baseline_run):
    _evaluator, trial, result = baseline_run
    assert result.status is TrialStatus.COMPLETED, result.error
    workdir = Path(trial.workdir)
    assert (workdir / "DONE").is_file()
    assert not (workdir / "FAILED").exists()


def test_every_step_left_its_artefacts(baseline_run):
    """One directory holds the whole trial: inputs, logs and outputs."""
    _evaluator, trial, _result = baseline_run
    workdir = Path(trial.workdir)
    for name in (
        "geometry.xml",  # what was simulated
        "params.yaml",  # and why
        "trial.json",
        "steering.py",
        "output_mu_pi.edm4hep.root",  # ddsim
        "events.edm4hep.root",  # job1
        "tracks.edm4hep.root",  # job4
        "ShipHits.root",  # job5
        "metrics.json",
    ):
        assert (workdir / name).is_file(), f"{name} is missing"


def test_the_metrics_are_complete_and_sane(baseline_run):
    _evaluator, trial, result = baseline_run
    metrics = result.metrics
    for key in (
        "nhits_sipad",
        "nhits_sitarget",
        "cost_proxy",
        "n_channels",
        "si_area_m2",
        "w_mass_kg",
        "cpu_time_s",
        "wall_time_s",
    ):
        assert key in metrics, f"{key} was not recorded"

    assert metrics["nhits_sipad"] > 0, (
        "no SiPad hits at all: the chain ran but saw nothing"
    )
    assert metrics["nevents"] == NEVENTS
    assert 0 < metrics["cpu_time_s"] <= metrics["wall_time_s"] * 8  # 8 cores of slack


def test_the_analytic_metrics_match_the_hand_computed_baseline(baseline_run):
    """Cross-check against test_metrics.py, through a completely different path.

    There the numbers come from a hand-written constants table; here they come
    from the compact file that ddsim actually loaded.
    """
    _evaluator, _trial, result = baseline_run
    assert result.metrics["cost_proxy"] == pytest.approx(1413.374, abs=1e-3)
    assert result.metrics["si_area_m2"] == pytest.approx(40.913312, abs=1e-5)
    assert result.metrics["w_mass_kg"] == pytest.approx(1859.748, abs=1e-3)
    assert result.metrics["n_channels"] == 3889920
    assert result.metrics["sipad_nlayers"] == 10


def test_the_same_seed_gives_the_same_physics(baseline_run, tmp_path_factory):
    """Reproducibility, the reason run_sim.py grew a --seed."""
    _evaluator, _trial, first = baseline_run
    _e2, _t2, again = run_payload(tmp_path_factory.mktemp("same_seed"), SEED)
    assert again.status is TrialStatus.COMPLETED, again.error
    assert again.metrics["nhits_sipad"] == first.metrics["nhits_sipad"]
    assert again.metrics["nhits_sitarget"] == first.metrics["nhits_sitarget"]


def test_a_different_seed_gives_different_physics(baseline_run, tmp_path_factory):
    """The other half: the seed must actually reach Geant4."""
    _evaluator, _trial, first = baseline_run
    _e2, _t2, other = run_payload(tmp_path_factory.mktemp("other_seed"), SEED + 1)
    assert other.status is TrialStatus.COMPLETED, other.error
    assert other.metrics["nhits_sipad"] != first.metrics["nhits_sipad"]
    # ... but nothing about the geometry may move with it.
    assert other.metrics["cost_proxy"] == pytest.approx(first.metrics["cost_proxy"])


def test_the_steering_carries_the_seed(baseline_run):
    _evaluator, trial, _result = baseline_run
    steering = (Path(trial.workdir) / "steering.py").read_text()
    assert f"randomSeed     = {SEED}" in steering


def test_metrics_json_is_readable_without_the_driver(baseline_run):
    """A trial directory has to be self-explanatory six months later."""
    _evaluator, trial, result = baseline_run
    on_disk = json.loads((Path(trial.workdir) / "metrics.json").read_text())
    assert on_disk["nhits_sipad"] == result.metrics["nhits_sipad"]
