"""The SND evaluator: search-space parameters in, a runnable trial out.

This is the whole of the detector-specific surface the loop sees. It renders the
geometry, refuses the ones that would not build, writes the job's inputs, and
afterwards reads back what the payload measured.

Each trial owns exactly one directory, which holds everything about it: the
compact XML that was simulated, the parameters that produced it, the payload's
inputs and logs, the ROOT files and `metrics.json`. A trial is therefore
reproducible from its directory alone, without the store and without the driver.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from ..core.evaluator import Evaluator
from ..core.types import Result, Trial
from .geometry import Geometry, repo_root
from .metrics import CostModel

SIM_FILE = "output_mu_pi.edm4hep.root"
HITS_FILE = "ShipHits.root"
COMPACT_FILE = "geometry.xml"
TRIAL_SPEC = "trial.json"


@dataclass
class EvaluatorConfig:
    """Everything about *how* a trial is evaluated, none of it about *which*."""

    nevents: int = 100
    output_level: str = "INFO"
    python: str = "python3"
    # Where trial directories go. The CLI fills this in from the run directory.
    runs_dir: str | None = None
    cost_model: dict[str, float] = field(default_factory=dict)
    # Paths into the existing project. Defaults are derived from the repo root;
    # they exist in the config so an experiment can point at a modified chain
    # without touching code.
    geometry_dir: str | None = None
    analysis_dir: str | None = None
    run_sim: str | None = None
    parse_geometry: str | None = None

    @classmethod
    def from_config(cls, cfg: Any) -> EvaluatorConfig:
        fields = set(cls.__dataclass_fields__)
        kwargs = {k: v for k, v in dict(cfg).items() if k in fields}
        if "cost_model" in kwargs and kwargs["cost_model"]:
            kwargs["cost_model"] = {
                k: float(v) for k, v in dict(kwargs["cost_model"]).items()
            }
        return cls(**kwargs)


class SNDEvaluator(Evaluator):
    def __init__(self, cfg: EvaluatorConfig | None = None) -> None:
        self.cfg = cfg or EvaluatorConfig()
        repo = repo_root()
        self.analysis_dir = Path(
            self.cfg.analysis_dir or repo / "Optimization" / "Analysis"
        )
        self.run_sim = Path(
            self.cfg.run_sim
            or repo / "Optimization" / "Simulation" / "run_scripts" / "run_sim.py"
        )
        self.parse_geometry = Path(
            self.cfg.parse_geometry
            or repo / "simulation" / "geometry" / "parse_geometry.py"
        )
        self.geometry = Geometry(geo_dir=self.cfg.geometry_dir)
        self.cost_model = CostModel.from_config(self.cfg.cost_model)
        # Absolute, always: the payload runs with its own workdir as cwd, so a
        # relative path here would be resolved against itself on the far side.
        runs_dir = Path(
            self.cfg.runs_dir or (repo / "Optimization" / "mobo" / "runs" / "adhoc")
        )
        self.trials_dir = runs_dir.resolve() / "trials"

    # ── the Evaluator contract ───────────────────────────────────────────────

    def workdir(self, trial: Trial) -> Path:
        return self.trials_dir / trial.name

    def validate(self, trial: Trial) -> str | None:
        """Dry-run the real renderer. No CPU is spent on a geometry that fails."""
        return self.geometry.check(trial.params)

    def prepare(self, trial: Trial) -> tuple[list[str], Path]:
        workdir = self.workdir(trial)
        workdir.mkdir(parents=True, exist_ok=True)

        constants = self.geometry.constants_for(trial.params)
        compact = workdir / COMPACT_FILE
        ok, lines, errors = self.geometry.write(constants, compact)
        if not ok:
            # validate() should have caught this; if it did not, the two
            # disagree and that is worth an exception rather than a silent run.
            raise RuntimeError("geometry rejected at prepare(): " + "; ".join(errors))

        # A human-readable snapshot next to the machine-readable one: this is
        # what you read when a trial six weeks old looks wrong.
        (workdir / "params.yaml").write_text(
            yaml.safe_dump(
                {
                    "trial_id": trial.trial_id,
                    "name": trial.name,
                    "tag": trial.tag,
                    "seed": trial.seed,
                    "params": _plain(trial.params),
                    "constants": _plain(constants),
                    "geometry_summary": lines,
                },
                sort_keys=True,
                default_flow_style=False,
            )
        )

        (workdir / TRIAL_SPEC).write_text(
            json.dumps(
                {
                    "trial_id": trial.trial_id,
                    "name": trial.name,
                    "seed": trial.seed,
                    "nevents": int(self.cfg.nevents),
                    "compact": str(compact),
                    "sim_file": SIM_FILE,
                    "hits_file": HITS_FILE,
                    "analysis_dir": str(self.analysis_dir),
                    "run_sim": str(self.run_sim),
                    "parse_geometry": str(self.parse_geometry),
                    "output_level": self.cfg.output_level,
                    "cost_model": {
                        "si_per_m2": self.cost_model.si_per_m2,
                        "w_per_kg": self.cost_model.w_per_kg,
                        "rho_w_g_cm3": self.cost_model.rho_w_g_cm3,
                    },
                },
                indent=2,
                sort_keys=True,
            )
        )

        payload = Path(__file__).resolve().parent / "payload.py"
        return [self.cfg.python, str(payload), "--workdir", str(workdir)], workdir

    def collect(self, trial: Trial) -> Result:
        """Read the finished workdir. Never raises; a broken run is a Result."""
        workdir = Path(trial.workdir) if trial.workdir else self.workdir(trial)
        metrics_path = workdir / "metrics.json"
        if metrics_path.is_file():
            try:
                metrics = json.loads(metrics_path.read_text())
            except (OSError, ValueError) as exc:
                return Result.failed(f"unreadable metrics.json: {exc}")
            if "nhits_sipad" in metrics:
                return Result.ok({k: float(v) for k, v in metrics.items()})
            return Result.failed("metrics.json has no nhits_sipad", metrics)

        failed = workdir / "FAILED"
        if failed.is_file():
            text = failed.read_text().strip().splitlines()
            return Result.failed(text[0] if text else "the payload reported a failure")
        return Result.failed("no metrics.json and no FAILED sentinel")

    def baseline_params(self) -> dict[str, Any]:
        """The committed design, in search-space coordinates. Trial 0."""
        return self.geometry.baseline_params()

    def describe(self) -> dict[str, Any]:
        return {
            "nevents": self.cfg.nevents,
            "output_level": self.cfg.output_level,
            "geometry_dir": str(self.geometry.geo_dir),
            "analysis_dir": str(self.analysis_dir),
            "baseline_params_file": str(self.geometry.params_file),
            "cost_model": {
                "si_per_m2": self.cost_model.si_per_m2,
                "w_per_kg": self.cost_model.w_per_kg,
                "rho_w_g_cm3": self.cost_model.rho_w_g_cm3,
            },
            "git_sha": _git_sha(),
        }

    # ── environment ──────────────────────────────────────────────────────────

    def check_environment(self) -> list[str]:
        """What would stop a payload from running here. Empty means fine.

        Sourcing key4hep alone is not enough and fails in a way that costs an
        hour to diagnose: the geometry uses this project's own DD4hep plugins
        (the SiPad/SiTarget/MTC builders and the CartesianStripXStereo
        segmentation), which live in `install/` and only reach dd4hep through
        the LD_LIBRARY_PATH that `init_key4ship.sh` prepends. Without it ddsim
        dies with "FAILED to create segmentation ... [Missing factory]".
        """
        problems = [
            f"{tool} is not on PATH"
            for tool in ("ddsim", "k4run")
            if shutil.which(tool) is None
        ]
        install = repo_root() / "install"
        if install.is_dir():
            search = os.environ.get("LD_LIBRARY_PATH", "").split(":")
            if not any(entry.startswith(str(install)) for entry in search if entry):
                problems.append(
                    f"{install} is not on LD_LIBRARY_PATH (the project's DD4hep "
                    "plugins would not load)"
                )
        return problems


def _plain(value: Any) -> Any:
    """OmegaConf/numpy scalars -> something yaml.safe_dump will accept."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, bool | str) or value is None:
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return float(value)
    return str(value)


def _git_sha() -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_root()), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None
