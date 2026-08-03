"""Entry points: mobo-run, mobo-resume, mobo-status, mobo-report.

`mobo-run` is a Hydra application — every knob in `conf/` is overridable on the
command line, and the *resolved* config is written into the run directory. That
saved copy is what `mobo-resume` reads, so a resumed run cannot silently differ
from the one it continues (a changed default in `conf/` will not leak into it).

    mobo-run experiment=snd_proxy executor=local evaluator.nevents=10 \\
             experiment.loop.max_trials=4
    mobo-status runs/snd_proxy
    mobo-resume runs/snd_proxy
    mobo-report runs/snd_proxy
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any

from omegaconf import DictConfig, OmegaConf

from .core.loop import AsyncLoop, LoopConfig
from .core.optimizer import MOBOptimizer, OptimizerConfig
from .core.pareto import hypervolume, objective_matrix, pareto_mask, ref_point_tensor
from .core.search_space import SearchSpace
from .core.store import TrialStore
from .core.types import ObjectiveSpec, TrialStatus
from .paths import package_dir

PACKAGE_DIR = package_dir()
CONFIG_NAME = "config.yaml"

log = logging.getLogger("mobo")


def conf_dir() -> Path:
    return Path(os.environ.get("MOBO_CONF_DIR", str(PACKAGE_DIR / "conf")))


def default_runs_dir() -> Path:
    return Path(os.environ.get("MOBO_RUNS_DIR", str(PACKAGE_DIR / "runs")))


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, str(level).upper(), logging.INFO),
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        force=True,
    )


# ── assembling an experiment ─────────────────────────────────────────────────


def objectives_from_config(cfg: Any) -> list[ObjectiveSpec]:
    return [
        ObjectiveSpec(
            name=str(o["name"]),
            direction=str(o.get("direction", "max")),
            ref_point=(None if o.get("ref_point") is None else float(o["ref_point"])),
        )
        for o in cfg
    ]


def ref_point_rule(cfg: Any) -> dict[str, float]:
    """objective name -> multiple of the baseline that is still acceptable."""
    return {
        str(o["name"]): float(o["ref_from_baseline"])
        for o in cfg
        if o.get("ref_from_baseline") is not None
    }


def make_executor(cfg: Any) -> Any:
    kind = str(cfg.get("kind", "local")).lower()
    if kind == "local":
        from .exec.local import LocalExecutor

        return LocalExecutor(
            max_parallel=int(cfg.get("max_parallel", 2)),
            timeout_hours=cfg.get("timeout_hours"),
        )
    if kind == "htcondor":
        from .exec.htcondor import CondorConfig, CondorExecutor

        return CondorExecutor(CondorConfig.from_config(cfg))
    raise ValueError(f"unknown executor kind: {kind!r}")


def build(cfg: DictConfig, run_dir: Path) -> tuple[AsyncLoop, TrialStore]:
    """Turn a resolved config into a ready-to-run loop."""
    from .ship.evaluator import EvaluatorConfig, SNDEvaluator

    experiment = cfg.experiment
    space = SearchSpace.from_config(experiment)
    objectives = objectives_from_config(experiment.objectives)

    optimizer = MOBOptimizer(
        space, objectives, OptimizerConfig.from_config(cfg.optimizer)
    )

    ev_cfg = EvaluatorConfig.from_config(cfg.evaluator)
    ev_cfg.runs_dir = str(run_dir)
    evaluator = SNDEvaluator(ev_cfg)

    executor = make_executor(cfg.executor)

    loop_cfg = LoopConfig.from_config(experiment.loop)
    loop_cfg.ref_point_from_baseline = ref_point_rule(experiment.objectives)
    if "poll_interval" in cfg.executor and cfg.executor.get("poll_interval"):
        loop_cfg.poll_interval = float(cfg.executor.poll_interval)

    store = TrialStore(run_dir / "trials.db")
    loop = AsyncLoop(
        optimizer,
        evaluator,
        executor,
        store,
        loop_cfg,
        global_seed=int(cfg.seed),
        on_progress=_progress_callback(run_dir),
    )
    return loop, store


def _progress_callback(run_dir: Path) -> Any:
    """Refresh the report after every completed trial, cheaply and safely."""

    def callback(loop: AsyncLoop) -> None:
        from .viz.report import write_report

        write_report(run_dir, loop.optimizer.objectives)

    return callback


def record_provenance(store: TrialStore, cfg: DictConfig, evaluator: Any) -> None:
    """Everything needed to explain a number in this store six months from now."""
    import botorch
    import torch

    from . import __version__

    store.set_meta("config", OmegaConf.to_container(cfg, resolve=True))
    store.set_meta("evaluator", evaluator.describe())
    store.set_meta(
        "versions",
        {
            "mobo": __version__,
            "python": sys.version.split()[0],
            "torch": torch.__version__,
            "botorch": botorch.__version__,
        },
    )
    store.set_meta("seed", int(cfg.seed))


# ── mobo-run ─────────────────────────────────────────────────────────────────


def run() -> None:
    import hydra

    @hydra.main(version_base=None, config_path=str(conf_dir()), config_name="config")
    def _run(cfg: DictConfig) -> None:
        setup_logging(cfg.get("log_level", "INFO"))
        run_dir = (
            Path(cfg.run_dir)
            if cfg.get("run_dir")
            else default_runs_dir() / cfg.experiment.name
        )
        run_dir.mkdir(parents=True, exist_ok=True)

        # The resolved config, saved before anything runs: this is what
        # mobo-resume will read, so the run is pinned to what it started with.
        config_path = run_dir / CONFIG_NAME
        if not config_path.exists():
            OmegaConf.save(cfg, config_path, resolve=True)
        else:
            log.info("continuing the existing run in %s", run_dir)

        loop, store = build(cfg, run_dir)
        record_provenance(store, cfg, loop.evaluator)

        problems = loop.evaluator.check_environment()
        if problems:
            sys.exit(
                "the environment is not ready to run trials:\n  - "
                + "\n  - ".join(problems)
                + "\nRun `source init_key4ship.sh` from the repo root first "
                "(key4hep's own setup.sh alone is not enough: the geometry "
                "needs this project's DD4hep plugins from install/)."
            )

        log.info("run dir: %s", run_dir)
        log.info("search space: %s", loop.optimizer.space)
        loop.run()
        _write_report(run_dir)
        print_status(run_dir)

    _run()


# ── mobo-resume ──────────────────────────────────────────────────────────────


def load_run(run_dir: Path) -> DictConfig:
    config_path = Path(run_dir) / CONFIG_NAME
    if not config_path.is_file():
        sys.exit(f"{config_path} not found — is {run_dir} a mobo run directory?")
    return OmegaConf.load(config_path)  # type: ignore[return-value]


def resume() -> None:
    ap = argparse.ArgumentParser(description="Continue an experiment from its store.")
    ap.add_argument("run_dir", help="the run directory (holds trials.db and config.yaml)")
    ap.add_argument(
        "--max-trials",
        type=int,
        default=None,
        help="raise the trial budget while resuming",
    )
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    cfg = load_run(run_dir)
    setup_logging(cfg.get("log_level", "INFO"))
    if args.max_trials is not None:
        cfg.experiment.loop.max_trials = int(args.max_trials)
        OmegaConf.save(cfg, run_dir / CONFIG_NAME, resolve=True)

    loop, _store = build(cfg, run_dir)
    log.info("resuming %s", run_dir)
    loop.run()
    _write_report(run_dir)
    print_status(run_dir)


# ── mobo-status ──────────────────────────────────────────────────────────────


def status() -> None:
    ap = argparse.ArgumentParser(description="Show the state of an experiment.")
    ap.add_argument("run_dir")
    ap.add_argument(
        "--all", action="store_true", help="list every trial, not just the front"
    )
    args = ap.parse_args()
    print_status(Path(args.run_dir).resolve(), verbose=args.all)


def print_status(run_dir: Path, verbose: bool = False) -> None:
    cfg = load_run(run_dir)
    store = TrialStore(run_dir / "trials.db")
    objectives = objectives_from_config(cfg.experiment.objectives)
    stored_ref = store.get_meta("ref_point")
    if stored_ref:
        objectives = [
            ObjectiveSpec(o.name, o.direction, stored_ref.get(o.name, o.ref_point))
            for o in objectives
        ]

    trials = store.all()
    counts = store.counts()
    print(f"\n{run_dir}")
    summary = ", ".join(f"{k.lower()}={v}" for k, v in sorted(counts.items())) or "none"
    print(f"  trials: {len(trials)}  ({summary})")

    y, kept = objective_matrix(trials, objectives)
    if not kept:
        print("  no completed trials yet")
        return

    ref = ref_point_tensor(objectives)
    if ref is not None:
        print(f"  reference point: { ({o.name: o.ref_point for o in objectives}) }")
        print(f"  hypervolume:     {hypervolume(y, ref):.6g}")

    mask = pareto_mask(y).tolist()
    front = [t for t, keep in zip(kept, mask, strict=False) if keep]
    print(f"  pareto front:    {len(front)} point(s)")

    names = [o.name for o in objectives]
    header = ["trial", "tag", "status", *names]
    rows = []
    for t in trials if verbose else front:
        rows.append(
            [
                t.name,
                t.tag or "",
                t.status.value,
                *[_fmt(t.metrics.get(n) if t.metrics else None) for n in names],
            ]
        )
    if rows:
        widths = [
            max(len(str(r[i])) for r in [header, *rows]) for i in range(len(header))
        ]
        print()
        for row in [header, *rows]:
            print(
                "  "
                + "  ".join(str(c).ljust(w) for c, w in zip(row, widths, strict=False))
            )

    failed = [t for t in trials if t.status is TrialStatus.FAILED]
    if failed:
        print("\n  failures:")
        for t in failed[:10]:
            print(f"    {t.name}  {t.error}")
    print()


def _fmt(value: Any) -> str:
    return "-" if value is None else (f"{float(value):.6g}")


# ── mobo-report ──────────────────────────────────────────────────────────────


def report() -> None:
    ap = argparse.ArgumentParser(description="Write the HTML report of an experiment.")
    ap.add_argument("run_dir")
    args = ap.parse_args()
    path = _write_report(Path(args.run_dir).resolve())
    print(path)


def _write_report(run_dir: Path) -> Path:
    from .viz.report import write_report

    cfg = load_run(run_dir)
    objectives = objectives_from_config(cfg.experiment.objectives)
    store = TrialStore(run_dir / "trials.db")
    stored_ref = store.get_meta("ref_point")
    if stored_ref:
        objectives = [
            ObjectiveSpec(o.name, o.direction, stored_ref.get(o.name, o.ref_point))
            for o in objectives
        ]
    return write_report(run_dir, objectives)


if __name__ == "__main__":
    run()
