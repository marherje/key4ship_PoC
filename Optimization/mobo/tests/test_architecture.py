"""The layering, enforced rather than documented.

Two invariants, both of which are cheap to check statically and expensive to
discover the hard way:

1. `core/`, `exec/` and `viz/` never import `ship/`. That is the whole claim of
   portability: moving this loop to another detector means writing a new
   Evaluator, not editing the optimizer.
2. `ship/payload.py` and `ship/metrics.py` import nothing that a worker node
   might not have. They run under the bare key4hep python — no venv, no torch,
   no mobo package installed — and an accidental `import torch` at the top of
   either would turn every Condor job into a failure that only shows up on the
   cluster.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "mobo"

# What the payload and metrics are allowed to import at module level.
WORKER_ALLOWED = {
    "ROOT",  # from key4hep, and only inside a function
    "yaml",
    "mobo",  # its own package, when loaded by path
}


def module_files(package: str) -> list[Path]:
    return sorted((SRC / package).rglob("*.py"))


def imported_modules(path: Path) -> set[str]:
    """Top-level names of everything imported by a file.

    Relative imports are resolved first, so `from ..ship.geometry import x`
    inside `mobo/exec/` comes back as `mobo.ship.geometry` and can be checked
    like any other name.
    """
    tree = ast.parse(path.read_text())
    package = ["mobo", *path.relative_to(SRC).parts[:-1]]
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package[: len(package) - node.level + 1]
                names.add(".".join([*base, node.module] if node.module else base))
            elif node.module:
                names.add(node.module)
    return names


def mobo_subpackage(name: str) -> str | None:
    """`mobo.core.types` -> "core"; anything else -> None."""
    parts = name.split(".")
    return parts[1] if len(parts) > 1 and parts[0] == "mobo" else None


@pytest.mark.parametrize("package", ["core", "exec", "viz"])
def test_no_ship_imports_outside_the_ship_layer(package):
    offenders = [
        f"{path.relative_to(SRC)} imports {name}"
        for path in module_files(package)
        for name in imported_modules(path)
        if mobo_subpackage(name) == "ship"
    ]
    assert not offenders, "the detector layer leaked into the generic one:\n" + "\n".join(
        offenders
    )


@pytest.mark.parametrize("package", ["core", "exec", "viz"])
def test_generic_layers_import_at_most_the_layers_below_them(package):
    """core -> exec (the Executor ABC) is fine; nothing may reach into ship."""
    allowed = {
        "core": {"core", "exec", "paths", None},
        "exec": {"exec", "core", "paths", None},
        "viz": {"viz", "core", "paths", None},
    }[package]
    for path in module_files(package):
        for name in imported_modules(path):
            if name.startswith("mobo") or name.startswith("."):
                sub = mobo_subpackage(name)
                assert sub in allowed, f"{path.relative_to(SRC)} imports {name}"


def test_the_payload_runs_on_a_bare_key4hep_python():
    """Nothing heavier than the stdlib at import time."""
    for name in ("payload.py", "metrics.py"):
        path = SRC / "ship" / name
        for module in imported_modules(path):
            if module.startswith("mobo"):
                pytest.fail(f"{name} imports {module}; it must be loadable by path alone")
            assert (
                module.split(".")[0] in sys.stdlib_module_names
                or module.split(".")[0] in WORKER_ALLOWED
            ), f"{name} imports {module}, which a worker node may not have"


def test_the_payload_can_be_run_as_a_script():
    """It is executed as `python3 payload.py`, never imported as a module."""
    import subprocess

    proc = subprocess.run(
        [sys.executable, str(SRC / "ship" / "payload.py"), "--help"],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stderr
    assert "--workdir" in proc.stdout


def test_core_does_not_import_hydra_or_the_cli():
    """The ask/tell core has to be usable from a notebook, not just from Hydra."""
    for path in module_files("core"):
        modules = {name.split(".")[0] for name in imported_modules(path)}
        assert "hydra" not in modules, path
