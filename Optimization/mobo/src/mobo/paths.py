"""Where things are on disk.

Its own module because both the detector layer (which needs the geometry and the
analysis chain) and the Condor executor (which needs `init_key4ship.sh`) want
the repository root, and the executor is not allowed to import the detector
layer — that separation is what keeps the loop portable, and it is enforced by
tests/test_architecture.py.
"""

from __future__ import annotations

import os
from pathlib import Path

# .../Optimization/mobo/src/mobo/paths.py
_PACKAGE_DIR = Path(__file__).resolve().parents[2]  # Optimization/mobo
_REPO = _PACKAGE_DIR.parent.parent  # the key4ship checkout


def package_dir() -> Path:
    """Optimization/mobo — where conf/ and runs/ live by default."""
    return _PACKAGE_DIR


def repo_root() -> Path:
    """The key4ship checkout. `MOBO_REPO_ROOT` wins, for unusual installs."""
    return Path(os.environ.get("MOBO_REPO_ROOT", str(_REPO))).resolve()


def init_script() -> Path:
    """The environment setup a worker node has to source before anything else."""
    return repo_root() / "init_key4ship.sh"
