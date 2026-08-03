"""The one interface that separates the optimizer from the physics.

An `Evaluator` turns a parameter dict into a command to run and, later, a run
directory into metrics. Everything detector-specific — rendering a compact XML,
choosing an event generator, reading an RNTuple — lives behind these three
methods, which is what makes the rest of the package portable to another
detector (implement this, write a config, done).

The split into `validate` / `prepare` / `collect` exists because the three
happen at different times and in different places: `validate` is a cheap
feasibility gate on the driver, `prepare` writes the inputs a worker will read,
and `collect` runs on the driver again after the worker is gone.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from .types import Result, Trial


class Evaluator(ABC):
    """Parameters in, metrics out — with the run in between done elsewhere."""

    @abstractmethod
    def prepare(self, trial: Trial) -> tuple[list[str], Path]:
        """Create the trial's workdir and return the command to run in it.

        The command must be self-contained: an executor may run it on a worker
        node with nothing but the shared filesystem and the experiment's own
        environment setup.
        """

    @abstractmethod
    def collect(self, trial: Trial) -> Result:
        """Read the finished workdir. Never raises: a broken run is a Result.

        Classifying the failure (which step died, and why) is part of the job —
        `Result.error` is what ends up in the store and in the report.
        """

    def validate(self, trial: Trial) -> str | None:
        """Cheap feasibility check, before any CPU is spent. None == feasible.

        The returned string is the reason, recorded on an INFEASIBLE trial.
        """
        del trial
        return None

    def baseline_params(self) -> dict[str, Any] | None:
        """The reference design, if there is one. Trial 0 of a fresh run.

        Anchors the reference point and doubles as an end-to-end sanity check of
        the whole chain: if the baseline fails, nothing else is worth running.
        """
        return None

    def describe(self) -> dict[str, Any]:
        """Free-form provenance recorded in the store (versions, paths, …)."""
        return {}
