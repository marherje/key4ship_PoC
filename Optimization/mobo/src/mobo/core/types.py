"""The vocabulary the whole package shares: parameters, objectives, trials.

Everything here is a plain dataclass or enum, JSON-serializable, with no
dependency on torch or on any detector. `store.py` persists these, `optimizer.py`
turns them into tensors, `exec/` moves them around and `ship/` fills them in.

Sign convention, stated once for the whole package: **the optimizer maximizes**.
An objective declared `direction="min"` is negated at the boundary
(`ObjectiveSpec.signed`), and only there — nothing downstream of that call has
to remember which way a given metric points.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

# ── parameters ───────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ParameterSpec:
    """One free dimension of the search space.

    `kind="int"` is optimized as a continuous relaxation and rounded on the way
    back to physical units (standard BoTorch practice: the GP sees a smooth
    space, the evaluator sees an integer).
    """

    name: str
    low: float
    high: float
    kind: str = "float"  # "float" | "int"
    log: bool = False

    def __post_init__(self) -> None:
        if self.kind not in ("float", "int"):
            raise ValueError(
                f"{self.name}: kind must be 'float' or 'int', not {self.kind!r}"
            )
        if not self.high > self.low:
            raise ValueError(
                f"{self.name}: need low < high, got [{self.low}, {self.high}]"
            )
        if self.log and self.low <= 0:
            raise ValueError(f"{self.name}: log scale needs low > 0, got {self.low}")
        if self.kind == "int" and (
            self.low != int(self.low) or self.high != int(self.high)
        ):
            raise ValueError(f"{self.name}: int bounds must be whole numbers")


@dataclass(frozen=True)
class ObjectiveSpec:
    """One optimized metric.

    `ref_point` is the worst value still considered useful, in *physical* units
    and physical direction; it is the hypervolume reference. It may be None at
    construction time and filled in later from the baseline trial (see
    `loop.py`), which is why this dataclass does not require it.
    """

    name: str
    direction: str = "max"  # "max" | "min"
    ref_point: float | None = None

    def __post_init__(self) -> None:
        if self.direction not in ("max", "min"):
            raise ValueError(
                f"{self.name}: direction must be 'max' or 'min', not {self.direction!r}"
            )

    @property
    def sign(self) -> float:
        """+1 for a maximized objective, -1 for a minimized one."""
        return 1.0 if self.direction == "max" else -1.0

    def signed(self, value: float) -> float:
        """Physical value -> the maximization convention used internally."""
        return self.sign * float(value)

    def unsigned(self, value: float) -> float:
        """Internal value -> physical units. Inverse of `signed`."""
        return self.sign * float(value)


# ── trials ───────────────────────────────────────────────────────────────────


class TrialStatus(str, Enum):
    """PROPOSED -> SUBMITTED -> RUNNING -> COMPLETED | FAILED | INFEASIBLE.

    INFEASIBLE is reserved for a geometry rejected *before* any CPU is spent on
    it (the dry-run gate); FAILED is a trial that ran and did not produce
    metrics. Both are excluded from the GP, but only FAILED is worth retrying.
    """

    PROPOSED = "PROPOSED"
    SUBMITTED = "SUBMITTED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    INFEASIBLE = "INFEASIBLE"

    @property
    def is_terminal(self) -> bool:
        return self in (TrialStatus.COMPLETED, TrialStatus.FAILED, TrialStatus.INFEASIBLE)

    @property
    def is_in_flight(self) -> bool:
        return self in (TrialStatus.SUBMITTED, TrialStatus.RUNNING)


@dataclass
class Result:
    """What an evaluation produced. `metrics` is always the *full* dict.

    Which subset of it gets optimized is a config decision (`ObjectiveSpec`),
    never a property of the evaluator: every run records everything it measured,
    so an objective can be added later without re-running anything.
    """

    status: TrialStatus
    metrics: dict[str, float] = field(default_factory=dict)
    error: str | None = None

    @classmethod
    def ok(cls, metrics: dict[str, float]) -> Result:
        return cls(status=TrialStatus.COMPLETED, metrics=dict(metrics))

    @classmethod
    def failed(cls, error: str, metrics: dict[str, float] | None = None) -> Result:
        return cls(status=TrialStatus.FAILED, metrics=dict(metrics or {}), error=error)

    @classmethod
    def infeasible(cls, error: str) -> Result:
        return cls(status=TrialStatus.INFEASIBLE, error=error)


@dataclass
class Trial:
    """One point of the search space and everything that happened to it."""

    trial_id: int
    params: dict[str, Any]
    unit_x: list[float]
    seed: int
    status: TrialStatus = TrialStatus.PROPOSED
    metrics: dict[str, float] | None = None
    error: str | None = None
    workdir: str | None = None
    job_id: str | None = None
    attempt: int = 0
    tag: str | None = None  # e.g. "baseline", "sobol", "qlognehvi"
    created_at: float | None = None
    submitted_at: float | None = None
    finished_at: float | None = None

    @property
    def name(self) -> str:
        return f"t{self.trial_id:05d}"

    def apply(self, result: Result) -> None:
        """Fold an evaluation result into the trial."""
        self.status = result.status
        self.metrics = dict(result.metrics) if result.metrics else None
        self.error = result.error


def derive_seed(global_seed: int, trial_id: int) -> int:
    """Deterministic per-trial seed.

    A hash rather than `global_seed + trial_id` so that two experiments whose
    global seeds differ by one do not share most of their per-trial seeds.
    Truncated to 31 bits because that is what Geant4/ddsim accepts.
    """
    digest = hashlib.sha256(f"{int(global_seed)}:{int(trial_id)}".encode()).digest()
    return int.from_bytes(digest[:4], "big") & 0x7FFFFFFF
