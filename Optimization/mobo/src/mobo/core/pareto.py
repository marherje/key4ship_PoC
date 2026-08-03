"""Pareto bookkeeping shared by the optimizer, the CLI and the report.

Kept apart from `optimizer.py` so that showing the status of a run does not
require building a GP: `mobo-status` and `viz/` read the store and call these.

Everything takes and returns the maximization convention (`ObjectiveSpec.signed`).
"""

from __future__ import annotations

from collections.abc import Sequence

import torch
from botorch.utils.multi_objective.hypervolume import Hypervolume
from botorch.utils.multi_objective.pareto import is_non_dominated

from .models import DTYPE
from .types import ObjectiveSpec, Trial, TrialStatus


def objective_matrix(
    trials: Sequence[Trial], objectives: Sequence[ObjectiveSpec]
) -> tuple[torch.Tensor, list[Trial]]:
    """(n, m) tensor of signed objective values, plus the trials it came from.

    Trials that are not COMPLETED, or whose metrics are missing one of the
    objectives, are dropped — a half-measured point is not an observation.
    """
    rows, kept = [], []
    for t in trials:
        if t.status is not TrialStatus.COMPLETED or not t.metrics:
            continue
        try:
            row = [obj.signed(float(t.metrics[obj.name])) for obj in objectives]
        except (KeyError, TypeError, ValueError):
            continue
        rows.append(row)
        kept.append(t)
    if not rows:
        return torch.empty((0, len(objectives)), dtype=DTYPE), []
    return torch.tensor(rows, dtype=DTYPE), kept


def ref_point_tensor(objectives: Sequence[ObjectiveSpec]) -> torch.Tensor | None:
    """Signed reference point, or None if any objective still lacks one."""
    if any(obj.ref_point is None for obj in objectives):
        return None
    return torch.tensor(
        [obj.signed(float(obj.ref_point)) for obj in objectives],  # type: ignore[arg-type]
        dtype=DTYPE,
    )


def pareto_mask(y: torch.Tensor) -> torch.Tensor:
    """Boolean mask of the non-dominated rows of a signed objective matrix."""
    if y.numel() == 0:
        return torch.zeros(y.shape[0], dtype=torch.bool)
    return is_non_dominated(y)


def hypervolume(y: torch.Tensor, ref_point: torch.Tensor) -> float:
    """Dominated hypervolume w.r.t. `ref_point`; 0.0 if nothing dominates it."""
    if y.numel() == 0:
        return 0.0
    # Points worse than the reference in any objective contribute nothing and
    # make the box decomposition unhappy, so drop them first.
    better = (y > ref_point).all(dim=-1)
    if not better.any():
        return 0.0
    front = y[better][pareto_mask(y[better])]
    return float(Hypervolume(ref_point=ref_point).compute(front))


def hypervolume_trace(y: torch.Tensor, ref_point: torch.Tensor) -> list[float]:
    """Hypervolume after each observation, in the order given.

    Monotonically non-decreasing by construction — a useful invariant to assert
    on, since a drop can only mean the points were reordered or the reference
    moved (both bugs).
    """
    return [hypervolume(y[: i + 1], ref_point) for i in range(y.shape[0])]
