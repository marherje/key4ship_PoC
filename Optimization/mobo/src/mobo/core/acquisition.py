"""Where the next points come from: Sobol for the initial design, qLogNEHVI after.

`qLogNoisyExpectedHypervolumeImprovement` is the state of the art for noisy
multi-objective BO and, more to the point here, it is the one acquisition that
gives us asynchrony for free: the trials currently in flight are passed as
`X_pending`, so the batch it proposes is conditioned on them and never
re-proposes a point that is already running.

Both functions return points in the unit cube.
"""

from __future__ import annotations

from collections.abc import Sequence

import torch
from botorch.acquisition.multi_objective.logei import (
    qLogNoisyExpectedHypervolumeImprovement,
)
from botorch.models.model import Model
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from torch.quasirandom import SobolEngine

from .models import DTYPE, unit_bounds


def sobol_points(dim: int, n: int, seed: int, skip: int = 0) -> torch.Tensor:
    """`n` scrambled Sobol points, deterministic in (seed, skip).

    `skip` is how many points of this sequence were already handed out, so a
    resumed run continues the sequence instead of replaying its head — the
    low-discrepancy property is a property of the *sequence*, and restarting it
    would put duplicate points in the design.
    """
    engine = SobolEngine(dimension=dim, scramble=True, seed=seed)
    if skip:
        engine.fast_forward(skip)
    return engine.draw(n, dtype=DTYPE)


def propose_qlognehvi(
    model: Model,
    ref_point: Sequence[float],
    x_baseline: torch.Tensor,
    q: int,
    x_pending: torch.Tensor | None = None,
    num_restarts: int = 10,
    raw_samples: int = 512,
    mc_samples: int = 128,
    sequential: bool = True,
    prune_baseline: bool = True,
    seed: int | None = None,
) -> torch.Tensor:
    """Maximize qLogNEHVI over the unit cube; returns a (q, d) tensor.

    Everything here is in the maximization convention (see `types.py`): the
    caller has already flipped the sign of any minimized objective, ref point
    included.
    """
    if seed is not None:
        torch.manual_seed(int(seed))

    sampler = SobolQMCNormalSampler(sample_shape=torch.Size([mc_samples]), seed=seed)
    acqf = qLogNoisyExpectedHypervolumeImprovement(
        model=model,
        ref_point=list(float(r) for r in ref_point),
        X_baseline=x_baseline,
        sampler=sampler,
        prune_baseline=prune_baseline,
        X_pending=x_pending,
    )
    candidates, _ = optimize_acqf(
        acq_function=acqf,
        bounds=unit_bounds(x_baseline.shape[-1], x_baseline.dtype),
        q=q,
        num_restarts=num_restarts,
        raw_samples=raw_samples,
        # Greedy sequential selection: cheaper than a joint q-batch optimization
        # and, for q of a few units, essentially as good.
        sequential=sequential,
        options={"batch_limit": 5, "maxiter": 200},
    )
    return candidates.detach()
