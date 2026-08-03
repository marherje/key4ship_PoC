"""GP surrogates. A thin wrapper over BoTorch — no custom kernels, on purpose.

One independent `SingleTaskGP` per objective, collected in a `ModelListGP`.
Independent models are the right default here: the objectives (a simulated hit
count and an analytic cost) have no shared latent structure worth the extra
hyperparameters, and a ModelListGP lets one objective be noisy while another is
noiseless.

Noise is *inferred* by default. With a fixed per-trial event budget the MC noise
on a metric is roughly homoscedastic, which is exactly what the inferred
likelihood assumes; `train_Yvar` is wired through for the day the metrics come
with a per-trial standard error (tracking efficiency and resolution will).
"""

from __future__ import annotations

import torch
from botorch.fit import fit_gpytorch_mll
from botorch.models import ModelListGP, SingleTaskGP
from botorch.models.transforms.input import Normalize
from botorch.models.transforms.outcome import Standardize
from gpytorch.mlls import ExactMarginalLogLikelihood, SumMarginalLogLikelihood

DTYPE = torch.double


def unit_bounds(dim: int, dtype: torch.dtype = DTYPE) -> torch.Tensor:
    """The (2, d) bounds tensor of the unit cube, as optimize_acqf wants it."""
    return torch.stack([torch.zeros(dim, dtype=dtype), torch.ones(dim, dtype=dtype)])


def build_model(
    train_x: torch.Tensor,
    train_y: torch.Tensor,
    train_yvar: torch.Tensor | None = None,
) -> ModelListGP:
    """Fit one GP per objective column of `train_y`.

    `train_x` is already in the unit cube, so `Normalize` is pinned to those
    bounds instead of the data's: the default (normalize by observed min/max)
    would rescale the space every time a new corner point arrives, which moves
    the length scales under the model for no reason.
    """
    if train_x.ndim != 2:
        raise ValueError(f"train_x must be (n, d), got {tuple(train_x.shape)}")
    if train_y.ndim != 2:
        raise ValueError(f"train_y must be (n, m), got {tuple(train_y.shape)}")
    if train_x.shape[0] != train_y.shape[0]:
        raise ValueError("train_x and train_y disagree on the number of points")

    dim = train_x.shape[-1]
    models = []
    for i in range(train_y.shape[-1]):
        yvar_i = train_yvar[:, i : i + 1] if train_yvar is not None else None
        models.append(
            SingleTaskGP(
                train_X=train_x,
                train_Y=train_y[:, i : i + 1],
                train_Yvar=yvar_i,
                input_transform=Normalize(d=dim, bounds=unit_bounds(dim, train_x.dtype)),
                outcome_transform=Standardize(m=1),
            )
        )

    model = ModelListGP(*models)
    mll = (
        SumMarginalLogLikelihood(model.likelihood, model)
        if len(models) > 1
        else ExactMarginalLogLikelihood(models[0].likelihood, models[0])
    )
    fit_gpytorch_mll(mll)
    return model
