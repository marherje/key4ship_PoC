"""Smoke test of the environment: the venv has what the loop needs."""

import pytest


def test_optimizer_stack_importable():
    import botorch
    import gpytorch
    import torch

    assert torch.__version__
    assert botorch.__version__
    assert gpytorch.__version__


def test_config_and_io_stack_importable():
    import hydra
    import omegaconf
    import pandas
    import yaml

    assert hydra.__version__
    assert omegaconf.__version__
    assert pandas.__version__
    assert yaml.safe_load("a: 1") == {"a": 1}


def test_package_importable():
    import mobo

    assert mobo.__version__


@pytest.mark.key4hep
def test_root_importable():
    """Only in a shell with the key4hep stack sourced (worker nodes, lxplus)."""
    import ROOT

    assert ROOT.gROOT is not None
