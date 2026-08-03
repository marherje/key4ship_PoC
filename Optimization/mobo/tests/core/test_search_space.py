"""The unit-cube map has to be exactly invertible where it claims to be."""

from __future__ import annotations

import math
import random

import pytest

from mobo.core.search_space import SearchSpace
from mobo.core.types import ParameterSpec


def test_shape_and_names(mixed_space):
    assert mixed_space.dim == 3
    assert mixed_space.names == ["thickness", "nlayers", "tol"]


def test_fixed_parameters_are_passed_through(mixed_space):
    params = mixed_space.to_params([0.5, 0.5, 0.5])
    assert params["frame_gap"] == 0.1
    assert params["mode"] == "auto"
    assert set(params) == {"thickness", "nlayers", "tol", "frame_gap", "mode"}


def test_corners_hit_the_bounds(mixed_space):
    lo = mixed_space.to_params([0.0, 0.0, 0.0])
    hi = mixed_space.to_params([1.0, 1.0, 1.0])
    assert lo["thickness"] == pytest.approx(5.0)
    assert hi["thickness"] == pytest.approx(15.0)
    assert lo["nlayers"] == 4
    assert hi["nlayers"] == 40
    assert lo["tol"] == pytest.approx(1e-4)
    assert hi["tol"] == pytest.approx(1e-1)


def test_roundtrip_params_unit_params(mixed_space):
    """params -> unit -> params is the identity for every attainable point.

    The other direction cannot be (rounding an integer is many-to-one), and
    this is the direction that matters: it is what reloading a trial from the
    store does.
    """
    rng = random.Random(1234)
    for _ in range(500):
        u = [rng.random() for _ in range(mixed_space.dim)]
        params = mixed_space.to_params(u)
        again = mixed_space.to_params(mixed_space.to_unit(params))
        assert again["thickness"] == pytest.approx(params["thickness"], rel=1e-12)
        assert again["nlayers"] == params["nlayers"]
        assert again["tol"] == pytest.approx(params["tol"], rel=1e-12)


def test_integers_stay_in_range_and_are_ints(mixed_space):
    rng = random.Random(7)
    seen = set()
    for _ in range(2000):
        n = mixed_space.to_params([rng.random() for _ in range(3)])["nlayers"]
        assert isinstance(n, int)
        assert 4 <= n <= 40
        seen.add(n)
    # The map must be onto: every layer count has to be reachable, or part of
    # the design space is silently unavailable to the optimizer.
    assert seen == set(range(4, 41))


def test_log_scale_is_monotone_and_geometric(mixed_space):
    values = [mixed_space.to_params([0.0, 0.0, u / 20.0])["tol"] for u in range(21)]
    assert values == sorted(values)
    ratios = [b / a for a, b in zip(values, values[1:], strict=False)]
    assert all(math.isclose(r, ratios[0], rel_tol=1e-9) for r in ratios)


def test_out_of_cube_is_clipped(space2d):
    assert space2d.to_params([-0.5, 1.5]) == {"x0": 0.0, "x1": 1.0}


def test_to_unit_rejects_missing_parameters(mixed_space):
    with pytest.raises(ValueError, match="missing parameter"):
        mixed_space.to_unit({"thickness": 10.0})


def test_bad_specs_are_rejected():
    with pytest.raises(ValueError, match="low < high"):
        ParameterSpec("a", 1.0, 1.0)
    with pytest.raises(ValueError, match="log scale"):
        ParameterSpec("a", 0.0, 1.0, log=True)
    with pytest.raises(ValueError, match="whole numbers"):
        ParameterSpec("a", 0.5, 4.0, kind="int")
    with pytest.raises(ValueError, match="duplicate"):
        SearchSpace([ParameterSpec("a", 0, 1), ParameterSpec("a", 0, 2)])
    with pytest.raises(ValueError, match="both free and fixed"):
        SearchSpace([ParameterSpec("a", 0, 1)], fixed={"a": 3})


def test_from_config():
    space = SearchSpace.from_config(
        {
            "parameters": {
                "w": {"low": 5, "high": 15},
                "n": {"low": 1, "high": 10, "kind": "int"},
                "eps": {"low": 1e-3, "high": 1.0, "log": True},
            },
            "fixed": {"gap": 0.1},
        }
    )
    assert space.dim == 3
    assert space.fixed == {"gap": 0.1}
    assert [s.kind for s in space.specs] == ["float", "int", "float"]
    assert [s.log for s in space.specs] == [False, False, True]
