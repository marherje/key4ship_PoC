"""The geometry bridge: limits, reparametrization, and the feasibility gate.

The claim under test is that the unit cube is feasible *by construction* — that
`derived limit x fill fraction` never produces a geometry the real renderer would
refuse. Everything here is checked against `config.py` itself rather than against
a copy of its formulas, so a change to `SND_compact_template.xml` shows up as a
failure here instead of as a run full of INFEASIBLE trials.

No ROOT and no key4hep needed: `config.py` and `parse_geometry.py` are pure
python, which is what makes this gate cheap enough to run per trial.
"""

from __future__ import annotations

import math
import random

import pytest

from mobo.core.search_space import SearchSpace
from mobo.ship.geometry import Geometry, _count_to_fill, _fill_to_count

# The experiment's search space, kept in step with conf/experiment/snd_proxy.yaml.
SPACE = {
    "parameters": {
        "SiPad_WThickness": {"low": 5.0, "high": 15.0},
        "SiPad_dim_z": {"low": 250, "high": 500},
        "sipad_fill": {"low": 0.3, "high": 1.0},
        "SiTarget_WThickness": {"low": 2.0, "high": 5.0},
        "SiTarget_spacing": {"low": 8.0, "high": 15.0},
        "sitarget_fill": {"low": 0.5, "high": 1.0},
        "xy_gap_frac": {"low": 0.1, "high": 1.0},
    },
    "fixed": {
        "SiPad_frame_gap": 0.1,
        "SiTarget_module_offset": 1,
        "SiPad_layer_gap": "auto",
    },
}


@pytest.fixture(scope="module")
def geometry() -> Geometry:
    return Geometry()


@pytest.fixture(scope="module")
def space() -> SearchSpace:
    return SearchSpace.from_config(SPACE)


# ── the baseline ─────────────────────────────────────────────────────────────


def test_baseline_is_inside_the_search_box(geometry, space):
    """Trial 0 is the baseline, so the baseline must be a point of the cube."""
    params = geometry.baseline_params()
    unit = space.to_unit(params)
    for spec, u in zip(space.specs, unit, strict=False):
        assert 0.0 <= u <= 1.0, (
            f"{spec.name} = {params[spec.name]} is outside [{spec.low}, {spec.high}]"
        )
        # And not *at* a clipped edge by accident: to_unit clamps, so a value
        # outside the range would silently look like a corner.
        assert spec.low <= float(params[spec.name]) <= spec.high


def test_baseline_params_reproduce_the_baseline_constants(geometry):
    """The reparametrization must be exact on the design it is anchored to.

    Compared as rendered, not as raw values: constrained lengths are passed to
    the template as verbatim strings (see `_exact_mm`), so `4.9` and `"4.9*mm"`
    are the same constant expressed two ways.
    """
    constants = geometry.constants_for(geometry.baseline_params())
    assert set(constants) == set(geometry.base_constants)
    for key, value in constants.items():
        assert geometry.config.format_value(key, value) == geometry.config.format_value(
            key, geometry.base_constants[key]
        ), key


def test_regression_baseline_xml_is_unchanged(geometry, tmp_path):
    """Rendering through the optimizer must equal rendering the old way.

    The old way is `make_variants.py` / `config.py`: base constants straight
    from parameters_template.yaml. Same output directory in both cases, since
    the include refs are rewritten relative to it.
    """
    direct = tmp_path / "direct.xml"
    ok, _lines, errors = geometry.config.write_variant(
        dict(geometry.base_constants), direct
    )
    assert ok, errors

    through_search_space = tmp_path / "reparametrized.xml"
    ok, _lines, errors = geometry.write(
        geometry.constants_for(geometry.baseline_params()), through_search_space
    )
    assert ok, errors

    assert direct.read_text() == through_search_space.read_text()


# ── derived limits vs the real renderer ──────────────────────────────────────


VARIANTS = [
    ("baseline", {}),
    ("thin_W", {"SiPad_WThickness": 5}),
    ("short_sipad", {"SiPad_dim_z": 200}),
    ("long_sipad", {"SiPad_dim_z": 500}),
    ("thick_W", {"SiPad_WThickness": 15}),
    ("wide_spacing", {"SiTarget_spacing": 15, "SiTarget_WThickness": 2}),
    ("tight_spacing", {"SiTarget_spacing": 8, "SiTarget_WThickness": 5}),
]


@pytest.mark.parametrize("name,overrides", VARIANTS, ids=[v[0] for v in VARIANTS])
def test_limits_are_exactly_the_boundary_of_feasibility(geometry, name, overrides):
    """N_max must build and N_max + 1 must not. That is what "max" has to mean.

    This is the cross-check that keeps `limits()` honest about the `auto` layer
    gap, whose rule lives in config.resolve_auto and not in the template — the
    one formula this module reproduces rather than reads.
    """
    params = dict(geometry.base_constants, **overrides)
    limits = geometry.limits(params)
    assert limits.sipad_nlayers_max >= 1
    assert limits.sitarget_nlayers_max >= 1

    at_max = dict(
        params,
        SiPad_NLayers=limits.sipad_nlayers_max,
        SiTarget_NLayers=limits.sitarget_nlayers_max,
        SiTarget_XY_plane_gap=limits.sitarget_xy_gap_max,
    )
    assert geometry.check_constants(at_max) is None, (
        f"{name}: the limit itself is infeasible"
    )

    over_sipad = dict(at_max, SiPad_NLayers=limits.sipad_nlayers_max + 1)
    assert geometry.check_constants(over_sipad) is not None, (
        f"{name}: SiPad limit is too low"
    )

    over_sitarget = dict(at_max, SiTarget_NLayers=limits.sitarget_nlayers_max + 1)
    assert geometry.check_constants(over_sitarget) is not None, (
        f"{name}: SiTarget limit is too low"
    )

    over_gap = dict(at_max, SiTarget_XY_plane_gap=limits.sitarget_xy_gap_max * 1.01)
    assert geometry.check_constants(over_gap) is not None, (
        f"{name}: XY gap limit is too low"
    )


@pytest.mark.parametrize("name,overrides", VARIANTS, ids=[v[0] for v in VARIANTS])
def test_limits_match_the_parser(geometry, name, overrides):
    """The template's own `*_max` constants, read back off a rendered variant.

    With a *fixed* layer gap the SiPad limit is exactly the template's
    floor(dim_z / layer_thickness), so all three can be compared directly.
    """
    params = dict(geometry.base_constants, **overrides, SiPad_layer_gap=1.0)
    limits = geometry.limits(params)

    resolved = geometry.resolve(
        dict(
            params,
            SiPad_NLayers=limits.sipad_nlayers_max,
            SiTarget_NLayers=limits.sitarget_nlayers_max,
            SiTarget_XY_plane_gap=limits.sitarget_xy_gap_max,
        )
    )
    assert limits.sipad_nlayers_max == int(resolved["SiPad_NLayers_max"]), name
    assert limits.sitarget_nlayers_max == int(resolved["SiTarget_NLayers_max"]), name
    assert limits.sitarget_xy_gap_max == pytest.approx(
        resolved["SiTarget_XY_plane_gap_max"]
    ), name


def test_the_coupling_between_the_two_detectors_is_respected(geometry):
    """SiTarget_dim_z = 1700 mm - SiPad_dim_z: growing one shrinks the other.

    The trap the README warns about — changing SiPad_dim_z and forgetting to
    recount SiTarget layers — cannot happen here, because the layer count is a
    fraction of a limit that already knows about the coupling.
    """
    short = geometry.limits(dict(geometry.base_constants, SiPad_dim_z=250))
    long = geometry.limits(dict(geometry.base_constants, SiPad_dim_z=500))
    assert short.sitarget_nlayers_max > long.sitarget_nlayers_max
    assert short.sipad_nlayers_max < long.sipad_nlayers_max


# ── the cube is feasible by construction ─────────────────────────────────────


def test_every_corner_of_the_cube_builds(geometry, space):
    """2^7 corners: the extremes are where a bad parametrization breaks."""
    import itertools

    for corner in itertools.product([0.0, 1.0], repeat=space.dim):
        params = space.to_params(corner)
        assert geometry.check(params) is None, (
            f"corner {corner} -> {geometry.check(params)}"
        )


def test_a_sample_of_the_cube_builds(geometry, space):
    import random

    rng = random.Random(20260801)
    for _ in range(60):
        params = space.to_params([rng.random() for _ in range(space.dim)])
        assert geometry.check(params) is None, params


def test_fills_land_inside_the_limits(geometry, space):
    import random

    rng = random.Random(7)
    for _ in range(40):
        params = space.to_params([rng.random() for _ in range(space.dim)])
        limits = geometry.limits(params)
        constants = geometry.constants_for(params)
        assert 1 <= constants["SiPad_NLayers"] <= limits.sipad_nlayers_max
        assert 1 <= constants["SiTarget_NLayers"] <= limits.sitarget_nlayers_max
        gap = float(str(constants["SiTarget_XY_plane_gap"]).replace("*mm", ""))
        assert 0 <= gap <= limits.sitarget_xy_gap_max


# ── the gate rejects what it should ──────────────────────────────────────────


def test_infeasible_geometry_is_rejected_without_raising(geometry):
    """The classic mistake: shrink SiPad and keep the old SiTarget layer count.

    SiTarget_dim_z = 1700 - 200 = 1500 mm, which at 11 mm spacing holds 136
    layers, not the 120+ the baseline would imply... and SiPad at 200 mm holds
    12 layers, not 22.
    """
    bad = dict(geometry.base_constants, SiPad_dim_z=200)  # keeps SiPad_NLayers=22
    reason = geometry.check_constants(bad)
    assert reason is not None
    assert "SiPad" in reason


def test_touching_envelopes_are_not_an_overlap(geometry):
    """SiTarget ends exactly where SiPad begins — by construction, not by luck.

    Those two z values come out of different expressions that are algebraically
    equal, so for a SiPad_dim_z that is not a round number they differ in the
    last bits. A strict inequality in the overlap check turns that into a
    spurious rejection, which an optimizer proposing arbitrary floats hits
    constantly (config.OVERLAP_TOL exists for this).
    """
    awkward = dict(geometry.base_constants, SiPad_dim_z=485.2903311418265)
    limits = geometry.limits(awkward)
    constants = dict(
        awkward,
        SiPad_NLayers=limits.sipad_nlayers_max,
        SiTarget_NLayers=limits.sitarget_nlayers_max,
    )
    assert geometry.check_constants(constants) is None

    # The tolerance that makes this work must stay far below anything physical,
    # so it can never hide a real overlap: 1 nm against detectors of 100s of mm.
    assert 0 < geometry.config.OVERLAP_TOL <= 1e-6


def test_unknown_parameter_names_are_a_hard_error(geometry):
    with pytest.raises(RuntimeError, match="neither a template constant"):
        geometry.constants_for({"SiPad_WThicknes": 10})  # typo


def test_check_reports_the_reason_rather_than_raising(geometry):
    """A geometry so short that `auto` cannot size its slack at all."""
    reason = geometry.check_constants(
        dict(geometry.base_constants, SiPad_dim_z=50, SiPad_NLayers=22)
    )
    assert reason is not None and reason


# ── the fill <-> count map ───────────────────────────────────────────────────


def test_fill_to_count_is_bounded_and_onto():
    maximum = 23
    counts = {_fill_to_count(f / 1000.0, maximum) for f in range(0, 1001)}
    assert min(counts) == 1 and max(counts) == maximum
    assert counts == set(range(1, maximum + 1))


def test_count_to_fill_round_trips_every_count():
    for maximum in (1, 5, 23, 120, 137):
        for count in range(1, maximum + 1):
            assert _fill_to_count(_count_to_fill(count, maximum), maximum) == count


def test_count_to_fill_lands_mid_interval():
    """Mid-interval, not on the edge: an edge value is one ulp from the wrong count."""
    assert _count_to_fill(22, 23) == pytest.approx(22.5 / 23)
    assert _count_to_fill(23, 23) == 1.0  # clipped, the top of the range
    assert math.isclose(_fill_to_count(1.0, 23), 23)


# ── the boundary of the box, where the acquisition likes to sit ──────────────


@pytest.mark.parametrize("free", ["xy_gap_frac", "sitarget_fill", "sipad_fill"])
def test_a_fill_pinned_at_its_maximum_always_builds(geometry, space, free):
    """The baseline sits at the top of two of these ranges, so the scan runs
    downward from a boundary — and boundaries are exactly where `optimize_acqf`
    puts its candidates. A limit that is only reachable in theory is not a
    limit, it is a third of the budget spent on INFEASIBLE trials.

    This caught a real one: `SiTarget_XY_plane_gap` is compared against its own
    derived maximum inside the template, and rendering it through
    `format_value`'s six significant digits rounded it *up* past that maximum
    for about a third of parameter combinations.
    """
    index = space.names.index(free)
    rng = random.Random(90210)
    for _ in range(40):
        unit = [rng.random() for _ in range(space.dim)]
        unit[index] = 1.0
        params = space.to_params(unit)
        assert geometry.check(params) is None, (
            f"{free}=1.0 rejected: {geometry.check(params)}"
        )


def test_the_gap_at_its_maximum_survives_rendering(geometry):
    """The rendered XML must resolve to a gap that is <= the resolved maximum.

    Checked against the parser rather than against our own arithmetic, since
    the failure mode was precisely a disagreement between the two.
    """
    rng = random.Random(7)
    for _ in range(25):
        params = dict(
            geometry.base_constants,
            SiTarget_spacing=rng.uniform(8.0, 15.0),
            SiTarget_WThickness=rng.uniform(2.0, 5.0),
        )
        limits = geometry.limits(params)
        constants = geometry.constants_for(dict(params, xy_gap_frac=1.0))
        resolved = geometry.resolve(
            dict(constants, SiTarget_NLayers=limits.sitarget_nlayers_max)
        )
        assert resolved["SiTarget_XY_plane_gap"] <= resolved["SiTarget_XY_plane_gap_max"]
