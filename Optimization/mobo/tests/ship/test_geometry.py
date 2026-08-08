"""The geometry bridge: derived limits and the feasibility gate.

The claim under test is that the search box is feasible *by construction* —
that no point of it produces a geometry the real renderer would refuse for any
reason other than the layers genuinely not fitting. That is what
`SiTarget_XY_plane_gap_frac` buys: a fraction of a budget the template derives,
so the whole [0, 1] interval is legal whatever the pitch turns out to be.
Everything here is checked against `config.py` itself rather than against a copy
of its formulas, so a change to `SND_compact_template.xml` shows up as a failure
here instead of as a run full of INFEASIBLE trials.

No ROOT and no key4hep needed: `config.py` and `parse_geometry.py` are pure
python, which is what makes this gate cheap enough to run per trial.
"""

from __future__ import annotations

import random

import pytest

from mobo.core.search_space import SearchSpace
from mobo.ship.geometry import Geometry

# The experiment's search space, kept in step with conf/experiment/snd_proxy.yaml.
SPACE = {
    "parameters": {
        "SiPad_WThickness": {"low": 10.0, "high": 30.0},
        "SiPad_dim_z": {"low": 220, "high": 1000},
        "SiPad_NLayers": {"low": 10, "high": 20, "kind": "int"},
        "SiTarget_WThickness": {"low": 2.0, "high": 5.0},
        "SiTarget_NLayers": {"low": 80, "high": 120, "kind": "int"},
        "SiTarget_XY_plane_gap_frac": {"low": 0.0, "high": 1.0},
    },
    "fixed": {
        "SiPad_frame_gap": 0.1,
        "SiTarget_module_offset": 1,
        "SiPad_layer_gap": "auto",
        "SiTarget_layer_gap": "auto",
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
    """Trial 0 must be the baseline exactly, not approximately.

    A round trip through the search space is now the identity — every parameter
    is a template constant — so this is cheap to assert and there is no longer
    any inversion that could drift. Still compared as *rendered*, since that is
    what actually reaches the XML.
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
    ("thin_W", {"SiPad_WThickness": 10}),
    ("short_sipad", {"SiPad_dim_z": 220}),
    ("long_sipad", {"SiPad_dim_z": 1000}),
    ("thick_W", {"SiPad_WThickness": 30, "SiPad_dim_z": 1000}),
    ("few_sitarget_layers", {"SiTarget_NLayers": 80, "SiTarget_WThickness": 2}),
    ("thick_sitarget_W", {"SiTarget_NLayers": 120, "SiTarget_WThickness": 5}),
]


@pytest.mark.parametrize("name,overrides", VARIANTS, ids=[v[0] for v in VARIANTS])
def test_limits_are_exactly_the_boundary_of_feasibility(geometry, name, overrides):
    """N_max must build and N_max + 1 must not. That is what "max" has to mean.

    This is the cross-check that keeps `limits()` honest about the `auto` layer
    gap, whose rule lives in config.resolve_auto and not in the template — the
    one formula this module reproduces rather than reads. Both detectors are on
    `auto` now, so the claim is made for both.

    The two bounds cannot be probed at the same point, and that is not a defect
    of the test. Since the pitch became `dim_z / NLayers`, the XY gap budget is
    what the pitch leaves over, so the layer count and the gap trade against each
    other: `sitarget_nlayers_max` is the count that fits with *no* gap, and
    `sitarget_xy_gap_max` is the gap that fits at the count actually asked for.
    Pinning both at once asks for a layer thinner than its own contents.
    """
    params = dict(geometry.base_constants, **overrides)
    limits = geometry.limits(params)
    assert limits.sipad_nlayers_max >= 1
    assert limits.sitarget_nlayers_max >= 1

    # -- the layer counts, at the zero gap their bound assumes ----------------
    at_max = dict(
        params,
        SiPad_NLayers=limits.sipad_nlayers_max,
        SiTarget_NLayers=limits.sitarget_nlayers_max,
        SiTarget_XY_plane_gap_frac=0,
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

    # -- the XY gap, at the layer count its bound was computed for ------------
    # The fraction *is* the bound now: 1 must build (that is what "the Y plane
    # can sit against the next absorber" means) and anything past it must not.
    # The gap in mm is read back off the rendered geometry rather than asserted
    # here, so the template stays the only place that owns the arithmetic.
    at_gap = dict(params, SiTarget_XY_plane_gap_frac=1.0)
    assert geometry.check_constants(at_gap) is None, (
        f"{name}: the XY gap limit is itself infeasible"
    )
    resolved = geometry.resolve(at_gap)
    assert resolved["SiTarget_XY_plane_gap"] == resolved["SiTarget_XY_plane_gap_max"], (
        f"{name}: frac = 1 must land exactly on the maximum, not near it"
    )

    over_gap = dict(params, SiTarget_XY_plane_gap_frac=1.01)
    assert geometry.check_constants(over_gap) is not None, (
        f"{name}: a fraction above 1 must be rejected"
    )


@pytest.mark.parametrize("name,overrides", VARIANTS, ids=[v[0] for v in VARIANTS])
def test_limits_match_the_parser(geometry, name, overrides):
    """The template's own `*_max` constants, read back off a rendered variant.

    With *fixed* layer gaps both limits are exactly the template's
    floor(dim_z / layer_thickness), so all three can be compared directly. The
    `auto` case cannot be checked this way — that is what the feasibility test
    above is for — because `auto` sizes the gap so the layers span the envelope,
    which makes floor(dim_z / thickness) return the count it was handed.
    """
    params = dict(
        geometry.base_constants,
        **overrides,
        SiPad_layer_gap=1.0,
        SiTarget_layer_gap=1.0,
    )
    limits = geometry.limits(params)

    resolved = geometry.resolve(
        dict(
            params,
            SiPad_NLayers=limits.sipad_nlayers_max,
            SiTarget_NLayers=limits.sitarget_nlayers_max,
            SiTarget_XY_plane_gap_frac=0,
        )
    )
    assert limits.sipad_nlayers_max == int(resolved["SiPad_NLayers_max"]), name
    assert limits.sitarget_nlayers_max == int(resolved["SiTarget_NLayers_max"]), name


def test_the_coupling_between_the_two_detectors_is_respected(geometry):
    """SiTarget_dim_z = SiDetector_total_dim_z - SiPad_dim_z: growing one shrinks
    the other.

    The trap the README warns about — changing SiPad_dim_z and forgetting to
    check that the SiTarget still fits — cannot happen here, because the limit
    already knows about the coupling: a longer SiPad leaves room for fewer
    SiTarget layers, and for more of its own.
    """
    short = geometry.limits(dict(geometry.base_constants, SiPad_dim_z=220))
    long = geometry.limits(dict(geometry.base_constants, SiPad_dim_z=1000))
    assert short.sitarget_nlayers_max > long.sitarget_nlayers_max
    assert short.sipad_nlayers_max < long.sipad_nlayers_max


# ── the cube is feasible by construction ─────────────────────────────────────


def test_every_feasible_corner_of_the_cube_builds(geometry, space):
    """The extremes are where a bad parametrization breaks.

    Not every corner is meant to build any more: with the layer counts searched
    directly, "many thick SiPad layers in a short envelope" is a corner that
    genuinely does not exist, and the gate is supposed to say so. What is tested
    is the implication — whenever the rigid content fits, the renderer accepts
    it — so a rejection can only ever come from the geometry not fitting, never
    from the parametrization mangling a feasible point.
    """
    import itertools

    built = 0
    for corner in itertools.product([0.0, 1.0], repeat=space.dim):
        params = space.to_params(corner)
        reason = geometry.check(params)
        if reason is None:
            built += 1
        else:
            # `no caben` is resolve_auto's way of saying the rigid content of the
            # layers exceeds the envelope. Anything else means the
            # parametrization produced a nonsensical geometry, which is the bug
            # this test exists to catch.
            assert "no caben" in reason, f"corner {corner} rejected for: {reason}"
    assert built >= 2 ** (space.dim - 1), (
        f"only {built}/{2**space.dim} corners build; the box is mostly empty"
    )


def test_most_of_the_cube_builds(geometry, space):
    """The interior has to be usable, not merely non-crashing.

    The only geometries the gate may refuse are the ones that do not fit; if
    that were most of the cube the optimizer would spend its budget discovering
    the parametrization instead of the physics, which is the whole reason the
    fill fractions existed. Two thirds is a floor, not a target — the measured
    figure at the ranges of snd_proxy is around 80%.
    """
    import random

    rng = random.Random(20260801)
    ok = 0
    trials = 60
    for _ in range(trials):
        params = space.to_params([rng.random() for _ in range(space.dim)])
        reason = geometry.check(params)
        if reason is None:
            ok += 1
        else:
            assert "no caben" in reason, (
                f"rejected for something other than not fitting: {reason}"
            )
    assert ok >= 2 * trials // 3, f"only {ok}/{trials} of the cube builds"


def test_the_xy_gap_lands_inside_its_limit(geometry, space):
    """Anywhere in the box, the derived gap stays a legal slice within budget.

    Read off the *rendered* geometry: the fraction is what this layer passes on,
    and the whole claim is that the template turns it into a length that is both
    positive (Geant4 rejects a slice of null extent) and no larger than the
    budget it is a fraction of.
    """
    import random

    rng = random.Random(7)
    for _ in range(40):
        params = space.to_params([rng.random() for _ in range(space.dim)])
        resolved = geometry.resolve(geometry.constants_for(params))
        gap, gap_max = (
            resolved["SiTarget_XY_plane_gap"],
            resolved["SiTarget_XY_plane_gap_max"],
        )
        assert 0 < gap <= gap_max


# ── the gate rejects what it should ──────────────────────────────────────────


def test_infeasible_geometry_is_rejected_without_raising(geometry):
    """The classic mistake: keep a layer count a shrunken envelope cannot hold.

    20 SiPad layers of 30 mm of tungsten are 711 mm of rigid content, and the
    envelope offered is 250 mm. `auto` has no slack to give and says so instead
    of raising.
    """
    bad = dict(
        geometry.base_constants,
        SiPad_dim_z=250,
        SiPad_WThickness=30,
        SiPad_NLayers=20,
    )
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
        SiTarget_XY_plane_gap_frac=0,  # the count bound is the one with no gap
    )
    assert geometry.check_constants(constants) is None

    # The tolerance that makes this work must stay far below anything physical,
    # so it can never hide a real overlap: 1 nm against detectors of 100s of mm.
    assert 0 < geometry.config.OVERLAP_TOL <= 1e-6


def test_unknown_parameter_names_are_a_hard_error(geometry):
    with pytest.raises(RuntimeError, match="not a template constant"):
        geometry.constants_for({"SiPad_WThicknes": 10})  # typo


def test_check_reports_the_reason_rather_than_raising(geometry):
    """A geometry so short that `auto` cannot size its slack at all."""
    reason = geometry.check_constants(
        dict(geometry.base_constants, SiPad_dim_z=50, SiPad_NLayers=22)
    )
    assert reason is not None and reason


# ── the boundary of the box, where the acquisition likes to sit ──────────────


def test_the_xy_gap_pinned_at_its_maximum_always_builds(geometry, space):
    """Boundaries are exactly where `optimize_acqf` puts its candidates, so a
    limit that is only reachable in theory is not a limit — it is a third of the
    budget spent on INFEASIBLE trials.

    This caught a real one, back when this layer turned the fraction into a
    length itself: the gap is compared against its own derived maximum inside
    the template, and rendering it through `format_value`'s six significant
    digits rounded it *up* past that maximum for about a third of parameter
    combinations. The interpolation lives in the template now and the fraction
    is written at full precision, so the failure mode is gone by construction —
    this test is what keeps it that way.

    Only the geometries that fit are asserted on: a fraction of 1.0 cannot
    rescue a layer count the envelope was never going to hold.
    """
    index = space.names.index("SiTarget_XY_plane_gap_frac")
    rng = random.Random(90210)
    checked = 0
    for _ in range(40):
        unit = [rng.random() for _ in range(space.dim)]
        unit[index] = 1.0
        params = space.to_params(unit)
        limits = geometry.limits(params)
        constants = geometry.constants_for(params)
        if (
            constants["SiPad_NLayers"] > limits.sipad_nlayers_max
            or constants["SiTarget_NLayers"] > limits.sitarget_nlayers_max
        ):
            continue
        assert geometry.check(params) is None, (
            f"SiTarget_XY_plane_gap_frac=1.0 rejected: {geometry.check(params)}"
        )
        checked += 1
    assert checked, "no feasible sample to check the gap maximum on"


def test_the_gap_at_its_maximum_survives_rendering(geometry):
    """The rendered XML must resolve to a gap that is <= the resolved maximum.

    Checked against the parser rather than against our own arithmetic, since
    the failure mode was precisely a disagreement between the two.
    """
    rng = random.Random(7)
    for _ in range(25):
        params = dict(
            geometry.base_constants,
            SiTarget_NLayers=rng.randint(80, 120),
            SiTarget_WThickness=rng.uniform(2.0, 5.0),
        )
        constants = geometry.constants_for(dict(params, SiTarget_XY_plane_gap_frac=1.0))
        resolved = geometry.resolve(constants)
        assert resolved["SiTarget_XY_plane_gap"] <= resolved["SiTarget_XY_plane_gap_max"]
