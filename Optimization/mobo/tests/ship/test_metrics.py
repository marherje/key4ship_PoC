"""The analytic metrics, against numbers worked out by hand.

`cost_proxy` is what half the optimization is steering by, so it gets checked
digit by digit rather than "looks about right": every factor in it (two silicon
planes per SiTarget layer, sensor area vs envelope area, mm³ -> cm³ -> kg) is a
place where a silent factor of two would produce a plausible-looking and
completely wrong Pareto front.
"""

from __future__ import annotations

import pytest

from mobo.ship.metrics import (
    CostModel,
    analytic_metrics,
    channel_counts,
    silicon_areas_m2,
    tungsten_masses_kg,
)

# The baseline geometry, as the template resolves it. Written out explicitly so
# this file tests the arithmetic and not the parser.
BASELINE = {
    # SiPad: 3x3 ASUs of 2x2 wafers; a wafer is 16 pads of 5.53 mm plus a
    # 0.61 mm rim per side -> 88.48 + 1.22 = 89.70 mm.
    "SiPad_NLayers": 22,
    "SiPad_NASUsX": 3,
    "SiPad_NASUsY": 3,
    "SiPad_NWafersX": 2,
    "SiPad_NWafersY": 2,
    "SiPad_WaferSizeX": 89.70,
    "SiPad_WaferSizeY": 89.70,
    "SiPad_WThickness": 10.0,
    "SiPad_dim_x": 540.0,
    "SiPad_dim_y": 540.0,
    "SiPad_NcellsX": 96,
    "SiPad_NcellsY": 96,
    # SiTarget: 4x2 sensors per plane, two planes (x and y) per layer.
    "SiTarget_NLayers": 120,
    "SiTarget_ncols": 4,
    "SiTarget_nrows": 2,
    "SiTarget_sensor_width": 99.25,
    "SiTarget_sensor_height": 199.5,
    "SiTarget_WThickness": 3.5,
    "SiTarget_env_width": 400.0,
    "SiTarget_env_height": 400.0,
    "SiTarget_strip_pitch": 0.0755,
}


def test_silicon_area_by_hand():
    # SiPad: one plane is (3*2*89.70)^2 = 538.2^2 mm^2, 22 of them.
    plane = (3 * 2 * 89.70) ** 2 / 1e6
    assert plane == pytest.approx(0.28965924)
    # SiTarget: (4*99.25) x (2*199.5) = 397 x 399 mm^2, TWO planes per layer.
    sitarget_plane = (4 * 99.25) * (2 * 199.5) / 1e6
    assert sitarget_plane == pytest.approx(0.158403)

    areas = silicon_areas_m2(BASELINE)
    assert areas["sipad"] == pytest.approx(22 * plane)
    assert areas["sitarget"] == pytest.approx(2 * 120 * sitarget_plane)
    assert areas["sitarget"] == pytest.approx(38.01672)


def test_tungsten_mass_by_hand():
    # SiPad: 22 plates of 10 mm over 540x540 mm -> 64152 cm^3 of tungsten.
    masses = tungsten_masses_kg(BASELINE)
    assert masses["sipad"] == pytest.approx(64152.0 * 19.3 / 1000.0)
    assert masses["sipad"] == pytest.approx(1238.1336)
    # SiTarget: 120 plates of 3.5 mm over 400x400 mm.
    assert masses["sitarget"] == pytest.approx(67200.0 * 19.3 / 1000.0)
    assert masses["sitarget"] == pytest.approx(1296.96)


def test_cost_proxy_by_hand():
    metrics = analytic_metrics(BASELINE, CostModel(si_per_m2=30.0, w_per_kg=0.1))
    si = 22 * (538.2**2) / 1e6 + 2 * 120 * (397 * 399) / 1e6
    w = (22 * 10 * 540 * 540 + 120 * 3.5 * 400 * 400) / 1e3 * 19.3 / 1e3
    assert metrics["si_area_m2"] == pytest.approx(si)
    assert metrics["w_mass_kg"] == pytest.approx(w)
    assert metrics["cost_proxy"] == pytest.approx(30.0 * si + 0.1 * w)
    assert metrics["cost_proxy"] == pytest.approx(1585.186, abs=1e-3)


def test_the_two_silicon_planes_per_sitarget_layer_are_not_forgotten():
    """The single likeliest factor-of-two in the whole cost model."""
    one_plane = (4 * 99.25) * (2 * 199.5) / 1e6
    assert silicon_areas_m2(BASELINE)["sitarget"] / (120 * one_plane) == pytest.approx(
        2.0
    )


def test_channel_counts_by_hand():
    counts = channel_counts(BASELINE)
    assert counts["sipad"] == 22 * 96 * 96
    # Per sensor: 99.25/0.0755 = 1314 strips across, 199.5/0.0755 = 2642 down;
    # 8 sensors, both planes, 120 layers.
    assert counts["sitarget"] == 120 * 8 * (1314 + 2642)


def test_coefficients_come_from_the_config():
    cheap = analytic_metrics(BASELINE, CostModel(si_per_m2=1.0, w_per_kg=0.0))
    assert cheap["cost_proxy"] == pytest.approx(cheap["si_area_m2"])
    only_w = analytic_metrics(BASELINE, CostModel(si_per_m2=0.0, w_per_kg=1.0))
    assert only_w["cost_proxy"] == pytest.approx(only_w["w_mass_kg"])


def test_cost_model_from_config_ignores_unknown_keys():
    model = CostModel.from_config({"si_per_m2": 12.0, "nonsense": 1})
    assert model.si_per_m2 == 12.0
    assert model.w_per_kg == CostModel().w_per_kg


# ── monotonicity: the properties the optimizer will exploit ──────────────────


@pytest.mark.parametrize(
    "key",
    ["SiPad_NLayers", "SiTarget_NLayers", "SiPad_WThickness", "SiTarget_WThickness"],
)
def test_cost_grows_with_every_knob(key):
    less = analytic_metrics(dict(BASELINE, **{key: BASELINE[key] * 0.5}))
    more = analytic_metrics(dict(BASELINE, **{key: BASELINE[key] * 2.0}))
    assert more["cost_proxy"] > less["cost_proxy"]


@pytest.mark.parametrize("key", ["SiPad_NLayers", "SiTarget_NLayers"])
def test_channels_grow_with_the_layer_count(key):
    less = analytic_metrics(dict(BASELINE, **{key: BASELINE[key] // 2}))
    more = analytic_metrics(dict(BASELINE, **{key: BASELINE[key] * 2}))
    assert more["n_channels"] > less["n_channels"]


def test_absorber_thickness_does_not_change_the_channel_count():
    """Thicker tungsten costs money and buys showers, not readout channels."""
    thick = analytic_metrics(dict(BASELINE, SiPad_WThickness=15.0))
    thin = analytic_metrics(dict(BASELINE, SiPad_WThickness=5.0))
    assert thick["n_channels"] == thin["n_channels"]
    assert thick["cost_proxy"] > thin["cost_proxy"]


# ── against the real geometry ────────────────────────────────────────────────


def test_matches_the_constants_the_parser_resolves():
    """The hand-written table above must really be the baseline geometry."""
    from mobo.ship.geometry import Geometry

    geometry = Geometry()
    resolved = geometry.resolve(geometry.base_constants)
    for key, value in BASELINE.items():
        assert float(resolved[key]) == pytest.approx(float(value)), key

    assert analytic_metrics(resolved)["cost_proxy"] == pytest.approx(
        analytic_metrics(BASELINE)["cost_proxy"]
    )
