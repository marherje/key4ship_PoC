"""Bridge to the existing DD4hep variant machinery, and the reparametrization.

Two jobs.

**1. Reparametrization.** The natural design parameters are not a good search
space. `SiPad_NLayers` is bounded by `floor(SiPad_dim_z / layer_thickness)`,
which depends on `SiPad_WThickness` and `SiPad_dim_z`; `SiTarget_XY_plane_gap`
is bounded by a budget built from `SiTarget_spacing`, `SiTarget_WThickness` and
the module offset. A box in those coordinates is mostly infeasible, and an
optimizer that spends its budget discovering that learns nothing about physics.
So the search space uses *fill fractions* instead — "how much of the space that
is available do we use" — which makes the unit cube feasible by construction:

    SiPad_NLayers          = floor(sipad_fill    x SiPad_NLayers_max(params))
    SiTarget_NLayers       = floor(sitarget_fill x SiTarget_NLayers_max(params))
    SiTarget_XY_plane_gap  =       xy_gap_frac   x SiTarget_XY_plane_gap_max(params)

**2. Not duplicating the template.** Those `*_max` are computed by rendering the
template with a probe value and resolving its `<define>` block with the very
parser the rest of the project uses — not by copying the formulas here. A change
to `SND_compact_template.xml` propagates on its own; only the `auto` slack rule
(which lives in `config.resolve_auto`, not in the template) is reproduced, and
`tests/ship/test_geometry.py` pins it against the real renderer.

`config.write_variant(..., dry_run=True)` is then the final gate: cheap, exact,
and it costs no CPU to reject a geometry that would not have built.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..paths import repo_root  # re-exported: the ship layer's callers expect it here


def geometry_dir() -> Path:
    return repo_root() / "Optimization" / "Simulation" / "geometry"


def load_module(name: str, path: Path) -> Any:
    """Import a file by path, exactly as `controller.py` does and for the reason.

    `config` and `make_variants` are generic module names, and key4hep's
    PYTHONPATH is crowded enough that a plain import can pick up something else
    entirely.
    """
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Free parameters this module derives rather than passes straight through.
FILL_PARAMS = {
    "sipad_fill": "SiPad_NLayers",
    "sitarget_fill": "SiTarget_NLayers",
    "xy_gap_frac": "SiTarget_XY_plane_gap",
}


@dataclass(frozen=True)
class Limits:
    """The derived upper bounds the fill fractions are fractions *of*."""

    sipad_nlayers_max: int
    sitarget_nlayers_max: int
    sitarget_xy_gap_max: float

    def as_dict(self) -> dict[str, float]:
        return {
            "SiPad_NLayers_max": float(self.sipad_nlayers_max),
            "SiTarget_NLayers_max": float(self.sitarget_nlayers_max),
            "SiTarget_XY_plane_gap_max": self.sitarget_xy_gap_max,
        }


class Geometry:
    """Renders and validates SND variants, and knows their derived limits."""

    def __init__(
        self,
        geo_dir: str | Path | None = None,
        params_file: str | Path | None = None,
        template: str | Path | None = None,
    ) -> None:
        self.geo_dir = Path(geo_dir) if geo_dir else geometry_dir()
        self.config = load_module("mobo_geo_config", self.geo_dir / "config.py")
        self.template = Path(template) if template else Path(self.config.TEMPLATE)
        self.params_file = Path(params_file) if params_file else Path(self.config.PARAMS)
        self.base_name, self.base_constants = self.config.load_params(self.params_file)

    # ── resolving the template ───────────────────────────────────────────────

    @property
    def _template_text(self) -> str:
        return self.template.read_text()

    def resolve(self, constants: dict[str, Any]) -> dict[str, float]:
        """Every `<define>` constant of the rendered template, evaluated.

        Pure python (xml.etree + the project's own resolver), so this works
        without ROOT and costs milliseconds — cheap enough to call per trial.
        """
        text = self.config.render(self._template_text, constants, self.geo_dir)
        return self.config._resolve_constants(text, self.geo_dir)

    def limits(self, params: dict[str, Any]) -> Limits:
        """Derived bounds for a set of *direct* parameters.

        Computed on a probe geometry — one layer everywhere, no XY gap — because
        the bounds must not depend on the very quantities they bound.
        """
        constants = self.direct_constants(params)
        auto = (
            str(constants.get("SiPad_layer_gap", "")).strip().lower() == self.config.AUTO
        )
        probe = dict(
            constants,
            SiPad_NLayers=1,
            SiTarget_NLayers=1,
            SiTarget_XY_plane_gap=0,
        )
        if auto:
            probe["SiPad_layer_gap"] = self.config.AUTO_MIN
        c = self.resolve(probe)

        if auto:
            # `auto` sizes the air slice so that N layers span dim_z exactly, so
            # the binding constraint is not floor(dim_z/thickness) but the point
            # where the slack itself would have to go negative:
            #     slack = dim_z/N - rigid - AUTO_EPS >= AUTO_MIN
            # (`rigid` is the layer without its slack; the probe injected
            # AUTO_MIN as the gap, so take it back out.)
            rigid = c["SiPad_LayerThickness"] - self.config.AUTO_MIN
            budget = rigid + self.config.AUTO_EPS + self.config.AUTO_MIN
            sipad_max = int(math.floor(c["SiPad_dim_z"] / budget)) if budget > 0 else 0
        else:
            sipad_max = int(c["SiPad_NLayers_max"])

        return Limits(
            sipad_nlayers_max=max(0, sipad_max),
            sitarget_nlayers_max=max(0, int(c["SiTarget_NLayers_max"])),
            sitarget_xy_gap_max=float(c["SiTarget_XY_plane_gap_max"]),
        )

    # ── parameters -> constants ──────────────────────────────────────────────

    def direct_constants(self, params: dict[str, Any]) -> dict[str, Any]:
        """Base constants overridden by the parameters that *are* constants."""
        constants = dict(self.base_constants)
        unknown = [k for k in params if k not in constants and k not in FILL_PARAMS]
        if unknown:
            # Same typo guard as make_variants.py: a parameter that matches
            # nothing is a bug in the config, not something to ignore.
            raise RuntimeError(
                "parameter(s) that are neither a template constant nor a known "
                "derived parameter: " + ", ".join(sorted(unknown))
            )
        for key, value in params.items():
            if key in constants:
                constants[key] = value
        return constants

    def constants_for(self, params: dict[str, Any]) -> dict[str, Any]:
        """Search-space parameters -> the full constants dict `config.py` wants."""
        constants = self.direct_constants(params)
        if not any(k in params for k in FILL_PARAMS):
            return constants

        limits = self.limits(params)
        if "sipad_fill" in params:
            constants["SiPad_NLayers"] = _fill_to_count(
                params["sipad_fill"], limits.sipad_nlayers_max
            )
        if "sitarget_fill" in params:
            constants["SiTarget_NLayers"] = _fill_to_count(
                params["sitarget_fill"], limits.sitarget_nlayers_max
            )
        if "xy_gap_frac" in params:
            constants["SiTarget_XY_plane_gap"] = self._length_within(
                "SiTarget_XY_plane_gap",
                float(params["xy_gap_frac"]) * limits.sitarget_xy_gap_max,
                limits.sitarget_xy_gap_max,
            )
        return constants

    def _length_within(self, name: str, value: float, limit: float) -> Any:
        """A length that still renders to at most `limit`.

        `config.format_value` writes a float as `"%g*mm"` — six significant
        digits — and the template then compares this constant against its own
        derived maximum. Six digits round *up* about half the time, which turns
        `gap == gap_max` into `gap > gap_max` and gets the geometry rejected.
        The baseline never showed it (4.9 mm is exact in six digits); an
        optimizer proposing arbitrary floats, whose acquisition maxima like to
        sit on boundaries, lost about a third of the trials it put there.

        So: ask the renderer what it would actually write. If that lands at or
        below the limit, keep the plain float and the XML is unchanged — which
        is what keeps trial 0 byte-identical to the committed baseline. Only
        when it would overshoot is the value written as a full-precision
        verbatim string instead, which `format_value` passes through untouched
        (`config.resolve_auto` writes `SiPad_layer_gap` the same way).
        """
        value = min(float(value), float(limit))
        rendered = str(self.config.format_value(name, value)).replace("*mm", "")
        try:
            if float(rendered) <= limit:
                return value
        except ValueError:  # not a plain number: leave it to the caller
            return value
        # 17 significant digits round-trip a float exactly, so what the parser
        # reads back is bit-for-bit the value checked here.
        return f"{value:.17g}*mm"

    def baseline_params(self) -> dict[str, Any]:
        """The template's own design, expressed in search-space coordinates.

        Integer fills are inverted to the *midpoint* of the interval that maps
        back to the baseline count: the edge value is exact in real arithmetic
        but one ulp of floating point either way changes the floor, and trial 0
        has to reproduce the baseline geometry byte for byte.
        """
        params: dict[str, Any] = dict(self.base_constants)
        limits = self.limits(self.base_constants)

        params["sipad_fill"] = _count_to_fill(
            int(self.base_constants["SiPad_NLayers"]), limits.sipad_nlayers_max
        )
        params["sitarget_fill"] = _count_to_fill(
            int(self.base_constants["SiTarget_NLayers"]), limits.sitarget_nlayers_max
        )
        params["xy_gap_frac"] = (
            float(self.base_constants["SiTarget_XY_plane_gap"])
            / limits.sitarget_xy_gap_max
            if limits.sitarget_xy_gap_max > 0
            else 0.0
        )
        return params

    # ── the gate ─────────────────────────────────────────────────────────────

    def check(self, params: dict[str, Any]) -> str | None:
        """Would this geometry build? None if yes, the reason if not.

        A dry run of the real renderer: identical validation, nothing written.
        """
        try:
            constants = self.constants_for(params)
        except RuntimeError as exc:
            return str(exc)
        return self.check_constants(constants)

    def check_constants(self, constants: dict[str, Any]) -> str | None:
        try:
            ok, _lines, errors = self.config.write_variant(
                constants, self.geo_dir / "variants" / "_gate.xml", dry_run=True
            )
        except RuntimeError as exc:
            # Raised before validation even runs: an unsizable `auto`, an
            # unresolved placeholder. Same outcome, so report it the same way.
            return str(exc)
        return None if ok else "; ".join(errors)

    def write(
        self, constants: dict[str, Any], out_path: str | Path
    ) -> tuple[bool, list[str], list[str]]:
        """Render and write one compact file. Returns (ok, summary, errors)."""
        return self.config.write_variant(constants, Path(out_path), dry_run=False)


# ── the fill <-> count map ───────────────────────────────────────────────────


def _fill_to_count(fill: float, maximum: int) -> int:
    """Fraction of the available layers -> a layer count in [1, maximum]."""
    if maximum < 1:
        return 0
    return max(1, min(maximum, int(math.floor(float(fill) * maximum))))


def _count_to_fill(count: int, maximum: int) -> float:
    """Inverse of `_fill_to_count`, landing in the middle of the interval."""
    if maximum < 1:
        return 0.0
    return min(1.0, (count + 0.5) / maximum)
