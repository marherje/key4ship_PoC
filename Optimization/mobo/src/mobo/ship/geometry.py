"""Bridge to the existing DD4hep variant machinery, and the derived limits.

Two jobs.

**1. Knowing what fits.** Both detectors carry their layer gap on `auto`, so
each one's pitch is `dim_z / NLayers` and the layers come out equidistant across
the whole envelope. That is what lets the layer counts be searched *directly*,
as integers: a count is then bounded only by the rigid content of a layer, not
by another quantity being optimized at the same time. It is also why there is no
`SiTarget_spacing` parameter — a SiTarget layer *is* the W-to-W pitch, so the
spacing was never independent of the count, and all it did on its own was decide
how much dead air each layer carried.

There is no reparametrization left to do here. Every search parameter is a
template constant, passed straight through. The one quantity whose budget really
is derived — where the Y plane sits in the air a layer has left once the pitch is
fixed — is a fraction *in the template itself*
(`SiTarget_XY_plane_gap_frac`, 0 against the X plane, 1 against the next
absorber), so the interpolation happens where the budget is defined and this
module never has to know the formula. That also retired `_length_within`: a
fraction is written to the XML at full precision and turned into a length by the
template, so there is no longer a length rendered at six significant digits that
can round *up* past the maximum it is about to be checked against.

**2. Not duplicating the template.** The `*_max` are computed by rendering the
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


# NOTE: there is deliberately no table of derived parameters here any more.
# `sipad_fill`, `sitarget_fill` and `xy_gap_frac` all used to live in one, each
# mapping a unit-interval search coordinate onto a length. The first two went
# when both layer gaps moved to `auto` (the pitch became dim_z/NLayers, so the
# counts stopped being bounded by a quantity under optimization and could be
# searched directly as integers); the last one went into the template, as
# `SiTarget_XY_plane_gap_frac`. Every search parameter is now a template
# constant and `constants_for` is a pass-through.


@dataclass(frozen=True)
class Limits:
    """The derived upper bounds: what fits, and what the XY fraction is *of*."""

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

        The probe zeroes `SiTarget_XY_plane_gap_frac` and puts both `auto` gaps
        on their floor, so what comes back is the layer stripped down to its
        rigid content — the bound must not depend on the very quantities it
        bounds.

        `SiPad_NLayers` and `SiTarget_NLayers` are *not* probed away, unlike in
        the fill-fraction days when they were derived: they are direct search
        parameters now, and `SiTarget_XY_plane_gap_max` legitimately depends on
        the layer count (the pitch is `dim_z / NLayers`), so it has to be
        resolved at the counts actually asked for.
        """
        constants = self.direct_constants(params)
        probe = dict(constants, SiTarget_XY_plane_gap_frac=0)
        auto_keys = [
            k
            for k in ("SiPad_layer_gap", "SiTarget_layer_gap")
            if str(constants.get(k, "")).strip().lower() == self.config.AUTO
        ]
        for key in auto_keys:
            probe[key] = self.config.AUTO_MIN
        c = self.resolve(probe)

        def nlayers_max(gap_key: str, thick_key: str, env_key: str, max_key: str) -> int:
            """How many layers fit, honouring `auto` when it is in play."""
            if gap_key not in auto_keys:
                return int(c[max_key])
            # `auto` sizes the air slice so that N layers span dim_z exactly, so
            # the binding constraint is not floor(dim_z/thickness) but the point
            # where the slack itself would have to go negative:
            #     slack = dim_z/N - rigid - AUTO_EPS >= AUTO_MIN
            # (`rigid` is the layer without its slack; the probe injected
            # AUTO_MIN as the gap, so take it back out.)
            rigid = c[thick_key] - self.config.AUTO_MIN
            budget = rigid + self.config.AUTO_EPS + self.config.AUTO_MIN
            return int(math.floor(c[env_key] / budget)) if budget > 0 else 0

        return Limits(
            sipad_nlayers_max=max(
                0,
                nlayers_max(
                    "SiPad_layer_gap",
                    "SiPad_LayerThickness",
                    "SiPad_dim_z",
                    "SiPad_NLayers_max",
                ),
            ),
            sitarget_nlayers_max=max(
                0,
                nlayers_max(
                    "SiTarget_layer_gap",
                    "SiTarget_LayerThickness",
                    "SiTarget_dim_z",
                    "SiTarget_NLayers_max",
                ),
            ),
            sitarget_xy_gap_max=float(c["SiTarget_XY_plane_gap_max"]),
        )

    # ── parameters -> constants ──────────────────────────────────────────────

    def direct_constants(self, params: dict[str, Any]) -> dict[str, Any]:
        """Base constants overridden by the parameters that *are* constants."""
        constants = dict(self.base_constants)
        unknown = [k for k in params if k not in constants]
        if unknown:
            # Same typo guard as make_variants.py: a parameter that matches
            # nothing is a bug in the config, not something to ignore.
            raise RuntimeError(
                "parameter(s) that are not a template constant: "
                + ", ".join(sorted(unknown))
            )
        for key, value in params.items():
            if key in constants:
                constants[key] = value
        return constants

    def constants_for(self, params: dict[str, Any]) -> dict[str, Any]:
        """Search-space parameters -> the full constants dict `config.py` wants.

        A pass-through, and that is the point: every search parameter is a
        template constant, including `SiTarget_XY_plane_gap_frac`, so nothing
        here can disagree with the template about what a parameter means.
        """
        return self.direct_constants(params)

    def baseline_params(self) -> dict[str, Any]:
        """The template's own design, expressed in search-space coordinates.

        Now the identity. The layer counts are search parameters themselves
        (`SearchSpace` rounds them back to the same integers) and the XY gap
        fraction is one too, so there is nothing left to invert — which is the
        strongest form the "trial 0 reproduces the baseline" claim can take.
        """
        return dict(self.base_constants)

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
