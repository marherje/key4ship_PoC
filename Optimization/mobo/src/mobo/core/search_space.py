"""The unit cube [0,1]^d <-> physical parameters map.

The optimizer only ever sees the unit cube: BoTorch's defaults (Normalize,
Sobol, `optimize_acqf` bounds) are all written for it, and a search space whose
dimensions differ by orders of magnitude would otherwise need per-dimension
care in three different places.

Two asymmetries are deliberate:

* integers round on the way *out* only (unit -> physical). The GP keeps a
  continuous relaxation, which is what makes it able to interpolate at all;
* `fixed` parameters are handed to the evaluator but never optimized. Freeing
  one later is a one-line config change, not a code change.

Deliberately torch-free: this module is pure arithmetic, and the tests that
hammer it should not pay for importing torch.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any

from .types import ParameterSpec


class SearchSpace:
    """A list of ParameterSpecs plus the parameters held constant."""

    def __init__(
        self,
        specs: Iterable[ParameterSpec],
        fixed: dict[str, Any] | None = None,
    ) -> None:
        self.specs: list[ParameterSpec] = list(specs)
        if not self.specs:
            raise ValueError("a search space needs at least one free parameter")
        self.fixed: dict[str, Any] = dict(fixed or {})

        names = [s.name for s in self.specs]
        dupes = sorted({n for n in names if names.count(n) > 1})
        if dupes:
            raise ValueError(f"duplicate parameter name(s): {', '.join(dupes)}")
        clash = sorted(set(names) & set(self.fixed))
        if clash:
            raise ValueError(f"parameter(s) both free and fixed: {', '.join(clash)}")

    # ── shape ────────────────────────────────────────────────────────────────

    @property
    def dim(self) -> int:
        return len(self.specs)

    @property
    def names(self) -> list[str]:
        return [s.name for s in self.specs]

    def __len__(self) -> int:
        return self.dim

    def __repr__(self) -> str:
        return "SearchSpace({})".format(
            ", ".join(
                "{}:{}[{:g},{:g}]{}".format(
                    s.name, s.kind, s.low, s.high, " log" if s.log else ""
                )
                for s in self.specs
            )
        )

    # ── the map ──────────────────────────────────────────────────────────────

    @staticmethod
    def _to_physical(spec: ParameterSpec, u: float) -> float | int:
        u = min(1.0, max(0.0, float(u)))
        if spec.log:
            value = math.exp(
                math.log(spec.low) + u * (math.log(spec.high) - math.log(spec.low))
            )
        else:
            value = spec.low + u * (spec.high - spec.low)
        if spec.kind == "int":
            return int(min(spec.high, max(spec.low, round(value))))
        # The endpoints matter (they are the corners the acquisition likes to
        # sit on) and floating point can push them a hair outside.
        return float(min(spec.high, max(spec.low, value)))

    @staticmethod
    def _to_unit(spec: ParameterSpec, value: float) -> float:
        if spec.log:
            u = (math.log(float(value)) - math.log(spec.low)) / (
                math.log(spec.high) - math.log(spec.low)
            )
        else:
            u = (float(value) - spec.low) / (spec.high - spec.low)
        return min(1.0, max(0.0, u))

    def to_params(self, unit_x: Sequence[float]) -> dict[str, Any]:
        """Unit cube point -> the full parameter dict (free + fixed)."""
        if len(unit_x) != self.dim:
            raise ValueError(f"expected {self.dim} coordinates, got {len(unit_x)}")
        params: dict[str, Any] = {
            spec.name: self._to_physical(spec, u)
            for spec, u in zip(self.specs, unit_x, strict=False)
        }
        params.update(self.fixed)
        return params

    def to_unit(self, params: dict[str, Any]) -> list[float]:
        """Parameter dict -> unit cube point. Extra keys (fixed, derived) ignored.

        Not an exact inverse of `to_params` for integer dimensions — it cannot
        be, since rounding is many-to-one — but `to_params(to_unit(p)) == p`
        holds for every attainable p, which is the direction that matters
        (re-loading a trial from the store must reproduce its geometry).
        """
        missing = [s.name for s in self.specs if s.name not in params]
        if missing:
            raise ValueError(f"missing parameter(s): {', '.join(missing)}")
        return [self._to_unit(spec, params[spec.name]) for spec in self.specs]

    def clip(self, unit_x: Sequence[float]) -> list[float]:
        return [min(1.0, max(0.0, float(u))) for u in unit_x]

    # ── config ───────────────────────────────────────────────────────────────

    @classmethod
    def from_config(cls, cfg: Any) -> SearchSpace:
        """Build from the `parameters:`/`fixed:` block of an experiment config.

        Expected shape (OmegaConf or plain dicts):

            parameters:
              SiPad_WThickness: {low: 5, high: 15}
              SiPad_NLayers:    {low: 4, high: 40, kind: int}
              foo:              {low: 1e-3, high: 1, log: true}
            fixed:
              SiPad_frame_gap: 0.1
        """
        raw_params = cfg.get("parameters", {})
        specs = []
        for name, body in dict(raw_params).items():
            body = dict(body)
            specs.append(
                ParameterSpec(
                    name=str(name),
                    low=float(body["low"]),
                    high=float(body["high"]),
                    kind=str(body.get("kind", "float")),
                    log=bool(body.get("log", False)),
                )
            )
        raw_fixed = dict(cfg["fixed"]) if "fixed" in cfg and cfg["fixed"] else {}
        return cls(specs, raw_fixed)
