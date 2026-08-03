#!/usr/bin/env python3
"""Render SND_compact_template.xml into a compact file (geometry variants).

The <!-- Input parameters --> block of SND_compact_template.xml holds the free
design parameters as XnameX placeholders. This script substitutes them with the
values of parameters_template.yaml and writes the result to
variants/<name>.xml. SND_compact.xml (the committed baseline) is never
overwritten — each run produces a new compact file, so several geometries can
coexist. See make_variants.py to generate a whole set in one go.

The shared, committed geometry definitions (elements/materials/detector XMLs,
parse_geometry.py and the SND_compact.xml baseline) stay in the base project,
under simulation/geometry — COMMON below. Only the free design parameters and
the template live here.

Usage (from Optimization/Simulation/geometry/):
    python3 config.py                    # -> variants/<name>.xml
    python3 config.py --dry-run          # validate, write nothing
    python3 config.py --name thin_W      # override the output name
    python3 config.py --output foo.xml   # override the full output path
    python3 config.py --show [xml]       # only print a geometry summary

Safety: the rendered XML is validated on a temporary copy BEFORE being written
(constants resolvable, NLayers/gaps within their derived *_max limits, detector
envelopes ordered along z and non-overlapping).
"""

import argparse
import os
import re
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
# Common (shared, committed) geometry definitions: elements/materials/detector
# XMLs, parse_geometry.py and the SND_compact.xml baseline live here.
COMMON = (HERE / ".." / ".." / ".." / "simulation" / "geometry").resolve()
TEMPLATE = HERE / "SND_compact_template.xml"
PARAMS = HERE / "parameters_template.yaml"
VARIANTS = HERE / "variants"
# The committed baseline: read-only for this script, and the default of --show.
COMPACT = COMMON / "SND_compact.xml"

PLACEHOLDER_RE = re.compile(r"X([A-Za-z_][A-Za-z_0-9]*)X")
REF_RE = re.compile(r'(ref=")([^"/][^"]*)(")')

# Parameters that may be given as `auto` in the YAML: a slack thickness sized so
# that NLayers layers span dim_z exactly, i.e. the layers come out equidistant.
# {parameter: (thickness constant, envelope constant, layer count constant)}
AUTO = "auto"
AUTO_SLACK = {
    "SiPad_layer_gap": ("SiPad_LayerThickness", "SiPad_dim_z", "SiPad_NLayers"),
}
# Floor an `auto` parameter never goes below: the slice has to stay in the
# geometry (the index of the sensitive slice behind it depends on it) and Geant4
# rejects a box of null extent, so a slack slice that is not needed shrinks to a
# thickness no particle will ever step in rather than to zero.
AUTO_MIN = 1e-6          # mm (1 nm; Geant4's own tolerance is ~1e-9 mm)
# Shaved off the computed slack so the layers always come out SHORTER than the
# envelope. Without it, floating-point noise on a tight fit can make the
# template's floor(dim_z/LayerThickness) return NLayers-1 and the
# NLayers <= NLayers_max check fail on a geometry that does fit. Kept two orders
# of magnitude above AUTO_MIN so the margin absorbs the minimum gap as well.
AUTO_EPS = 1e-4          # mm

# Tolerance of the envelope-overlap check. Adjacent detectors are meant to
# touch: SiTarget ends exactly where SiPad begins, since SiPad spans
# [-SiPad_dim_z, 0] and SiTarget is placed right upstream of it. Those two z
# values are algebraically equal but computed by different expressions, so for a
# SiPad_dim_z that is not a round number they differ in the last bits and a
# strict `>` reports an overlap of ~1e-13 mm on a perfectly valid geometry.
# 1e-6 mm (1 nm) is below anything physical here and matches the scale already
# used for the minimum slice thickness.
OVERLAP_TOL = 1e-6       # mm

# (value, derived upper limit) pairs checked on the resolved geometry.
LIMIT_PAIRS = [
    ("SiPad_NLayers", "SiPad_NLayers_max"),
    ("SiTarget_NLayers", "SiTarget_NLayers_max"),
    ("SiTarget_XY_plane_gap", "SiTarget_XY_plane_gap_max"),
]


# ── rendering ───────────────────────────────────────────────────────────────

def discover_placeholders(text):
    """Names of the XnameX placeholders present in the template."""
    return set(PLACEHOLDER_RE.findall(text))


def format_value(name, val):
    """YAML value -> XML value string."""
    if isinstance(val, str):
        return val                      # verbatim: allows DD4hep expressions
    if name.endswith("NLayers"):
        return str(int(val))
    return "%g*mm" % float(val)


def _fix_refs(text, out_dir):
    """Point the <include>/<gdmlFile> refs of the template at COMMON.

    The template names its includes bare (ref="elements.xml") because they used
    to sit next to it; they now live in COMMON, and DD4hep resolves a ref
    relative to the file doing the including. So every ref that COMMON can
    satisfy is rewritten relative to the directory the variant is written to
    (from variants/ that comes out as ../../../../simulation/geometry/…).
    """
    def sub(m):
        target = COMMON / m.group(2)
        if not target.exists():
            return m.group(0)
        return m.group(1) + os.path.relpath(target, out_dir) + m.group(3)

    return REF_RE.sub(sub, text)


def render(template_text, constants, out_dir):
    """Substitute the placeholders and fix the relative refs.

    Raises RuntimeError if the placeholder set and the constants keys do not
    match exactly, so adding a placeholder to the template forces declaring it
    and a mistyped key never passes unnoticed.
    """
    placeholders = discover_placeholders(template_text)
    keys = set(constants)
    missing = sorted(placeholders - keys)
    extra = sorted(keys - placeholders)
    problems = []
    if missing:
        problems.append("no value given for template placeholder(s): "
                        + ", ".join(missing))
    if extra:
        problems.append("parameter(s) with no XnameX placeholder in %s: %s"
                        % (TEMPLATE.name, ", ".join(extra)))
    if problems:
        raise RuntimeError("; ".join(problems))

    text = PLACEHOLDER_RE.sub(
        lambda m: format_value(m.group(1), constants[m.group(1)]), template_text)
    left = discover_placeholders(text)
    if left:
        raise RuntimeError("unsubstituted placeholder(s): " + ", ".join(sorted(left)))
    return _fix_refs(text, out_dir)


# ── automatic slack ─────────────────────────────────────────────────────────

def _resolve_constants(text, out_dir):
    """<define> constants of a rendered compact, fully evaluated."""
    with tempfile.NamedTemporaryFile("w", suffix=".xml", delete=False,
                                     dir=str(out_dir if Path(out_dir).exists()
                                             else HERE)) as tmp:
        tmp.write(text)
        tmp_path = tmp.name
    try:
        return _geometry(tmp_path)._constants
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def resolve_auto(constants, template_text, out_dir):
    """Replace the `auto` parameters by the slack that fills their envelope.

    An `auto` parameter is an air slice whose only job is to take up whatever
    the rest of the layer leaves over, so that the layers tile the envelope
    evenly instead of being packed at its centre (the plugin stacks them
    contiguously; see SiPadDetector.cpp). Its value cannot be written down in
    the YAML because it depends on the slice thicknesses, which live in the
    template — so it is measured: the template is rendered once with the
    parameter at AUTO_MIN, the resolved layer thickness is read back, and the
    leftover is spread over the layers.

    Returns (constants, notes); `constants` is unchanged if nothing is `auto`.
    """
    asked = {k for k, v in constants.items()
             if str(v).strip().lower() == AUTO}
    unsupported = sorted(asked - set(AUTO_SLACK))
    if unsupported:
        raise RuntimeError("`auto` no esta soportado para: %s (solo para %s)"
                           % (", ".join(unsupported), ", ".join(sorted(AUTO_SLACK))))
    auto = sorted(asked)
    if not auto:
        return constants, []

    constants = dict(constants)
    probe = dict(constants, **{k: AUTO_MIN for k in auto})
    c = _resolve_constants(render(template_text, probe, out_dir), out_dir)

    notes = []
    for key in auto:
        thick_key, env_key, n_key = AUTO_SLACK[key]
        n = int(c[n_key])
        if n < 1:
            raise RuntimeError("%s: cannot size %s, %s = %d" % (key, key, n_key, n))
        # The probe thickness still carries the AUTO_MIN we injected.
        rigid = c[thick_key] - AUTO_MIN
        slack = c[env_key] / n - rigid - AUTO_EPS
        if slack < AUTO_MIN:
            raise RuntimeError(
                "%s = auto: %d capas de %.4f mm no caben en %s = %.4f mm "
                "(faltan %.4f mm); reduce %s o el grosor de la capa."
                % (key, n, rigid, env_key, c[env_key], n * rigid - c[env_key], n_key))
        # Verbatim string: %g would round to 6 significant digits and the layers
        # would no longer land on an exact multiple of the envelope.
        constants[key] = "%.10g*mm" % slack
        notes.append("  %-24s %s  (auto: %d capas equidistantes de %.4f mm en %s)"
                     % (key, constants[key], n, rigid + slack, env_key))
    return constants, notes


# ── validation ──────────────────────────────────────────────────────────────

def _geometry(compact_path):
    """Parse a compact file with a freshly reloaded parse_geometry."""
    sys.path.insert(0, str(COMMON))
    import importlib
    import parse_geometry
    importlib.reload(parse_geometry)
    return parse_geometry.SNDGeometry(compact_path)


def summarize(compact_path):
    """Resolve the geometry and return (summary lines, error lines)."""
    c = _geometry(compact_path)._constants

    lines, errors = [], []

    for value, limit in LIMIT_PAIRS:
        if value not in c or limit not in c:
            continue
        if c[value] > c[limit]:
            errors.append("%s = %g > %s = %g" % (value, c[value], limit, c[limit]))
        else:
            lines.append("  %-24s %9.3f  (max %.3f, margen %.3f)"
                         % (value, c[value], c[limit], c[limit] - c[value]))

    # (label, z centre constant, full length constant); a detector whose
    # constants are absent (e.g. the EmulsionTarget, dropped from the PoC
    # geometry) is simply skipped.
    envelopes = [
        ("EmulsionTarget", "EmTgt_z_position", "EmTgt_env_z", None),
        ("SiTarget",       "SiTarget_z_position", "SiTarget_dim_z", "SiTarget_NLayers"),
        ("SiPad",          "SiPad_z_position", "SiPad_dim_z", "SiPad_NLayers"),
        ("MTC40",          "MTC40_z_position", "MTC_total_length", "MTC_NLayers"),
        ("MTC50",          "MTC50_z_position", "MTC_total_length", "MTC_NLayers"),
        ("MTC60",          "MTC60_z_position", "MTC_total_length", "MTC_NLayers"),
    ]
    spans = []
    for label, z_key, len_key, n_key in envelopes:
        if z_key not in c or len_key not in c:
            continue
        if n_key in c:
            label = "%s (%d capas)" % (label, int(c[n_key]))
        half = c[len_key] / 2.0
        spans.append((label, c[z_key] - half, c[z_key] + half))
    spans.sort(key=lambda s: s[1])

    lines += ["  %-22s z = [%9.1f, %9.1f] mm" % s for s in spans]
    for (na, _, hia), (nb, lob, _) in zip(spans, spans[1:]):
        if hia > lob + OVERLAP_TOL:
            errors.append(f"SOLAPE: {na} termina en {hia:.1f} y {nb} empieza en {lob:.1f}")
    # NOTE: the default particle gun (multiplePG_base.py) sits at z = -1000 mm,
    # intentionally *inside* SiTarget's span; EmulsionTarget legitimately sits
    # upstream of that, so this is only a warning (printed via `lines`), not a
    # hard error that would block writing the XML.
    if spans and spans[0][1] <= -1000.0:
        lines.append(f"  AVISO: el primer detector empieza en {spans[0][1]:.1f} mm "
                     "<= -1000 mm (posicion por defecto del gun); esto es "
                     "esperado si EmulsionTarget/SiTarget estan aguas arriba "
                     "del gun por diseno.")
    return lines, errors


# ── public API (also used by make_variants.py) ──────────────────────────────

def write_variant(constants, out_path, dry_run=False, template=None):
    """Render, validate and (unless dry_run) write one compact file.

    Returns (ok, lines, errors); nothing is written when errors is non-empty.
    """
    out_path = Path(out_path).resolve()
    if out_path == COMPACT:
        raise RuntimeError(
            "%s is the committed baseline and is never generated; "
            "choose another output name." % COMPACT.name)

    template_text = Path(template or TEMPLATE).read_text()
    out_dir = out_path.parent
    constants, notes = resolve_auto(constants, template_text, out_dir)
    text = render(template_text, constants, out_dir)

    with tempfile.NamedTemporaryFile("w", suffix=".xml", delete=False,
                                     dir=str(out_dir if out_dir.exists() else HERE)) as tmp:
        tmp.write(text)
        tmp_path = tmp.name
    try:
        lines, errors = summarize(tmp_path)
    finally:
        Path(tmp_path).unlink(missing_ok=True)
    lines = notes + lines

    if errors or dry_run:
        return (not errors), lines, errors

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text)
    return True, lines, errors


def load_params(path):
    """Read a parameters YAML -> (name, constants dict)."""
    params = yaml.safe_load(Path(path).read_text()) or {}
    constants = params.get("constants")
    if not isinstance(constants, dict):
        raise RuntimeError(f"{path}: missing a `constants:` mapping")
    return params.get("name") or Path(path).stem, constants


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--params", default=str(PARAMS), help="parameters YAML path")
    ap.add_argument("--name", default=None,
                    help="output variant name (default: `name` in the YAML)")
    ap.add_argument("--output", default=None,
                    help="full output path (overrides --name and the variants/ dir)")
    ap.add_argument("--dry-run", action="store_true", help="validate, write nothing")
    ap.add_argument("--show", nargs="?", const=str(COMPACT), default=None,
                    metavar="XML", help="only print the summary of an existing compact file")
    args = ap.parse_args()

    if args.show:
        lines, errors = summarize(args.show)
        print("Geometria de %s:" % Path(args.show).name)
        print("\n".join(lines))
        for e in errors:
            print("  !!", e)
        return 1 if errors else 0

    name, constants = load_params(args.params)
    if args.name:
        name = args.name
    out_path = Path(args.output) if args.output else VARIANTS / f"{name}.xml"

    print("Variante '%s' -> %s" % (name, out_path))
    if out_path.exists():
        print("  (el fichero ya existe y sera sobrescrito)")
    print("Parametros de entrada:")
    for k in sorted(constants):
        print("  %-24s %s" % (k, format_value(k, constants[k])))

    try:
        ok, lines, errors = write_variant(constants, out_path, dry_run=args.dry_run)
    except RuntimeError as exc:
        print("  !! ERROR:", exc)
        print("No se ha escrito nada.")
        return 1

    print("Geometria resultante:")
    print("\n".join(lines))
    if not ok:
        for e in errors:
            print("  !! ERROR:", e)
        print("No se ha escrito nada.")
        return 1

    if args.dry_run:
        print("(dry-run: nada escrito)")
        return 0

    print("Escrito %s" % out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
