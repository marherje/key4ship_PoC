#!/usr/bin/env python3
"""Apply parameters.yaml to SND_compact.xml (semi-automatic geometry config).

Reads the layer counts and z positions from parameters.yaml and patches the
corresponding <constant> entries of SND_compact.xml. Everything downstream
(ddsim steering files, ACTSGeoSvc, parse_geometry-based job configs) reads
the XML, so nothing else needs to change.

Usage (from simulation/geometry/):
    python3 config.py               # apply parameters.yaml
    python3 config.py --dry-run     # show what would change, don't write
    python3 config.py --show        # only print the current geometry summary

Safety: the new XML is validated on a temporary copy (constants resolvable,
detector envelopes ordered along z and non-overlapping, and the default gun
position upstream of the first detector) BEFORE overwriting SND_compact.xml.
"""

import argparse
import re
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
COMPACT = HERE / "SND_compact.xml"
PARAMS = HERE / "parameters.yaml"

# Canonical derived expressions used when a z_position is set to `auto`
# (these are the original self-stacking definitions of SND_compact.xml).
AUTO_EXPR = {
    "SiTarget_z_position":
        "-1*(Ecal_dim_z/2.0 + env_safety + SiTarget_LayerThickness*SiTarget_NLayers/2.0 +10*mm)",
    "Ecal_z_position": "0*mm",
    "MTC40_z_position": "Ecal_dim_z/2.0 + 100.0*mm + MTC_total_length/2.0",
    "MTC50_z_position": "MTC40_z_position + MTC_total_length + MTC_inter_station_gap",
    "MTC60_z_position": "MTC50_z_position + MTC_total_length + MTC_inter_station_gap",
}


def _zvalue(name, zcfg):
    """YAML z entry -> XML value string ('auto' -> canonical expression)."""
    if zcfg is None or (isinstance(zcfg, str) and zcfg.lower() == "auto"):
        return AUTO_EXPR[name]
    return "%s*mm" % float(zcfg)


def build_constant_map(params):
    """Map <constant> name -> new value string from the YAML parameters."""
    consts = {
        "SiTarget_NLayers": str(int(params["sitarget"]["n_layers"])),
        "Ecal_NLayers":     str(int(params["sipad"]["n_layers"])),
        "MTC_NLayers":      str(int(params["mtc"]["n_layers"])),
        "SiTarget_z_position": _zvalue("SiTarget_z_position",
                                       params["sitarget"].get("z_position")),
        "Ecal_z_position":     _zvalue("Ecal_z_position",
                                       params["sipad"].get("z_position")),
    }
    stations = params["mtc"].get("stations", {})
    for st in ("MTC40", "MTC50", "MTC60"):
        consts[f"{st}_z_position"] = _zvalue(
            f"{st}_z_position", stations.get(st, {}).get("z_position"))
    return consts


def patch_xml(text, consts):
    """Replace the value="..." of each named <constant> in the XML text."""
    changes = []
    for name, value in consts.items():
        pattern = re.compile(
            r'(<constant\s+name="%s"\s+value=")([^"]*)("\s*/>)' % re.escape(name))
        m = pattern.search(text)
        if not m:
            raise RuntimeError(f"Constant '{name}' not found in {COMPACT}")
        if m.group(2) != value:
            changes.append((name, m.group(2), value))
            text = pattern.sub(lambda mo: mo.group(1) + value + mo.group(3),
                               text, count=1)
    return text, changes


def summarize(compact_path):
    """Resolve the geometry and return (summary lines, error lines)."""
    sys.path.insert(0, str(HERE))
    import importlib
    import parse_geometry
    importlib.reload(parse_geometry)
    g = parse_geometry.SNDGeometry(compact_path)
    c = g._constants

    spans = []
    st_half = c["SiTarget_LayerThickness"] * c["SiTarget_NLayers"] / 2.0
    spans.append(("SiTarget (%d capas)" % int(c["SiTarget_NLayers"]),
                  c["SiTarget_z_position"] - st_half,
                  c["SiTarget_z_position"] + st_half))
    spans.append(("SiPad/ECAL (%d capas)" % int(c["Ecal_NLayers"]),
                  c["Ecal_z_position"] - c["Ecal_dim_z"] / 2.0,
                  c["Ecal_z_position"] + c["Ecal_dim_z"] / 2.0))
    for st in ("MTC40", "MTC50", "MTC60"):
        spans.append(("%s (%d capas)" % (st, int(c["MTC_NLayers"])),
                      c[f"{st}_z_position"] - c["MTC_total_length"] / 2.0,
                      c[f"{st}_z_position"] + c["MTC_total_length"] / 2.0))

    lines = ["  %-22s z = [%9.1f, %9.1f] mm" % s for s in spans]
    errors = []
    for (na, _, hia), (nb, lob, _) in zip(spans, spans[1:]):
        if hia > lob:
            errors.append(f"SOLAPE: {na} termina en {hia:.1f} y {nb} empieza en {lob:.1f}")
    if spans[0][1] <= -1000.0:
        errors.append(f"El primer detector empieza en {spans[0][1]:.1f} mm <= -1000 mm "
                      "(posicion por defecto del gun); adelanta el gun o el detector.")
    return lines, errors


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--params", default=str(PARAMS), help="parameters.yaml path")
    ap.add_argument("--dry-run", action="store_true", help="don't write the XML")
    ap.add_argument("--show", action="store_true", help="only print current geometry")
    ap.add_argument("--output", default=None,
                    help="write the patched XML to this file instead of "
                         "overwriting SND_compact.xml (for coexisting "
                         "geometry variants; no .bak is made)")
    args = ap.parse_args()

    if args.show:
        lines, errors = summarize(COMPACT)
        print("Geometria actual (%s):" % COMPACT.name)
        print("\n".join(lines))
        for e in errors:
            print("  !!", e)
        return 0

    params = yaml.safe_load(open(args.params))
    consts = build_constant_map(params)

    target = Path(args.output) if args.output else COMPACT

    original = COMPACT.read_text()
    patched, changes = patch_xml(original, consts)

    if not changes and target == COMPACT:
        print("SND_compact.xml ya coincide con parameters.yaml — nada que hacer.")
        lines, errors = summarize(COMPACT)
        print("\n".join(lines))
        return 0

    print("Cambios:")
    for name, old, new in changes:
        print(f"  {name}: '{old}'  ->  '{new}'")

    # Validate on a temp copy before touching the real file.
    with tempfile.NamedTemporaryFile("w", suffix=".xml", delete=False) as tmp:
        tmp.write(patched)
        tmp_path = tmp.name
    try:
        lines, errors = summarize(tmp_path)
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    print("Geometria resultante:")
    print("\n".join(lines))
    if errors:
        for e in errors:
            print("  !! ERROR:", e)
        print("No se ha escrito SND_compact.xml.")
        return 1

    if args.dry_run:
        print("(dry-run: nada escrito)")
        return 0

    if target == COMPACT:
        shutil.copy(COMPACT, str(COMPACT) + ".bak")
        COMPACT.write_text(patched)
        print(f"Escrito {COMPACT.name} (backup en {COMPACT.name}.bak).")
        print("Recuerda re-simular: la geometria de los datos existentes ya no coincide.")
    else:
        target.write_text(patched)
        print(f"Escrito {target} (variante de geometria; SND_compact.xml intacto).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
