#!/usr/bin/env python3
"""Generate several compact files at once from SND_compact_template.xml.

Each entry of the variants YAML gives a name and the parameters that differ
from the base parameters file; everything else is inherited. One compact file
per entry is written to variants/<name>.xml (see config.py, which does the
rendering and validation).

Variants YAML format:

    base: parameters_template.yaml   # optional, this is the default
    variants:
      thin_W:      { SiPad_WThickness: 5 }
      short_sipad: { SiPad_dim_z: 200, SiPad_NLayers: 12 }

Usage (from Optimization/Simulation/geometry/):
    python3 make_variants.py variants.yaml
    python3 make_variants.py variants.yaml --dry-run
    python3 make_variants.py variants.yaml --only thin_W --out-dir /tmp/geo

An invalid variant is reported and skipped; the others are still written, and
the exit code is 1 if any of them failed.
"""

import argparse
import sys
from pathlib import Path

import yaml

import config

HERE = Path(__file__).resolve().parent


def load_variants(path):
    """Read the variants YAML -> (base params path, {name: overrides})."""
    spec = yaml.safe_load(Path(path).read_text()) or {}
    variants = spec.get("variants")
    if not isinstance(variants, dict) or not variants:
        raise RuntimeError(f"{path}: missing a non-empty `variants:` mapping")
    base = spec.get("base", str(config.PARAMS))
    base = Path(base)
    if not base.is_absolute():
        base = Path(path).resolve().parent / base
    return base, variants


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("spec", help="variants YAML")
    ap.add_argument("--out-dir", default=str(config.VARIANTS),
                    help="output directory (default: variants/)")
    ap.add_argument("--only", action="append", default=None, metavar="NAME",
                    help="only this variant (repeatable)")
    ap.add_argument("--dry-run", action="store_true", help="validate, write nothing")
    args = ap.parse_args()

    base_path, variants = load_variants(args.spec)
    _, base_constants = config.load_params(base_path)
    out_dir = Path(args.out_dir)

    if args.only:
        unknown = sorted(set(args.only) - set(variants))
        if unknown:
            print("Variantes desconocidas: %s" % ", ".join(unknown))
            return 1
        variants = {k: v for k, v in variants.items() if k in args.only}

    print("Base: %s" % base_path)
    print("Salida: %s%s" % (out_dir, " (dry-run)" if args.dry_run else ""))

    results = []
    for name, overrides in variants.items():
        overrides = overrides or {}
        print("\n=== %s ===" % name)
        unknown = sorted(set(overrides) - set(base_constants))
        if unknown:
            msg = "parametro(s) inexistente(s) en la base: " + ", ".join(unknown)
            print("  !! ERROR:", msg)
            results.append((name, False, msg))
            continue

        constants = dict(base_constants)
        constants.update(overrides)
        for k in sorted(overrides):
            print("  %-24s %s" % (k, config.format_value(k, overrides[k])))

        try:
            ok, lines, errors = config.write_variant(
                constants, out_dir / f"{name}.xml", dry_run=args.dry_run)
        except Exception as exc:                       # keep going on the rest
            print("  !! ERROR:", exc)
            results.append((name, False, str(exc)))
            continue

        print("\n".join(lines))
        if not ok:
            for e in errors:
                print("  !! ERROR:", e)
            results.append((name, False, "; ".join(errors)))
        else:
            results.append((name, True, "" if args.dry_run
                            else str(out_dir / f"{name}.xml")))

    print("\nResumen:")
    for name, ok, info in results:
        print("  %-20s %-6s %s" % (name, "OK" if ok else "ERROR", info))
    failed = sum(1 for _, ok, _ in results if not ok)
    if failed:
        print("%d variante(s) fallaron." % failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
