"""UQFF Star-Magic command-line interface.

After ``pip install uqff``, exposes the ``uqff`` command:

  uqff predict <closure_name>       fetch a UQFF closure value
  uqff list [--filter SUBSTR]       list all available closure names
  uqff status                       print the production status summary
  uqff surfaces                     list public calculate_* surfaces
  uqff version                      print version + key metrics
  uqff gate                         run the fidelity gate (if installed)

Pass ``--json`` to predict/status/surfaces for machine-readable output.
"""
import argparse
import json
import sys
from typing import Any

import uqff_pure_calculator as u


_VERSION = "5.27.1"


def _dump(obj: Any) -> str:
    """JSON-dump that tolerates non-serializable values."""
    return json.dumps(obj, indent=2, default=lambda v: repr(v))


def _cmd_predict(args: argparse.Namespace) -> int:
    name = args.name.lower().strip()
    result = u.calculate_paradox({"paradox": name})
    value = result.get("value") if isinstance(result, dict) else None
    if value is None:
        print(f"ERROR: closure '{name}' not found in dispatcher.", file=sys.stderr)
        print("Hint: use lowercase keys; try `uqff list --filter <substr>` to search.",
              file=sys.stderr)
        return 1
    if args.json:
        print(_dump(value))
    else:
        print(f"closure: {name}")
        print(_dump(value))
    return 0


def _cmd_list(args: argparse.Namespace) -> int:
    keys = sorted(u.PARADOX_TO_CLOSURE.keys())
    if args.filter:
        flt = args.filter.lower()
        keys = [k for k in keys if flt in k.lower()]
    if not keys:
        print(f"No closures match filter '{args.filter}'.", file=sys.stderr)
        return 1
    if args.json:
        print(_dump(keys))
    else:
        for k in keys:
            print(k)
        print(f"\nTotal: {len(keys)} closure name(s)")
    return 0


def _cmd_status(args: argparse.Namespace) -> int:
    r = u.calculate_status_report({})
    summary = r.get("value", {}).get("summary", {}) if isinstance(r, dict) else {}
    if not summary:
        print("ERROR: status_report returned no summary.", file=sys.stderr)
        return 1
    if args.json:
        print(_dump(summary))
    else:
        print("UQFF Star-Magic — Production Status Summary")
        print("=" * 50)
        for k, v in summary.items():
            if isinstance(v, dict):
                print(f"\n{k}:")
                for kk, vv in v.items():
                    print(f"  {kk}: {vv}")
            else:
                print(f"{k}: {v}")
    return 0


def _cmd_surfaces(args: argparse.Namespace) -> int:
    publics = sorted(n for n in dir(u) if n.startswith("calculate_"))
    if args.json:
        print(_dump(publics))
    else:
        for s in publics:
            print(s)
        print(f"\nTotal: {len(publics)} public calculate_* surfaces")
    return 0


def _cmd_version(args: argparse.Namespace) -> int:
    summary = {}
    try:
        r = u.calculate_status_report({})
        summary = r.get("value", {}).get("summary", {}) if isinstance(r, dict) else {}
    except Exception:
        pass
    closures = summary.get("total_closures", "unknown")
    primitives = summary.get("truly_independent_primitives", "unknown")
    derivative = summary.get("derivative_primitives", "unknown")
    print(f"uqff {_VERSION}")
    print(f"  closures: {closures}")
    print(f"  truly_independent_primitives: {primitives}")
    print(f"  derivative_primitives: {derivative}")
    print(f"  python: {sys.version.split()[0]}")
    return 0


def _cmd_gate(args: argparse.Namespace) -> int:
    import importlib.util
    spec = importlib.util.find_spec("uqff_fidelity_tests")
    if spec is None or spec.origin is None:
        print("ERROR: uqff_fidelity_tests not found in install.", file=sys.stderr)
        print("The fidelity gate ships as a separate file; install from source "
              "via `pip install -e .` to get it.", file=sys.stderr)
        return 1
    import subprocess
    return subprocess.call([sys.executable, spec.origin])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="uqff",
        description="UQFF Star-Magic command-line interface "
                    "(see https://github.com/Daniel8Murphy0007/Star-Magic)",
    )
    parser.add_argument("--version", action="version", version=f"uqff {_VERSION}")

    sub = parser.add_subparsers(dest="cmd", required=True, metavar="COMMAND")

    p_predict = sub.add_parser("predict", help="fetch a UQFF closure value by name")
    p_predict.add_argument("name", help="closure name (lowercase, e.g. 'hubble_tension')")
    p_predict.add_argument("--json", action="store_true", help="emit JSON")
    p_predict.set_defaults(func=_cmd_predict)

    p_list = sub.add_parser("list", help="list available closure names")
    p_list.add_argument("--filter", help="substring filter (case-insensitive)")
    p_list.add_argument("--json", action="store_true", help="emit JSON array")
    p_list.set_defaults(func=_cmd_list)

    p_status = sub.add_parser("status", help="print the production status summary")
    p_status.add_argument("--json", action="store_true", help="emit JSON")
    p_status.set_defaults(func=_cmd_status)

    p_surf = sub.add_parser("surfaces", help="list public calculate_* surfaces")
    p_surf.add_argument("--json", action="store_true", help="emit JSON")
    p_surf.set_defaults(func=_cmd_surfaces)

    p_ver = sub.add_parser("version", help="print version + key metrics")
    p_ver.set_defaults(func=_cmd_version)

    p_gate = sub.add_parser("gate", help="run the fidelity gate (source install only)")
    p_gate.set_defaults(func=_cmd_gate)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
