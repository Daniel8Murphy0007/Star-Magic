"""UQFF Star-Magic command-line interface.

After ``pip install uqff``, exposes the ``uqff`` command:

  uqff predict <closure_name>       fetch a UQFF closure value (searches all namespaces)
  uqff search <substring>           search ALL closure / observable / paradox names
  uqff list [--filter SUBSTR]       list dispatch keys (use --all to include all sources)
  uqff status                       print the production status summary
  uqff surfaces                     list public calculate_* surfaces
  uqff version                      print version + key metrics
  uqff gate                         run the fidelity gate (if installed)

Pass ``--json`` to commands for machine-readable output.

The dispatcher searches multiple closure name spaces in order:
  1. PARADOX_TO_CLOSURE (794 paradox dispatch keys)
  2. PARADOX_TO_MILLENNIUM (8 Clay Millennium prize problems + aliases)
  3. calculate_lenr_full sub-keys (10 reactor variants incl. holmlid_D_minus_1)
  4. calculate_nuclear_magic sub-keys (magic numbers, BE/A, deuteron, alpha)
  5. Bucket observable names from 9 bucket surfaces (cosmology, particle physics, etc.)
"""
import argparse
import json
import sys
from typing import Any, Iterable

import uqff_pure_calculator as u


_VERSION = "5.28.0"


# ---------------------------------------------------------------------------
# Multi-namespace closure registry
# ---------------------------------------------------------------------------

_BUCKET_SURFACES = (
    "calculate_cosmology",
    "calculate_particle_physics",
    "calculate_gw_events",
    "calculate_agn_jet",
    "calculate_astrophysics",
    "calculate_high_energy_astro",
    "calculate_qgp",
    "calculate_higgs_precision",
    "calculate_bsm_constraints",
)


def _all_paradox_keys() -> list[str]:
    return sorted(u.PARADOX_TO_CLOSURE.keys())


def _all_millennium_keys() -> list[str]:
    keys = set()
    try:
        keys.update(getattr(u, "PARADOX_TO_MILLENNIUM", {}).keys())
    except Exception:
        pass
    keys.update([
        "yang_mills_mass_gap", "riemann_hypothesis", "navier_stokes_smoothness",
        "hodge_conjecture", "poincare_conjecture", "p_vs_np",
        "bsd_conjecture", "info_paradox",
        "yang_mills", "riemann", "navier_stokes", "hodge", "poincare",
        "pvsnp", "bsd", "page_curve", "black_hole_info", "hawking_info",
    ])
    return sorted(keys)


def _all_lenr_keys() -> list[str]:
    try:
        report = u.calculate_lenr_full({})["value"]
        if isinstance(report, dict):
            return sorted(report.keys())
    except Exception:
        pass
    return []


def _all_nuclear_keys() -> list[str]:
    try:
        report = u.calculate_nuclear_magic({})["value"]
        if isinstance(report, dict):
            return sorted(report.keys())
    except Exception:
        pass
    return []


def _all_bucket_observables() -> dict[str, list[str]]:
    """Returns {surface_name: [observable_names...]}."""
    out: dict[str, list[str]] = {}
    for surf in _BUCKET_SURFACES:
        fn = getattr(u, surf, None)
        if fn is None:
            continue
        try:
            value = fn({}).get("value", {})
            obs_list = value.get("observables", []) if isinstance(value, dict) else []
            names = [o.get("observable", "") for o in obs_list if isinstance(o, dict)]
            names = [n for n in names if n]
            if names:
                out[surf] = names
        except Exception:
            continue
    return out


def _dump(obj: Any) -> str:
    return json.dumps(obj, indent=2, default=lambda v: repr(v))


# ---------------------------------------------------------------------------
# Predict (with multi-source fallback)
# ---------------------------------------------------------------------------

def _try_paradox(name: str):
    """Try the standard PARADOX_TO_CLOSURE / PARADOX_TO_MILLENNIUM dispatcher."""
    result = u.calculate_paradox({"paradox": name})
    value = result.get("value") if isinstance(result, dict) else None
    return value


def _try_lenr_full(name: str):
    """Try calculate_lenr_full sub-keys (holmlid_D_minus_1, parkhomov_NiH, ...).
    Case-insensitive lookup because sub-keys have mixed case."""
    try:
        report = u.calculate_lenr_full({})["value"]
        if isinstance(report, dict):
            name_lc = name.lower()
            for k, v in report.items():
                if k.lower() == name_lc:
                    return {k: v}
    except Exception:
        pass
    return None


def _try_nuclear(name: str):
    """Try calculate_nuclear_magic sub-keys (magic_numbers, be_per_a_peak_MeV, ...).
    Case-insensitive lookup because sub-keys have mixed case."""
    try:
        report = u.calculate_nuclear_magic({})["value"]
        if isinstance(report, dict):
            name_lc = name.lower()
            for k, v in report.items():
                if k.lower() == name_lc:
                    return {k: v}
    except Exception:
        pass
    return None


def _try_bucket_observable(name: str):
    """Search every bucket surface's observables list for a matching name."""
    name_lc = name.lower().strip()
    for surf in _BUCKET_SURFACES:
        fn = getattr(u, surf, None)
        if fn is None:
            continue
        try:
            value = fn({}).get("value", {})
            obs_list = value.get("observables", []) if isinstance(value, dict) else []
            for o in obs_list:
                if isinstance(o, dict):
                    label = o.get("observable", "")
                    if label and (name_lc == label.lower() or name_lc in label.lower()):
                        return {
                            "source_surface": surf,
                            "observable": label,
                            **{k: o.get(k) for k in ("uqff_derived", "anchor", "residual_pct")},
                        }
        except Exception:
            continue
    return None


def _cmd_predict(args: argparse.Namespace) -> int:
    name = args.name.lower().strip()

    # Walk through each namespace in order.
    # Order matters: paradox keys are lowercase-canonical; LENR/nuclear sub-keys
    # have mixed case (holmlid_D_minus_1, magic_numbers, etc.).
    for source_name, fn in [
        ("PARADOX_TO_CLOSURE", _try_paradox),
        ("calculate_lenr_full", _try_lenr_full),
        ("calculate_nuclear_magic", _try_nuclear),
        ("bucket_observables", _try_bucket_observable),
    ]:
        value = fn(name)
        if value is not None:
            if args.json:
                print(_dump({"source": source_name, "name": name, "value": value}))
            else:
                print(f"closure: {name}  (source: {source_name})")
                print(_dump(value))
            return 0

    print(f"ERROR: closure '{name}' not found in any namespace.", file=sys.stderr)
    print("Hint: use `uqff search <substr>` to find candidates across all sources.",
          file=sys.stderr)
    return 1


def _cmd_search(args: argparse.Namespace) -> int:
    needle = args.substr.lower().strip()
    if not needle:
        print("ERROR: search substring required", file=sys.stderr)
        return 1

    results = {}

    paradox_hits = [k for k in _all_paradox_keys() if needle in k.lower()]
    if paradox_hits:
        results["PARADOX_TO_CLOSURE"] = paradox_hits

    millennium_hits = [k for k in _all_millennium_keys() if needle in k.lower()]
    if millennium_hits:
        results["PARADOX_TO_MILLENNIUM"] = millennium_hits

    lenr_hits = [k for k in _all_lenr_keys() if needle in k.lower()]
    if lenr_hits:
        results["calculate_lenr_full"] = lenr_hits

    nuclear_hits = [k for k in _all_nuclear_keys() if needle in k.lower()]
    if nuclear_hits:
        results["calculate_nuclear_magic"] = nuclear_hits

    bucket_obs = _all_bucket_observables()
    for surf, names in bucket_obs.items():
        hits = [n for n in names if needle in n.lower()]
        if hits:
            results[surf] = hits

    total = sum(len(v) for v in results.values())
    if args.json:
        print(_dump({"substring": needle, "total_matches": total, "by_source": results}))
        return 0 if total else 1

    if total == 0:
        print(f"No matches for '{needle}' across any closure namespace.")
        print(f"Searched: {len(_all_paradox_keys())} paradox + {len(_all_millennium_keys())} millennium "
              f"+ {len(_all_lenr_keys())} LENR + {len(_all_nuclear_keys())} nuclear "
              f"+ {sum(len(v) for v in bucket_obs.values())} bucket observables")
        return 1

    print(f"=== Matches for '{needle}' across {len(results)} namespace(s), {total} total hit(s) ===\n")
    for source, hits in results.items():
        print(f"[{source}] {len(hits)} match(es):")
        for h in hits:
            print(f"  {h}")
        print()

    print(f"To inspect any match: uqff predict <name>")
    return 0


def _cmd_list(args: argparse.Namespace) -> int:
    all_groups = {}

    if args.all:
        all_groups["PARADOX_TO_CLOSURE"] = _all_paradox_keys()
        all_groups["PARADOX_TO_MILLENNIUM"] = _all_millennium_keys()
        all_groups["calculate_lenr_full"] = _all_lenr_keys()
        all_groups["calculate_nuclear_magic"] = _all_nuclear_keys()
        all_groups.update(_all_bucket_observables())
    else:
        all_groups["PARADOX_TO_CLOSURE"] = _all_paradox_keys()

    if args.filter:
        flt = args.filter.lower()
        for k in list(all_groups.keys()):
            all_groups[k] = [n for n in all_groups[k] if flt in n.lower()]
            if not all_groups[k]:
                del all_groups[k]

    if args.json:
        print(_dump(all_groups))
        return 0

    total = sum(len(v) for v in all_groups.values())
    if not all_groups:
        flt_msg = f" matching filter '{args.filter}'" if args.filter else ""
        scope = "ALL namespaces" if args.all else "PARADOX_TO_CLOSURE"
        print(f"No closures in {scope}{flt_msg}.")
        if not args.all:
            print("Hint: try `uqff list --all` to search across all namespaces.")
        return 1

    for source, hits in all_groups.items():
        print(f"\n[{source}] {len(hits)}:")
        for h in hits:
            print(f"  {h}")
    print(f"\nTotal: {total} name(s) across {len(all_groups)} namespace(s)")
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
        print("UQFF Star-Magic - Production Status Summary")
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
    print(f"  closures (PARADOX_TO_CLOSURE): {closures}")
    print(f"  truly_independent_primitives: {primitives}")
    print(f"  derivative_primitives: {derivative}")
    print(f"  python: {sys.version.split()[0]}")
    print(f"\nTo see ALL closure namespaces:  uqff list --all")
    print(f"To search across all sources:   uqff search <substring>")
    return 0


def _cmd_gate(args: argparse.Namespace) -> int:
    import importlib.util
    spec = importlib.util.find_spec("uqff_fidelity_tests")
    if spec is None or spec.origin is None:
        print("ERROR: uqff_fidelity_tests not found.", file=sys.stderr)
        return 1
    import subprocess
    return subprocess.call([sys.executable, spec.origin])


def _cmd_serve(args):
    """Launch the REST API server (requires uqff[api] extras)."""
    try:
        import uqff_api
    except ImportError:
        print("ERROR: uqff_api requires FastAPI. Install with:", file=sys.stderr)
        print("  pip install 'uqff[api]'", file=sys.stderr)
        return 1
    uqff_api.run(host=args.host, port=args.port, reload=args.reload)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="uqff",
        description="UQFF Star-Magic command-line interface")
    parser.add_argument("--version", action="version", version=f"uqff {_VERSION}")

    sub = parser.add_subparsers(dest="cmd", required=True, metavar="COMMAND")

    p_predict = sub.add_parser("predict", help="fetch a UQFF closure value (searches all namespaces)")
    p_predict.add_argument("name", help="closure name (lowercase)")
    p_predict.add_argument("--json", action="store_true", help="emit JSON")
    p_predict.set_defaults(func=_cmd_predict)

    p_search = sub.add_parser("search", help="search ALL namespaces for closure names")
    p_search.add_argument("substr", help="substring (case-insensitive)")
    p_search.add_argument("--json", action="store_true", help="emit JSON")
    p_search.set_defaults(func=_cmd_search)

    p_list = sub.add_parser("list", help="list closure names")
    p_list.add_argument("--filter", help="substring filter")
    p_list.add_argument("--all", action="store_true", help="include all 5 namespaces")
    p_list.add_argument("--json", action="store_true", help="emit JSON")
    p_list.set_defaults(func=_cmd_list)

    p_status = sub.add_parser("status", help="production status summary")
    p_status.add_argument("--json", action="store_true", help="emit JSON")
    p_status.set_defaults(func=_cmd_status)

    p_surf = sub.add_parser("surfaces", help="list public calculate_* surfaces")
    p_surf.add_argument("--json", action="store_true", help="emit JSON")
    p_surf.set_defaults(func=_cmd_surfaces)

    p_ver = sub.add_parser("version", help="print version + key metrics")
    p_ver.set_defaults(func=_cmd_version)

    p_gate = sub.add_parser("gate", help="run the fidelity gate")
    p_gate.set_defaults(func=_cmd_gate)

    p_serve = sub.add_parser("serve",
        help="launch the REST API server (requires `pip install \'uqff[api]\'`)")
    p_serve.add_argument("--host", default="127.0.0.1")
    p_serve.add_argument("--port", type=int, default=8000)
    p_serve.add_argument("--reload", action="store_true",
        help="auto-reload on code change")
    p_serve.set_defaults(func=_cmd_serve)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
