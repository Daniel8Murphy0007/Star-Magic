"""UQFF Star-Magic command-line interface.

After ``pip install uqff``, exposes the ``uqff`` command:

  CLOSURE QUERIES
    uqff predict <closure_name>     fetch a UQFF closure value (all namespaces)
    uqff search <substring>         search ALL closure / observable / paradox names
    uqff list [--filter SUBSTR]     list dispatch keys (use --all for every namespace)

  PROOF CORPUS ACCESS (v5.29.0 — bundled with the package)
    uqff proofs list                list 1,994+ whitepapers
    uqff proofs show <NAME>         print a specific whitepaper
    uqff proofs path                print install dir for bundled proofs
    uqff millennium                 run run_millennium_proofs.py (7 Clay proofs)
    uqff axioms                     print AXIOMS_AND_THEOREMS.md
    uqff manuscript                 print path to the compiled manuscript PDF
    uqff lean                       print location of the Lean 4 scaffold
    uqff atlas                      print CLOSURE_ATLAS.md
    uqff gold-standard              print Gold_Standard_Pure_UQFF.md
    uqff grok-archives              list bundled Grok proof-conversation archives

  STATUS / TESTING
    uqff status                     print production status summary
    uqff surfaces                   list public calculate_* surfaces
    uqff version                    print version + key metrics
    uqff gate                       run the fidelity gate
    uqff serve [--port N]           launch REST API (requires uqff[api])

Pass ``--json`` to query commands for machine-readable output.

Closure-search order:
  1. PARADOX_TO_CLOSURE (794 paradox dispatch keys)
  2. PARADOX_TO_MILLENNIUM (8 Clay Millennium prize problems + aliases)
  3. calculate_lenr_full sub-keys (10 reactor variants incl. holmlid_D_minus_1)
  4. calculate_nuclear_magic sub-keys (magic numbers, BE/A, deuteron, alpha)
  5. Bucket observable names from 9 bucket surfaces
"""
import argparse
import json
import os
import sys
from typing import Any

import uqff_pure_calculator as u


_VERSION = "5.60.0"


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
    result = u.calculate_paradox({"paradox": name})
    value = result.get("value") if isinstance(result, dict) else None
    return value


def _try_lenr_full(name: str):
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


def _try_assimilation_dispatch(name: str):
    """Try the Phase E/F/G dispatch (assimilation_dispatch.DISPATCH) for `name`.
    Case-insensitive lookup since dispatch keys are mixed-case (e.g. LCDM_BAO_...).
    Returns the calculate_analytic_closures qcalcgeom_solve result (decomposed view)
    or None if the observable is not in the dispatch or sympy is unavailable.
    """
    try:
        import assimilation_dispatch as _ad
        import uqff_pure_calculator as _u
    except Exception:
        return None
    canonical = None
    target_lc = (name or "").lower()
    for key in _ad.DISPATCH:
        if key.lower() == target_lc:
            canonical = key
            break
    if canonical is None:
        return None
    try:
        r = _u.calculate_analytic_closures({
            "qcalcgeom_solve": {"observable": canonical, "decompose": True}})
    except Exception:
        return None
    val = r.get("value") if isinstance(r, dict) else None
    if not isinstance(val, dict) or val.get("value") is None:
        return None
    return val


def _cmd_assimilate(args: argparse.Namespace) -> int:
    """Route an observable through the qcalcgeom_solver bus with geometry/numeric controls."""
    try:
        import assimilation_dispatch as _ad
        import uqff_pure_calculator as _u
    except Exception as _e:
        print(f"ERROR: assimilation_dispatch / uqff_pure_calculator import failed: {_e}", file=sys.stderr)
        return 1
    canonical = args.name
    if canonical not in _ad.DISPATCH:
        target_lc = canonical.lower()
        for k in _ad.DISPATCH:
            if k.lower() == target_lc:
                canonical = k
                break
    payload = {"observable": canonical}
    if args.geometry and args.geometry != "auto":
        payload["geometry"] = args.geometry
    if args.numeric and args.numeric != "numerical":
        payload["numeric"] = args.numeric
    if args.decompose:
        payload["decompose"] = True
    res = _u.calculate_analytic_closures({"qcalcgeom_solve": payload})
    val = res.get("value") if isinstance(res, dict) else None
    if val is None:
        print(f"ERROR: observable '{args.name}' not in dispatch or solver bus returned None.",
              file=sys.stderr)
        print(f"Hint: `uqff list --dispatch` to enumerate the 114 known observables.",
              file=sys.stderr)
        return 1
    if args.json:
        print(_dump(val))
    else:
        if isinstance(val, dict):
            print(f"observable: {args.name}")
            for k in ("value", "target", "residual_pct", "geometry_used",
                      "numeric_system", "overdetermination_N",
                      "assimilation_status"):
                if k in val:
                    print(f"  {k}: {val[k]}")
        else:
            print(f"value: {val}")
    return 0


def _list_dispatch_observables():
    """Return sorted list of all observables in assimilation_dispatch.DISPATCH, or [] if unavailable."""
    try:
        import assimilation_dispatch as _ad
        return sorted(_ad.DISPATCH.keys())
    except Exception:
        return []


def _list_dispatch_domains():
    """Return list of dispatch domains, or [] if unavailable."""
    try:
        import assimilation_dispatch as _ad
        return _ad.domains()
    except Exception:
        return []


def _cmd_predict(args: argparse.Namespace) -> int:
    name = args.name.lower().strip()
    for source_name, fn in [
        ("PARADOX_TO_CLOSURE", _try_paradox),
        ("calculate_lenr_full", _try_lenr_full),
        ("calculate_nuclear_magic", _try_nuclear),
        ("bucket_observables", _try_bucket_observable),
        ("assimilation_dispatch", _try_assimilation_dispatch),
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
        return 1
    print(f"=== Matches for '{needle}' across {len(results)} namespace(s), {total} total hit(s) ===\n")
    for source, hits in results.items():
        print(f"[{source}] {len(hits)} match(es):")
        for h in hits:
            print(f"  {h}")
        print()
    print("To inspect any match: uqff predict <name>")
    return 0


def _cmd_list(args: argparse.Namespace) -> int:
    if getattr(args, 'domain', None) or getattr(args, 'dispatch', False):
        try:
            import assimilation_dispatch as _ad
        except Exception as _e:
            print(f'ERROR: assimilation_dispatch unavailable: {_e}', file=sys.stderr)
            return 1
        dom = getattr(args, 'domain', None)
        if dom:
            names = sorted(n for n, r in _ad.DISPATCH.items() if r.get('domain') == dom)
        else:
            names = sorted(_ad.DISPATCH.keys())
        flt = getattr(args, 'filter', None)
        if flt:
            names = [n for n in names if flt.lower() in n.lower()]
        if args.json:
            print(_dump(names))
        else:
            for n in names:
                rec = _ad.DISPATCH[n]
                print(f"{n:<42s}  {rec['domain']:<6s}  owner={rec['owner_geometry']:<10s}  resid={rec.get('residual_pct')}%")
            print(f"\n{len(names)} observable(s).")
        return 0
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
        print("No matching closures.")
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
    whitepapers = summary.get("whitepapers_bundled", "unknown")
    print(f"uqff {_VERSION}")
    print(f"  closures (PARADOX_TO_CLOSURE): {closures}")
    print(f"  truly_independent_primitives: {primitives}")
    print(f"  derivative_primitives: {derivative}")
    print(f"  whitepapers_bundled: {whitepapers}")
    print(f"  python: {sys.version.split()[0]}")
    print("\nProof corpus commands:")
    print("  uqff proofs list      uqff axioms        uqff lean")
    print("  uqff millennium       uqff manuscript    uqff atlas")
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
    try:
        import uqff_api
    except ImportError:
        print("ERROR: uqff_api requires FastAPI. Install with:", file=sys.stderr)
        print("  pip install 'uqff[api]'", file=sys.stderr)
        return 1
    uqff_api.run(host=args.host, port=args.port, reload=args.reload)
    return 0


# ---------------------------------------------------------------------------
# v5.29.0 — Proof corpus access commands
# ---------------------------------------------------------------------------

def _pkg_install_dir() -> str:
    return os.path.dirname(os.path.abspath(u.__file__))


def _find_corpus_path(*candidates: str) -> str | None:
    """Search for a file/dir across plausible install locations."""
    pkg_dir = _pkg_install_dir()
    search_roots = [
        pkg_dir,
        os.path.join(pkg_dir, ".."),
        os.path.join(sys.prefix, "share", "uqff"),
        os.path.join(sys.prefix, "share", "uqff", "proofs"),
        os.path.join(sys.prefix, "share", "uqff", "gold_standard"),
        os.path.join(sys.prefix, "share", "uqff", "grok_archives"),
        os.path.join(sys.prefix, "share", "uqff", "submission_prep"),
        os.path.join(sys.prefix, "share", "uqff", "audit_reports"),
        os.getcwd(),
    ]
    for root in search_roots:
        for cand in candidates:
            p = os.path.join(root, cand)
            if os.path.exists(p):
                return os.path.abspath(p)
    return None


def _cmd_proofs(args: argparse.Namespace) -> int:
    pkg_dir = _pkg_install_dir()
    wp_dir = os.path.join(pkg_dir, "whitepapers")
    if not os.path.isdir(wp_dir):
        wp_dir = _find_corpus_path("whitepapers")
    if not wp_dir or not os.path.isdir(wp_dir):
        print("ERROR: whitepapers/ not found in install. Reinstall uqff>=5.29.0.", file=sys.stderr)
        return 1

    action = args.action
    if action == "path":
        print(wp_dir)
        return 0
    if action == "list":
        files = sorted(f for f in os.listdir(wp_dir)
                       if f.endswith((".md", ".tex")))
        if args.filter:
            flt = args.filter.lower()
            files = [f for f in files if flt in f.lower()]
        if args.json:
            print(_dump(files))
        else:
            for f in files:
                print(f)
            print(f"\nTotal: {len(files)} whitepaper(s) in {wp_dir}")
        return 0
    if action == "show":
        if not args.name:
            print("ERROR: `uqff proofs show <NAME>` requires a whitepaper name.", file=sys.stderr)
            return 1
        # accept partial / case-insensitive match
        target = args.name.lower()
        candidates = [f for f in os.listdir(wp_dir)
                      if target in f.lower() and f.endswith((".md", ".tex"))]
        if not candidates:
            print(f"ERROR: no whitepaper matches '{args.name}'.", file=sys.stderr)
            return 1
        # prefer exact match
        exact = [f for f in candidates if f.lower() == target or f.lower() == target + ".md"]
        chosen = exact[0] if exact else candidates[0]
        if len(candidates) > 1 and not exact:
            print(f"# Multiple matches; showing {chosen}. Others:", file=sys.stderr)
            for c in candidates[1:5]:
                print(f"#   {c}", file=sys.stderr)
        path = os.path.join(wp_dir, chosen)
        with open(path, encoding="utf-8") as f:
            sys.stdout.write(f.read())
        return 0
    print(f"ERROR: unknown proofs action '{action}'.", file=sys.stderr)
    return 1


def _cmd_millennium(args: argparse.Namespace) -> int:
    """Run run_millennium_proofs.py — the 7-Clay-problem proof runner."""
    try:
        import run_millennium_proofs as rmp
        rmp.main()
        return 0
    except ImportError:
        # fallback: invoke directly if installed as a script
        path = _find_corpus_path("run_millennium_proofs.py")
        if not path:
            print("ERROR: run_millennium_proofs.py not found.", file=sys.stderr)
            return 1
        import subprocess
        return subprocess.call([sys.executable, path])


def _print_file(filename: str, label: str) -> int:
    """Find filename in the corpus and print its contents."""
    path = _find_corpus_path(filename)
    if not path:
        print(f"ERROR: {label} not found ({filename}).", file=sys.stderr)
        return 1
    with open(path, encoding="utf-8") as f:
        sys.stdout.write(f.read())
    return 0


def _cmd_axioms(args: argparse.Namespace) -> int:
    return _print_file("AXIOMS_AND_THEOREMS.md", "axiom-theorem inventory")


def _cmd_atlas(args: argparse.Namespace) -> int:
    return _print_file("CLOSURE_ATLAS.md", "closure atlas")


def _cmd_gold_standard(args: argparse.Namespace) -> int:
    return _print_file("Gold_Standard_Pure_UQFF.md", "Gold Standard reference")


def _cmd_manuscript(args: argparse.Namespace) -> int:
    """Print the path to the bundled manuscript PDF."""
    candidates = [
        os.path.join("Manuscript_1_12Feb2026", "uqff_production_arxiv.pdf"),
        os.path.join("Manuscript 1_12Feb2026", "uqff_production_arxiv.pdf"),
    ]
    for c in candidates:
        path = _find_corpus_path(c)
        if path:
            print(path)
            return 0
    print("ERROR: manuscript PDF not found.", file=sys.stderr)
    return 1


def _cmd_lean(args: argparse.Namespace) -> int:
    """Print the path to the Lean 4 formal-verification scaffold."""
    pkg_dir = _pkg_install_dir()
    formal = os.path.join(pkg_dir, "formal")
    if not os.path.isdir(formal):
        formal = _find_corpus_path("formal")
    if not formal or not os.path.isdir(formal):
        print("ERROR: formal/ Lean scaffold not found.", file=sys.stderr)
        return 1
    print(f"Lean 4 scaffold at: {formal}")
    print("\nFiles:")
    for root, dirs, files in os.walk(formal):
        for f in sorted(files):
            rel = os.path.relpath(os.path.join(root, f), formal)
            print(f"  {rel}")
    print("\nTo build (requires Lean 4 + Lake):")
    print(f"  cd {formal} && lake build")
    return 0


def _cmd_grok_archives(args: argparse.Namespace) -> int:
    archives = [
        "grok._b9afa8b6_3b85_31May2026.md",
        "grok_8461fe4e_c903.md",
        "grok_b8e305e6_1f29.md",
    ]
    found = []
    for a in archives:
        p = _find_corpus_path(a)
        if p:
            size_mb = os.path.getsize(p) / 1024 / 1024
            found.append((a, p, size_mb))
    if not found:
        print("ERROR: no grok archives found.", file=sys.stderr)
        return 1
    print("Bundled Grok proof-conversation archives:\n")
    for name, path, size_mb in found:
        print(f"  {name}  ({size_mb:.2f} MB)")
        print(f"    {path}")
    print(f"\nTotal: {len(found)} archive(s)")
    print("\nTo view: cat <path>  or  less <path>")
    return 0


# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="uqff",
        description="UQFF Star-Magic command-line interface")
    parser.add_argument("--version", action="version", version=f"uqff {_VERSION}")

    sub = parser.add_subparsers(dest="cmd", required=True, metavar="COMMAND")

    # ---- Closure queries ----
    p_predict = sub.add_parser("predict", help="fetch a UQFF closure value (all namespaces)")
    p_predict.add_argument("name")
    p_predict.add_argument("--json", action="store_true")
    p_predict.set_defaults(func=_cmd_predict)

    p_search = sub.add_parser("search", help="search ALL namespaces for closure names")
    p_search.add_argument("substr")
    p_search.add_argument("--json", action="store_true")
    p_search.set_defaults(func=_cmd_search)

    p_list = sub.add_parser("list", help="list closure names")
    p_list.add_argument("--filter")
    p_list.add_argument("--all", action="store_true")
    p_list.add_argument("--dispatch", action="store_true",
        help="list only the 114 Phase E/F/G dispatch observables")
    p_list.add_argument("--domain", choices=["SI","SM","LCDM","astro","GR","chem","CM","bio","geo","KK"],
        help="filter dispatch observables by domain (implies --dispatch)")
    p_list.add_argument("--json", action="store_true")
    p_list.set_defaults(func=_cmd_list)

    # ---- Status / testing ----
    p_status = sub.add_parser("status", help="production status summary")
    p_status.add_argument("--json", action="store_true")
    p_status.set_defaults(func=_cmd_status)

    p_surf = sub.add_parser("surfaces", help="list public calculate_* surfaces")
    p_surf.add_argument("--json", action="store_true")
    p_surf.set_defaults(func=_cmd_surfaces)

    p_ver = sub.add_parser("version", help="print version + key metrics")
    p_ver.set_defaults(func=_cmd_version)

    p_gate = sub.add_parser("gate", help="run the fidelity gate")
    p_gate.set_defaults(func=_cmd_gate)

    p_assim = sub.add_parser("assimilate",
        help="route an observable through the qcalcgeom_solver 4x3 dispatch matrix (Phase E/F/G)")
    p_assim.add_argument("name", help="observable name (see `uqff list --dispatch`)")
    p_assim.add_argument("--geometry", default="auto",
        choices=["auto", "qcalcgeom", "bsfg", "dpm", "d26"],
        help="geometry backend (default: auto = owner geometry)")
    p_assim.add_argument("--numeric", default="numerical",
        choices=["symbolic", "numerical", "discrete", "all"],
        help="numeric backend (default: numerical)")
    p_assim.add_argument("--decompose", action="store_true",
        help="return the full 8-field solver-bus result dict")
    p_assim.add_argument("--json", action="store_true")
    p_assim.set_defaults(func=_cmd_assimilate)

    p_serve = sub.add_parser("serve", help="launch the REST API (requires uqff[api])")
    p_serve.add_argument("--host", default="127.0.0.1")
    p_serve.add_argument("--port", type=int, default=8000)
    p_serve.add_argument("--reload", action="store_true")
    p_serve.set_defaults(func=_cmd_serve)

    # ---- v5.29.0: Proof-corpus access ----
    p_proofs = sub.add_parser("proofs", help="access bundled whitepaper proof corpus")
    p_proofs.add_argument("action", choices=["list", "show", "path"],
                          help="list = enumerate; show = print one; path = install dir")
    p_proofs.add_argument("name", nargs="?", help="whitepaper name for `show`")
    p_proofs.add_argument("--filter", help="substring filter for `list`")
    p_proofs.add_argument("--json", action="store_true")
    p_proofs.set_defaults(func=_cmd_proofs)

    p_mill = sub.add_parser("millennium", help="run the 7-Clay-Millennium proof runner")
    p_mill.set_defaults(func=_cmd_millennium)

    p_ax = sub.add_parser("axioms", help="print AXIOMS_AND_THEOREMS.md")
    p_ax.set_defaults(func=_cmd_axioms)

    p_atlas = sub.add_parser("atlas", help="print CLOSURE_ATLAS.md")
    p_atlas.set_defaults(func=_cmd_atlas)

    p_gs = sub.add_parser("gold-standard", help="print Gold_Standard_Pure_UQFF.md")
    p_gs.set_defaults(func=_cmd_gold_standard)

    p_man = sub.add_parser("manuscript", help="print path to compiled manuscript PDF")
    p_man.set_defaults(func=_cmd_manuscript)

    p_lean = sub.add_parser("lean", help="print Lean 4 formal-verification scaffold location")
    p_lean.set_defaults(func=_cmd_lean)

    p_grok = sub.add_parser("grok-archives", help="list bundled Grok proof-conversation archives")
    p_grok.set_defaults(func=_cmd_grok_archives)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
