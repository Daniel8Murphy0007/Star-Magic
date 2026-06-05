#!/usr/bin/env python3
"""
UQFF Pure Calculator — Test / Harness companion (Plan sec wire/ship/hook).

This is the ONLY runtime entry-point for the pure calculator (per Plan sec 12:
the calculator itself contains no __main__, no I/O, no timestamps). Everything
that needs side effects — printing, JSON dumps, CI smoke gates, the runtime
introspection report — lives here.

Usage:
    python uqff_pure_calculator_Test.py                # full smoke pass
    python uqff_pure_calculator_Test.py --inventory    # dump _dispatch_keys()
    python uqff_pure_calculator_Test.py --report       # diff% honesty report
    python uqff_pure_calculator_Test.py --quick        # quick gate (CI)

Exit code 0 on pass; non-zero on failure. NOT REPLACEMENT.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, Dict, List, Tuple

import uqff_pure_calculator as u


# === Public surface assertions ===================================================

_EXPECTED_PUBLIC = [
    "calculate_resonant_adpm",
    "calculate_scm",
    "calculate_f_u_bi",
    "calculate_f_u_bi_i",
    "calculate_triadic_g",
    "calculate_vacuum_ledger",
    "calculate_analytic_closures",
]


def _check(name: str, ok: bool, detail: str = "") -> Tuple[str, bool, str]:
    return (name, bool(ok), detail)


def _assert_shape(out: Dict[str, Any]) -> bool:
    return isinstance(out, dict) and "value" in out and "provenance" in out \
        and isinstance(out["provenance"], str) and "NOT REPLACEMENT" in out["provenance"]


# === Test cases ==================================================================

def test_public_surface() -> List[Tuple[str, bool, str]]:
    """The 7 calculate_* functions exist, are callable, and return the
    {value, provenance} contract with NOT REPLACEMENT tagged."""
    out: List[Tuple[str, bool, str]] = []
    for name in _EXPECTED_PUBLIC:
        fn = getattr(u, name, None)
        if not callable(fn):
            out.append(_check(f"public::{name}", False, "missing or not callable"))
            continue
        res = fn({})
        out.append(_check(f"public::{name}::shape", _assert_shape(res),
                          f"keys={sorted(res.keys()) if isinstance(res, dict) else type(res).__name__}"))
    return out


def test_version_and_inventory() -> List[Tuple[str, bool, str]]:
    """__version__ + _dispatch_keys() inventory present (Fix 6)."""
    out: List[Tuple[str, bool, str]] = []
    out.append(_check("version::present",
                      isinstance(u.__version__, str) and len(u.__version__) > 0,
                      f"version={u.__version__!r}"))
    inv = u._dispatch_keys()
    out.append(_check("inventory::shape", isinstance(inv, dict),
                      f"type={type(inv).__name__}"))
    out.append(_check("inventory::publics_match",
                      sorted(inv.get("public_functions", [])) == sorted(_EXPECTED_PUBLIC),
                      f"got={sorted(inv.get('public_functions', []))}"))
    out.append(_check("inventory::lagrangian_17",
                      len(inv.get("lagrangian_sectors", [])) == 17,
                      f"count={len(inv.get('lagrangian_sectors', []))}"))
    out.append(_check("inventory::millennium_8",
                      len(inv.get("millennium_targets", [])) == 8,
                      f"count={len(inv.get('millennium_targets', []))}"))
    out.append(_check("inventory::field_leaves_14",
                      len(inv.get("universal_field_leaves", [])) == 14,
                      f"count={len(inv.get('universal_field_leaves', []))}"))
    return out


def test_precision_honesty() -> List[Tuple[str, bool, str]]:
    """Fixes 1 + 2 + 5: provenance diff% is honest — not a baked-in 0.000%.

    The 4-term vacuum ledger has a real ~0.117% residual vs Planck Lambda;
    Yang-Mills closure has a real ~2300% residual vs lattice anchor. Both must
    surface in the provenance, not be silently rounded to '0.000%'.
    """
    out: List[Tuple[str, bool, str]] = []

    vac = u.calculate_vacuum_ledger({})
    vac_prov = vac["provenance"]
    out.append(_check("honesty::vac_ledger_diff_present",
                      "diff=" in vac_prov, "missing diff= token"))
    # The vacuum ledger residual is computed from the actual analytic chain
    # and must not be silently '0.0000%' unless the floats are bit-equal.
    vd = vac["value"]
    if vd["total"] != vd["planck_target"]:
        out.append(_check("honesty::vac_ledger_no_false_zero",
                          "diff=0.0000%" not in vac_prov,
                          f"total={vd['total']} target={vd['planck_target']}"))

    # Yang-Mills must not falsely claim 0.000%.
    ym = u._millennium("yang_mills")
    if ym:
        _, ym_prov = ym
        out.append(_check("honesty::ym_diff_present",
                          "diff=" in ym_prov, "missing diff= token"))
        out.append(_check("honesty::ym_no_false_zero",
                          "diff=0.000%" not in ym_prov and "diff=0.0000%" not in ym_prov,
                          "ym claims false 0% match vs lattice anchor"))
        out.append(_check("honesty::ym_stationarity_routed",
                          "_stationarity_residual" in ym_prov,
                          "Fix 5: Millennium provenance must cite the single closure operator"))
    return out


def test_dse_parallel_paths() -> List[Tuple[str, bool, str]]:
    """Fix 3 + DSE wrapper: CP1 returns numeric scalar (not None) when the
    resolver returns a dict, and the public DSE result complies with the
    {value, provenance} contract."""
    out: List[Tuple[str, bool, str]] = []
    res = u.calculate_analytic_closures({"dse": True, "system": "sgr_a_star", "r": 1.0e10})
    out.append(_check("dse::shape", _assert_shape(res), "contract violation"))
    ch = res["value"]["channels"]
    out.append(_check("dse::cp1_numeric",
                      isinstance(ch.get("cp1"), (int, float)) and ch["cp1"] is not None,
                      f"cp1={ch.get('cp1')!r}"))
    out.append(_check("dse::cp2_numeric",
                      isinstance(ch.get("cp2"), (int, float)),
                      f"cp2={ch.get('cp2')!r}"))
    out.append(_check("dse::cp3_numeric",
                      isinstance(ch.get("cp3"), (int, float)),
                      f"cp3={ch.get('cp3')!r}"))
    out.append(_check("dse::cp4_numeric",
                      isinstance(ch.get("cp4"), (int, float)),
                      f"cp4={ch.get('cp4')!r}"))
    return out


def test_spinor_and_predictions() -> List[Tuple[str, bool, str]]:
    """Fix 4: spinor closure and P1-P14 predictions are reachable via
    dataset['spinor'] and dataset['prediction'] keys."""
    out: List[Tuple[str, bool, str]] = []
    sp = u.calculate_analytic_closures({"spinor": True})
    out.append(_check("spinor::shape", _assert_shape(sp), "contract violation"))
    out.append(_check("spinor::has_locks",
                      isinstance(sp["value"], dict)
                      and "anchor_lock_1" in sp["value"]
                      and "anchor_lock_2_natural" in sp["value"],
                      f"keys={list(sp['value'].keys()) if isinstance(sp['value'], dict) else 'N/A'}"))

    for pid in ("p1", "p6", "p7", "p11", "p12", "p14", "kk", "xi_test", "ledger"):
        r = u.calculate_analytic_closures({"prediction": pid})
        out.append(_check(f"prediction::{pid}::shape",
                          isinstance(r, dict) and "value" in r and "provenance" in r,
                          f"got={r if not isinstance(r, dict) else 'dict'}"))
    return out


def test_lagrangian_sectors() -> List[Tuple[str, bool, str]]:
    """Slice 2: all 17 lagrangian sectors route through the single
    stationarity primitive and report residual=0 on empty terms."""
    out: List[Tuple[str, bool, str]] = []
    inv = u._dispatch_keys()
    for paper in inv["lagrangian_sectors"]:
        r = u.calculate_analytic_closures({"lagrangian_sector": paper})
        info = r["value"]
        out.append(_check(f"lagrangian::{paper}::closure",
                          info.get("residual") == 0.0,
                          f"residual={info.get('residual')}"))
    return out


def test_slice3_leaves() -> List[Tuple[str, bool, str]]:
    """Slice 3: all 14 universal-field leaves are reachable and numeric."""
    out: List[Tuple[str, bool, str]] = []
    for key in ("ug1", "ug2", "ug3", "ug4", "u_i", "u_m", "u_b",
                "f_env_layer27", "h_res", "a_res", "f_res", "u_dp",
                "k_nuc", "s_shell"):
        r = u.calculate_analytic_closures({key: True, "r": 1.0e9})
        ok = isinstance(r, dict) and isinstance(r.get("value"), (int, float))
        out.append(_check(f"leaf::{key}::numeric", ok,
                          f"value={r.get('value') if isinstance(r, dict) else r!r}"))
    return out


# === Report builders =============================================================

def _run_all() -> List[Tuple[str, bool, str]]:
    results: List[Tuple[str, bool, str]] = []
    results.extend(test_public_surface())
    results.extend(test_version_and_inventory())
    results.extend(test_precision_honesty())
    results.extend(test_dse_parallel_paths())
    results.extend(test_spinor_and_predictions())
    results.extend(test_lagrangian_sectors())
    results.extend(test_slice3_leaves())
    return results


def _print_results(results: List[Tuple[str, bool, str]]) -> int:
    passed = 0
    failed = 0
    for name, ok, detail in results:
        tag = "PASS" if ok else "FAIL"
        line = f"[{tag}] {name}"
        if detail:
            line += f"  ({detail})"
        print(line)
        if ok:
            passed += 1
        else:
            failed += 1
    print()
    print(f"Total: {passed + failed} | Passed: {passed} | Failed: {failed}")
    print(f"uqff_pure_calculator version: {u.__version__}")
    return 0 if failed == 0 else 1


def _print_inventory() -> int:
    print(json.dumps(u._dispatch_keys(), indent=2, sort_keys=True))
    return 0


def _print_honesty_report() -> int:
    """Surface every diff% currently produced by the public functions so a
    reviewer can confirm there are no false 0.000% claims."""
    print(f"# UQFF Pure Calculator — diff% honesty report")
    print(f"# version: {u.__version__}\n")

    vac = u.calculate_vacuum_ledger({})
    print("## 4-term vacuum ledger (calculate_vacuum_ledger):")
    print(f"   total       = {vac['value']['total']:.6e}")
    print(f"   planck_tgt  = {vac['value']['planck_target']:.6e}")
    print(f"   provenance  = {vac['provenance'][-160:]}\n")

    print("## 8 Millennium closures (analytic-vs-anchor residual):")
    for key in sorted(u.MILLENNIUM_TARGETS.keys()):
        out = u._millennium(key)
        if out is None:
            continue
        val, prov = out
        # Pull the diff= ... token out of provenance.
        diff_tok = next((t for t in prov.split("|") if "diff=" in t), "").strip()
        ref = u.MILLENNIUM_TARGETS[key][0]
        print(f"   {key:<18s} UQFF={val!r:<22s} REF={ref!r:<12s} {diff_tok}")
    print()

    print("## Spinor closure (Map sec 9 row 9):")
    sp = u.calculate_analytic_closures({"spinor": True})
    if isinstance(sp.get("value"), dict):
        for k, v in sp["value"].items():
            print(f"   {k:<28s} = {v!r}")
    return 0


# === Entry-point =================================================================

def main(argv: List[str]) -> int:
    p = argparse.ArgumentParser(description="UQFF pure-calculator test harness")
    p.add_argument("--inventory", action="store_true",
                   help="dump _dispatch_keys() JSON and exit")
    p.add_argument("--report", action="store_true",
                   help="print diff%% honesty report and exit")
    p.add_argument("--quick", action="store_true",
                   help="quick CI gate: public surface + version only")
    args = p.parse_args(argv)

    if args.inventory:
        return _print_inventory()
    if args.report:
        return _print_honesty_report()
    if args.quick:
        results = test_public_surface() + test_version_and_inventory()
        return _print_results(results)

    return _print_results(_run_all())


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
