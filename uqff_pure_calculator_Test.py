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


def test_f_ubii_proofs_complete() -> List[Tuple[str, bool, str]]:
    """Gap 1 (04Jun2026 audit): 17 named F_UBii proofs all callable, finite, non-zero."""
    import math as _m
    out: List[Tuple[str, bool, str]] = []
    expected = {"virx", "termv", "upar", "coup", "orbdec", "kn", "fermi", "kne",
                "whim", "ps", "sfe", "hawk", "bd", "roche", "ent", "dec", "lobe"}
    inv = u._f_ubii_inventory()
    proof_names = set(inv.get("proof_names", [])) if isinstance(inv, dict) else set()
    out.append(_check("f_ubii::inventory_size_17",
                      inv.get("total_proofs") == 17 and expected.issubset(proof_names),
                      f"total={inv.get('total_proofs')!r}, missing={sorted(expected - proof_names)}"))
    for name in sorted(expected):
        try:
            res = u._f_ubii_proof(name)
            val = res["value"] if isinstance(res, dict) else res
            v = val[0] if isinstance(val, tuple) else val
            ok = isinstance(v, (int, float)) and _m.isfinite(v) and v != 0.0
            out.append(_check(f"f_ubii::{name}::finite_nonzero", ok, f"v={v!r}"))
        except Exception as e:
            out.append(_check(f"f_ubii::{name}::callable", False, f"exc={e!r}"))
    return out


def test_quantum_chain_ladder() -> List[Tuple[str, bool, str]]:
    """Gap 2 (04Jun2026 audit): 26-rung Quantum Chain ladder + derive_from_quantum_chain."""
    import math as _m
    out: List[Tuple[str, bool, str]] = []
    out.append(_check("quantum_chain::n_levels_26",
                      u.N_QUANTUM_CHAIN_LEVELS == 26,
                      f"got {u.N_QUANTUM_CHAIN_LEVELS}"))
    rungs_ok = True
    for n in range(1, u.N_QUANTUM_CHAIN_LEVELS + 1):
        e = u._quantum_chain_E_n(n)
        if not (isinstance(e, (int, float)) and _m.isfinite(e) and e > 0.0):
            rungs_ok = False
            break
    out.append(_check("quantum_chain::26_rungs_positive_finite", rungs_ok,
                      f"E_1={u._quantum_chain_E_n(1):.3e} E_26={u._quantum_chain_E_n(26):.3e}"))
    inv = u._quantum_chain_inventory()
    out.append(_check("quantum_chain::inventory_shape",
                      isinstance(inv, dict) and "n_levels" in inv and "E_total_J" in inv,
                      f"keys={list(inv.keys()) if isinstance(inv,dict) else type(inv).__name__}"))
    deriv = u._derive_from_quantum_chain()
    needed = {"kg_chain", "rho_mass_kg_m3", "rho_E_J_m3"}
    out.append(_check("quantum_chain::derive_three_quantities",
                      isinstance(deriv, dict) and needed.issubset(deriv.keys()),
                      f"missing={sorted(needed - set(deriv.keys())) if isinstance(deriv,dict) else 'not-dict'}"))
    return out


def test_paradox_proofs_dispatched() -> List[Tuple[str, bool, str]]:
    """Gap 5 (04Jun2026 audit): 8 paradox proofs route through Millennium ports;
    info_paradox provenance must surface the F_UBii_ent operator value."""
    import math as _m
    out: List[Tuple[str, bool, str]] = []
    out.append(_check("paradox::registry_8_entries",
                      len(u.PARADOX_TO_MILLENNIUM) == 8,
                      f"got {len(u.PARADOX_TO_MILLENNIUM)}"))
    for paradox_name in sorted(u.PARADOX_TO_MILLENNIUM.keys()):
        try:
            res = u._paradox_proof(paradox_name)
            val = res["value"] if isinstance(res, dict) else res
            v = val[0] if isinstance(val, tuple) else val
            ok = isinstance(v, (int, float)) and _m.isfinite(v)
            out.append(_check(f"paradox::{paradox_name}::finite", ok, f"v={v!r}"))
        except Exception as e:
            out.append(_check(f"paradox::{paradox_name}::callable", False, f"exc={e!r}"))
    info = u._paradox_proof("info_paradox")
    prov = info[1] if isinstance(info, tuple) else (info.get("provenance", "") if isinstance(info, dict) else "")
    out.append(_check("paradox::info_paradox::cites_F_UBii_ent",
                      "F_UBii_ent" in prov,
                      f"prov[:120]={prov[:120]!r}"))
    return out


def test_new_inflation_primitives() -> List[Tuple[str, bool, str]]:
    """Gaps 3+4 (04Jun2026 audit): 12 new primitives (10 inflation/cosmology + Lambda_QCD + f_b)."""
    import math as _m
    out: List[Tuple[str, bool, str]] = []
    new_keys = ("omega_gw_h2", "r_tensor_scalar", "dn_s_dlnk",
                "f_nl_equil", "f_nl_orth", "epsilon_slow_roll",
                "eta_slow_roll", "n_efolds", "t_reh_gev",
                "h_inflation_gev", "lambda_qcd_mev", "f_baryon")
    out.append(_check("primitives::ledger_has_12_new_keys",
                      all(k in u._LEDGER_PRIMITIVE for k in new_keys),
                      f"missing={[k for k in new_keys if k not in u._LEDGER_PRIMITIVE]}"))
    out.append(_check("primitives::ledger_count_at_least_160",
                      len(u._LEDGER_PRIMITIVE) >= 160,
                      f"got {len(u._LEDGER_PRIMITIVE)}"))
    for k in new_keys:
        try:
            v = u._LEDGER_PRIMITIVE[k]()
            ok = isinstance(v, (int, float)) and _m.isfinite(v) and v != 0.0
            out.append(_check(f"primitive::{k}::finite_nonzero", ok, f"v={v!r}"))
        except Exception as e:
            out.append(_check(f"primitive::{k}::callable", False, f"exc={e!r}"))
    return out


def test_dispatch_inventory_counts() -> List[Tuple[str, bool, str]]:
    """Gap 1-5: dispatch_keys exposes complete inventory surfaces with correct counts."""
    out: List[Tuple[str, bool, str]] = []
    d = u._dispatch_keys()
    out.append(_check("dispatch::f_ubii_proofs_17",
                      "f_ubii_proofs" in d and len(d["f_ubii_proofs"]) == 17,
                      f"got {len(d.get('f_ubii_proofs', []))}"))
    out.append(_check("dispatch::paradox_proofs_8",
                      "paradox_proofs" in d and len(d["paradox_proofs"]) == 8,
                      f"got {len(d.get('paradox_proofs', []))}"))
    out.append(_check("dispatch::quantum_chain_n_levels_26",
                      isinstance(d.get("quantum_chain"), dict)
                      and d["quantum_chain"].get("n_levels") == 26,
                      f"got {d.get('quantum_chain')!r}"))
    out.append(_check("dispatch::ledger_primitive_keys_count",
                      "ledger_primitive_keys" in d and len(d["ledger_primitive_keys"]) >= 160,
                      f"got {len(d.get('ledger_primitive_keys', []))}"))
    out.append(_check("dispatch::regime_corpus_1018",
                      isinstance(d.get("regime_corpus"), dict)
                      and d["regime_corpus"].get("total_variants") == 1018,
                      f"got {d.get('regime_corpus')!r}"))
    return out


def test_constant_closure_report() -> List[Tuple[str, bool, str]]:
    """Session 261 polish: closure report mechanically grades each fundamental.

    Verifies the report exposes 16 named constants + _summary, that h/c/k_B/e/N_A
    are correctly tagged as 'identity', that sigma_SB / Lambda / hbar are
    'derived' (closure error < 1%), that Delta_SCm is 'hardcoded', and that
    the documented broken-formula constants (G, eps_0, mu_0, Rydberg, Bohr_a0,
    Compton, Wien_b, mu_B) are tagged 'broken' pending Layer 45 repair.
    """
    out: List[Tuple[str, bool, str]] = []
    r = u._constant_closure_report()
    expected_constants = {
        "h", "hbar", "c", "k_B", "e", "N_A", "G", "sigma_SB",
        "eps_0", "mu_0", "Rydberg", "Bohr_a0", "Compton", "Wien_b",
        "mu_B", "Lambda", "Delta_SCm",
    }
    out.append(_check("closure_report::all_17_constants_present",
                      expected_constants.issubset(set(r.keys())),
                      f"missing: {expected_constants - set(r.keys())}"))
    out.append(_check("closure_report::summary_present",
                      "_summary" in r and isinstance(r["_summary"], dict),
                      f"summary keys: {list(r.get('_summary', {}).keys())}"))
    # Identity assertions (SI base anchors are not UQFF-derived; they are definitional)
    for name in ("h", "c", "k_B", "e", "N_A"):
        out.append(_check(f"closure_report::{name}_is_identity",
                          r[name]["status"] == "identity",
                          f"status={r[name]['status']}"))
    # Genuine derivations (closure error < 1%)
    for name in ("hbar", "sigma_SB", "Lambda"):
        out.append(_check(f"closure_report::{name}_is_derived",
                          r[name]["status"] == "derived" and r[name]["err_pct"] < 1.0,
                          f"status={r[name]['status']} err={r[name]['err_pct']:.4f}%"))
    # Delta_SCm: hardcoded primitive, geometric derivation pending
    out.append(_check("closure_report::delta_scm_is_hardcoded",
                      r["Delta_SCm"]["status"] == "hardcoded",
                      f"status={r['Delta_SCm']['status']}"))
    # Honest 'broken' tags (UQFF saturation formulas pending Layer 45 repair)
    for name in ("G", "eps_0", "mu_0", "Rydberg", "Bohr_a0", "Compton", "Wien_b", "mu_B"):
        out.append(_check(f"closure_report::{name}_flagged_broken",
                          r[name]["status"] == "broken",
                          f"status={r[name]['status']} err={r[name]['err_pct']:.4g}%"))
    # Summary counters must match expected polish-session totals
    s = r["_summary"]
    out.append(_check("closure_report::summary_total_17",
                      s.get("total") == 17, f"total={s.get('total')}"))
    out.append(_check("closure_report::summary_derived_3",
                      s.get("derived") == 3, f"derived={s.get('derived')}"))
    out.append(_check("closure_report::summary_identity_5",
                      s.get("identity") == 5, f"identity={s.get('identity')}"))
    out.append(_check("closure_report::summary_broken_8",
                      s.get("broken") == 8, f"broken={s.get('broken')}"))
    out.append(_check("closure_report::summary_hardcoded_1",
                      s.get("hardcoded") == 1, f"hardcoded={s.get('hardcoded')}"))
    # Dispatch routes
    out.append(_check("closure_report::derive_constant_route",
                      isinstance(u._derive_constant("closure_report"), dict),
                      "_derive_constant('closure_report') returns dict"))
    out.append(_check("closure_report::dispatch_keys_summary",
                      isinstance(u._dispatch_keys().get("constant_closure_report"), dict)
                      and u._dispatch_keys()["constant_closure_report"].get("total") == 17,
                      "dispatch_keys exposes summary"))
    return out


def test_io_ports_session262() -> List[Tuple[str, bool, str]]:
    """Session 262: IPData/OPData symbolic IO wiring.

    The pure calculator must orchestrate all 7 calculate_* surfaces from a
    single symbolic input (IPData.InputParameters or dict) and persist the
    result via OPData. Validates the offline-safe lazy-import pattern, the
    98-field IPData schema exposure, and the spontaneous-answer contract.
    """
    out: List[Tuple[str, bool, str]] = []

    # --- io_surface metadata is published via _dispatch_keys ---
    d = u._dispatch_keys()
    out.append(_check("io::dispatch_io_surface_present",
                      "io_surface" in d and isinstance(d["io_surface"], dict),
                      f"keys: {list(d.keys())[-5:]}"))
    io = d.get("io_surface", {})
    out.append(_check("io::ipdata_available",
                      io.get("ipdata_available") is True,
                      f"got {io.get('ipdata_available')}"))
    out.append(_check("io::opdata_available",
                      io.get("opdata_available") is True,
                      f"got {io.get('opdata_available')}"))
    out.append(_check("io::ipdata_schema_keys_at_least_90",
                      len(io.get("ipdata_schema_keys", [])) >= 90,
                      f"got {len(io.get('ipdata_schema_keys', []))}"))
    out.append(_check("io::calculate_surfaces_7",
                      len(io.get("calculate_surfaces", [])) == 7,
                      f"got {len(io.get('calculate_surfaces', []))}"))

    # --- _solve_symbolic with bare kwargs ---
    res = u._solve_symbolic(M=1.989e30, r=1.496e11, t_n=0.0, omega=1.25e12)
    out.append(_check("io::solve_symbolic_returns_dict",
                      isinstance(res, dict) and "available_equations" in res,
                      "missing available_equations"))
    out.append(_check("io::solve_symbolic_no_errors",
                      res.get("errors") == {},
                      f"errors: {res.get('errors')}"))
    out.append(_check("io::solve_symbolic_all_7_fired",
                      len(res.get("available_equations", [])) == 7,
                      f"available: {res.get('available_equations')}"))
    out.append(_check("io::solve_symbolic_param_count_4",
                      res.get("input_param_count") == 4,
                      f"got {res.get('input_param_count')}"))

    # --- _solve_from_input with IPData.InputParameters ---
    import IPData as ipd
    p = ipd.InputParameters(
        query_name="harness_smoke",
        M=1.989e30, r=1.496e11, T=5778.0, B=1e-3, omega=1.25e12,
    )
    res = u._solve_from_input(p, query_name="harness_smoke_test")
    out.append(_check("io::ipdata_solve_query_id_present",
                      isinstance(res.get("query_id"), str) and len(res["query_id"]) > 0,
                      f"qid: {res.get('query_id')}"))
    out.append(_check("io::ipdata_solve_no_errors",
                      res.get("errors") == {},
                      f"errors: {res.get('errors')}"))
    out.append(_check("io::ipdata_solve_long_form_7",
                      len(res.get("long_form_equations", [])) == 7,
                      f"got {len(res.get('long_form_equations', []))}"))
    # Each long_form entry must have a provenance string
    lf = res.get("long_form_equations", [])
    out.append(_check("io::ipdata_solve_long_form_has_provenance",
                      all(isinstance(e.get("provenance"), str) and len(e["provenance"]) > 0
                          for e in lf),
                      f"provenance missing on {sum(1 for e in lf if not e.get('provenance'))} entries"))

    # --- OPData recall round-trip ---
    back = u._recall(res["query_id"])
    out.append(_check("io::opdata_recall_round_trip",
                      back is not None,
                      f"recall returned None for {res['query_id']}"))
    if back:
        out.append(_check("io::opdata_recall_query_name",
                          back.get("input_params", {}).get("query_name") == "harness_smoke",
                          f"got {back.get('input_params', {}).get('query_name')}"))

    # --- _input_to_dataset coercion contract ---
    out.append(_check("io::coerce_none_to_empty",
                      u._input_to_dataset(None) == {},
                      "None coercion failed"))
    out.append(_check("io::coerce_dict_passthrough",
                      u._input_to_dataset({"M": 1.0, "r": 2.0}) == {"M": 1.0, "r": 2.0},
                      "dict coercion failed"))
    coerced = u._input_to_dataset(p)
    out.append(_check("io::coerce_ipdata_to_dict",
                      isinstance(coerced, dict) and coerced.get("M") == 1.989e30,
                      f"coerced.M = {coerced.get('M')}"))

    # --- offline-safe contract: surfaces subset works ---
    res = u._solve_from_input({"M": 1e30, "r": 1e10},
                              surfaces=["calculate_f_u_bi"], store=False)
    out.append(_check("io::surfaces_subset_one_surface",
                      res.get("available_equations") == ["calculate_f_u_bi"],
                      f"got {res.get('available_equations')}"))

    return out


def test_new_polish_primitives() -> List[Tuple[str, bool, str]]:
    """Session 261 polish: hbar, k_b, delta_scm_j exposed as ledger primitives + dispatch."""
    out: List[Tuple[str, bool, str]] = []
    import math as _m
    hbar_val = u._hbar_primitive_sat()
    out.append(_check("polish::hbar_value",
                      abs(hbar_val - u.PLANCK_H / (2.0 * _m.pi)) < 1e-45,
                      f"got {hbar_val}"))
    out.append(_check("polish::k_b_primitive",
                      u._k_b_primitive_sat() == u.K_B, "identity"))
    out.append(_check("polish::delta_scm_primitive",
                      u._delta_scm_primitive_sat() == u.DELTA_SCM_J, "matches DELTA_SCM_J"))
    # Ledger primitive registration
    out.append(_check("polish::hbar_in_ledger",
                      "hbar" in u._LEDGER_PRIMITIVE, "registered"))
    out.append(_check("polish::k_b_in_ledger",
                      "k_b" in u._LEDGER_PRIMITIVE, "registered"))
    out.append(_check("polish::delta_scm_j_in_ledger",
                      "delta_scm_j" in u._LEDGER_PRIMITIVE, "registered"))
    out.append(_check("polish::ledger_count_at_least_163",
                      len(u._LEDGER_PRIMITIVE) >= 163, f"got {len(u._LEDGER_PRIMITIVE)}"))
    # Dispatch routes
    out.append(_check("polish::derive_hbar",
                      isinstance(u._derive_constant("hbar"), float),
                      "_derive_constant('hbar') returns float"))
    out.append(_check("polish::derive_delta_scm",
                      isinstance(u._derive_constant("delta_scm"), float),
                      "_derive_constant('delta_scm') returns float"))
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
    results.extend(test_f_ubii_proofs_complete())
    results.extend(test_quantum_chain_ladder())
    results.extend(test_paradox_proofs_dispatched())
    results.extend(test_new_inflation_primitives())
    results.extend(test_dispatch_inventory_counts())
    results.extend(test_constant_closure_report())
    results.extend(test_new_polish_primitives())
    results.extend(test_io_ports_session262())
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
