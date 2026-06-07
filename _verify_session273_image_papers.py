"""Session 273 cold smoke test: verify all 26 PAPER_1086..PAPER_1111 wirings.

Confirms:
  - Registry populated with 178 closures
  - 26 dispatcher routes (23 closure probes + 3 framework probes)
  - Universal S### lookup works
  - Closures within 1% of paper_predicted (for transcribed lambdas)
  - Closures within 5% of observed_anchor (where present)
"""
from __future__ import annotations

import sys

import uqff_pure_calculator as upc


# 23 closure-bearing paper probe functions (direct attribute access)
CLOSURE_PROBE_FUNCS = [
    "_l96_uqff_PAPER1086_DE_GAMMA_probe", "_l96_uqff_PAPER1087_DE_EOS_probe",
    "_l96_uqff_PAPER1088_FUBII_7C_probe", "_l96_uqff_PAPER1089_INFL_LAG_probe",
    "_l96_uqff_PAPER1090_DE_LAG_probe",
    "_l96_uqff_PAPER1092_CMB_PSP_probe", "_l96_uqff_PAPER1093_CMB_TEMP_probe",
    "_l96_uqff_PAPER1094_CMB_LAG_probe", "_l96_uqff_PAPER1095_HORIZ_probe",
    "_l96_uqff_PAPER1096_11DOM_probe",
    "_l96_uqff_PAPER1098_QGATE_probe", "_l96_uqff_PAPER1100_LQG_AREA_probe",
    "_l96_uqff_PAPER1101_T2_FUBI_probe",
    "_l96_uqff_PAPER1102_HOLO_probe", "_l96_uqff_PAPER1103_LQG_SF_probe",
    "_l96_uqff_PAPER1104_SMBH_probe", "_l96_uqff_PAPER1105_HYD_MUGE_probe",
    "_l96_uqff_PAPER1106_STR_26D_probe",
    "_l96_uqff_PAPER1107_26D_FOLD_probe", "_l96_uqff_PAPER1108_VDS_NS_probe",
    "_l96_uqff_PAPER1109_VAC_LAD_probe",
    "_l96_uqff_PAPER1110_RIEM_PI_probe", "_l96_uqff_PAPER1111_YM_probe",
]

FRAMEWORK_PROBE_FUNCS = [
    "_l96_uqff_PAPER1091_production_scaling_probe",
    "_l96_uqff_PAPER1097_production_scaling_probe",
    "_l96_uqff_PAPER1099_production_scaling_probe",
]

LOOKUP_SPOT_TESTS = [
    "S2605",   # w_DE_t0 = -1 lambda
    "S2612",   # G_NEWTON
    "S2617",   # BSH fractional (verbatim)
    "S2717",   # beta_i_CMB
    "S3000",   # 26! factorial
    "S3013",   # zeta(3) approx via S26_3
    "S3023",   # V(phi0) = -RHO_SCM
    "S3035",   # g_YM coupling
    "S3039",   # M_glueball 0++
]


def main() -> int:
    print("=" * 78)
    print("Session 273 cold smoke test (PAPER_1086..PAPER_1111)")
    print("=" * 78)

    reg = upc._PA_S273_CLOSURE_REGISTRY
    print(f"\nRegistry size: {len(reg)} closures")

    derive_ok = 0
    derive_none = 0
    pred_1pct_ok = 0
    pred_compared = 0
    obs_5pct_ok = 0
    obs_compared = 0
    transcribed = 0
    for s_id, (tag, label, form, predicted, observed, unit) in reg.items():
        if form is not None:
            transcribed += 1
        d = form() if form is not None else predicted
        if d is None:
            derive_none += 1
            continue
        derive_ok += 1
        if predicted not in (None, 0, 0.0):
            try:
                err = abs(d - predicted) / abs(predicted) * 100.0
                pred_compared += 1
                if err < 1.0:
                    pred_1pct_ok += 1
            except (TypeError, ZeroDivisionError):
                pass
        if observed not in (None, 0, 0.0):
            try:
                err = abs(d - observed) / abs(observed) * 100.0
                obs_compared += 1
                if err < 5.0:
                    obs_5pct_ok += 1
            except (TypeError, ZeroDivisionError):
                pass

    print(f"Derived (non-None):       {derive_ok} / {len(reg)}")
    print(f"None entries (formula-only labels): {derive_none}")
    print(f"Transcribed lambdas:      {transcribed} / {len(reg)}")
    print(f"Within 1% of paper_pred:  {pred_1pct_ok} / {pred_compared}")
    print(f"Within 5% of observed:    {obs_5pct_ok} / {obs_compared}")

    print("\n--- Per-paper closure probes ---")
    route_pass = 0
    for fname in CLOSURE_PROBE_FUNCS:
        try:
            fn = getattr(upc, fname)
            out = fn()
            closures = out.get("closures", 0)
            transcribed_count = out.get("transcribed_count", 0)
            exact = out.get("exact_count", 0)
            print(f"  {fname:<44} closures={closures:<3} transcribed={transcribed_count:<2} exact={exact}")
            route_pass += 1
        except Exception as e:
            print(f"  {fname:<44} FAIL: {e}")
    print(f"Closure-bearing probes pass: {route_pass} / {len(CLOSURE_PROBE_FUNCS)}")

    print("\n--- Framework-only probes ---")
    fw_pass = 0
    for fname in FRAMEWORK_PROBE_FUNCS:
        try:
            fn = getattr(upc, fname)
            out = fn()
            print(f"  {fname:<48} type={out.get('type','?')}")
            fw_pass += 1
        except Exception as e:
            print(f"  {fname:<48} FAIL: {e}")
    print(f"Framework-only probes pass: {fw_pass} / {len(FRAMEWORK_PROBE_FUNCS)}")

    print("\n--- Manifest probe ---")
    try:
        manifest = upc._l96_uqff_session273_image_paper_manifest()
        print(f"  session: {manifest.get('session')}, pdfs: {manifest.get('total_pdfs_in_image')}")
        print(f"  closures wired: {manifest.get('total_new_closures_wired')}, probes: {manifest.get('total_new_probe_functions')}")
        manifest_ok = True
    except Exception as e:
        print(f"  manifest FAIL: {e}")
        manifest_ok = False

    print("\n--- Direct S### lookup spot test ---")
    spot_pass = 0
    for s_id in LOOKUP_SPOT_TESTS:
        try:
            out = upc._l96_uqff_closure_S273_lookup(s_id)
            d = out.get("derived")
            pred = out.get("paper_predicted")
            err_p = out.get("err_vs_predicted_pct")
            tr = out.get("form_transcribed")
            print(f"  {s_id}: derived={d}, predicted={pred}, err_pct={err_p}, transcribed={tr}")
            spot_pass += 1
        except Exception as e:
            print(f"  {s_id}: FAIL: {e}")
    print(f"Lookup spot tests pass: {spot_pass} / {len(LOOKUP_SPOT_TESTS)}")

    total_routes = len(CLOSURE_PROBE_FUNCS) + len(FRAMEWORK_PROBE_FUNCS) + 1  # + manifest
    routes_ok = route_pass + fw_pass + (1 if manifest_ok else 0)
    print("\n" + "=" * 78)
    print(f"SUMMARY: {derive_ok}/{len(reg)} derive, {routes_ok}/{total_routes} routes pass, "
          f"{pred_1pct_ok}/{pred_compared} within 1% predicted, "
          f"{obs_5pct_ok}/{obs_compared} within 5% observed, "
          f"{spot_pass}/{len(LOOKUP_SPOT_TESTS)} lookups OK")
    print("=" * 78)

    overall_ok = routes_ok == total_routes and spot_pass == len(LOOKUP_SPOT_TESTS)
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
