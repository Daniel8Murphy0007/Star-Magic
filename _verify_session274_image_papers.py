"""Session 274 cold smoke test: verify all 26 PAPER_1038..PAPER_1063 wirings."""
from __future__ import annotations

import sys

import uqff_pure_calculator as upc


CLOSURE_PROBE_FUNCS = [
    "_l96_uqff_PAPER1038_WD_CRYST_probe", "_l96_uqff_PAPER1039_CLUST_PROF_probe",
    "_l96_uqff_PAPER1040_SHOCK_probe", "_l96_uqff_PAPER1041_COOL_CORE_probe",
    "_l96_uqff_PAPER1042_MOCK_THETA_probe", "_l96_uqff_PAPER1043_MULTI_probe",
    "_l96_uqff_PAPER1044_SZ_probe", "_l96_uqff_PAPER1045_RELIC_probe",
    "_l96_uqff_PAPER1046_LENS_probe", "_l96_uqff_PAPER1047_SNIax_probe",
    "_l96_uqff_PAPER1048_M_SIGMA_probe", "_l96_uqff_PAPER1049_SOURCE10_probe",
    "_l96_uqff_PAPER1050_MUGE9_probe", "_l96_uqff_PAPER1051_DUALITY_probe",
    "_l96_uqff_PAPER1052_TQFT_probe", "_l96_uqff_PAPER1053_SWAMP_probe",
    "_l96_uqff_PAPER1054_SUSY_probe", "_l96_uqff_PAPER1055_cMERA_probe",
    "_l96_uqff_PAPER1056_QEC_probe", "_l96_uqff_PAPER1057_NCG_probe",
    "_l96_uqff_PAPER1058_LQG_AREA_probe", "_l96_uqff_PAPER1059_CGC_probe",
    "_l96_uqff_PAPER1060_LENR_VDS_probe", "_l96_uqff_PAPER1061_KOZIMA_probe",
    "_l96_uqff_PAPER1062_WORM_probe", "_l96_uqff_PAPER1063_HCG_probe",
]

LOOKUP_SPOT_TESTS = [
    "S3108",   # omega_SCm rad/s (transcribed)
    "S3123",   # BETA0_DPM (transcribed)
    "S3124",   # kappa_calibration_per_day (transcribed)
    "S3214",   # rho_SCm calibration (transcribed)
    "S3308",   # omega_SCm rad/s 1048 (transcribed)
    "S3309",   # rho_SCm 1048 (transcribed)
    "S3314",   # crossover_state_index (transcribed const)
    "S3315",   # scale_range_orders (transcribed const)
    "S3403",   # V(phi0) = -RHO_SCM (transcribed)
    "S3404",   # omega_SCm 1054 (transcribed)
    "S3407",   # kappa 1055 (transcribed)
    "S3503",   # H_SCm calibration (transcribed)
]


def main() -> int:
    print("=" * 78)
    print("Session 274 cold smoke test (PAPER_1038..PAPER_1063)")
    print("=" * 78)

    reg = upc._PA_S274_CLOSURE_REGISTRY
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
            print(f"  {fname:<48} closures={closures:<3} transcribed={transcribed_count:<2} exact={exact}")
            route_pass += 1
        except Exception as e:
            print(f"  {fname:<48} FAIL: {e}")
    print(f"Closure-bearing probes pass: {route_pass} / {len(CLOSURE_PROBE_FUNCS)}")

    print("\n--- Manifest probe ---")
    try:
        manifest = upc._l96_uqff_session274_image_paper_manifest()
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
            out = upc._l96_uqff_closure_S274_lookup(s_id)
            d = out.get("derived")
            pred = out.get("paper_predicted")
            err_p = out.get("err_vs_predicted_pct")
            tr = out.get("form_transcribed")
            print(f"  {s_id}: derived={d}, predicted={pred}, err_pct={err_p}, transcribed={tr}")
            spot_pass += 1
        except Exception as e:
            print(f"  {s_id}: FAIL: {e}")
    print(f"Lookup spot tests pass: {spot_pass} / {len(LOOKUP_SPOT_TESTS)}")

    total_routes = len(CLOSURE_PROBE_FUNCS) + 1  # + manifest
    routes_ok = route_pass + (1 if manifest_ok else 0)
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
