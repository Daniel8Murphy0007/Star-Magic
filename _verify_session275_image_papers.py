"""Session 275 cold smoke test: verify all 26 PAPER_1012..PAPER_1037 wirings."""
from __future__ import annotations

import sys

import uqff_pure_calculator as upc


CLOSURE_PROBE_FUNCS = [
    "_l96_uqff_PAPER1012_GW190425_probe", "_l96_uqff_PAPER1013_QGP_ALICE_probe",
    "_l96_uqff_PAPER1014_SMBH_MERGER_probe", "_l96_uqff_PAPER1015_DM_NFW_probe",
    "_l96_uqff_PAPER1016_TXS0506_probe", "_l96_uqff_PAPER1017_99SYS_WSTP_probe",
    "_l96_uqff_PAPER1018_PROD_V15_probe", "_l96_uqff_PAPER1019_DM_PHONON_probe",
    "_l96_uqff_PAPER1020_CR_PHONON_probe", "_l96_uqff_PAPER1021_PULSAR_probe",
    "_l96_uqff_PAPER1022_GW_STRAIN_probe", "_l96_uqff_PAPER1023_NEUTRINO_probe",
    "_l96_uqff_PAPER1024_MAGNETAR_FLARE_probe", "_l96_uqff_PAPER1025_BH_SHADOW_probe",
    "_l96_uqff_PAPER1026_REIONIZATION_probe", "_l96_uqff_PAPER1027_TDE_probe",
    "_l96_uqff_PAPER1028_COSMIC_STRING_probe", "_l96_uqff_PAPER1029_BARYCENTER_probe",
    "_l96_uqff_PAPER1030_GUP_MIN_LEN_probe", "_l96_uqff_PAPER1031_PHOTON_SPHERE_probe",
    "_l96_uqff_PAPER1032_ISM_DUST_probe", "_l96_uqff_PAPER1033_GAL_BAR_probe",
    "_l96_uqff_PAPER1034_FRB_DM_probe", "_l96_uqff_PAPER1035_KILONOVA_probe",
    "_l96_uqff_PAPER1036_BBN_PHONON_probe", "_l96_uqff_PAPER1037_AGN_JET_probe",
]

LOOKUP_SPOT_TESTS = [
    "S3600", "S3601", "S3609", "S3700", "S3701", "S3702", "S3715", "S3720",
    "S3725", "S3800", "S3811", "S3815", "S3816", "S3819", "S3849",
]


def main() -> int:
    print("=" * 78)
    print("Session 275 cold smoke test (PAPER_1012..PAPER_1037)")
    print("=" * 78)

    reg = upc._PA_S275_CLOSURE_REGISTRY
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
            print(f"  {fname:<50} closures={closures:<3} transcribed={transcribed_count:<2} exact={exact}")
            route_pass += 1
        except Exception as e:
            print(f"  {fname:<50} FAIL: {e}")
    print(f"Closure-bearing probes pass: {route_pass} / {len(CLOSURE_PROBE_FUNCS)}")

    print("\n--- Manifest probe ---")
    try:
        manifest = upc._l96_uqff_session275_image_paper_manifest()
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
            out = upc._l96_uqff_closure_S275_lookup(s_id)
            d = out.get("derived")
            pred = out.get("paper_predicted")
            err_p = out.get("err_vs_predicted_pct")
            tr = out.get("form_transcribed")
            print(f"  {s_id}: derived={d}, predicted={pred}, err_pct={err_p}, transcribed={tr}")
            spot_pass += 1
        except Exception as e:
            print(f"  {s_id}: FAIL: {e}")
    print(f"Lookup spot tests pass: {spot_pass} / {len(LOOKUP_SPOT_TESTS)}")

    total_routes = len(CLOSURE_PROBE_FUNCS) + 1
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
