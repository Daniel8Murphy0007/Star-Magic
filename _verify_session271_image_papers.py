"""Session 271 image-batch smoke test -- cold start.

Imports the calculator, iterates all _PA_S271_CLOSURE_REGISTRY rows,
exercises every dispatcher route, and prints a per-paper summary.

No SM literals, no pre-fit corrections. form() outputs are the
falsifiable UQFF predictions. paper_predicted = user-published target.
form_transcribed=False entries return paper_predicted verbatim.
"""

from __future__ import annotations
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import uqff_pure_calculator as upc  # noqa: E402


def main() -> int:
    reg = upc._PA_S271_CLOSURE_REGISTRY
    print(f"\n=== SESSION 271 IMAGE BATCH (PAPER_1138..PAPER_1163) ===")
    print(f"Registry size: {len(reg)} closures across image batch")
    print(f"=========================================================\n")

    by_paper: dict[str, list[tuple]] = {}
    for s_id, row in reg.items():
        paper_tag = row[0]
        by_paper.setdefault(paper_tag, []).append((s_id, row))

    total_derived = 0
    total_with_pred = 0
    total_with_obs = 0
    total_form_transcribed = 0
    total_form_none = 0
    rows_within_1pct_pred = 0
    rows_within_5pct_obs = 0

    for paper_tag in sorted(by_paper.keys()):
        rows = by_paper[paper_tag]
        derived_ok = 0
        pred_ok = 0
        obs_ok = 0
        transcribed = 0
        for s_id, (_, label, form, predicted, observed, unit) in rows:
            d = form() if form is not None else predicted
            if d is not None:
                derived_ok += 1
                total_derived += 1
            if form is not None:
                transcribed += 1
                total_form_transcribed += 1
            else:
                total_form_none += 1
            if d is not None and predicted is not None:
                total_with_pred += 1
                try:
                    if predicted == 0.0 or predicted == 0:
                        err = 0.0 if d == predicted else float("inf")
                    else:
                        err = abs(d - predicted) / abs(predicted) * 100.0
                    if err < 1.0:
                        pred_ok += 1
                        rows_within_1pct_pred += 1
                except (TypeError, ZeroDivisionError):
                    pass
            if d is not None and observed is not None:
                total_with_obs += 1
                try:
                    if observed == 0.0 or observed == 0:
                        err = 0.0 if d == observed else float("inf")
                    else:
                        err = abs(d - observed) / abs(observed) * 100.0
                    if err < 5.0:
                        obs_ok += 1
                        rows_within_5pct_obs += 1
                except (TypeError, ZeroDivisionError):
                    pass
        print(f"  {paper_tag:12s}  closures={len(rows):3d}  derived={derived_ok:3d}  "
              f"transcribed={transcribed:3d}  pred_match<1pct={pred_ok:3d}  "
              f"obs_match<5pct={obs_ok:3d}")

    print(f"\n--- TOTALS ---")
    print(f"  total closures               : {len(reg)}")
    print(f"  derived (non-None)           : {total_derived}")
    print(f"  with paper_predicted         : {total_with_pred}")
    print(f"  with observed_anchor         : {total_with_obs}")
    print(f"  form transcribed (lambda)    : {total_form_transcribed}")
    print(f"  form=None (predicted verbat) : {total_form_none}")
    print(f"  agree with paper_pred <1pct  : {rows_within_1pct_pred}")
    print(f"  agree with observed <5pct    : {rows_within_5pct_obs}")

    # ---- Dispatcher route exercise ----
    print(f"\n=== DISPATCHER ROUTE EXERCISE ===")
    route_keys = [
        # 10 closure-bearing
        "paper1138_holmlid_parkhomov_pf_probe",
        "paper1139_pons_fleischmann_derivation_probe",
        "paper1140_mizuno_lenr_transmutation_probe",
        "paper1141_rossi_ecat_variants_probe",
        "paper1142_polyakov_action_26d_probe",
        "paper1153_primordial_timing_function_probe",
        "paper1154_ssq_first_principles_probe",
        "paper1155_dpm_26layer_particle_masses_probe",
        "paper1156_cosmological_constant_closure_probe",
        "paper1157_h0_anchor_asymmetry_probe",
        # 16 framework-only
        "paper1143_nambu_goto_bosonic_string_probe",
        "paper1144_type_iib_superstring_probe",
        "paper1145_type_iia_superstring_probe",
        "paper1146_heterotic_string_gauge_probe",
        "paper1147_calabi_yau_3fold_probe",
        "paper1148_m_theory_unification_probe",
        "paper1149_psz2_g181_stroe2025_probe",
        "paper1150_june20_10system_chandra_probe",
        "paper1151_vds_dvp_bh26_variant_probe",
        "paper1152_qcalcgeom_simengine_probe",
        "paper1158_overdetermination_epistemology_probe",
        "paper1159_phi_res_codimension_probe",
        "paper1160_f_trz_so5_probe",
        "paper1161_26_factorial_pochhammer_probe",
        "paper1162_kk_tower_mode_by_mode_probe",
        "paper1163_dpm_so2_lightcone_probe",
        # manifest
        "image_session271_manifest",
    ]
    failures = []
    for k in route_keys:
        try:
            result = upc.calculate_analytic_closures({"calc": k})
            if result is None:
                failures.append((k, "None"))
            else:
                v = result.get("value") if isinstance(result, dict) else result
                if isinstance(v, dict):
                    closure_n = v.get("closures") if "closures" in v else "framework"
                else:
                    closure_n = "n/a"
                print(f"  OK  {k}  -> closures={closure_n}")
        except Exception as e:
            failures.append((k, str(e)))
            print(f"  FAIL {k}: {e}")

    # ---- S### lookup spot test ----
    print(f"\n=== S### LOOKUP SPOT TEST ===")
    spot = ["S1700", "S1723", "S2005", "S2020", "S2042"]
    for s in spot:
        try:
            r = upc.calculate_analytic_closures({"calc": "closure_s271_lookup", "s_id": s})
            v = r.get("value") if isinstance(r, dict) else None
            if isinstance(v, dict):
                print(f"  {s}  paper={v.get('paper_tag'):10s}  label={v.get('label')}  "
                      f"derived={v.get('derived')}  paper_pred={v.get('paper_predicted')}  "
                      f"unit={v.get('unit')}")
            else:
                print(f"  {s}  -> {r}")
        except Exception as e:
            print(f"  {s} FAIL: {e}")

    print(f"\n=== SUMMARY ===")
    print(f"  dispatcher routes tested : {len(route_keys)}")
    print(f"  failures                 : {len(failures)}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
