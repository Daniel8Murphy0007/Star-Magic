"""Session 270 image-batch smoke test — cold start.

Imports the calculator, iterates all _PA_S270_CLOSURE_REGISTRY rows,
exercises every dispatcher route, and prints a per-paper summary.

No SM literals, no pre-fit corrections. form() outputs are the
falsifiable UQFF predictions. paper_predicted = user-published target.
form_transcribed=False entries return paper_predicted verbatim.
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import uqff_pure_calculator as upc  # noqa: E402


def main() -> int:
    reg = upc._PA_S270_CLOSURE_REGISTRY
    print(f"\n=== SESSION 270 IMAGE BATCH (PAPER_1164..PAPER_1185) ===")
    print(f"Registry size: {len(reg)} closures across image batch")
    print(f"=========================================================\n")

    # Group by paper_tag
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
    rows_within_1pct_obs = 0

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
                    if predicted == 0.0:
                        err = 0.0 if d == 0.0 else float("inf")
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
                    if observed == 0.0:
                        err = 0.0 if d == 0.0 else float("inf")
                    else:
                        err = abs(d - observed) / abs(observed) * 100.0
                    if err < 5.0:
                        obs_ok += 1
                        rows_within_1pct_obs += 1
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
    print(f"  agree with observed <5pct    : {rows_within_1pct_obs}")

    # ---- Dispatcher route exercise ----
    print(f"\n=== DISPATCHER ROUTE EXERCISE ===")
    route_keys = [
        "paper1164_t22_moduli_stabilization_probe",
        "paper1165_beta_i_triangular_probe",
        "paper1166_v_ua_polynomial_probe",
        "paper1168_falsifiable_predictions_probe",
        "paper1169_numerical_confrontation_probe",
        "paper1170_vacuum_energy_ledger_probe",
        "paper1171_kk_regulator_probe",
        "paper1172_R26_gauss_bonnet_probe",
        "paper1173_kk_tower_hbar_probe",
        "paper1174_closed_ledger_p6_p10_probe",
        "paper1175_ligo_o5_ringdown_probe",
        "paper1176_euclid_sigma8_probe",
        "paper1178_desi_y5_w_probe",
        "paper1180_cmb_s4_mu_distortion_probe",
        "paper1181_grand_unification_probe",
        "paper1182_millennium_prize_probe",
        "paper1184_chandra_bridge_probe",
        "paper1184_open_problems_probe",
        "paper1185_neutrino_gw_probe",
        "paper1185_quantum_gravity_probe",
        "paper1167_lagrangian_master_synthesis_probe",
        "paper1177_joint_falsifier_triple_probe",
        "paper1179_quadruple_falsifier_probe",
        "paper1181_gap_verification_jobb_probe",
        "paper1183_first_principles_variational_probe",
        "paper1183_aggressive_paradox_probe",
        "image_session270_manifest",
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
    spot = ["S1100", "S271", "S301", "S1520", "S1540"]
    for s in spot:
        try:
            r = upc.calculate_analytic_closures({"calc": "closure_s270_lookup", "s_id": s})
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
