"""Session 272 image-batch smoke test -- cold start.

Imports the calculator, iterates all _PA_S272_CLOSURE_REGISTRY rows,
exercises every dispatcher route, and prints a per-paper summary.

No SM literals, no pre-fit corrections. form() outputs are the
falsifiable UQFF predictions. paper_predicted = user-published target.
form_transcribed=False entries return paper_predicted verbatim
(used when paper formula references primitives outside the locked
vacuum-ledger, including paper-stated [SSq]=0.57 vs locked SSQ=0.505).
"""

from __future__ import annotations
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import uqff_pure_calculator as upc  # noqa: E402


def main() -> int:
    reg = upc._PA_S272_CLOSURE_REGISTRY
    print(f"\n=== SESSION 272 IMAGE BATCH (PAPER_1112..PAPER_1137) ===")
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
        print(f"  {paper_tag:14s}  closures={len(rows):3d}  derived={derived_ok:3d}  "
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
        # 25 closure-bearing
        "paper1113_cms_higgs_kappa_probe",
        "paper1114_atlas_higgs_width_probe",
        "paper1115_scs_21cm_probe",
        "paper1116_electroweak_axion_probe",
        "paper1117_scs_radio_frb_probe",
        "paper1118_chiral_graphene_probe",
        "paper1119_lorentz_heaviside_probe",
        "paper1120_higgs_decay_modes_probe",
        "paper1121_interstellar_shock_probe",
        "paper1122_bow_shock_ism_probe",
        "paper1123_h2o_maser_probe",
        "paper1124_cgm_dwarf_probe",
        "paper1125_agn_msigma_probe",
        "paper1126_psr_j0030_lenr_probe",
        "paper1127_lqg_holonomy_probe",
        "paper1128_string_theory_26d_probe",
        "paper1129_vds_dvp_bh_probe",
        "paper1130_26d_folding_probe",
        "paper1131_vacuum_manifold_primordial_probe",
        "paper1132_primordial_split_ladder_probe",
        "paper1133_holmlid_rydberg_probe",
        "paper1134_riemann_closure_probe",
        "paper1135_hub_reactor_probe",
        "paper1136_ker_reactor_probe",
        "paper1137_rossi_parkhomov_probe",
        # 1 framework-only
        "paper1112_production_scaling_probe",
        # manifest
        "image_session272_manifest",
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

    # ---- Direct S### lookup spot test (bypasses dispatcher) ----
    print(f"\n=== DIRECT S### LOOKUP SPOT TEST ===")
    spot = ["S2107", "S2131", "S2306", "S2401", "S2522", "S2509"]
    for s in spot:
        try:
            v = upc._l96_uqff_closure_S272_lookup(s)
            print(f"  {s}  paper={v.get('paper_tag'):12s}  label={v.get('label'):44s}  "
                  f"derived={v.get('derived')}  paper_pred={v.get('paper_predicted')}  "
                  f"unit={v.get('unit')}")
        except Exception as e:
            print(f"  {s} FAIL: {e}")

    print(f"\n=== SUMMARY ===")
    print(f"  dispatcher routes tested : {len(route_keys)}")
    print(f"  failures                 : {len(failures)}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
