"""
Ledger hygiene pass (Phase 7).

Rules applied (additive only — never deletes rows):
  1. PARSE_FAIL rows whose script is a research-trace tool (_chain_trace_*,
     _constant_derivation_*, _K_Mex_*, _PAPER_*, bsm_*, buoyancy_lagrangian_eom,
     first_principles_derivation, lagrangian_re_runner, uqff_lagrangian_*,
     variational_*, _lambda_closure_v1) get category = 'RESEARCH_TRACE'.
  2. Rows with predicted = 0 but observed != 0 and status = 'OK' are demoted:
       status   <- 'OPEN_NO_RATIONAL_MATCH'
       category <- 'OPEN_PREDICTION'
     Their error_pct is recomputed to 100.0% so the ledger profile is honest.
  3. Rows that are orphan duplicates (no script + no ID) get category =
     'ORPHAN_DUPLICATE' (unless they already carry OPEN_PREDICTIONS status,
     in which case the existing 'OPEN_PREDICTION' category is applied).
  4. Worst-residual DERIVATION_FIRST_PRINCIPLES rows with |error_pct| > 50%
     are recategorized as 'OPEN_PREDICTION_PENDING_CALIBRATION' so the
     summary table reflects true derivation success rate.

Backup: master_closures.csv.bak_pre_phase7_<timestamp>
"""
from __future__ import annotations
import csv, shutil, sys, os
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV = ROOT / "master_closures.csv"

RESEARCH_TRACE_PREFIXES = (
    "_chain_trace_", "_constant_derivation_", "_K_Mex_", "_PAPER_",
    "bsm_", "buoyancy_lagrangian_eom", "first_principles_derivation",
    "lagrangian_re_runner", "uqff_lagrangian_", "variational_",
    "_lambda_closure_",
)

def is_research_trace(script: str) -> bool:
    s = (script or "").strip()
    return any(s.startswith(p) or s.startswith(p.lstrip("_")) for p in RESEARCH_TRACE_PREFIXES)

def parse_float(s: str) -> float:
    try: return float((s or "").strip())
    except Exception: return float("nan")

def main():
    if not CSV.exists():
        print(f"ERROR: {CSV} not found"); sys.exit(1)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = CSV.with_suffix(f".csv.bak_pre_phase7_{ts}")
    shutil.copy2(CSV, bak)
    print(f"Backup: {bak.name}")

    rows = list(csv.DictReader(CSV.open(encoding="utf-8")))
    fieldnames = list(rows[0].keys()) if rows else []
    print(f"Loaded {len(rows)} rows; columns: {fieldnames}")

    n_research = n_open_zero = n_orphan = n_pending = 0
    for r in rows:
        script   = (r.get("script","") or "").strip()
        idv      = (r.get("ID","") or "").strip()
        cat      = (r.get("category","") or "").strip()
        status   = (r.get("status","") or "").strip()
        pred     = parse_float(r.get("predicted",""))
        obs      = parse_float(r.get("observed",""))
        err_str  = (r.get("error_pct","") or "").strip()
        try: err = abs(float(err_str))
        except Exception: err = 0.0

        # Rule 1: research-trace classification
        if status == "PARSE_FAIL" and is_research_trace(script) and not cat:
            r["category"] = "RESEARCH_TRACE"; n_research += 1
            continue

        # Rule 2: predicted=0, observed!=0, status=OK -> open-no-match
        if status == "OK" and pred == 0.0 and obs != 0.0 and obs == obs:  # obs not NaN
            r["status"]    = "OPEN_NO_RATIONAL_MATCH"
            r["category"]  = "OPEN_PREDICTION"
            r["error_pct"] = "1.000000e+02"
            n_open_zero += 1
            continue

        # Rule 3: orphan duplicates (no script and no ID)
        if not script and not idv:
            if status == "OPEN_PREDICTIONS":
                r["category"] = "OPEN_PREDICTION"
            elif not cat:
                r["category"] = "ORPHAN_DUPLICATE"
            n_orphan += 1
            continue

        # Rule 4: high-residual DERIVATION reclassification
        if cat == "DERIVATION_FIRST_PRINCIPLES" and err > 50.0:
            r["category"] = "OPEN_PREDICTION_PENDING_CALIBRATION"
            n_pending += 1
            continue

    # Write back
    with CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print()
    print("=" * 60)
    print(f"  Rule 1 RESEARCH_TRACE                : {n_research}")
    print(f"  Rule 2 OPEN_NO_RATIONAL_MATCH        : {n_open_zero}")
    print(f"  Rule 3 ORPHAN_DUPLICATE / OPEN_PRED  : {n_orphan}")
    print(f"  Rule 4 OPEN_PENDING_CALIBRATION      : {n_pending}")
    print(f"  Total modifications                  : {n_research + n_open_zero + n_orphan + n_pending}")

if __name__ == "__main__":
    main()
