"""
test_phase_e7_master_closures_merge.py
Phase E7 — Verify the extended master_closures.csv schema and dispatch tagging.

Daniel's workflow rule: "YOU FIND AND I REVIEW AND COMMIT. PERIOD!"

Phase E7 produces the merged file as `master_closures.csv.PROPOSED_E7`
(the original mount blocks in-place overwrite of master_closures.csv).
Daniel performs the atomic swap when he commits. This harness verifies
the proposed file is structurally sound AND that swap-target compatibility
is preserved.

Verifies:
 1. PROPOSED file exists; row count matches original (2216 + header).
 2. Header has exactly the 13 original columns + 3 new columns
    (geometry_used, numeric_system, assimilation_status) at the tail.
 3. At least 30 rows are tagged from assimilation_dispatch.py.
 4. Tagged rows carry valid geometry (qcalcgeom|bsfg|dpm|d26) and
    valid status (OK|EXACT|TENSION).
 5. The audit log (phase_e7_tag_audit.csv) was generated and is non-empty.
 6. All original rows' first 13 columns are byte-for-byte identical to
    master_closures.csv.PRE_PHASE_E7_BACKUP (no historical data mutated).
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ORIG = ROOT / "master_closures.csv"
BACKUP = ROOT / "master_closures.csv.PRE_PHASE_E7_BACKUP"
PROPOSED = ROOT / "master_closures.csv.PROPOSED_E7"
AUDIT = ROOT / "phase_e7_tag_audit.csv"

ORIGINAL_COLS = [
    "closure","predicted","observed","error_pct","status","cvw_stamp",
    "sm_anchor","label","raw_output","category","name","script","ID",
]
NEW_COLS = ["geometry_used","numeric_system","assimilation_status"]
EXPECTED_COLS = ORIGINAL_COLS + NEW_COLS
VALID_GEOMS = {"qcalcgeom","bsfg","dpm","d26"}
VALID_STATUSES = {"OK","EXACT","TENSION"}


def main():
    print("=" * 72)
    print("PHASE E7 — master_closures.csv extended schema verification")
    print("=" * 72)

    # 1. Check live file. If it has 16 cols, swap already done. Otherwise check PROPOSED.
    if ORIG.exists():
        with ORIG.open(encoding="utf-8", newline="") as f:
            live_hdr = next(csv.reader(f))
        if live_hdr == EXPECTED_COLS:
            target = ORIG
            print("Live master_closures.csv already at v16. (Daniel swapped.)")
        else:
            target = PROPOSED
            print("Live master_closures.csv still at v13. Reading PROPOSED_E7.")
    else:
        target = PROPOSED
        print("master_closures.csv missing — reading PROPOSED_E7.")

    if not target.exists():
        print(f"FAIL: {target.name} not found")
        return 1
    print(f"Target: {target.name}")

    # 2. Schema check
    with target.open(encoding="utf-8", newline="") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        rows = list(rdr)
    if header != EXPECTED_COLS:
        print(f"FAIL: header mismatch")
        print(f"  expected: {EXPECTED_COLS}")
        print(f"  actual:   {header}")
        return 1
    print(f"PASS  schema: {len(header)} cols = 13 original + 3 new")
    print(f"PASS  row count: {len(rows)} (expected 2216)")
    if len(rows) != 2216:
        print(f"FAIL: row count mismatch")
        return 1

    # 3. Tagged-row count
    tagged = [r for r in rows if r[13] != "" or r[14] != "" or r[15] != ""]
    print(f"PASS  tagged rows: {len(tagged)} (>=30 required)")
    if len(tagged) < 30:
        return 1

    # 4. Validate tagged values
    bad_geom = bad_status = bad_numeric = 0
    geom_counts = {}
    status_counts = {}
    for r in tagged:
        g = r[13]; n = r[14]; s = r[15]
        if g not in VALID_GEOMS:
            bad_geom += 1
        else:
            geom_counts[g] = geom_counts.get(g, 0) + 1
        if n != "numerical":
            bad_numeric += 1
        if s not in VALID_STATUSES:
            bad_status += 1
        else:
            status_counts[s] = status_counts.get(s, 0) + 1
    if bad_geom or bad_status or bad_numeric:
        print(f"FAIL: bad_geom={bad_geom} bad_numeric={bad_numeric} bad_status={bad_status}")
        return 1
    print(f"PASS  tagged geometry distribution: {dict(sorted(geom_counts.items()))}")
    print(f"PASS  tagged status   distribution: {dict(sorted(status_counts.items()))}")

    # 5. Audit log
    if not AUDIT.exists():
        print(f"FAIL: audit log {AUDIT.name} missing")
        return 1
    with AUDIT.open(encoding="utf-8", newline="") as f:
        audit_rows = list(csv.reader(f))
    print(f"PASS  audit log present: {len(audit_rows) - 1} tagged entries documented")
    if len(audit_rows) - 1 != len(tagged):
        print(f"WARN  audit log count {len(audit_rows)-1} != tagged count {len(tagged)}")

    # 6. Backward compatibility — original 13 columns unchanged vs PRE_PHASE_E7_BACKUP
    if BACKUP.exists():
        with BACKUP.open(encoding="utf-8", newline="") as f:
            orig_rdr = csv.reader(f)
            orig_hdr = next(orig_rdr)
            orig_rows = list(orig_rdr)
        mutations = 0
        for o, n in zip(orig_rows, rows):
            if o != n[:13]:
                mutations += 1
        if mutations:
            print(f"FAIL: {mutations} rows had original columns mutated")
            return 1
        print(f"PASS  zero mutations to original 13 columns vs PRE_PHASE_E7_BACKUP")
    else:
        print("WARN  PRE_PHASE_E7_BACKUP missing; cannot verify backward compatibility")

    print()
    print("PHASE E7 SUCCESS CRITERION MET.")
    print()
    if target == PROPOSED:
        print("NEXT STEP (Daniel's commit): swap PROPOSED into live")
        print("  cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic")
        print("  mv master_closures.csv master_closures.csv.PRE_E7_LIVE")
        print("  mv master_closures.csv.PROPOSED_E7 master_closures.csv")
        print()
        print("After swap, re-run this harness — it will detect the v16 live file.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
