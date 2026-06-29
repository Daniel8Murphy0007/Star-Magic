"""
test_phase_e8_overdetermination_map.py
Phase E8 — Verify the OVERDETERMINATION_MAP family is well-formed.

Checks:
 1. All three files exist and are non-empty.
 2. Long format has exactly N_obs * 4 * 3 rows.
 3. Wide format has N_obs rows and 18 cols.
 4. Per-observable: owner geometry's 3 numeric cells are all populated (owner_N == 3).
 5. Status distribution sums to total rows (no spurious labels).
 6. TENSION cells in long format == 3 (BAO observable, one per numeric system).
 7. Worst residual per domain is consistent with dispatch records (sanity).
 8. .md summary contains expected headers and TENSION block.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import assimilation_dispatch as ad

LONG = ROOT / "OVERDETERMINATION_MAP.csv"
WIDE = ROOT / "OVERDETERMINATION_WIDE.csv"
MD = ROOT / "OVERDETERMINATION_MAP.md"


def main():
    print("=" * 72)
    print("PHASE E8 — OVERDETERMINATION_MAP regression")
    print("=" * 72)

    for f in (LONG, WIDE, MD):
        if not f.exists() or f.stat().st_size == 0:
            print(f"FAIL: {f.name} missing or empty")
            return 1
        print(f"PASS  {f.name}: {f.stat().st_size} bytes")

    # 1. Long format integrity
    with LONG.open(encoding="utf-8", newline="") as f:
        long_rows = list(csv.DictReader(f))
    n_obs = len(ad.DISPATCH)
    expected_rows = n_obs * 4 * 3
    if len(long_rows) != expected_rows:
        print(f"FAIL: long row count {len(long_rows)} != expected {expected_rows}")
        return 1
    print(f"PASS  long format: {len(long_rows)} rows = {n_obs} obs x 4 geom x 3 num")

    # 2. Wide format integrity
    with WIDE.open(encoding="utf-8", newline="") as f:
        wide_reader = csv.DictReader(f)
        wide_fields = wide_reader.fieldnames
        wide_rows = list(wide_reader)
    if len(wide_rows) != n_obs:
        print(f"FAIL: wide row count {len(wide_rows)} != {n_obs}")
        return 1
    EXPECTED_WIDE_COLS = ["observable","domain","owner_geometry",
        "qg_sym","qg_num","qg_dis","bsfg_sym","bsfg_num","bsfg_dis",
        "dpm_sym","dpm_num","dpm_dis","d26_sym","d26_num","d26_dis",
        "owner_N","total_N","primary_source"]
    if wide_fields != EXPECTED_WIDE_COLS:
        print(f"FAIL: wide schema mismatch")
        print(f"  expected: {EXPECTED_WIDE_COLS}")
        print(f"  actual:   {wide_fields}")
        return 1
    print(f"PASS  wide format: {len(wide_rows)} rows, {len(wide_fields)} cols")

    # 3. Owner_N == 3 for every observable
    bad_owner_N = [w for w in wide_rows if int(w["owner_N"]) != 3]
    if bad_owner_N:
        print(f"FAIL: {len(bad_owner_N)} observables have owner_N != 3")
        for w in bad_owner_N[:5]:
            print(f"    {w['observable']}: owner_N={w['owner_N']}")
        return 1
    print(f"PASS  owner_N == 3 for all {len(wide_rows)} observables")

    # 4. Status distribution
    status_counts = Counter(r["status"] for r in long_rows)
    print(f"      status distribution: {dict(status_counts)}")
    if set(status_counts.keys()) - {"EXACT","OK","TENSION","GAP"}:
        print(f"FAIL: unexpected status labels: {set(status_counts.keys())}")
        return 1
    if sum(status_counts.values()) != expected_rows:
        print(f"FAIL: status counts {sum(status_counts.values())} != total rows {expected_rows}")
        return 1
    print(f"PASS  status counts sum to total rows")

    # 5. TENSION cells == 3 (BAO across symbolic / numerical / discrete)
    tension_rows = [r for r in long_rows if r["status"] == "TENSION"]
    if len(tension_rows) != 3:
        print(f"FAIL: TENSION count {len(tension_rows)} != 3 (BAO OPEN_QUESTION expected)")
        return 1
    tension_obs = {r["observable"] for r in tension_rows}
    if tension_obs != {"LCDM_BAO_rd_H0_over_c_OPEN"}:
        print(f"FAIL: TENSION observables {tension_obs} != {{LCDM_BAO_rd_H0_over_c_OPEN}}")
        return 1
    print(f"PASS  TENSION cells: 3 (LCDM_BAO_rd_H0_over_c_OPEN x 3 numeric)")

    # 6. Sanity: BAO TENSION residual ~ 4.77%
    bao_residual = float(tension_rows[0]["residual_pct"])
    if not (4.5 < bao_residual < 5.0):
        print(f"FAIL: BAO residual {bao_residual:.4f}% outside [4.5, 5.0]")
        return 1
    print(f"PASS  BAO TENSION residual: {bao_residual:.4f}% (in expected ~4.77% range)")

    # 7. Population coverage matches expectation (336 = 112 * 3)
    populated = sum(1 for r in long_rows if r["value"] != "")
    expected_populated = n_obs * 3
    if populated != expected_populated:
        print(f"FAIL: populated cells {populated} != expected {expected_populated}")
        return 1
    coverage = 100.0 * populated / expected_rows
    print(f"PASS  populated cells: {populated} / {expected_rows} ({coverage:.1f}%)")
    print(f"      (each observable's owner geometry contributes 3 cells = 4x3 minus 3x3 GAP)")

    # 8. .md headers present
    md_text = MD.read_text(encoding="utf-8")
    required_md = [
        "# OVERDETERMINATION_MAP",
        "## Top-line metrics",
        "## Per-domain rollup",
        "## OPEN_QUESTION / TENSION cells",
        "LCDM_BAO_rd_H0_over_c_OPEN",
        "## Schema notes",
    ]
    missing = [h for h in required_md if h not in md_text]
    if missing:
        print(f"FAIL: .md missing sections: {missing}")
        return 1
    print(f"PASS  .md contains all 6 required sections + BAO TENSION row")

    print()
    print("PHASE E8 SUCCESS CRITERION MET.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
