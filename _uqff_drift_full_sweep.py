#!/usr/bin/env python3
"""
_uqff_drift_full_sweep.py  --  CI drift sweep (Session 256)

Verifies that all Session 246-256 closures still produce in-tolerance results.
Exits 0 on full pass, 1 on any failure. Designed for unattended CI invocation.

Coverage:
  - uqff_closed_constants assertions (K_MEXICAN_HAT, PHI_RES, F_TRZ, ...)
  - CP4 calculators #253-#257 smoke tests
  - QCalcGeom test runner (_qcalcgeom_tests.exe) -> all tests PASS expected (version-agnostic)
"""
from __future__ import annotations

import os
import subprocess
import sys
import traceback


FAILED: list[str] = []


def _ok(label: str) -> None:
    print(f"  [PASS] {label}")


def _fail(label: str, detail: str = "") -> None:
    FAILED.append(label)
    print(f"  [FAIL] {label}  {detail}")


def check_constants() -> None:
    print("\n[1/3] uqff_closed_constants drift check")
    try:
        import uqff_closed_constants as U  # type: ignore
        expect = {
            "K_MEXICAN_HAT":   25.0 / 12.0,
            "PHI_RES":         5.0 / 6.0,
            "F_TRZ":           1.0 / 10.0,
            "D_CRIT":          26,
            "D_PHYS":          4,
            "D_BSFG":          6,
            "DIM_SO5":         10,
            "RHO_SCM_DEFAULT": 7.09e-37,
            "V_UA_DEFAULT":    1.0e8,
        }
        for k, v in expect.items():
            actual = getattr(U, k, None)
            if actual is None:
                _fail(f"const::{k}", "missing")
            elif abs(float(actual) - float(v)) > 1e-12 * max(1.0, abs(float(v))):
                _fail(f"const::{k}", f"got {actual} expected {v}")
            else:
                _ok(f"const::{k} = {actual}")
    except Exception:
        _fail("import uqff_closed_constants", traceback.format_exc().splitlines()[-1])


def check_cp4_calculators() -> None:
    print("\n[2/3] CP4 #253-#257 smoke tests")
    cases = [
        ("UQFFVacuumEnergyLedgerCalculator", "rho_Lambda_closed", 5.96e-10, 0.005),
        ("UQFFKKTowerRegulatorCalculator",   "rho_KK",            5.86e-10, 0.005),
    ]
    try:
        import CondensedPhysics4 as CP4  # type: ignore
    except Exception:
        _fail("import CondensedPhysics4", traceback.format_exc().splitlines()[-1])
        return

    for cls_name, key, expect, tol in cases:
        try:
            cls = getattr(CP4, cls_name)
            r   = cls().compute()
            v   = r.get(key)
            if v is None:
                _fail(f"{cls_name}.{key}", "missing key")
                continue
            if abs(v - expect) / expect <= tol:
                _ok(f"{cls_name}.{key} = {v:.4e} (target {expect:.4e}, tol {tol*100:.2f}%)")
            else:
                _fail(f"{cls_name}.{key}", f"{v:.4e} vs {expect:.4e}")
        except Exception:
            _fail(cls_name, traceback.format_exc().splitlines()[-1])


def check_qcalcgeom() -> None:
    print("\n[3/3] QCalcGeom test runner (_qcalcgeom_tests.exe)")
    exe = os.path.join(os.path.dirname(__file__), "_qcalcgeom_tests.exe")
    if not os.path.exists(exe):
        _fail("_qcalcgeom_tests.exe", "not built (run g++ build first)")
        return
    try:
        out = subprocess.run([exe], capture_output=True, text=True, timeout=120)
        tail = out.stdout.splitlines()[-12:]
        for line in tail:
            print(f"    {line}")
        # Version-agnostic: accept any N/N pass count, fail only on any FAILED>0
        import re as _re
        m_total = _re.search(r"Total tests:\s+(\d+)", out.stdout)
        m_pass  = _re.search(r"PASSED:\s+(\d+)", out.stdout)
        m_fail  = _re.search(r"FAILED:\s+(\d+)", out.stdout)
        if m_total and m_pass and m_fail and m_fail.group(1) == "0" and m_pass.group(1) == m_total.group(1):
            _ok(f"QCalcGeom {m_pass.group(1)}/{m_total.group(1)} PASS")
        else:
            _fail("QCalcGeom", f"unexpected output (pass={m_pass and m_pass.group(1)}, fail={m_fail and m_fail.group(1)}, total={m_total and m_total.group(1)})")
    except Exception:
        _fail("QCalcGeom run", traceback.format_exc().splitlines()[-1])


def main() -> int:
    print("=" * 70)
    print(" UQFF Drift Full Sweep  --  Session 256")
    print("=" * 70)
    check_constants()
    check_cp4_calculators()
    check_qcalcgeom()
    print("\n" + "=" * 70)
    if FAILED:
        print(f" RESULT: {len(FAILED)} FAILURE(S)")
        for f in FAILED:
            print(f"   - {f}")
        return 1
    print(" RESULT: ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
