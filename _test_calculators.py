#!/usr/bin/env python3
"""
_test_calculators.py — Comprehensive test runner for ALL CP calculators.

Tests every calculator class in CP1 / CP2 / CP3 / CP4 by:
  1. Extracting class names via AST (no import required)
  2. Importing the module
  3. Instantiating each class (no-arg or with dataset fallback)
  4. Calling compute(dataset=TEST_DATASET)
  5. Reporting pass/fail counts + writing _test_results.json

Usage:
    python _test_calculators.py

Output:
    _test_results.json  — full pass/fail report
    Console             — per-file + grand-total summary
"""

import sys
import os
import ast
import inspect
import importlib
import json
import traceback
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

# ── ensure the repo root is on the path ──────────────────────────────────────
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# ── Standard UQFF test dataset (covers the widest range of calculator inputs) ─
TEST_DATASET: Dict[str, Any] = {
    # Core body parameters
    "M":           1.989e30,      # mass [kg] (solar mass)
    "mass":        1.989e30,
    "r":           6.96e8,        # radius [m] (solar radius)
    "radius":      6.96e8,
    "R":           6.96e8,
    "distance":    1.496e11,      # 1 AU [m]
    "d":           1.496e11,
    # Electromagnetic
    "B0":          1e-4,          # surface magnetic field [T]
    "B":           1e-4,
    "omega_s":     2.5e-6,        # angular velocity [rad/s]
    "omega":       2.5e-6,
    "omega0":      2.5e-6,
    "v":           1e5,           # velocity [m/s]
    "v_sw":        4e5,           # solar wind velocity [m/s]
    # Thermal
    "T":           5778.0,        # temperature [K]
    "temperature": 5778.0,
    "T_eff":       5778.0,
    # Observational
    "luminosity":  3.828e26,      # [W]
    "L":           3.828e26,
    "z":           0.0,           # redshift
    "SFR":         1.0,           # star formation rate [M_sun/yr]
    "flux":        1e-10,         # flux [W/m²]
    "alpha":       1.7,           # spectral index
    # Density / pressure
    "rho":         1e-20,         # density [kg/m³]
    "rho_vac":     7.09e-37,      # SCm vacuum density [J/m³]
    "rho_sw":      1e-20,         # solar wind density
    "n_e":         1e6,           # electron number density [m⁻³]
    "P":           1e5,           # pressure [Pa]
    "E":           1e40,          # energy [J]
    "E_react":     1e46,          # reactor efficiency [W/m³]
    # Timing / angles
    "t":           0.0,           # time [s]
    "tn":          0.0,           # normalized time
    "theta":       0.1,           # angle [rad]
    "phi":         0.0,
    # UQFF calibrated constants
    "kappa":       0.0005,        # E_react decay constant [day⁻¹]
    "beta_i":      0.603,         # buoyancy coupling
    "SSq":         0.57,          # self-similar quotient
    "H_SCm":       0.99,
    "k_eta":       1e-113,
    "U_UA":        1e-4,
    # Galactic parameters
    "Mbh":         4.15e6 * 1.989e30,   # Sgr A* mass [kg]
    "M_bh":        4.15e6 * 1.989e30,
    "Omega_g":     7.3e-16,       # galactic spin [rad/s]
    "d_g":         2.55e20,       # galactic distance [m]
    # Quantum / atomic
    "n":           1,             # principal quantum number
    "l":           0,             # orbital quantum number
    "frequency":   1e6,           # [Hz]
    "wavelength":  3e2,           # [m]
    "Z":           1,             # atomic number
    # Additional fields many calculators use
    "gamma":       1.33,          # adiabatic index / Lorentz factor
    "Gamma":       1.5,           # Lorentz factor
    "chi":         0.5,
    "epsilon":     0.1,
    "delta":       0.01,
    "sigma":       5.67e-8,       # Stefan-Boltzmann [W/m²K⁴]
    "metallicity": 0.014,
    "age":         4.6e9,         # stellar age [yr]
    "spin":        0.5,           # dimensionless spin
    "name":        "Sun",
    # Extended UQFF fields
    "Ug1":         1e-10,
    "Ug2":         1e-10,
    "Ug3":         1e-10,
    "Ug4":         1e-10,
    "Ubi":         1e-10,
    "Um":          1e-10,
    "F_U":         1e-10,
    "C_concentration": 1.0,
    "delta_def":   0.01,
    "k1":          1e-12,
    "k2":          1e-12,
    "k3":          1e-12,
    "k4":          1e-12,
    "QA":          1.0,
    "delta_sw":    0.01,
    "HSCm":        0.99,
    "rho_A":       7.09e-37,
    # DPM
    "I_dpm":       1e6,           # DPM current [A]
    "A_dpm":       1e-4,          # DPM area [m²]
    "omega1":      2.5e-6,
    "omega2":      2.4e-6,
    "f_dpm":       1e6,           # DPM frequency [Hz]
    "V_sys":       1e27,          # system volume [m³]
    "E_vac_neb":   7.09e-37,
    # Misc
    "mu":          1e-6,
    "c":           2.998e8,
    "G":           6.674e-11,
    "h":           6.626e-34,
    "hbar":        1.055e-34,
    "k_B":         1.381e-23,
    "N":           1,             # count
    "x":           1.0,
    "y":           0.0,
    "p":           1e-10,
    "q":           1.0,
}

# ── CP file definitions ────────────────────────────────────────────────────────
CP_FILES = [
    ("CP1", "CondensedPhysics",  "CondensedPhysics.py"),
    ("CP2", "CondensedPhysics2", "CondensedPhysics2.py"),
    ("CP3", "CondensedPhysics3", "CondensedPhysics3.py"),
    ("CP4", "CondensedPhysics4", "CondensedPhysics4.py"),
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_class_names_ast(filepath: str) -> List[str]:
    """Return all top-level class names defined in *filepath* via AST (no import)."""
    with open(filepath, encoding="utf-8-sig", errors="ignore") as fh:
        content = fh.read()
    tree = ast.parse(content)
    return [
        node.name
        for node in ast.walk(tree)
        if isinstance(node, ast.ClassDef) and node.col_offset == 0
    ]


def try_instantiate(cls: type) -> Tuple[Optional[object], Optional[str]]:
    """Try to create an instance of *cls*.  Returns (obj, error_msg)."""
    # Attempt 1: no-arg constructor
    try:
        return cls(), None
    except TypeError:
        pass
    except Exception as e:
        return None, f"__init__() raised {type(e).__name__}: {str(e)[:200]}"

    # Attempt 2: pass TEST_DATASET as first positional arg
    try:
        return cls(TEST_DATASET), None
    except Exception as e:
        return None, f"__init__(dataset) raised {type(e).__name__}: {str(e)[:200]}"


def test_class(cls: type) -> Tuple[bool, str]:
    """Instantiate *cls* and call compute(dataset=TEST_DATASET).
    Returns (passed: bool, error_message: str)."""

    obj, err = try_instantiate(cls)
    if obj is None:
        return False, err or "instantiation failed (unknown)"

    if not hasattr(obj, "compute"):
        return False, "no compute() method"

    # Attempt 1: compute(dataset=TEST_DATASET)
    try:
        obj.compute(dataset=TEST_DATASET)
        return True, ""
    except TypeError:
        pass
    except Exception as e:
        return False, f"compute(dataset=...) raised {type(e).__name__}: {str(e)[:200]}"

    # Attempt 2: compute(TEST_DATASET)  (positional)
    try:
        obj.compute(TEST_DATASET)
        return True, ""
    except TypeError:
        pass
    except Exception as e:
        return False, f"compute(dataset) raised {type(e).__name__}: {str(e)[:200]}"

    # Attempt 3: compute()  (no args, some calculators ignore dataset)
    try:
        obj.compute()
        return True, ""
    except Exception as e:
        return False, f"compute() raised {type(e).__name__}: {str(e)[:200]}"


# ── Main test runner ──────────────────────────────────────────────────────────

def run_tests() -> Tuple[int, int]:
    start_time = datetime.now()
    print(f"\n{'='*70}")
    print(f"  UQFF Calculator Test Suite  —  {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}")

    per_file_results: Dict[str, Any] = {}
    grand_pass = 0
    grand_fail = 0
    grand_skip = 0

    for tag, mod_name, filename in CP_FILES:
        filepath = os.path.join(REPO_ROOT, filename)

        print(f"\n{'─'*70}")
        print(f"  {tag}  ({filename})")
        print(f"{'─'*70}")

        # ── Step 1: AST class discovery ───────────────────────────────────────
        print(f"  [1/3] AST scan ... ", end="", flush=True)
        try:
            ast_names = get_class_names_ast(filepath)
        except Exception as e:
            print(f"FAILED: {e}")
            per_file_results[tag] = {"error": f"AST failed: {e}"}
            grand_fail += 1
            continue
        print(f"{len(ast_names)} class names found")

        # ── Step 2: Module import ─────────────────────────────────────────────
        print(f"  [2/3] Import module ... ", end="", flush=True)
        try:
            mod = importlib.import_module(mod_name)
            print(f"OK")
        except Exception as e:
            msg = f"import failed: {type(e).__name__}: {str(e)[:300]}"
            print(f"FAILED\n       {msg}")
            per_file_results[tag] = {
                "error": msg,
                "ast_count": len(ast_names),
            }
            grand_fail += len(ast_names)
            continue

        # ── Step 3: Resolve classes ───────────────────────────────────────────
        resolved: Dict[str, type] = {}
        skipped: List[str] = []
        for name in ast_names:
            obj = getattr(mod, name, None)
            if obj is not None and inspect.isclass(obj):
                resolved[name] = obj
            else:
                skipped.append(name)

        print(f"  [3/3] Resolved {len(resolved)}/{len(ast_names)} class objects "
              f"({len(skipped)} not found in module namespace)")

        # ── Step 4: Test each class ───────────────────────────────────────────
        print(f"  Testing {len(resolved)} classes ...", flush=True)
        file_pass = 0
        file_fail = 0
        failures: List[Dict[str, str]] = []

        total = len(resolved)
        items = list(resolved.items())
        for idx, (cls_name, cls) in enumerate(items, 1):
            # Progress every 50 classes
            if idx % 50 == 0:
                pct = idx * 100 // total
                print(f"    [{idx:>4}/{total}  {pct:>3}%]  pass={file_pass}  fail={file_fail}", flush=True)

            ok, err = test_class(cls)
            if ok:
                file_pass += 1
            else:
                file_fail += 1
                failures.append({"class": cls_name, "error": err})

        # Print final progress line
        print(f"    [{total:>4}/{total}  100%]  pass={file_pass}  fail={file_fail}")

        # ── Summary for this file ─────────────────────────────────────────────
        pass_rate = file_pass * 100 // max(total, 1)
        print(f"\n  RESULT  {tag}: {file_pass} PASS  {file_fail} FAIL  "
              f"({pass_rate}% pass rate)  {len(skipped)} skipped (not in namespace)")

        if failures:
            print(f"\n  FAILURES ({len(failures)} total):")
            for f in failures[:30]:
                print(f"    ✗  {f['class']}")
                print(f"       {f['error']}")
            if len(failures) > 30:
                print(f"    ... and {len(failures) - 30} more — see _test_results.json")

        per_file_results[tag] = {
            "ast_count": len(ast_names),
            "resolved": len(resolved),
            "skipped_from_namespace": skipped,
            "pass": file_pass,
            "fail": file_fail,
            "pass_rate_pct": pass_rate,
            "failures": failures,
        }
        grand_pass += file_pass
        grand_fail += file_fail
        grand_skip += len(skipped)

    # ── Grand total ───────────────────────────────────────────────────────────
    elapsed = (datetime.now() - start_time).total_seconds()
    total_tested = grand_pass + grand_fail
    grand_pct = grand_pass * 100 // max(total_tested, 1)

    print(f"\n{'='*70}")
    print(f"  GRAND TOTAL:")
    print(f"    Tested  : {total_tested}")
    print(f"    PASS    : {grand_pass}  ({grand_pct}%)")
    print(f"    FAIL    : {grand_fail}")
    print(f"    Skipped : {grand_skip}  (classes in AST but not in module namespace)")
    print(f"    Elapsed : {elapsed:.1f}s")
    print(f"{'='*70}\n")

    # ── Write JSON report ─────────────────────────────────────────────────────
    report = {
        "timestamp":    start_time.isoformat(),
        "elapsed_sec":  elapsed,
        "grand_total":  {"pass": grand_pass, "fail": grand_fail,
                         "tested": total_tested, "pass_rate_pct": grand_pct},
        "per_file":     per_file_results,
    }
    report_path = os.path.join(REPO_ROOT, "_test_results.json")
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(f"  Full report → {report_path}\n")

    return grand_pass, grand_fail


if __name__ == "__main__":
    p, f = run_tests()
    sys.exit(0 if f == 0 else 1)
