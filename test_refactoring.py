#!/usr/bin/env python3
"""Quick test of QCalc.py refactoring - verify violations removed and errors raised correctly."""

from QCalc import CONSTANTS, ComputeParams, EnhancedBuoyancyCalculator
from QCalc_validation import ReferenceSystemLibrary

print("=" * 80)
print("QCalc.py Refactoring Verification Test")
print("=" * 80)

# Test 1: Verify violations removed from CONSTANTS
print("\n[TEST 1] Verify violations removed from CONSTANTS")
print("-" * 80)
try:
    assert 'M_bh_SgrA' not in CONSTANTS, "ERROR: M_bh_SgrA still in CONSTANTS"
    assert 'd_g_SunSgrA' not in CONSTANTS, "ERROR: d_g_SunSgrA still in CONSTANTS"
    print(f"✓ CONSTANTS dictionary clean ({len(CONSTANTS)} entries, violations removed)")
except AssertionError as e:
    print(f"✗ {e}")

# Test 2: Verify ReferenceSystemLibrary exists and has correct values
print("\n[TEST 2] Verify ReferenceSystemLibrary created correctly")
print("-" * 80)
try:
    sgr = ReferenceSystemLibrary.SGR_A_STAR
    assert sgr.M_bh == 8.15e36, f"ERROR: Sgr A* M_bh incorrect: {sgr.M_bh}"
    assert sgr.d_g == 2.44e20, f"ERROR: Sgr A* d_g incorrect: {sgr.d_g}"
    print(f"✓ Sgr A*: M_bh = {sgr.M_bh:.3e} kg, d_g = {sgr.d_g:.3e} m")
    print(f"✓ Source: {sgr.source}")
    
    sun = ReferenceSystemLibrary.SUN
    assert sun.M == 1.989e30, f"ERROR: Sun M incorrect: {sun.M}"
    print(f"✓ Sun: M = {sun.M:.3e} kg")
    
    sgr1745 = ReferenceSystemLibrary.SGR_1745
    assert sgr1745.B == 4.4e13, f"ERROR: SGR 1745 B incorrect: {sgr1745.B}"
    print(f"✓ SGR 1745: B = {sgr1745.B:.3e} T")
    
except Exception as e:
    print(f"✗ {e}")

# Test 3: Verify EnhancedBuoyancyCalculator raises errors when M_bh/d_g missing
print("\n[TEST 3] Verify proper errors raised when M_bh/d_g missing")
print("-" * 80)
try:
    calc = EnhancedBuoyancyCalculator()
    params = ComputeParams(M=1e30, r=1e10, t=0)
    
    # Should raise ValueError
    try:
        calc.compute_results(params, {'Ug1': 1e-9})
        print("✗ ERROR: Should have raised ValueError for missing M_bh")
    except ValueError as e:
        if 'M_bh required' in str(e):
            print(f"✓ Correctly raised ValueError for missing M_bh")
            print(f"  Message: {str(e)[:100]}...")
        else:
            print(f"✗ Wrong error message: {e}")
            
except Exception as e:
    print(f"✗ Unexpected error: {e}")

# Test 4: Verify calculator works WITH proper M_bh/d_g params
print("\n[TEST 4] Verify calculator works with explicit M_bh/d_g")
print("-" * 80)
try:
    calc = EnhancedBuoyancyCalculator()
    params = ComputeParams(
        M=1e30, 
        r=1e10, 
        t=0,
        M_bh=ReferenceSystemLibrary.SGR_A_STAR.M_bh,
        d_g=ReferenceSystemLibrary.SGR_A_STAR.d_g
    )
    
    results = calc.compute_results(params, {'Ug1': 1e-9, 'Ug2': 1e-10, 'Ug3': 1e-11, 'Ug4': 1e-12})
    print(f"✓ Calculation succeeded with {len(results)} results")
    print(f"  Example: {results[0].name} = {results[0].result:.3e} {results[0].unit}")
    
except Exception as e:
    print(f"✗ Calculation failed: {e}")

# Test 5: Verify backward compatibility - old test files need updates
print("\n[TEST 5] Backward compatibility check")
print("-" * 80)
print("⚠️  WARNING: test_phase3.py imports CONSTANTS['M_bh_SgrA'] and will need updates")
print("   Update pattern:")
print("     OLD: M_bh_SgrA = CONSTANTS['M_bh_SgrA']")
print("     NEW: from QCalc_validation import ReferenceSystemLibrary")
print("          M_bh_SgrA = ReferenceSystemLibrary.SGR_A_STAR.M_bh")

print("\n" + "=" * 80)
print("Refactoring verification complete!")
print("=" * 80)
