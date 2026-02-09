#!/usr/bin/env python3
"""
Test Document 10: ENSDF (NNDC 2025) Pb-206 n=8 Bindings Verification
"""

import sys
import numpy as np
sys.path.insert(0, '.')
from CondensedPhysics import NUCLEAR_SHELL_MODEL

def test_nuclear_binding_shell_levels():
    """Run all tests for Document 10: ENSDF NNDC 2025."""
    print("=" * 70)
    print("ENSDF (NNDC 2025) Pb-206 n=8 BINDINGS - VERIFICATION TESTS")
    print("=" * 70)
    
    tests_passed = 0
    tests_total = 0
    
    # TEST 1: ENSDF NNDC 2025 Pb-206 verification
    tests_total += 1
    try:
        result = NUCLEAR_SHELL_MODEL.verify_ENSDF_NNDC_2025_Pb206()
        if result['UQFF_verified']:
            tests_passed += 1
            print(f"TEST 1: ENSDF NNDC 2025 Pb-206 - n_solved={result['verification']['n_solved']:.0f} PASSED")
        else:
            print(f"TEST 1: ENSDF NNDC 2025 Pb-206 - FAILED")
    except Exception as e:
        print(f"TEST 1: ENSDF NNDC 2025 - ERROR: {e}")
    
    # TEST 2: 29 ENSDF levels loaded
    tests_total += 1
    try:
        result = NUCLEAR_SHELL_MODEL.verify_ENSDF_NNDC_2025_Pb206()
        n_levels = result['ENSDF_data']['n_levels']
        if n_levels == 29:
            tests_passed += 1
            print(f"TEST 2: ENSDF levels count - {n_levels} levels PASSED")
        else:
            print(f"TEST 2: ENSDF levels count - {n_levels} levels (expected 29) FAILED")
    except Exception as e:
        print(f"TEST 2: ENSDF levels - ERROR: {e}")
    
    # TEST 3: n=8 -> 10^{-12} J
    tests_total += 1
    try:
        result = NUCLEAR_SHELL_MODEL.verify_ENSDF_NNDC_2025_Pb206()
        E_8 = result['verification']['E_n_predicted_J']
        expected = 1e-12
        if 0.9e-12 <= E_8 <= 1.1e-12:
            tests_passed += 1
            print(f"TEST 3: E_8 = {E_8:.2e} J (expected ~10^-12) PASSED")
        else:
            print(f"TEST 3: E_8 = {E_8:.2e} J (expected ~10^-12) FAILED")
    except Exception as e:
        print(f"TEST 3: E_8 calculation - ERROR: {e}")
    
    # TEST 4: Total binding B ≈ 1.71 GeV
    tests_total += 1
    try:
        result = NUCLEAR_SHELL_MODEL.verify_ENSDF_NNDC_2025_Pb206()
        B_GeV = result['binding']['B_total_GeV']
        if 1.5 <= B_GeV <= 1.9:
            tests_passed += 1
            print(f"TEST 4: Total binding B = {B_GeV:.2f} GeV PASSED")
        else:
            print(f"TEST 4: Total binding B = {B_GeV:.2f} GeV (expected ~1.71) FAILED")
    except Exception as e:
        print(f"TEST 4: Total binding - ERROR: {e}")
    
    # TEST 5: Polynomial fit R² > 0.95
    tests_total += 1
    try:
        result = NUCLEAR_SHELL_MODEL.verify_ENSDF_NNDC_2025_Pb206()
        R2_8 = result['polynomial_fit']['R_squared_8']
        if R2_8 > 0.95:
            tests_passed += 1
            print(f"TEST 5: Polynomial fit deg=8 R²={R2_8:.4f} PASSED")
        else:
            print(f"TEST 5: Polynomial fit deg=8 R²={R2_8:.4f} (expected >0.95) FAILED")
    except Exception as e:
        print(f"TEST 5: Polynomial fit - ERROR: {e}")
    
    # TEST 6: Maximum level ~10 MeV
    tests_total += 1
    try:
        result = NUCLEAR_SHELL_MODEL.verify_ENSDF_NNDC_2025_Pb206()
        max_MeV = result['verification']['max_level_MeV']
        if 9.5 <= max_MeV <= 10.5:
            tests_passed += 1
            print(f"TEST 6: Max level = {max_MeV:.1f} MeV PASSED")
        else:
            print(f"TEST 6: Max level = {max_MeV:.1f} MeV (expected ~10) FAILED")
    except Exception as e:
        print(f"TEST 6: Max level - ERROR: {e}")
    
    # TEST 7: Long-form ENSDF proof
    tests_total += 1
    try:
        proof = NUCLEAR_SHELL_MODEL.long_form_ENSDF_proof()
        required_sections = ['INTRODUCTION', 'UQFF MODEL', 'MATHEMATICAL', 'EMPIRICAL', 'CONCLUSION']
        has_all = all(s in proof for s in required_sections)
        if len(proof) > 5000 and has_all:
            tests_passed += 1
            print(f"TEST 7: Long-form ENSDF proof - {len(proof)} chars, all sections PASSED")
        else:
            print(f"TEST 7: Long-form ENSDF proof - {len(proof)} chars, has_sections={has_all} FAILED")
    except Exception as e:
        print(f"TEST 7: Long-form proof - ERROR: {e}")
    
    # TEST 8: Built-in run_tests
    tests_total += 1
    try:
        builtin_result = NUCLEAR_SHELL_MODEL.run_tests()
        passed = builtin_result['tests_passed']
        total = builtin_result['tests_total']
        if builtin_result['all_passed']:
            tests_passed += 1
            print(f"TEST 8: Built-in tests - {passed}/{total} passed PASSED")
        else:
            print(f"TEST 8: Built-in tests - {passed}/{total} passed FAILED")
    except Exception as e:
        print(f"TEST 8: Built-in tests - ERROR: {e}")
    
    print("=" * 70)
    print(f"RESULTS: {tests_passed}/{tests_total} TESTS PASSED")
    if tests_passed == tests_total:
        print("✓ ALL TESTS PASSED - Document 10 implementation verified")
    else:
        print(f"✗ {tests_total - tests_passed} TEST(S) FAILED")
    print("=" * 70)
    
    return tests_passed == tests_total

if __name__ == "__main__":
    success = test_nuclear_binding_shell_levels()
    sys.exit(0 if success else 1)
