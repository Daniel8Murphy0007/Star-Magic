#!/usr/bin/env python3
"""
Test Fermi LAT 4LAC Blazar Model - E_react Verification
Tests enhanced methods for HEASARC verification
"""
import sys
import numpy as np
sys.path.insert(0, r'c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic')

from CondensedPhysics import FERMI_4LAC_BLAZAR_MODEL

def test_E_react_at_t0():
    """Test E_react = 10^46 at t=0."""
    result = FERMI_4LAC_BLAZAR_MODEL.compute_E_react(0.0)
    
    assert result['E_react'] == 1e46, f"Expected 10^46, got {result['E_react']}"
    assert result['decay_factor'] == 1.0, f"Expected decay=1.0, got {result['decay_factor']}"
    
    print(f"TEST 1: E_react(t=0) = {result['E_react']:.2e} PASSED")
    return True

def test_E_react_decay():
    """Test E_react half-life and 1/e decay."""
    model = FERMI_4LAC_BLAZAR_MODEL
    
    # At half-life, should be 50%
    result_half = model.compute_E_react(model.t_half)
    expected_half = 0.5 * model.E_react_0
    assert abs(result_half['E_react'] - expected_half) / expected_half < 0.01
    
    # At τ, should be ~36.8% (1/e)
    result_tau = model.compute_E_react(model.tau)
    expected_tau = model.E_react_0 / np.e
    assert abs(result_tau['E_react'] - expected_tau) / expected_tau < 0.01
    
    print(f"TEST 2: E_react decay - t_half={model.t_half:.0f}d → 50%, τ={model.tau:.0f}d → 36.8% PASSED")
    return True

def test_luminosity_in_4LAC_range():
    """Test luminosity falls within 4LAC observed range."""
    result = FERMI_4LAC_BLAZAR_MODEL.compute_blazar_luminosity(t_days=0.0)
    
    assert result['in_4LAC_range'], f"L_γ = {result['L_gamma']:.2e} outside 4LAC range"
    assert 39 <= result['L_gamma_log10'] <= 47, f"log10(L) = {result['L_gamma_log10']:.2f} outside range"
    
    print(f"TEST 3: L_γ = 10^{result['L_gamma_log10']:.2f} W in 4LAC range PASSED")
    return True

def test_verify_4LAC_range():
    """Test comprehensive 4LAC range verification."""
    result = FERMI_4LAC_BLAZAR_MODEL.verify_4LAC_range()
    
    assert result['all_verified'], "Not all test cases verified"
    assert len(result['test_cases']) == 5, f"Expected 5 test cases, got {len(result['test_cases'])}"
    
    print(f"TEST 4: 4LAC Range - {len(result['test_cases'])} test cases, all_verified={result['all_verified']} PASSED")
    return True

def test_verify_HEASARC_4LAC():
    """Test HEASARC-specific verification."""
    result = FERMI_4LAC_BLAZAR_MODEL.verify_HEASARC_4LAC()
    
    # Check HEASARC catalog data
    assert result['HEASARC_catalog']['total_AGNs'] == 3407
    assert result['HEASARC_catalog']['blazars_fraction'] == 0.98
    
    # Check all luminosity tests pass
    assert result['all_verified'], "Not all luminosity tests verified"
    
    print(f"TEST 5: HEASARC 4LAC - {result['HEASARC_catalog']['total_AGNs']} AGNs, all_verified={result['all_verified']} PASSED")
    return True

def test_long_form_solution():
    """Test long form solution output."""
    result = FERMI_4LAC_BLAZAR_MODEL.long_form_solution()
    
    # Check key sections
    assert 'FERMI LAT 4LAC' in result, "Missing Fermi LAT 4LAC"
    assert 'E_react' in result, "Missing E_react"
    assert 'HEASARC' in result, "Missing HEASARC reference"
    assert 'STEP 1' in result, "Missing Step 1"
    assert 'STEP 7' in result, "Missing Step 7"
    assert 'CONCLUSION' in result, "Missing conclusion"
    
    print(f"TEST 6: Long Form - {len(result)} chars, all sections present PASSED")
    return True

def test_run_tests():
    """Test the built-in test suite."""
    all_passed, tests, summary = FERMI_4LAC_BLAZAR_MODEL.run_tests()
    
    assert all_passed, f"Built-in tests failed: {summary}"
    assert len(tests) == 5, f"Expected 5 tests, got {len(tests)}"
    
    print(f"TEST 7: Built-in Tests - {summary} PASSED")
    return True

if __name__ == '__main__':
    import numpy as np
    
    print("=" * 70)
    print("FERMI LAT 4LAC BLAZAR MODEL - E_REACT VERIFICATION TESTS")
    print("=" * 70)
    
    tests = [
        test_E_react_at_t0,
        test_E_react_decay,
        test_luminosity_in_4LAC_range,
        test_verify_4LAC_range,
        test_verify_HEASARC_4LAC,
        test_long_form_solution,
        test_run_tests,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"{test.__name__}: FAILED - {e}")
            failed += 1
    
    print("=" * 70)
    print(f"RESULTS: {passed}/{len(tests)} TESTS PASSED")
    if failed == 0:
        print("✓ ALL TESTS PASSED - Document 9 implementation verified")
    else:
        print(f"✗ {failed} tests failed")
    print("=" * 70)
