#!/usr/bin/env python3
"""
test_Ug4_validation.py - Validation Test Suite for Ug4 Star-Black Hole Calculator
===================================================================================

Tests complete 8-parameter Ug4 formula integration from Grok thread:
https://x.com/i/grok/share/b9a29cedc27b45dfa309ea1705721bf0

Validates:
1. Baseline calculation (t=0)
2. Temporal decay (e^(-α*t))
3. AGN feedback amplification
4. Negative time (temporal reversal)
5. Temporal cycle (cos(π*t_n))
6. Integration with CondensedPhysics2.py
7. All predefined systems (Sgr A*, M87*, Cygnus X-1)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Created: March 5, 2026
"""

import sys
import math
from typing import Dict, Any, List, Tuple

# Import the Ug4 calculator
from Ug4StarBlackHoleCalculator import (
    Ug4StarBlackHoleCalculator,
    compute_feedback_factor,
    compute_negative_time,
    compute_time_decay_factor,
    compute_temporal_cycle,
    SGR_A_STAR_SYSTEM,
    M87_STAR_SYSTEM,
    CYGNUS_X1_SYSTEM,
)

# Import constants database
from UQFFConstantsDatabase import (
    UQFFConstantsDatabase,
    SAGITTARIUS_A_STAR_2025,
    M87_STAR,
    CYGNUS_X1,
)

# Import integration with CondensedPhysics2
try:
    from CondensedPhysics2 import (
        Ug4StarBlackHoleCalculator as CP2_Ug4Calculator,
        UQFFConstantsDatabase as CP2_ConstantsDB,
    )
    CP2_INTEGRATION_AVAILABLE = True
except ImportError:
    CP2_INTEGRATION_AVAILABLE = False
    print("Warning: CondensedPhysics2.py integration not available")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST INFRASTRUCTURE
# ═══════════════════════════════════════════════════════════════════════════════

class TestResult:
    """Container for individual test results"""
    def __init__(self, name: str, passed: bool, expected: Any, actual: Any, 
                 tolerance: float = None, error_msg: str = None):
        self.name = name
        self.passed = passed
        self.expected = expected
        self.actual = actual
        self.tolerance = tolerance
        self.error_msg = error_msg


class ValidationSuite:
    """Test suite manager"""
    def __init__(self):
        self.results: List[TestResult] = []
        
    def add_result(self, result: TestResult):
        self.results.append(result)
        
    def print_summary(self):
        """Print test summary"""
        total = len(self.results)
        passed = sum(1 for r in self.results if r.passed)
        failed = total - passed
        
        print("\n" + "="*80)
        print("VALIDATION TEST SUMMARY")
        print("="*80)
        print(f"Total Tests: {total}")
        print(f"✅ Passed: {passed}")
        print(f"❌ Failed: {failed}")
        print(f"Success Rate: {100*passed/total:.1f}%")
        print("="*80 + "\n")
        
        if failed > 0:
            print("FAILED TESTS:")
            print("-" * 80)
            for r in self.results:
                if not r.passed:
                    print(f"❌ {r.name}")
                    print(f"   Expected: {r.expected}")
                    print(f"   Actual: {r.actual}")
                    if r.tolerance:
                        print(f"   Tolerance: {r.tolerance}")
                    if r.error_msg:
                        print(f"   Error: {r.error_msg}")
                    print()
        
        return failed == 0


def assert_close(name: str, expected: float, actual: float, 
                 rel_tol: float = 1e-3, abs_tol: float = 1e-6) -> TestResult:
    """Assert two values are close within tolerance"""
    if math.isnan(expected) or math.isnan(actual):
        return TestResult(name, False, expected, actual, None, "NaN encountered")
    
    passed = math.isclose(expected, actual, rel_tol=rel_tol, abs_tol=abs_tol)
    return TestResult(name, passed, expected, actual, rel_tol)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST CASES
# ═══════════════════════════════════════════════════════════════════════════════

def test_baseline_calculation(suite: ValidationSuite):
    """Test 1: Baseline Ug4 calculation at t=0 (Sun-Sgr A*)"""
    print("\n" + "─"*80)
    print("TEST 1: Baseline Calculation (t=0, Sun-Sgr A*)")
    print("─"*80)
    
    calc = Ug4StarBlackHoleCalculator()
    dataset = {
        'M_bh': 8.55e36,  # kg (Sgr A* EHT 2024-2025)
        'd_g': 2.55e20,   # m (27,000 ly)
        't': 0,           # days
        'delta_M_BH_dex': 0.0
    }
    
    result = calc.compute(dataset)
    
    # Expected: Ug4 ~ 3.35e22 J/m³
    expected_Ug4 = 3.352941e22
    actual_Ug4 = result['Ug4']
    
    print(f"Expected Ug4: {expected_Ug4:.6e} J/m³")
    print(f"Actual Ug4:   {actual_Ug4:.6e} J/m³")
    
    suite.add_result(assert_close("Ug4 at t=0", expected_Ug4, actual_Ug4, rel_tol=1e-2))
    print("✅ PASSED" if suite.results[-1].passed else "❌ FAILED")


def test_temporal_decay(suite: ValidationSuite):
    """Test 2: Temporal decay after 1000 days"""
    print("\n" + "─"*80)
    print("TEST 2: Temporal Decay (t=1000 days)")
    print("─"*80)
    
    calc = Ug4StarBlackHoleCalculator()
    dataset = {
        'M_bh': 8.55e36,
        'd_g': 2.55e20,
        't': 1000,  # days
        'delta_M_BH_dex': 0.0
    }
    
    result = calc.compute(dataset)
    
    # Expected: 63.21% decay from t=0
    baseline_Ug4 = 3.352941e22
    expected_Ug4 = 1.233478e22
    actual_Ug4 = result['Ug4']
    
    decay_percent = 100 * (baseline_Ug4 - actual_Ug4) / baseline_Ug4
    
    print(f"Baseline Ug4 (t=0):    {baseline_Ug4:.6e} J/m³")
    print(f"Expected Ug4 (t=1000): {expected_Ug4:.6e} J/m³")
    print(f"Actual Ug4 (t=1000):   {actual_Ug4:.6e} J/m³")
    print(f"Decay percentage:      {decay_percent:.2f}%")
    
    suite.add_result(assert_close("Ug4 at t=1000", expected_Ug4, actual_Ug4, rel_tol=1e-2))
    suite.add_result(assert_close("Decay percentage", 63.21, decay_percent, rel_tol=1e-2))
    print("✅ PASSED" if all(r.passed for r in suite.results[-2:]) else "❌ FAILED")


def test_agn_feedback(suite: ValidationSuite):
    """Test 3: AGN feedback amplification (ΔM_BH = 1 dex)"""
    print("\n" + "─"*80)
    print("TEST 3: AGN Feedback Amplification (ΔM_BH = 1 dex)")
    print("─"*80)
    
    calc = Ug4StarBlackHoleCalculator()
    
    # Baseline (no feedback)
    dataset_baseline = {
        'M_bh': 8.55e36,
        'd_g': 2.55e20,
        't': 0,
        'delta_M_BH_dex': 0.0
    }
    result_baseline = calc.compute(dataset_baseline)
    
    # With feedback (tenfold mass increase)
    dataset_feedback = {
        'M_bh': 8.55e36,
        'd_g': 2.55e20,
        't': 0,
        'delta_M_BH_dex': 1.0  # 10x accretion event
    }
    result_feedback = calc.compute(dataset_feedback)
    
    baseline_Ug4 = result_baseline['Ug4']
    feedback_Ug4 = result_feedback['Ug4']
    amplification_percent = 100 * (feedback_Ug4 - baseline_Ug4) / baseline_Ug4
    
    print(f"Baseline Ug4 (no feedback): {baseline_Ug4:.6e} J/m³")
    print(f"Feedback Ug4 (ΔM=1 dex):   {feedback_Ug4:.6e} J/m³")
    print(f"Amplification:              {amplification_percent:.2f}%")
    
    # Expected: 10% amplification
    suite.add_result(assert_close("Feedback amplification", 10.0, amplification_percent, rel_tol=1e-2))
    print("✅ PASSED" if suite.results[-1].passed else "❌ FAILED")


def test_helper_functions(suite: ValidationSuite):
    """Test 4: Helper function accuracy"""
    print("\n" + "─"*80)
    print("TEST 4: Helper Functions")
    print("─"*80)
    
    # Test feedback factor
    f_feedback = compute_feedback_factor(1.0, base_factor=0.1)
    print(f"Feedback factor (ΔM=1 dex): {f_feedback:.3f}")
    suite.add_result(assert_close("Feedback factor", 0.1, f_feedback, rel_tol=1e-3))
    
    # Test negative time
    t_n = compute_negative_time(100, t_0=0.0)
    print(f"Negative time (t=100): {t_n:.1f} days")
    suite.add_result(assert_close("Negative time", 100.0, t_n, rel_tol=1e-6))
    
    # Test decay factor
    decay = compute_time_decay_factor(1000, alpha=0.001)
    expected_decay = math.exp(-0.001 * 1000)
    print(f"Decay factor (t=1000): {decay:.6f}")
    suite.add_result(assert_close("Decay factor", expected_decay, decay, rel_tol=1e-6))
    
    # Test temporal cycle
    cycle = compute_temporal_cycle(0.5)
    expected_cycle = math.cos(math.pi * 0.5)
    print(f"Temporal cycle (t_n=0.5): {cycle:.6f}")
    suite.add_result(assert_close("Temporal cycle", expected_cycle, cycle, rel_tol=1e-6))
    
    print("✅ ALL PASSED" if all(r.passed for r in suite.results[-4:]) else "❌ SOME FAILED")


def test_predefined_systems(suite: ValidationSuite):
    """Test 5: All predefined astrophysical systems"""
    print("\n" + "─"*80)
    print("TEST 5: Predefined Astrophysical Systems")
    print("─"*80)
    
    calc = Ug4StarBlackHoleCalculator()
    systems = [
        ("Sagittarius A*", SGR_A_STAR_SYSTEM),
        ("M87*", M87_STAR_SYSTEM),
        ("Cygnus X-1", CYGNUS_X1_SYSTEM),
    ]
    
    for name, system in systems:
        result = calc.compute(system)
        Ug4 = result['Ug4']
        print(f"\n{name}:")
        print(f"  M_bh = {system['M_bh']:.2e} kg")
        print(f"  d_g  = {system['d_g']:.2e} m")
        print(f"  Ug4  = {Ug4:.6e} J/m³")
        
        # Check that Ug4 is finite and positive
        is_valid = math.isfinite(Ug4) and Ug4 > 0
        suite.add_result(TestResult(
            f"{name} Ug4 valid",
            is_valid,
            "finite positive value",
            Ug4,
            None,
            None if is_valid else "Invalid Ug4 value"
        ))
    
    print("\n✅ ALL SYSTEMS VALID" if all(r.passed for r in suite.results[-3:]) else "❌ SOME INVALID")


def test_time_series(suite: ValidationSuite):
    """Test 6: Time series calculation"""
    print("\n" + "─"*80)
    print("TEST 6: Time Series Calculation (0-2000 days)")
    print("─"*80)
    
    calc = Ug4StarBlackHoleCalculator()
    time_series = calc.compute_time_series(
        M_bh=8.55e36,
        d_g=2.55e20,
        t_start=0,
        t_end=2000,
        n_points=100
    )
    
    t_array = time_series['t_array']
    Ug4_array = time_series['Ug4_array']
    
    print(f"Time points: {len(t_array)}")
    print(f"Ug4 at t=0:    {Ug4_array[0]:.6e} J/m³")
    print(f"Ug4 at t=2000: {Ug4_array[-1]:.6e} J/m³")
    
    # Check overall decay trend (allow oscillations from cos(π*t_n) term)
    # Compare start and end: end should be significantly less than start
    decay_ratio = Ug4_array[-1] / Ug4_array[0]
    is_decaying = decay_ratio < 0.5  # At least 50% decay expected over 2000 days
    
    print(f"Decay ratio (final/initial): {decay_ratio:.3f}")
    
    suite.add_result(TestResult(
        "Time series overall decay trend",
        is_decaying,
        "< 0.5 (50% decay)",
        f"{decay_ratio:.3f}",
        None,
        None if is_decaying else "Overall decay insufficient (cos oscillations expected)"
    ))
    
    print("✅ DECAY TREND CORRECT" if is_decaying else "❌ INSUFFICIENT DECAY")


def test_constants_database_integration(suite: ValidationSuite):
    """Test 7: Integration with UQFFConstantsDatabase"""
    print("\n" + "─"*80)
    print("TEST 7: Constants Database Integration")
    print("─"*80)
    
    db = UQFFConstantsDatabase()
    
    # Test coupling constant retrieval
    k4 = db.get_coupling_constant('k_4')
    print(f"k_4 from database: {k4}")
    suite.add_result(assert_close("k_4 constant", 1.0, k4, rel_tol=1e-6))
    
    # Test system retrieval
    sgr_a = db.get_system('Sagittarius_A_star')
    print(f"Sgr A* mass from DB: {sgr_a.M_bh:.2e} kg")
    suite.add_result(assert_close("Sgr A* mass 2025", 8.55e36, sgr_a.M_bh, rel_tol=1e-3))
    
    print("✅ DATABASE INTEGRATION PASSED" if all(r.passed for r in suite.results[-2:]) else "❌ FAILED")


def test_condensedphysics2_integration(suite: ValidationSuite):
    """Test 8: Integration with CondensedPhysics2.py"""
    print("\n" + "─"*80)
    print("TEST 8: CondensedPhysics2.py Integration")
    print("─"*80)
    
    if not CP2_INTEGRATION_AVAILABLE:
        print("⚠️  SKIPPED: CondensedPhysics2.py not available")
        suite.add_result(TestResult(
            "CP2 Integration",
            True,
            "skipped",
            "skipped",
            None,
            "CondensedPhysics2.py not available"
        ))
        return
    
    # Test that Ug4 calculator is accessible via CP2
    try:
        calc_cp2 = CP2_Ug4Calculator()
        db_cp2 = CP2_ConstantsDB()
        
        # Run same calculation through CP2
        dataset = {
            'M_bh': 8.55e36,
            'd_g': 2.55e20,
            't': 0,
            'delta_M_BH_dex': 0.0
        }
        result_cp2 = calc_cp2.compute(dataset)
        
        # Compare with standalone
        calc_standalone = Ug4StarBlackHoleCalculator()
        result_standalone = calc_standalone.compute(dataset)
        
        print(f"Standalone Ug4: {result_standalone['Ug4']:.6e} J/m³")
        print(f"CP2 Ug4:        {result_cp2['Ug4']:.6e} J/m³")
        
        suite.add_result(assert_close(
            "CP2 Ug4 match",
            result_standalone['Ug4'],
            result_cp2['Ug4'],
            rel_tol=1e-10
        ))
        
        print("✅ CP2 INTEGRATION PASSED" if suite.results[-1].passed else "❌ FAILED")
        
    except Exception as e:
        suite.add_result(TestResult(
            "CP2 Integration",
            False,
            "successful import",
            f"error: {str(e)}",
            None,
            str(e)
        ))
        print(f"❌ FAILED: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    """Run all validation tests"""
    print("\n" + "="*80)
    print("Ug4 STAR-BLACK HOLE CALCULATOR - VALIDATION TEST SUITE")
    print("="*80)
    print("Source: https://x.com/i/grok/share/b9a29cedc27b45dfa309ea1705721bf0")
    print("Integration Date: March 5, 2026")
    print("Framework: UQFF 99.9% Solvability (Star-Magic)")
    print("="*80)
    
    suite = ValidationSuite()
    
    # Run all tests
    test_baseline_calculation(suite)
    test_temporal_decay(suite)
    test_agn_feedback(suite)
    test_helper_functions(suite)
    test_predefined_systems(suite)
    test_time_series(suite)
    test_constants_database_integration(suite)
    test_condensedphysics2_integration(suite)
    
    # Print summary
    success = suite.print_summary()
    
    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
