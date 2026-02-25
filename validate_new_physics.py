#!/usr/bin/env python3
"""
Validation Script for New Physics Calculators (Feb 25, 2026)

Tests:
1. PulsarTimingArrayUQFFCalculator
2. CosmicRayPropagationUQFFCalculator
3. GravitationalLensingUQFFCalculator
4. StringTheoryCompactificationCalculator
5. SelfExpandingMixin (set_learning_rate)
6. ValidationCoverageFramework
"""

import sys
import os

# Add parent dir to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from CondensedPhysics import (
    # New calculators
    PulsarTimingArrayUQFFCalculator,
    CosmicRayPropagationUQFFCalculator,
    GravitationalLensingUQFFCalculator,
    StringTheoryCompactificationCalculator,
    ValidationCoverageFramework,
    SelfExpandingMixin,
    
    # Global instances
    PTA_UQFF_CALC,
    PTA_NANOGRAV,
    CR_PROPAGATION_CALC,
    GRAV_LENSING_CALC,
    STRING_COMPACT_CALC,
    VALIDATION_FRAMEWORK,
    
    # Convenience functions
    compute_pta_spectrum,
    compute_gzk_horizon,
    compute_einstein_radius,
    compute_kk_spectrum,
    compute_uqff_26d_projection,
)


def test_pta_calculator():
    """Test Pulsar Timing Array Calculator."""
    print("\n" + "="*70)
    print("TEST: PulsarTimingArrayUQFFCalculator")
    print("="*70)
    
    # Run built-in tests
    results = PulsarTimingArrayUQFFCalculator.run_tests()
    
    for test in results['tests']:
        status = "✓" if test['passed'] else "✗"
        print(f"  {status} {test['name']}")
    
    # Additional tests
    print("\n  Additional tests:")
    
    # Test convenience function
    spec = compute_pta_spectrum(n_freq=10)
    assert 'frequencies' in spec, "Spectrum computation failed"
    print(f"  ✓ compute_pta_spectrum: {len(spec['frequencies'])} frequency bins")
    
    # Test global instance
    assert PTA_NANOGRAV.n_pulsars == 47, "NANOGrav config wrong"
    print(f"  ✓ PTA_NANOGRAV: {PTA_NANOGRAV.n_pulsars} pulsars, {PTA_NANOGRAV.T_obs / PTA_NANOGRAV.year_to_s:.0f} years")
    
    return results['passed'] == results['total']


def test_cosmic_ray_calculator():
    """Test Cosmic Ray Propagation Calculator."""
    print("\n" + "="*70)
    print("TEST: CosmicRayPropagationUQFFCalculator")
    print("="*70)
    
    results = CosmicRayPropagationUQFFCalculator.run_tests()
    
    for test in results['tests']:
        status = "✓" if test['passed'] else "✗"
        print(f"  {status} {test['name']}")
    
    # Additional tests
    print("\n  Additional tests:")
    
    # GZK horizon
    horizon = compute_gzk_horizon(1e20)
    assert 'lambda_Mpc' in horizon, "GZK horizon failed"
    print(f"  ✓ GZK horizon at 10^20 eV: {horizon['lambda_Mpc']:.1f} Mpc")
    
    # Spectrum features
    features = CR_PROPAGATION_CALC.spectrum_features()
    assert 'GZK_cutoff' in features, "Features missing"
    print(f"  ✓ GZK cutoff: {features['GZK_cutoff']['energy_EeV']:.0f} EeV")
    print(f"  ✓ Ankle: {features['ankle']['energy_EeV']:.0f} EeV")
    print(f"  ✓ Knee: {features['knee']['energy_PeV']:.0f} PeV")
    
    return results['passed'] == results['total']


def test_gravitational_lensing_calculator():
    """Test Gravitational Lensing Calculator."""
    print("\n" + "="*70)
    print("TEST: GravitationalLensingUQFFCalculator")
    print("="*70)
    
    results = GravitationalLensingUQFFCalculator.run_tests()
    
    for test in results['tests']:
        status = "✓" if test['passed'] else "✗"
        print(f"  {status} {test['name']}")
    
    # Additional tests
    print("\n  Additional tests:")
    
    # Einstein radius for galaxy cluster
    er = compute_einstein_radius(M_Msun=1e12, D_l_Mpc=1000, D_s_Mpc=2000)
    assert 'theta_E_arcsec' in er, "Einstein radius failed"
    print(f"  ✓ Einstein radius (10^12 M☉): {er['theta_E_arcsec']:.2f} arcsec")
    
    # Pre-defined systems
    for name in ['SgrA', 'M87', 'Abell2218']:
        sys_result = GRAV_LENSING_CALC.compute(mode='system', name=name)
        print(f"  ✓ {name}: θ_E = {sys_result['theta_E_arcsec']:.4f} arcsec")
    
    return results['passed'] == results['total']


def test_string_theory_calculator():
    """Test String Theory Compactification Calculator."""
    print("\n" + "="*70)
    print("TEST: StringTheoryCompactificationCalculator")
    print("="*70)
    
    results = StringTheoryCompactificationCalculator.run_tests()
    
    for test in results['tests']:
        status = "✓" if test['passed'] else "✗"
        print(f"  {status} {test['name']}")
    
    # Additional tests
    print("\n  Additional tests:")
    
    # KK tower
    kk = compute_kk_spectrum(R_m=1e-19, n_max=5)
    assert len(kk['masses_GeV']) == 5, "KK spectrum wrong"
    print(f"  ✓ KK spectrum (R=10^-19 m): E_1 = {kk['masses_GeV'][0]:.2e} GeV")
    
    # Calabi-Yau quintic
    cy = STRING_COMPACT_CALC.compute(mode='calabi_yau', h11=1, h21=101)
    assert cy['euler'] == -200, "Calabi-Yau Euler wrong"
    print(f"  ✓ Quintic CY: χ = {cy['euler']}, moduli = {cy['total_moduli']}")
    
    # UQFF 26D projection
    proj = compute_uqff_26d_projection(r=1e10)
    assert proj['total_layers'] == 26, "26D layers wrong"
    print(f"  ✓ UQFF 26D: {proj['extended_fraction']*100:.1f}% extended, {proj['compact_fraction']*100:.1f}% compact")
    
    return results['passed'] == results['total']


def test_self_expanding_mixin():
    """Test SelfExpandingMixin methods."""
    print("\n" + "="*70)
    print("TEST: SelfExpandingMixin (set_learning_rate)")
    print("="*70)
    
    class TestClass(SelfExpandingMixin):
        def __init__(self):
            self.init_self_expanding()
    
    instance = TestClass()
    
    # Test set_learning_rate
    instance.set_learning_rate(0.01)
    assert instance.get_learning_rate() == 0.01, "set_learning_rate failed"
    print("  ✓ set_learning_rate(0.01)")
    
    # Test set_enable_logging
    instance.set_enable_logging(True)
    assert instance.enable_logging == True, "set_enable_logging failed"
    print("  ✓ set_enable_logging(True)")
    
    # Test set_enable_dynamic_terms
    instance.set_enable_dynamic_terms(True)
    assert instance.enable_dynamic_terms == True, "set_enable_dynamic_terms failed"
    print("  ✓ set_enable_dynamic_terms(True)")
    
    # Test invalid learning rate
    try:
        instance.set_learning_rate(-0.1)
        print("  ✗ Should reject negative learning rate")
        return False
    except ValueError:
        print("  ✓ Correctly rejected negative learning rate")
    
    return True


def test_validation_framework():
    """Test ValidationCoverageFramework."""
    print("\n" + "="*70)
    print("TEST: ValidationCoverageFramework")
    print("="*70)
    
    # Discover calculators (limited scope for speed)
    test_calculators = [
        PulsarTimingArrayUQFFCalculator,
        CosmicRayPropagationUQFFCalculator,
        GravitationalLensingUQFFCalculator,
        StringTheoryCompactificationCalculator,
    ]
    
    results = VALIDATION_FRAMEWORK.run_all_tests(test_calculators)
    
    print(f"  ✓ Discovered {results['total_calculators']} calculators")
    print(f"  ✓ {results['with_run_tests']} have run_tests()")
    print(f"  ✓ Coverage: {results['coverage_percent']:.0f}%")
    print(f"  ✓ Tests: {results['total_passed']}/{results['total_tests']} passed ({results['pass_rate']:.0f}%)")
    
    # Generate report (just test it works)
    report = VALIDATION_FRAMEWORK.generate_report(results)
    assert 'VALIDATION COVERAGE REPORT' in report, "Report generation failed"
    print(f"  ✓ Report generated ({len(report)} chars)")
    
    return results['pass_rate'] == 100


def main():
    """Run all validation tests."""
    print("="*70)
    print("NEW PHYSICS CALCULATORS VALIDATION SUITE")
    print("Feb 25, 2026")
    print("="*70)
    
    tests = [
        ("PulsarTimingArrayUQFFCalculator", test_pta_calculator),
        ("CosmicRayPropagationUQFFCalculator", test_cosmic_ray_calculator),
        ("GravitationalLensingUQFFCalculator", test_gravitational_lensing_calculator),
        ("StringTheoryCompactificationCalculator", test_string_theory_calculator),
        ("SelfExpandingMixin", test_self_expanding_mixin),
        ("ValidationCoverageFramework", test_validation_framework),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            result = test_func()
            if result:
                passed += 1
                print(f"\n  → PASSED: {name}")
            else:
                failed += 1
                print(f"\n  → FAILED: {name}")
        except Exception as e:
            failed += 1
            print(f"\n  → ERROR: {name}")
            print(f"    Exception: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print(f"VALIDATION SUMMARY: {passed}/{len(tests)} tests passed")
    print("="*70)
    
    if failed == 0:
        print("\n✓ ALL TESTS PASSED - New calculators validated!")
        print("  Ready for commit and push.")
        return 0
    else:
        print(f"\n✗ {failed} test(s) failed - review errors above")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
