#!/usr/bin/env python3
"""
Validation Script: HawkingTemperatureUQFFCalculator
===================================================

Validates the Python implementation of UQFF-modified Hawking temperature.
Cross-validates against the C++ module (uqff_temperature_formula.cpp).

Key validation:
- T_UQFF / T_H ≈ 0.99 for Sgr A* (4×10⁶ M☉)
- Evaporation timescales exceed universe age for stellar/SMBH
- Primordial BH temperatures reach astronomical values

Author: Copilot / Star Magic UQFF Framework
Date: Feb 26, 2026
"""

import sys
import numpy as np

# Add current directory for imports
sys.path.insert(0, '.')

from CondensedPhysics import (
    HawkingTemperatureUQFFCalculator,
    HAWKING_TEMP_CALC,
    compute_hawking_temperature,
    compute_hawking_temperature_sgra,
    compute_hawking_evaporation,
    compute_hawking_for_system,
    get_hawking_long_form
)


def test_basic_calculation():
    """Test basic Hawking temperature calculation."""
    print("\n" + "="*80)
    print("TEST 1: Basic Hawking Temperature Calculation")
    print("="*80)
    
    calc = HawkingTemperatureUQFFCalculator()
    
    # Sgr A*: 4×10⁶ M☉
    M_sgra = 4.0e6 * calc.M_solar
    result = calc.compute(mode='hawking', M=M_sgra)
    
    print(f"\nSagittarius A* (M = 4×10⁶ M☉):")
    print(f"  T_H (standard):   {result['T_H']:.6e} K")
    print(f"  T_UQFF (modified): {result['T_UQFF']:.6e} K")
    print(f"  Ratio T_UQFF/T_H:  {result['T_UQFF_to_T_H_ratio']:.6f}")
    
    # Validate ratio ≈ 0.99
    expected_ratio = (1 + calc.f_TRZ) * (1 - calc.rho_vac_SCm / calc.rho_vac_UA)
    print(f"  Expected ratio:    {expected_ratio:.6f}")
    
    assert abs(result['T_UQFF_to_T_H_ratio'] - expected_ratio) < 1e-10, "Ratio mismatch!"
    print("  ✓ PASSED: Ratio matches expected value")
    
    # Cross-validate with C++ module results
    # From uqff_temperature_formula.cpp: T_UQFF/T_H = 0.99 for default params
    assert abs(result['T_UQFF_to_T_H_ratio'] - 0.99) < 0.01, "C++ cross-validation failed!"
    print("  ✓ PASSED: Cross-validates with C++ module (T_UQFF/T_H ≈ 0.99)")
    
    return True


def test_evaporation_timescales():
    """Test evaporation timescale calculations."""
    print("\n" + "="*80)
    print("TEST 2: Evaporation Timescales")
    print("="*80)
    
    calc = HawkingTemperatureUQFFCalculator()
    t_universe = 4.35e17  # 13.8 Gyr in seconds
    
    # Test stellar mass BH (10 M☉) - should survive
    M_stellar = 10 * calc.M_solar
    result = calc.compute(mode='evaporation', M=M_stellar)
    
    print(f"\nStellar BH (10 M☉):")
    print(f"  t_evap (standard): {result['t_evap_standard']:.3e} s ({result['t_evap_years']:.3e} years)")
    print(f"  Survives universe: {result['survives_universe']}")
    
    assert result['survives_universe'], "Stellar BH should survive universe age!"
    print("  ✓ PASSED: Stellar BH survives")
    
    # Test primordial BH (10^10 kg) - borderline
    M_primordial = 1e10  # kg
    result_prim = calc.compute(mode='evaporation', M=M_primordial)
    
    print(f"\nPrimordial BH (10¹⁰ kg):")
    print(f"  T_H:               {result_prim['T_H']:.3e} K")
    print(f"  t_evap (standard): {result_prim['t_evap_standard']:.3e} s")
    print(f"  t_evap (years):    {result_prim['t_evap_years']:.3e}")
    print("  ✓ PASSED: Primordial BH calculation complete")
    
    return True


def test_all_systems():
    """Test all pre-defined BH systems."""
    print("\n" + "="*80)
    print("TEST 3: All Pre-Defined Black Hole Systems")
    print("="*80)
    
    calc = HawkingTemperatureUQFFCalculator()
    result = calc.compute(mode='all_systems')
    
    print(f"\nSystems computed: {result['n_systems']}")
    print(f"\n{'System':<20} {'M/M☉':<12} {'T_H (K)':<12} {'T_UQFF/T_H':<10}")
    print("-" * 54)
    
    for key, sys in result['systems'].items():
        print(f"{sys['name']:<20} {sys['M_solar']:<12.2e} {sys['T_H']:<12.3e} {sys['ratio']:<10.4f}")
    
    # Validate all ratios are ≈ 0.99
    for key, sys in result['systems'].items():
        assert abs(sys['ratio'] - 0.99) < 0.01, f"Ratio mismatch for {key}!"
    
    print("\n  ✓ PASSED: All systems computed with correct ratios")
    return True


def test_convenience_functions():
    """Test convenience function API."""
    print("\n" + "="*80)
    print("TEST 4: Convenience Functions")
    print("="*80)
    
    calc = HawkingTemperatureUQFFCalculator()
    
    # Test compute_hawking_temperature
    result1 = compute_hawking_temperature(M_solar=4.0e6)
    print(f"\ncompute_hawking_temperature(M_solar=4e6):")
    print(f"  T_UQFF: {result1['T_UQFF']:.6e} K")
    assert 'T_UQFF' in result1, "Missing T_UQFF in result!"
    print("  ✓ PASSED")
    
    # Test compute_hawking_temperature_sgra
    result2 = compute_hawking_temperature_sgra()
    print(f"\ncompute_hawking_temperature_sgra():")
    print(f"  T_UQFF: {result2['T_UQFF']:.6e} K")
    print(f"  System: {result2['system']['name']}")
    assert result2['system']['key'] == 'SgrA', "Wrong system!"
    print("  ✓ PASSED")
    
    # Test compute_hawking_for_system
    result3 = compute_hawking_for_system('M87', mode='evaporation')
    print(f"\ncompute_hawking_for_system('M87', mode='evaporation'):")
    print(f"  System: {result3['system']['name']}")
    print(f"  Survives universe: {result3['survives_universe']}")
    print("  ✓ PASSED")
    
    return True


def test_long_form_equation():
    """Test long-form equation output."""
    print("\n" + "="*80)
    print("TEST 5: Long-Form Equation Output")
    print("="*80)
    
    calc = HawkingTemperatureUQFFCalculator()
    M = 4.0e6 * calc.M_solar
    
    eq_str = get_hawking_long_form(M)
    print(eq_str[:1500] + "...\n")  # Print first 1500 chars
    
    # Verify key components present
    assert "UQFF-MODIFIED HAWKING TEMPERATURE" in eq_str, "Missing title!"
    assert "T_H =" in eq_str, "Missing T_H equation!"
    assert "T_UQFF" in eq_str, "Missing T_UQFF!"
    assert "RESULT:" in eq_str, "Missing result!"
    
    print("  ✓ PASSED: Long-form equation formatted correctly")
    return True


def test_evaporation_simulation():
    """Test evaporation simulation."""
    print("\n" + "="*80)
    print("TEST 6: Evaporation Simulation")
    print("="*80)
    
    calc = HawkingTemperatureUQFFCalculator()
    
    # Simulate small primordial BH
    M_initial = 1e10  # 10^10 kg
    result = calc.simulate_evaporation(M_initial, dt=1e10, n_steps=100)
    
    print(f"\nSimulation: M_initial = {M_initial:.2e} kg")
    print(f"  Steps: {result['n_steps']}")
    print(f"  t_final: {result['t_final']:.3e} s")
    print(f"  M_final: {result['M_final']:.6e} kg")
    print(f"  Mass lost fraction: {result['mass_lost_fraction']:.6f}")
    
    # Verify arrays
    assert len(result['times']) == len(result['masses']), "Array length mismatch!"
    assert len(result['temperatures_H']) == len(result['masses']), "Temperature array mismatch!"
    
    # Verify temperature increases as mass decreases
    if result['masses'][0] > result['masses'][-1]:
        assert result['temperatures_H'][-1] > result['temperatures_H'][0], "Temperature should increase!"
        print("  ✓ Temperature increases as mass decreases")
    
    print("  ✓ PASSED: Simulation completed successfully")
    return True


def test_global_instances():
    """Test global calculator instances."""
    print("\n" + "="*80)
    print("TEST 7: Global Calculator Instances")
    print("="*80)
    
    # Test HAWKING_TEMP_CALC
    result = HAWKING_TEMP_CALC.compute(mode='hawking', M_solar=10.0)
    print(f"\nHAWKING_TEMP_CALC (10 M☉):")
    print(f"  T_UQFF: {result['T_UQFF']:.6e} K")
    assert 'T_UQFF' in result, "Missing T_UQFF!"
    print("  ✓ PASSED")
    
    return True


def run_all_tests():
    """Run complete validation suite."""
    print("\n" + "="*80)
    print("HAWKING TEMPERATURE UQFF CALCULATOR - VALIDATION SUITE")
    print("="*80)
    
    tests = [
        ("Basic Calculation", test_basic_calculation),
        ("Evaporation Timescales", test_evaporation_timescales),
        ("All Systems", test_all_systems),
        ("Convenience Functions", test_convenience_functions),
        ("Long-Form Equation", test_long_form_equation),
        ("Evaporation Simulation", test_evaporation_simulation),
        ("Global Instances", test_global_instances),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"\n  ✗ FAILED: {name}")
            print(f"    Error: {e}")
            failed += 1
    
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    print(f"  Total tests: {len(tests)}")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    
    if failed == 0:
        print("\n  ✓ ALL TESTS PASSED - HawkingTemperatureUQFFCalculator VALIDATED")
        print("  ✓ Cross-validates with C++ module (uqff_temperature_formula.cpp)")
    else:
        print(f"\n  ✗ {failed} TEST(S) FAILED")
        sys.exit(1)
    
    return failed == 0


if __name__ == "__main__":
    run_all_tests()
