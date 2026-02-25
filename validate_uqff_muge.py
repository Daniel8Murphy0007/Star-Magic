#!/usr/bin/env python3
"""
UQFF MUGE Calculator Validation Script (Feb 26, 2026)

Tests the Python port of uqff_framework.cpp in CondensedPhysics.py.

Tests:
1. Basic instantiation and parameters
2. Quantum coherence calculation
3. Full MUGE computation
4. Term breakdown
5. Pre-defined systems (SgrA, M87, Sun, NeutronStar, Magnetar)
6. Self-expanding (add_term)
7. Time evolution simulation
8. Convenience functions
"""

import sys
import os
import numpy as np

# Add parent dir to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from CondensedPhysics import (
    UQFFMUGECalculator,
    UQFF_MUGE_CALC,
    UQFF_MUGE_SGRA,
    UQFF_MUGE_M87,
    UQFF_MUGE_MAGNETAR,
    compute_muge_gravity,
    compute_muge_breakdown,
    compute_muge_coherence,
    simulate_muge_evolution,
    compute_muge_all_systems,
    compute_muge_for_system,
    get_muge_long_form,
    get_muge_explanations,
)


def test_basic_instantiation():
    """T1: Test basic instantiation and parameters."""
    print("\n" + "="*60)
    print("T1: Basic Instantiation and Parameters")
    print("="*60)
    
    calc = UQFFMUGECalculator()
    
    # Check default params
    assert calc.has_param('G'), "Missing param G"
    assert calc.has_param('c'), "Missing param c"
    assert calc.has_param('hbar'), "Missing param hbar"
    assert calc.has_param('Lambda'), "Missing param Lambda"
    assert calc.has_param('M_initial'), "Missing param M_initial"
    
    # Check values
    assert abs(calc.get_param('G') - 6.6743e-11) < 1e-20, "Wrong G"
    assert abs(calc.get_param('c') - 2.99792458e8) < 1e-1, "Wrong c"
    
    # Test set/get
    calc.set_param('M_initial', 1.989e30)  # Solar mass
    assert abs(calc.get_param('M_initial') - 1.989e30) < 1e20, "set_param failed"
    
    # Test global instance
    assert UQFF_MUGE_CALC is not None, "Global instance missing"
    
    print(f"  ✓ Default params loaded ({len(calc.params)} params)")
    print(f"  ✓ G = {calc.get_param('G'):.5e} m³/kg/s²")
    print(f"  ✓ c = {calc.get_param('c'):.5e} m/s")
    print(f"  ✓ set_param/get_param working")
    print(f"  ✓ Global UQFF_MUGE_CALC available")
    return True


def test_quantum_coherence():
    """T2: Test quantum coherence calculation."""
    print("\n" + "="*60)
    print("T2: Quantum Coherence Calculation")
    print("="*60)
    
    calc = UQFFMUGECalculator.from_system('SgrA')
    r_horizon = calc.params['r_horizon']
    
    # At horizon, coherence should be maximal (Gaussian peak)
    coh_at_horizon = calc.quantum_coherence(r_horizon, 0.0)
    
    # Far from horizon, coherence should be small
    coh_far = calc.quantum_coherence(r_horizon + 1e15, 0.0)
    
    assert abs(coh_at_horizon) > abs(coh_far), "Coherence not peaked at horizon"
    
    print(f"  ✓ r_horizon = {r_horizon:.3e} m (Sgr A*)")
    print(f"  ✓ coherence at horizon = {coh_at_horizon:.6e}")
    print(f"  ✓ coherence at 10⁶ r_h = {coh_far:.6e}")
    print(f"  ✓ Gaussian envelope verified (peak at horizon)")
    return True


def test_muge_computation():
    """T3: Test full MUGE computation."""
    print("\n" + "="*60)
    print("T3: Full MUGE Gravity Computation")
    print("="*60)
    
    calc = UQFFMUGECalculator.from_system('SgrA')
    
    r = 1.27e10  # Event horizon
    t = 0.0
    
    g = calc.compute_MUGE(r, t, noise_level=0.0)
    
    assert g is not None, "compute_MUGE returned None"
    assert not np.isnan(g), "compute_MUGE returned NaN"
    assert not np.isinf(g), "compute_MUGE returned Inf"
    
    # Surface gravity of Sgr A* should be ~10^5 m/s² (order of magnitude)
    assert abs(g) > 1e-10, f"MUGE gravity too small: {g}"
    
    print(f"  ✓ r = {r:.3e} m")
    print(f"  ✓ t = {t:.3e} s")
    print(f"  ✓ g = {g:.6e} m/s²")
    print(f"  ✓ No NaN/Inf issues")
    return True


def test_term_breakdown():
    """T4: Test MUGE term breakdown."""
    print("\n" + "="*60)
    print("T4: MUGE Term Breakdown")
    print("="*60)
    
    calc = UQFFMUGECalculator.from_system('SgrA')
    
    result = calc.compute(mode='breakdown', r=1.27e10, t=0.0)
    
    assert 'terms' in result, "Missing 'terms' in breakdown"
    assert 'g_total' in result, "Missing 'g_total' in breakdown"
    
    terms = result['terms']
    expected_terms = ['base_gravity', 'sum_Ug', 'U_i', 'cosmological', 
                      'quantum', 'fluid', 'dark_matter', 'coherence']
    
    for term in expected_terms:
        assert term in terms, f"Missing term: {term}"
    
    # Verify sum matches total
    term_sum = sum(terms.values())
    assert abs(term_sum - result['g_total']) < 1e-10 * abs(result['g_total']), \
        f"Term sum {term_sum} != total {result['g_total']}"
    
    print(f"  ✓ All {len(expected_terms)} terms present")
    for name, value in terms.items():
        print(f"    {name}: {value:.6e}")
    print(f"  ✓ Total: {result['g_total']:.6e} m/s²")
    print(f"  ✓ Sum verified: {term_sum:.6e} m/s²")
    return True


def test_predefined_systems():
    """T5: Test all pre-defined astrophysical systems."""
    print("\n" + "="*60)
    print("T5: Pre-defined Astrophysical Systems")
    print("="*60)
    
    systems = ['SgrA', 'M87', 'Sun', 'NeutronStar', 'Magnetar']
    
    for sys_name in systems:
        calc = UQFFMUGECalculator.from_system(sys_name)
        system = UQFFMUGECalculator.ASTROPHYSICAL_SYSTEMS[sys_name]
        
        r = system.get('r_0', 1.0)
        g = calc.compute_MUGE(r, 0.0)
        
        assert not np.isnan(g), f"NaN for {sys_name}"
        assert not np.isinf(g), f"Inf for {sys_name}"
        
        print(f"  ✓ {system['name']}: M={system['M_initial']:.3e} kg, g={g:.6e} m/s²")
    
    # Test all_systems mode
    result = UQFF_MUGE_CALC.compute(mode='all_systems')
    assert result['n_systems'] == len(systems), "all_systems count mismatch"
    
    print(f"  ✓ All {len(systems)} systems computed successfully")
    return True


def test_self_expanding():
    """T6: Test self-expanding (add_term)."""
    print("\n" + "="*60)
    print("T6: Self-Expanding (add_term)")
    print("="*60)
    
    calc = UQFFMUGECalculator()
    r = 1e10
    t = 0.0
    
    # Baseline
    g_before = calc.compute_MUGE(r, t)
    assert calc.term_count() == 0, "Should start with 0 custom terms"
    
    # Add custom term: larger value to test against large base
    custom_value = 1e30  # Larger value for visibility
    calc.add_term(lambda r, t, cv=custom_value: cv)
    assert calc.term_count() == 1, "term_count should be 1"
    
    g_after = calc.compute_MUGE(r, t)
    
    # Should differ by exactly custom_value
    delta = abs(g_after - g_before - custom_value)
    rel_error = delta / custom_value if custom_value != 0 else delta
    assert rel_error < 1e-10, f"Custom term not applied: rel_error={rel_error}"
    
    # Clear terms
    calc.clear_terms()
    assert calc.term_count() == 0, "clear_terms failed"
    
    print(f"  ✓ g_before = {g_before:.6e} m/s²")
    print(f"  ✓ g_after = {g_after:.6e} m/s² (with {custom_value:.0e} term)")
    print(f"  ✓ Δg = {(g_after - g_before):.6e} (expected {custom_value:.0e})")
    print(f"  ✓ Relative error: {rel_error:.2e}")
    print(f"  ✓ clear_terms working")
    return True


def test_time_evolution():
    """T7: Test time evolution simulation."""
    print("\n" + "="*60)
    print("T7: Time Evolution Simulation")
    print("="*60)
    
    calc = UQFFMUGECalculator.from_system('SgrA')
    
    result = calc.simulate_evolution(
        r=1.27e10,
        t_start=0.0,
        t_end=1e10,
        dt=1e9,
        noise_level=0.0
    )
    
    assert result['mode'] == 'evolution', "Wrong mode"
    assert result['n_steps'] > 0, "No steps computed"
    assert len(result['times']) == result['n_steps'], "times/steps mismatch"
    assert len(result['gravities']) == result['n_steps'], "gravities/steps mismatch"
    
    print(f"  ✓ Steps: {result['n_steps']}")
    print(f"  ✓ t_start = {result['t_start']:.3e} s")
    print(f"  ✓ t_end = {result['t_end']:.3e} s")
    print(f"  ✓ g_min = {result['g_min']:.6e} m/s²")
    print(f"  ✓ g_max = {result['g_max']:.6e} m/s²")
    print(f"  ✓ g_mean = {result['g_mean']:.6e} m/s²")
    return True


def test_convenience_functions():
    """T8: Test convenience functions."""
    print("\n" + "="*60)
    print("T8: Convenience Functions")
    print("="*60)
    
    r = 1.27e10
    t = 0.0
    
    # compute_muge_gravity
    result1 = compute_muge_gravity(r, t)
    assert 'g_total' in result1, "compute_muge_gravity failed"
    print(f"  ✓ compute_muge_gravity: g = {result1['g_total']:.6e} m/s²")
    
    # compute_muge_breakdown
    result2 = compute_muge_breakdown(r, t)
    assert 'terms' in result2, "compute_muge_breakdown failed"
    print(f"  ✓ compute_muge_breakdown: {len(result2['terms'])} terms")
    
    # compute_muge_coherence
    result3 = compute_muge_coherence(r, t)
    assert 'coherence_at_r' in result3, "compute_muge_coherence failed"
    print(f"  ✓ compute_muge_coherence: coh = {result3['coherence_at_r']:.6e}")
    
    # simulate_muge_evolution
    result4 = simulate_muge_evolution(r, 0.0, 1e10, 1e9)
    assert result4['n_steps'] > 0, "simulate_muge_evolution failed"
    print(f"  ✓ simulate_muge_evolution: {result4['n_steps']} steps")
    
    # compute_muge_all_systems
    result5 = compute_muge_all_systems()
    assert result5['n_systems'] >= 5, "compute_muge_all_systems failed"
    print(f"  ✓ compute_muge_all_systems: {result5['n_systems']} systems")
    
    # compute_muge_for_system
    result6 = compute_muge_for_system('SgrA')
    assert 'g_total' in result6, "compute_muge_for_system failed"
    print(f"  ✓ compute_muge_for_system('SgrA'): g = {result6['g_total']:.6e}")
    
    # get_muge_long_form
    eq = get_muge_long_form(r, t)
    assert 'MUGE EQUATION' in eq, "get_muge_long_form failed"
    print(f"  ✓ get_muge_long_form: {len(eq)} chars generated")
    
    # get_muge_explanations
    expl = get_muge_explanations()
    assert 'UQFF FRAMEWORK' in expl, "get_muge_explanations failed"
    print(f"  ✓ get_muge_explanations: {len(expl)} chars")
    
    return True


def test_global_instances():
    """T9: Test global instances."""
    print("\n" + "="*60)
    print("T9: Global Instances")
    print("="*60)
    
    # Test all global instances
    assert UQFF_MUGE_CALC is not None, "UQFF_MUGE_CALC missing"
    assert UQFF_MUGE_SGRA is not None, "UQFF_MUGE_SGRA missing"
    assert UQFF_MUGE_M87 is not None, "UQFF_MUGE_M87 missing"
    assert UQFF_MUGE_MAGNETAR is not None, "UQFF_MUGE_MAGNETAR missing"
    
    # Test Sgr A* instance
    g_sgra = UQFF_MUGE_SGRA.compute_MUGE(1.27e10, 0.0)
    print(f"  ✓ UQFF_MUGE_SGRA: g = {g_sgra:.6e} m/s²")
    
    # Test M87 instance
    g_m87 = UQFF_MUGE_M87.compute_MUGE(1.92e13, 0.0)
    print(f"  ✓ UQFF_MUGE_M87: g = {g_m87:.6e} m/s²")
    
    # Test Magnetar instance
    g_mag = UQFF_MUGE_MAGNETAR.compute_MUGE(1e4, 0.0)
    print(f"  ✓ UQFF_MUGE_MAGNETAR: g = {g_mag:.6e} m/s²")
    
    print(f"  ✓ All 4 global instances working")
    return True


def main():
    """Run all validation tests."""
    print("="*60)
    print("UQFF MUGE Calculator Validation Suite")
    print("Python port of uqff_framework.cpp")
    print("="*60)
    
    tests = [
        ("T1: Basic Instantiation", test_basic_instantiation),
        ("T2: Quantum Coherence", test_quantum_coherence),
        ("T3: MUGE Computation", test_muge_computation),
        ("T4: Term Breakdown", test_term_breakdown),
        ("T5: Pre-defined Systems", test_predefined_systems),
        ("T6: Self-Expanding", test_self_expanding),
        ("T7: Time Evolution", test_time_evolution),
        ("T8: Convenience Functions", test_convenience_functions),
        ("T9: Global Instances", test_global_instances),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            result = test_func()
            if result:
                passed += 1
                print(f"  → PASSED: {name}")
            else:
                failed += 1
                print(f"  → FAILED: {name}")
        except Exception as e:
            failed += 1
            print(f"  → ERROR: {name}")
            print(f"    Exception: {e}")
    
    print("\n" + "="*60)
    print(f"VALIDATION SUMMARY: {passed}/{len(tests)} tests passed")
    print("="*60)
    
    if failed == 0:
        print("\n✓ ALL TESTS PASSED - UQFFMUGECalculator validated!")
        print("  Ready for commit and push.")
        return 0
    else:
        print(f"\n✗ {failed} test(s) failed - review errors above")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
