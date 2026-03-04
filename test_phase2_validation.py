"""
================================================================================
PHASE 2 VALIDATION TEST SUITE - QUANTUM LEVEL 26 + DPM COSMOLOGY
================================================================================

**Module**: test_phase2_validation.py
**Created**: March 4, 2026
**Purpose**: Comprehensive validation for Phase 2 integration

Tests:
1. QuantumLevel26Framework.py (630 lines)
   - 26-level energy density calculations
   - Universal Inertia per level
   - Phase transitions (solid/liquid/gas/plasma)
   - Cross-scale quantum coupling
   
2. DPMCosmologyModule.py (565 lines)
   - 26-center pre-Big Bang configuration
   - Universal inflation force at t=0
   - Universal Nuclear Core Model
   - Belly Button resonance
   - Inflation dynamics

3. CondensedPhysics2.py integration
   - Import verification
   - Export verification
   - Full pipeline test

**Expected Results**: 25/25 tests pass (100%)

---
©2026 Daniel T. Murphy, daniel.murphy00@gmail.com – All Rights Reserved
================================================================================
"""

import sys
import math
from typing import Dict, List, Tuple

# Test infrastructure setup
class TestResult:
    """Store result of single test."""
    def __init__(self, name: str, passed: bool, details: str = ""):
        self.name = name
        self.passed = passed
        self.details = details

class ValidationSuite:
    """Manage test suite execution and reporting."""
    def __init__(self, name: str):
        self.name = name
        self.tests: List[TestResult] = []
    
    def add_test(self, name: str, passed: bool, details: str = ""):
        """Add test result."""
        self.tests.append(TestResult(name, passed, details))
        status = "[PASS]" if passed else "[FAIL]"
        print(f"  {status}: {name}")
        if details and not passed:
            print(f"    Details: {details}")
    
    def summary(self) -> Tuple[int, int]:
        """Return (passed_count, total_count)."""
        passed = sum(1 for t in self.tests if t.passed)
        total = len(self.tests)
        return passed, total

def assert_close(actual: float, expected: float, rel_tol: float = 0.01, name: str = "") -> bool:
    """
    Assert two floats are close within relative tolerance.
    
    Args:
        actual: Actual value
        expected: Expected value
        rel_tol: Relative tolerance (default 1%)
        name: Test name for error reporting
    
    Returns:
        True if close, False otherwise
    """
    if expected == 0:
        return abs(actual) < abs(rel_tol)
    
    rel_error = abs((actual - expected) / expected)
    if rel_error >= rel_tol:
        print(f"      Expected: {expected:.6e}, Got: {actual:.6e}, Rel Error: {rel_error*100:.2f}%")
        return False
    return True


# ============================================================================
# TEST SUITE 1: QUANTUM LEVEL 26 FRAMEWORK
# ============================================================================

def test_quantum_level_26_framework(suite: ValidationSuite):
    """Test QuantumLevel26Framework.py module."""
    
    print("\n" + "="*80)
    print("TEST SUITE 1: QUANTUM LEVEL 26 FRAMEWORK")
    print("="*80)
    
    try:
        from QuantumLevel26Framework import (
            QuantumLevel26Calculator,
            PhaseTransitionCalculator,
            CrossScaleCouplingCalculator,
            QuantumLevel,
            QUANTUM_LEVELS_26
        )
    except ImportError as e:
        suite.add_test("Import QuantumLevel26Framework", False, str(e))
        return
    
    suite.add_test("Import QuantumLevel26Framework", True)
    
    # Test 1.1: Energy density for matter phases (Levels 10-13)
    calc = QuantumLevel26Calculator()
    
    # Expected: E_i = ρ_vac,[SCm] × level² = 1e-8 × level²
    E_10 = calc.compute_energy_density(10)
    expected_E_10 = 1e-8 * (10**2)  # 1e-6 J/m³
    suite.add_test(
        "Level 10 (Solids) Energy Density",
        assert_close(E_10, expected_E_10, name="E_10"),
        f"E_10 = {E_10:.2e} J/m³"
    )
    
    E_13 = calc.compute_energy_density(13)
    expected_E_13 = 1e-8 * (13**2)  # 1.69e-6 J/m³
    suite.add_test(
        "Level 13 (Plasma) Energy Density",
        assert_close(E_13, expected_E_13, name="E_13"),
        f"E_13 = {E_13:.2e} J/m³"
    )
    
    # Test 1.2: All 26 levels
    all_levels = calc.compute_all_levels()
    suite.add_test(
        "All 26 Levels Calculated",
        len(all_levels) == 26 and all(E > 0 for E in all_levels.values()),
        f"{len(all_levels)} levels computed"
    )
    
    # Test 1.3: Universal Inertia at t=0
    Ui_10 = calc.compute_universal_inertia(10, t_n=0.0, f_TRZ=0.01)
    # Expected: non-zero positive value
    suite.add_test(
        "Universal Inertia Level 10",
        Ui_10 > 0 and Ui_10 < 1e15,  # Reasonable range
        f"Ui_10 = {Ui_10:.2e} J/m³"
    )
    
    # Test 1.4: Total field energy
    E_total = calc.compute_total_field_energy(t_n=0.0)
    suite.add_test(
        "Total Field Energy (All 26 Levels)",
        E_total > 0,
        f"E_total = {E_total:.2e} J/m³"
    )
    
    # Test 1.5: Level lookup by scale
    level_nano = calc.get_level_by_scale(1e-9)  # nanometer → Level 10
    suite.add_test(
        "Level Lookup by Scale (nanometer)",
        level_nano == 10,
        f"1e-9 m → Level {level_nano}"
    )
    
    # Test 1.6: Phase transitions
    phase_calc = PhaseTransitionCalculator()
    transitions = phase_calc.compute_matter_phase_transitions()
    
    delta_E_solid_liquid = transitions['solid_to_liquid']
    expected_delta = 1e-8 * (11**2 - 10**2)  # 21e-8 = 2.1e-7 J/m³
    suite.add_test(
        "Solid → Liquid Transition Energy",
        assert_close(delta_E_solid_liquid, expected_delta, name="ΔE_10→11"),
        f"ΔE = {delta_E_solid_liquid:.2e} J/m³"
    )
    
    # Test 1.7: Cross-scale coupling
    coupling_calc = CrossScaleCouplingCalculator(coherence_length=5.0)
    coupling_10_11 = coupling_calc.compute_coupling_strength(10, 11)
    suite.add_test(
        "Adjacent Level Coupling (10↔11)",
        0.5 < coupling_10_11 < 1.0,  # Strong coupling for adjacent levels
        f"ψ_10↔11 = {coupling_10_11:.4f}"
    )
    
    coupling_10_26 = coupling_calc.compute_coupling_strength(10, 26)
    suite.add_test(
        "Distant Level Coupling (10↔26)",
        0.0 < coupling_10_26 < 0.1,  # Weak coupling for distant levels
        f"ψ_10↔26 = {coupling_10_26:.4f}"
    )
    
    # Test 1.8: Quantum level data structure
    level_info = calc.get_level_info(10)
    suite.add_test(
        "Level 10 Info Complete",
        isinstance(level_info, QuantumLevel) and 
        "SOLIDS" in level_info.state_description,
        f"State: {level_info.state_description}"
    )


# ============================================================================
# TEST SUITE 2: DPM COSMOLOGY MODULE
# ============================================================================

def test_dpm_cosmology_module(suite: ValidationSuite):
    """Test DPMCosmologyModule.py module."""
    
    print("\n" + "="*80)
    print("TEST SUITE 2: DPM COSMOLOGY MODULE")
    print("="*80)
    
    try:
        from DPMCosmologyModule import (
            DPMCosmologyCalculator,
            UniversalNuclearCoreModel,
            InflationDynamicsCalculator,
            DPMCenter,
            generate_26_centers
        )
    except ImportError as e:
        suite.add_test("Import DPMCosmologyModule", False, str(e))
        return
    
    suite.add_test("Import DPMCosmologyModule", True)
    
    # Test 2.1: 26-center generation
    centers = generate_26_centers()
    suite.add_test(
        "Generate 26 DPM Centers",
        len(centers) == 26 and all(isinstance(c, DPMCenter) for c in centers.values()),
        f"{len(centers)} centers created"
    )
    
    # Test 2.2: DPM cosmology calculator
    dpm = DPMCosmologyCalculator()
    
    # Test 2.3: Center energy
    E_center_1 = dpm.compute_center_energy(1)
    suite.add_test(
        "DPM Center 1 Energy",
        E_center_1 > 0,
        f"E_1 = {E_center_1:.2e} J"
    )
    
    # Test 2.4: Total pre-inflationary energy
    E_preinflationary = dpm.compute_total_preinflationary_energy()
    suite.add_test(
        "Total Pre-Inflationary Energy",
        E_preinflationary > 0,
        f"E_total = {E_preinflationary:.2e} J"
    )
    
    # Test 2.5: Inflation force at t=0
    F_U = dpm.compute_inflation_force_at_t0()
    suite.add_test(
        "Inflation Force at t=0",
        F_U > 1e10,  # Should be >> F_core ≈ 10^10 N
        f"F_U(t=0) = {F_U:.2e} N"
    )
    
    # Test 2.6: Center separations
    d_1_2 = dpm.compute_center_separation(1, 2)
    d_1_26 = dpm.compute_center_separation(1, 26)
    suite.add_test(
        "Center Separation (adjacent vs distant)",
        d_1_2 < d_1_26,  # Adjacent centers closer than distant
        f"d_1↔2 = {d_1_2:.2e} m < d_1↔26 = {d_1_26:.2e} m"
    )
    
    # Test 2.7: Universal Nuclear Core Model
    nuclear = UniversalNuclearCoreModel()
    
    # Test 2.8: UA-SCm coupling
    g_Fe56 = nuclear.compute_ua_scm_coupling(56)  # Iron-56 (reference)
    expected_g_Fe56 = (1e-8 / 1e-11) * ((56/56)**(1/3))  # = 1000
    suite.add_test(
        "UA-SCm Coupling for Fe-56",
        assert_close(g_Fe56, expected_g_Fe56, rel_tol=0.01),
        f"g_coupling = {g_Fe56:.1f}"
    )
    
    # Test 2.9: Nuclear binding energy
    B_Fe56 = nuclear.compute_nuclear_binding_energy(56, 26)
    suite.add_test(
        "Fe-56 Binding Energy",
        450 < B_Fe56 < 550,  # Experimental ~492 MeV
        f"B(Fe-56) = {B_Fe56:.1f} MeV"
    )
    
    # Test 2.10: Belly Button resonance
    f_bb_0 = nuclear.compute_belly_button_resonance(0)  # t=0
    f_bb_1yr = nuclear.compute_belly_button_resonance(31536000)  # 1 year
    suite.add_test(
        "Belly Button Resonance (time decay)",
        abs(f_bb_0) > abs(f_bb_1yr),  # Decays over time
        f"f_bb(0) = {f_bb_0:.4f} > f_bb(1yr) = {f_bb_1yr:.4f}"
    )
    
    # Test 2.11: Inflation dynamics
    inflation = InflationDynamicsCalculator()
    
    # Test 2.12: Scale factor (careful with overflow)
    try:
        a_0 = inflation.compute_scale_factor(0)
        suite.add_test(
            "Scale Factor at t=0",
            assert_close(a_0, 1.0, rel_tol=0.01),
            f"a(0) = {a_0:.4f}"
        )
    except (OverflowError, ValueError) as e:
        suite.add_test("Scale Factor at t=0", False, f"Numerical overflow: {e}")
    
    # Test 2.13: Center mixing entropy
    S_mixing = inflation.compute_center_mixing_entropy()
    suite.add_test(
        "26-Center Mixing Entropy",
        S_mixing > 0,
        f"S = {S_mixing:.2e} J/K"
    )
    
    # Test 2.14: Quantum level formation times
    t_form_1 = inflation.compute_quantum_level_formation_time(1)
    t_form_26 = inflation.compute_quantum_level_formation_time(26)
    suite.add_test(
        "Level Formation Time Progression",
        t_form_1 < t_form_26,  # Higher levels form later
        f"t_1 = {t_form_1:.2e} s < t_26 = {t_form_26:.2e} s"
    )


# ============================================================================
# TEST SUITE 3: CONDENSEDPHYSICS2 INTEGRATION
# ============================================================================

def test_condensedphysics2_integration(suite: ValidationSuite):
    """Test CondensedPhysics2.py integration."""
    
    print("\n" + "="*80)
    print("TEST SUITE 3: CONDENSEDPHYSICS2 INTEGRATION")
    print("="*80)
    
    try:
        from CondensedPhysics2 import (
            # Phase 2: Quantum Level 26 Framework (with Phase2_ prefix)
            Phase2_QuantumLevel26Calculator,
            Phase2_PhaseTransitionCalculator,
            Phase2_CrossScaleCouplingCalculator,
            Phase2_QuantumLevel,
            Phase2_QUANTUM_LEVELS_26,
            
            # Phase 2: DPM Cosmology Module (with Phase2_ prefix)
            Phase2_DPMCosmologyCalculator,
            Phase2_UniversalNuclearCoreModel,
            Phase2_InflationDynamicsCalculator,
            Phase2_DPMCenter,
            phase2_generate_26_centers
        )
    except ImportError as e:
        suite.add_test("Import Phase 2 from CP2", False, str(e))
        return
    
    suite.add_test("Import Phase 2 from CP2", True)
    
    # Test 3.1: Quantum Level Calculator instantiation
    ql_calc = Phase2_QuantumLevel26Calculator()
    E_10_cp2 = ql_calc.compute_energy_density(10)
    suite.add_test(
        "CP2: Quantum Level Calculator Works",
        E_10_cp2 > 0,
        f"E_10 = {E_10_cp2:.2e} J/m³"
    )
    
    # Test 3.2: DPM Cosmology Calculator instantiation
    dpm_calc = Phase2_DPMCosmologyCalculator()
    F_U_cp2 = dpm_calc.compute_inflation_force_at_t0()
    suite.add_test(
        "CP2: DPM Cosmology Calculator Works",
        F_U_cp2 > 0,
        f"F_U = {F_U_cp2:.2e} N"
    )
    
    # Test 3.3: Cross-module consistency (cross import from direct module)
    from QuantumLevel26Framework import QuantumLevel26Calculator as QLC_direct
    calc_direct = QLC_direct()
    E_10_direct = calc_direct.compute_energy_density(10)
    
    suite.add_test(
        "CP2 vs Direct Import Consistency",
        assert_close(E_10_cp2, E_10_direct, rel_tol=1e-10),
        f"CP2: {E_10_cp2:.2e}, Direct: {E_10_direct:.2e}"
    )


# ============================================================================
# MAIN TEST EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("="*80)
    print("PHASE 2 VALIDATION TEST SUITE - COMPLETE FRAMEWORK")
    print("="*80)
    print("\nTesting:")
    print("  • QuantumLevel26Framework.py (630 lines)")
    print("  • DPMCosmologyModule.py (565 lines)")
    print("  • CondensedPhysics2.py integration")
    print("="*80)
    
    # Initialize suite
    suite = ValidationSuite("Phase 2 Complete Validation")
    
    # Run all test suites
    test_quantum_level_26_framework(suite)
    test_dpm_cosmology_module(suite)
    test_condensedphysics2_integration(suite)
    
    # Final summary
    print("\n" + "="*80)
    print("VALIDATION TEST SUMMARY")
    print("="*80)
    
    passed, total = suite.summary()
    success_rate = (passed / total) * 100 if total > 0 else 0
    
    print(f"\nTotal Tests: {total}")
    print(f"[+] Passed: {passed}")
    print(f"[-] Failed: {total - passed}")
    print(f"Success Rate: {success_rate:.1f}%")
    
    print("\n" + "="*80)
    
    if passed == total:
        print("[SUCCESS] ALL TESTS PASSED - PHASE 2 INTEGRATION VERIFIED")
        print("="*80)
        print("\nPhase 2 Complete:")
        print("  * 26-Quantum Level Framework: OPERATIONAL")
        print("  * DPM Cosmology Module: OPERATIONAL")
        print("  * CondensedPhysics2 Integration: VERIFIED")
        print("  * Total New Code: 1,195 lines")
        print("  * Project Completion: 86% (Phase 1: 55% + Phase 2: 31%)")
        print("="*80)
        sys.exit(0)
    else:
        print("[ERROR] SOME TESTS FAILED - CHECK DETAILS ABOVE")
        print("="*80)
        sys.exit(1)
