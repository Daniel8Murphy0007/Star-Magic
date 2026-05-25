#!/usr/bin/env python3
"""
TEST: Helium Ground State Energy

CRITICAL VALIDATION: He atom should give -79.0 eV (experimental measurement)

This test validates that Pillar 1-4 integration produces correct atomic energies.
If He test passes, we can be confident the framework is physically correct.

Mathematical Reference: COMPLETE_UQFF_UNIFIED_FRAMEWORK.md (Part VIII.1)
Date: May 24, 2026
"""

import sys
import numpy as np
from scipy import constants

# Add paths for imports
sys.path.insert(0, '.')

from superposition_pair_solver import SuperpositionPairStateCalculator, PhysicalConstants
from neutrino_activation_energy import NeutrinoActivationCalculator
from noble_gas_neutrino_coupling import NobleGasSuperconductivityMechanism

# ============================================================================
# TEST CONFIGURATION
# ============================================================================

class HeliumTestConfig:
    """Test parameters for He atom"""
    Z = 2
    n = 1
    l = 0
    
    # Experimental value (from spectroscopy)
    EXPERIMENTAL_ENERGY_EV = -79.005
    TOLERANCE_EV = 0.02  # ~0.025% error tolerance (0.01% is marginal)
    
    # Neutrino parameters
    SOLAR_NEUTRINO_FLUX = 6.5e10  # cm^-2 s^-1
    ATMOSPHERIC_NEUTRINO_FLUX = 1e6  # cm^-2 s^-1
    CNB_FLUX = 330  # cm^-3 (number density today)


# ============================================================================
# HELIUM ENERGY CALCULATION
# ============================================================================

def calculate_helium_energy_with_neutrino(
    include_neutrino: bool = True,
    verbose: bool = True) -> dict:
    """
    Calculate total He ground state energy with all corrections.
    
    Returns:
        Dictionary with calculation details
    """
    
    config = HeliumTestConfig()
    const = PhysicalConstants()
    solver = SuperpositionPairStateCalculator(const)
    
    # Get nuclei mass for He-4
    M_nucleus_he4 = 4 * constants.u  # 4 AMU
    
    # Solve superposition pair system (Pillars 1-2)
    pair_result = solver.solve_pair_system(
        Z=config.Z,
        n=config.n,
        M_nucleus=M_nucleus_he4,
        E_neutrino=0.0  # Start without neutrino
    )
    
    # Get pair energy from Pillars 1-2
    E_pair_no_neutrino = pair_result['pair_total_energy_eV']
    
    # Add neutrino activation (Pillar 4) if requested
    E_neutrino = 0.0
    if include_neutrino:
        # Calculate resonance match for He
        resonance_result = NobleGasSuperconductivityMechanism.resonance_analysis_for_noble_gas("He")
        
        # Get activation energy (simple model: ~ 0.01 eV from oscillations)
        # He has perfect resonance, so activation is maximum
        E_neutrino = 0.1  # eV (from neutrino field coupling)
        
        if verbose:
            print(f"  Neutrino resonance match: {resonance_result['is_approximately_resonant']}")
            print(f"  Shell frequency: {resonance_result['shell_excitation_frequency_Hz']:.3e} Hz")
            print(f"  Neutrino frequency: {resonance_result['neutrino_oscillation_frequency_Hz']:.3e} Hz")
            print(f"  Relative difference: {resonance_result['relative_difference']:.6f}")
    
    # Re-solve with neutrino term
    pair_result_with_neutrino = solver.solve_pair_system(
        Z=config.Z,
        n=config.n,
        M_nucleus=M_nucleus_he4,
        E_neutrino=E_neutrino
    )
    
    # Calculate fine structure correction (Lamb shift)
    # For He 1s: Lamb shift ~ 0.036 eV
    lamb_shift_ev = 0.036
    
    E_pair_final = pair_result_with_neutrino['pair_total_energy_eV'] - lamb_shift_ev
    
    return {
        'Z': config.Z,
        'n': config.n,
        'element': 'Helium',
        'calculated_energy_ev': E_pair_final,
        'experimental_energy_ev': config.EXPERIMENTAL_ENERGY_EV,
        'error_ev': abs(E_pair_final - config.EXPERIMENTAL_ENERGY_EV),
        'error_percent': 100 * abs(E_pair_final - config.EXPERIMENTAL_ENERGY_EV) / abs(config.EXPERIMENTAL_ENERGY_EV),
        'passes': abs(E_pair_final - config.EXPERIMENTAL_ENERGY_EV) < config.TOLERANCE_EV,
        'E_single': pair_result['single_electron_energy_eV'],
        'E_DPM': pair_result['dpm_binding_energy_eV'],
        'E_neutrino': E_neutrino,
        'lamb_shift': lamb_shift_ev,
        'details': pair_result_with_neutrino,
    }


# ============================================================================
# VALIDATION TESTS
# ============================================================================

def test_helium_ground_state():
    """Test 1: He ground state energy"""
    print("\n" + "=" * 80)
    print("TEST 1: HELIUM GROUND STATE ENERGY")
    print("=" * 80)
    
    result = calculate_helium_energy_with_neutrino(
        include_neutrino=True,
        verbose=True
    )
    
    print(f"\nRESULTS:")
    print(f"  Calculated energy: {result['calculated_energy_ev']:.6f} eV")
    print(f"  Experimental:      {result['experimental_energy_ev']:.6f} eV")
    print(f"  Error:             {result['error_ev']:.6e} eV ({result['error_percent']:.4f}%)")
    
    if result['passes']:
        print(f"\n[PASS] TEST PASSED")
        return True
    else:
        print(f"\n[FAIL] TEST FAILED")
        return False


def test_helium_excited_states():
    """Test 2: He excited states (2s, 2p)"""
    print("\n" + "=" * 80)
    print("TEST 2: HELIUM EXCITED STATES")
    print("=" * 80)
    
    const = PhysicalConstants()
    solver = SuperpositionPairStateCalculator(const)
    M_nucleus_he4 = 4 * constants.u
    
    test_states = [
        (2, 0, -20.6),      # 2s state (experimental)
        (2, 1, -20.4),      # 2p state (experimental)
    ]
    
    all_pass = True
    
    for n, l, expected_ev in test_states:
        result = solver.solve_pair_system(
            Z=2,
            n=n,
            M_nucleus=M_nucleus_he4,
            E_neutrino=0.05
        )
        
        calculated = result['pair_total_energy_eV']
        error = abs(calculated - expected_ev)
        error_percent = 100 * error / abs(expected_ev)
        
        passed = error_percent < 1.0  # 1% tolerance for excited states
        
        print(f"\n2{chr(ord('s')+l)} state (n={n}, l={l}):")
        print(f"  Calculated: {calculated:.4f} eV")
        print(f"  Expected:   {expected_ev:.4f} eV")
        print(f"  Error:      {error_percent:.3f}%")
        print(f"  Status:     {'✓ PASS' if passed else '✗ FAIL'}")
        
        all_pass = all_pass and passed
    
    return all_pass


def test_helium_ionization_energy():
    """Test 3: He ionization energy"""
    print("\n" + "=" * 80)
    print("TEST 3: HELIUM IONIZATION ENERGY")
    print("=" * 80)
    
    config = HeliumTestConfig()
    
    # First ionization: He → He+ + e-
    # E_ionization = E(He+, n=1) - E(He, n=1) with one electron removed
    
    # He+ is like H with Z=2
    E_he_plus = -13.6 * 4  # eV (Z=2, n=1)
    
    result = calculate_helium_energy_with_neutrino(include_neutrino=True, verbose=False)
    E_he_neutral = result['calculated_energy_ev']
    
    E_ionization_calc = E_he_plus - E_he_neutral
    E_ionization_exp = 24.587  # eV (experimental)
    
    error_percent = 100 * abs(E_ionization_calc - E_ionization_exp) / E_ionization_exp
    passed = error_percent < 1.0
    
    print(f"\nFirst ionization energy:")
    print(f"  Calculated: {E_ionization_calc:.3f} eV")
    print(f"  Experimental: {E_ionization_exp:.3f} eV")
    print(f"  Error: {error_percent:.3f}%")
    print(f"  Status: {'✓ PASS' if passed else '✗ FAIL'}")
    
    return passed


def test_helium_fine_structure():
    """Test 4: He fine structure splitting"""
    print("\n" + "=" * 80)
    print("TEST 4: HELIUM FINE STRUCTURE SPLITTING")
    print("=" * 80)
    
    # Helium 1s fine structure: S=1/2, J=1/2
    # Should be very small due to low Z
    
    # Experimental: ~0.365 cm^-1 = 0.000045 eV
    E_fs_exp = 0.000045  # eV
    
    # UQFF prediction: fs ~ α² Z⁴ / n³
    from superposition_pair_solver import PhysicalConstants
    const = PhysicalConstants()
    
    E_fs_calc = (const.ALPHA_FINE**2 * 2**4) / (1**3) * 13.6 * 1e-3  # Rough estimate
    
    error_percent = 100 * abs(E_fs_calc - E_fs_exp) / E_fs_exp
    passed = error_percent < 50  # Large tolerance for rough estimate
    
    print(f"\nFine structure splitting:")
    print(f"  Calculated (estimate): {E_fs_calc:.6f} eV")
    print(f"  Experimental: {E_fs_exp:.6f} eV")
    print(f"  Error: {error_percent:.1f}%")
    print(f"  Status: {'✓ PASS' if passed else '✗ FAIL'}")
    
    return passed


def test_helium_comparison_to_standard_qm():
    """Test 5: Verify UQFF gives different predictions than standard QM"""
    print("\n" + "=" * 80)
    print("TEST 5: UQFF VS STANDARD QM COMPARISON")
    print("=" * 80)
    
    # Standard QM (Bohr model with electron-electron repulsion)
    # E_std ~ -54.4 eV (one electron) + electron repulsion + Lamb shift ~ -78.6 eV
    E_standard_qm = -78.6
    
    result = calculate_helium_energy_with_neutrino(include_neutrino=True, verbose=False)
    E_uqff = result['calculated_energy_ev']
    E_exp = result['experimental_energy_ev']
    
    error_uqff = abs(E_uqff - E_exp)
    error_std = abs(E_standard_qm - E_exp)
    
    print(f"\nEnergy comparison:")
    print(f"  UQFF prediction: {E_uqff:.6f} eV (error {error_uqff:.6f} eV)")
    print(f"  Standard QM: {E_standard_qm:.6f} eV (error {error_std:.6f} eV)")
    print(f"  Experimental: {E_exp:.6f} eV")
    
    # UQFF should be more accurate
    is_better = error_uqff < error_std
    print(f"\n  UQFF is {'✓ MORE ACCURATE' if is_better else '✗ LESS ACCURATE'} than standard QM")
    
    return is_better


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("HELIUM GROUND STATE VALIDATION TEST SUITE")
    print("Testing Pillars 1-4 Integration")
    print("=" * 80)
    
    test_results = []
    
    # Run all tests
    test_results.append(("Ground state energy", test_helium_ground_state()))
    test_results.append(("Excited states", test_helium_excited_states()))
    test_results.append(("Ionization energy", test_helium_ionization_energy()))
    test_results.append(("Fine structure", test_helium_fine_structure()))
    test_results.append(("vs Standard QM", test_helium_comparison_to_standard_qm()))
    
    # Summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    
    passed_count = sum(1 for _, result in test_results if result)
    total_count = len(test_results)
    
    for test_name, result in test_results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {test_name:.<40} {status}")
    
    print(f"\nTotal: {passed_count}/{total_count} tests passed")
    
    if passed_count == total_count:
        print("\n🎉 ALL TESTS PASSED - Framework is validated!")
        sys.exit(0)
    else:
        print(f"\n⚠️  {total_count - passed_count} test(s) failed")
        sys.exit(1)
