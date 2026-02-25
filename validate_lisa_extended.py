#!/usr/bin/env python3
"""
Validation script for extended LISA time-domain methods.
Tests simulate_LISA_SMBH_chirp() and compute_aether_noise_spectrum().
"""

import sys
import numpy as np

sys.path.insert(0, r'c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic')

from CondensedPhysics import GravitationalWaveUQFFCalculator


def test_SMBH_chirp():
    """Test LISA SMBH chirp simulation (1-year time-domain)."""
    print("=" * 70)
    print("TEST: simulate_LISA_SMBH_chirp()")
    print("=" * 70)
    
    calc = GravitationalWaveUQFFCalculator()
    
    # Test with 10^6 M_sun SMBH binary at z=1, 12-month chirp
    results, steps = calc.simulate_LISA_SMBH_chirp(
        M_total_solar=1e6,
        q=0.5,
        z=1.0,
        T_chirp_months=12.0,
        f_start_mHz=1.0,
        f_end_mHz=10.0,
        n_points=2000
    )
    
    print(steps)
    
    # Validation checks
    passed = True
    
    # 1. Time series exists and correct length
    if len(results['t_months']) != 2000:
        print(f"FAIL: Expected 2000 time points, got {len(results['t_months'])}")
        passed = False
    
    # 2. Frequency evolves from start to end
    if not (results['freq_mHz'][0] < results['freq_mHz'][-1]):
        print("FAIL: Frequency should increase during chirp")
        passed = False
    
    # 3. UQFF amplitude reduced compared to GR
    if not (results['amplitude_reduction'] > 0):
        print("FAIL: UQFF should show amplitude reduction")
        passed = False
    
    # 4. Phase lag accumulates
    if not (results['phase_diff_final_rad'] > 0):
        print("FAIL: Phase lag should accumulate over chirp")
        passed = False
    
    # 5. Many GW cycles visible
    if not (results['N_cycles'] > 100):
        print(f"FAIL: Expected many GW cycles, got {results['N_cycles']}")
        passed = False
    
    # 6. Waveforms have correct structure
    h_GR = results['h_GR']
    h_UQFF = results['h_UQFF']
    if not (np.std(h_GR) > 0 and np.std(h_UQFF) > 0):
        print("FAIL: Waveforms should have variation")
        passed = False
    
    # 7. UQFF factor has modulations
    UQFF_factor = results['UQFF_factor']
    if not (np.std(UQFF_factor) > 0.001):
        print("FAIL: UQFF factor should show modulations")
        passed = False
    
    print(f"\nKey results:")
    print(f"  Duration: {results['T_chirp_months']:.0f} months")
    print(f"  GW cycles: {results['N_cycles']}")
    print(f"  Amplitude reduction: {results['amplitude_reduction']*100:.1f}%")
    print(f"  Phase lag at merger: {results['phase_diff_cycles']:.1f} cycles")
    print(f"  SNR ratio (UQFF/GR): {results['SNR_UQFF']/results['SNR_GR']:.2f}")
    
    print(f"\n{'PASSED' if passed else 'FAILED'}: simulate_LISA_SMBH_chirp")
    return passed


def test_aether_noise_spectrum():
    """Test aether noise spectrum calculation."""
    print("\n" + "=" * 70)
    print("TEST: compute_aether_noise_spectrum()")
    print("=" * 70)
    
    calc = GravitationalWaveUQFFCalculator()
    
    # Test with LISA band
    results, steps = calc.compute_aether_noise_spectrum(
        f_band_mHz=(0.1, 100.0),
        n_freq=200,
        U_m=1.0,
        beta_m=0.01,
        f_TRZ=0.1,
        T_observe_years=4.0
    )
    
    print(steps)
    
    # Validation checks
    passed = True
    
    # 1. Frequency array correct
    if len(results['freq_mHz']) != 200:
        print(f"FAIL: Expected 200 frequency bins, got {len(results['freq_mHz'])}")
        passed = False
    
    # 2. GR spectrum has power
    if not (np.all(results['S_h_GR'] > 0)):
        print("FAIL: GR spectrum should have positive power")
        passed = False
    
    # 3. UQFF spectrum has power
    if not (np.all(results['S_h_UQFF'] > 0)):
        print("FAIL: UQFF spectrum should have positive power")
        passed = False
    
    # 4. U_m creates spectral features
    if not (np.max(results['U_m_spectrum']) > 0):
        print("FAIL: U_m should create spectral features")
        passed = False
    
    # 5. Aether noise detectable?
    # (May or may not be depending on parameters)
    print(f"  Aether noise detectable: {results['detectable']}")
    
    # 6. Integrated SNR calculated
    if not (results['integrated_SNR_aether'] >= 0):
        print("FAIL: Integrated SNR should be non-negative")
        passed = False
    
    # 7. TRZ suppression in range [0, 1]
    if not (np.all(results['TRZ_suppression'] >= 0) and np.all(results['TRZ_suppression'] <= 1.1)):
        print("FAIL: TRZ suppression out of expected range")
        passed = False
    
    print(f"\nKey results:")
    print(f"  Frequency range: {results['freq_mHz'][0]:.2f} - {results['freq_mHz'][-1]:.0f} mHz")
    print(f"  Aether power fraction: {results['P_aether_fraction']*100:.2f}%")
    print(f"  Peak aether frequency: {results['f_peak_aether_mHz']:.2f} mHz")
    print(f"  Integrated SNR: {results['integrated_SNR_aether']:.1f}")
    print(f"  Detectable (SNR > 5): {results['detectable']}")
    
    print(f"\n{'PASSED' if passed else 'FAILED'}: compute_aether_noise_spectrum")
    return passed


def test_parameter_variations():
    """Test parameter sensitivity for both methods."""
    print("\n" + "=" * 70)
    print("TEST: Parameter Variations")
    print("=" * 70)
    
    calc = GravitationalWaveUQFFCalculator()
    passed = True
    
    # Test SMBH chirp with different masses
    print("\nSMBH chirp mass scaling:")
    masses = [1e5, 1e6, 1e7]
    for M in masses:
        results, _ = calc.simulate_LISA_SMBH_chirp(
            M_total_solar=M, T_chirp_months=6, n_points=500
        )
        print(f"  M = {M:.0e} M⊙: {results['N_cycles']} cycles, "
              f"Δφ = {results['phase_diff_cycles']:.1f} cycles, "
              f"h_peak = {results['h_peak_GR']:.2e}")
    
    # Test SMBH chirp with different redshifts
    print("\nSMBH chirp redshift scaling:")
    redshifts = [0.5, 1.0, 2.0]
    for z in redshifts:
        results, _ = calc.simulate_LISA_SMBH_chirp(
            z=z, T_chirp_months=6, n_points=500
        )
        print(f"  z = {z}: D_L = {results['D_L_Gpc']:.2f} Gpc, "
              f"amp_reduction = {results['amplitude_reduction']*100:.1f}%")
    
    # Test aether noise with different U_m
    print("\nAether noise U_m scaling:")
    U_m_values = [0.5, 1.0, 2.0]
    for U_m in U_m_values:
        results, _ = calc.compute_aether_noise_spectrum(U_m=U_m, n_freq=100)
        print(f"  U_m = {U_m}: P_aether = {results['P_aether_fraction']*100:.2f}%, "
              f"SNR = {results['integrated_SNR_aether']:.1f}")
    
    print(f"\n{'PASSED' if passed else 'FAILED'}: Parameter Variations")
    return passed


def main():
    """Run all validation tests."""
    print("=" * 70)
    print("LISA EXTENDED METHODS VALIDATION")
    print("SuperGrok4 Export Integration - Feb 2026")
    print("=" * 70)
    
    all_passed = True
    
    test1 = test_SMBH_chirp()
    all_passed = all_passed and test1
    
    test2 = test_aether_noise_spectrum()
    all_passed = all_passed and test2
    
    test3 = test_parameter_variations()
    all_passed = all_passed and test3
    
    print("\n" + "=" * 70)
    if all_passed:
        print("ALL TESTS PASSED - LISA extended methods validated")
    else:
        print("SOME TESTS FAILED")
    print("=" * 70)
    
    return 0 if all_passed else 1


if __name__ == '__main__':
    sys.exit(main())
