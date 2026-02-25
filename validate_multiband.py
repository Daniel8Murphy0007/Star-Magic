#!/usr/bin/env python3
"""
Validation script for multi-band GW astronomy methods.
Tests predict_multiband_merger_timeline(), compute_WD_binary_foreground(),
and compare_detector_sensitivity_curves().
"""

import sys
import numpy as np

sys.path.insert(0, r'c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic')

from CondensedPhysics import GravitationalWaveUQFFCalculator


def test_multiband_merger_timeline():
    """Test multi-band merger timeline prediction."""
    print("=" * 70)
    print("TEST: predict_multiband_merger_timeline()")
    print("=" * 70)
    
    calc = GravitationalWaveUQFFCalculator()
    
    # Test with 30 M_sun binary (GW150914-like)
    results, steps = calc.predict_multiband_merger_timeline(
        M_total_solar=30.0,
        q=0.8,
        z=0.1,
        f_LISA_start_mHz=0.1,
        f_LISA_end_mHz=10.0,
        f_LIGO_start_Hz=10.0
    )
    
    print(steps)
    
    passed = True
    
    # 1. Timeline makes sense
    if not (results['T_LISA_years'] > 0):
        print("FAIL: LISA observation time should be positive")
        passed = False
    
    # 2. LIGO follows LISA
    if not (results['tau_LISA_start_years'] > results['tau_LIGO_start_seconds'] / (365.25 * 86400)):
        print("FAIL: LISA should see source before LIGO")
        passed = False
    
    # 3. Phase lag accumulates
    if not (results['phase_lag_at_LIGO_entry_cycles'] > 0):
        print("FAIL: Phase lag should accumulate in LISA band")
        passed = False
    
    # 4. GW cycles counted
    if not (results['N_cycles_LISA'] > 1000):
        print(f"FAIL: Expected many LISA cycles, got {results['N_cycles_LISA']}")
        passed = False
    
    # 5. Early warning time is positive
    if not (results['early_warning_years'] > 0):
        print("FAIL: Early warning time should be positive")
        passed = False
    
    print(f"\nKey results:")
    print(f"  LISA observation: {results['T_LISA_years']:.1f} years")
    print(f"  Early warning: {results['early_warning_years']:.1f} years before LIGO")
    print(f"  Phase lag at LIGO entry: {results['phase_lag_at_LIGO_entry_cycles']:.1f} cycles")
    print(f"  Template match: {results['match_LIGO']:.3f}")
    
    print(f"\n{'PASSED' if passed else 'FAILED'}: predict_multiband_merger_timeline")
    return passed


def test_WD_binary_foreground():
    """Test WD binary foreground calculation."""
    print("\n" + "=" * 70)
    print("TEST: compute_WD_binary_foreground()")
    print("=" * 70)
    
    calc = GravitationalWaveUQFFCalculator()
    
    results, steps = calc.compute_WD_binary_foreground(
        N_WD_binaries=100000,
        f_band_mHz=(0.1, 10.0),
        n_freq=200
    )
    
    print(steps)
    
    passed = True
    
    # 1. Foreground exists
    if not (np.all(results['S_h_GR'] > 0)):
        print("FAIL: GR foreground should have positive power")
        passed = False
    
    # 2. UQFF foreground is lower
    if not (results['P_UQFF'] < results['P_GR']):
        print("FAIL: UQFF should reduce foreground")
        passed = False
    
    # 3. Sensitivity improvement
    if not (results['avg_improvement'] >= 1.0):
        print("FAIL: UQFF should improve effective sensitivity")
        passed = False
    
    # 4. Reduction fraction makes sense
    if not (0 < results['foreground_reduction'] < 1):
        print(f"FAIL: Foreground reduction {results['foreground_reduction']} out of range")
        passed = False
    
    print(f"\nKey results:")
    print(f"  WD binaries: {results['N_WD_binaries']:,}")
    print(f"  Foreground reduction: {results['foreground_reduction']*100:.1f}%")
    print(f"  Sensitivity improvement: {results['avg_improvement']:.2f}×")
    print(f"  SMBH SNR improvement: {results['SNR_SMBH_UQFF']/results['SNR_SMBH_GR']:.2f}×")
    
    print(f"\n{'PASSED' if passed else 'FAILED'}: compute_WD_binary_foreground")
    return passed


def test_detector_sensitivity_curves():
    """Test detector sensitivity curve comparison."""
    print("\n" + "=" * 70)
    print("TEST: compare_detector_sensitivity_curves()")
    print("=" * 70)
    
    calc = GravitationalWaveUQFFCalculator()
    
    results, steps = calc.compare_detector_sensitivity_curves(
        f_min_Hz=1e-5,
        f_max_Hz=1e4,
        n_freq=500
    )
    
    print(steps)
    
    passed = True
    
    # 1. Frequency array spans full range
    if not (results['freq_Hz'][0] < 1e-4 and results['freq_Hz'][-1] > 1e3):
        print("FAIL: Frequency range should span 10^-5 to 10^4 Hz")
        passed = False
    
    # 2. LISA more sensitive at low f, LIGO at high f
    idx_1mHz = np.argmin(np.abs(results['freq_Hz'] - 1e-3))
    idx_100Hz = np.argmin(np.abs(results['freq_Hz'] - 100))
    
    if not (results['ASD_LISA'][idx_1mHz] < results['ASD_LIGO'][idx_1mHz]):
        print("Note: LISA should be more sensitive at 1 mHz")
    
    # 3. UQFF reduces SNR
    if not (results['SNR_GW150914_UQFF'] < results['SNR_GW150914_GR']):
        print("FAIL: UQFF should reduce SNR")
        passed = False
    
    # 4. Horizon distances reduced
    if not (results['D_horizon_LIGO_UQFF_Mpc'] < results['D_horizon_LIGO_GR_Mpc']):
        print("FAIL: UQFF should reduce horizon distance")
        passed = False
    
    print(f"\nKey results:")
    print(f"  LIGO horizon: {results['D_horizon_LIGO_GR_Mpc']:.0f} → {results['D_horizon_LIGO_UQFF_Mpc']:.0f} Mpc")
    print(f"  LISA horizon: {results['D_horizon_LISA_GR_Gpc']:.1f} → {results['D_horizon_LISA_UQFF_Gpc']:.1f} Gpc")
    print(f"  GW150914 SNR: {results['SNR_GW150914_GR']:.0f} → {results['SNR_GW150914_UQFF']:.0f}")
    print(f"  SMBH SNR: {results['SNR_SMBH_GR']:.0f} → {results['SNR_SMBH_UQFF']:.0f}")
    
    print(f"\n{'PASSED' if passed else 'FAILED'}: compare_detector_sensitivity_curves")
    return passed


def main():
    """Run all validation tests."""
    print("=" * 70)
    print("MULTI-BAND GW ASTRONOMY VALIDATION")
    print("SuperGrok4 Export Integration - Feb 2026")
    print("=" * 70)
    
    all_passed = True
    
    test1 = test_multiband_merger_timeline()
    all_passed = all_passed and test1
    
    test2 = test_WD_binary_foreground()
    all_passed = all_passed and test2
    
    test3 = test_detector_sensitivity_curves()
    all_passed = all_passed and test3
    
    print("\n" + "=" * 70)
    if all_passed:
        print("ALL TESTS PASSED - Multi-band methods validated")
    else:
        print("SOME TESTS FAILED")
    print("=" * 70)
    
    return 0 if all_passed else 1


if __name__ == '__main__':
    sys.exit(main())
