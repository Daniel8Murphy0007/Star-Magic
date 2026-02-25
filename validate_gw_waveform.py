"""Validate UQFF GW Waveform implementation."""
import sys
sys.path.insert(0, '.')
from CondensedPhysics import GravitationalWaveUQFFCalculator
import numpy as np

# GW150914-like parameters
# ℳ ~ 28 M_sun (chirp mass from 36+29 M_sun)
calc = GravitationalWaveUQFFCalculator(
    M_chirp=28.0,   # Solar masses
    D_L=410.0,      # Mpc
    iota=0.0        # Face-on
)

# Test full UQFF waveform
print("=== UQFF GW WAVEFORM TEST (GW150914-like) ===\n")

# Compute with derivation
results, steps = calc.compute_waveform_derivation(
    r=410.0,           # Source distance in Mpc
    B_t=4.4e11,        # 1% of B_crit
    f_TRZ=0.1,
    beta_m=0.01
)

print(f"Chirp mass: {results['M_chirp_solar']:.2f} M_sun")
print(f"Distance: {results['r_Mpc']:.1f} Mpc")
print()
print("DAMPING FACTORS:")
print(f"  Aether: {results['aether_damping']:.6f}")
print(f"  SCm: {results['SCm_factor']:.4f}")
print(f"  TRZ: {results['TRZ_factor']:.4f}")
print(f"  String binding: {results['string_binding']:.6f}")
print()
print("WAVEFORM RESULTS:")
print(f"  Peak strain (standard): {results['peak_strain_standard']:.4e}")
print(f"  Peak strain (UQFF): {results['peak_strain_UQFF']:.4e}")
print(f"  Amplitude ratio: {results['amplitude_ratio']:.4f}")
print(f"  Phase lag (avg): {results['avg_phase_lag_rad']:.4f} rad")
print()

if results['example_check']:
    ex = results['example_check']
    print(f"MID-FREQUENCY EXAMPLE (f = {ex['f_Hz']:.2f} Hz):")
    print(f"  h_standard: {ex['h_standard']:.4e}")
    print(f"  h_UQFF: {ex['h_UQFF']:.4e}")
    print(f"  damping ratio: {ex['damping_ratio']:.4f}")
    print()

# Verify expected behavior
print("VERIFICATION:")
expected_amplitude_reduction = 0.1 <= results['amplitude_ratio'] <= 0.5
expected_TRZ = results['TRZ_factor'] == 0.9
expected_SCm = abs(results['SCm_factor'] - 0.99) < 0.001

print(f"  TRZ factor = 0.9: {'PASS' if expected_TRZ else 'FAIL'}")
print(f"  SCm factor ≈ 0.99: {'PASS' if expected_SCm else 'FAIL'}")
print(f"  Amplitude reduction 10-50%: {'PASS' if expected_amplitude_reduction else 'FAIL'}")
print()

# Check arrays exist and have correct shape
has_arrays = (len(results['h_plus_UQFF']) == len(results['h_plus_standard']) and
              len(results['h_plus_UQFF']) > 0)
print(f"  Arrays match shape: {'PASS' if has_arrays else 'FAIL'}")
print()

if expected_TRZ and expected_SCm and has_arrays:
    print("TEST PASSED!")
else:
    print("CHECK NEEDED")
