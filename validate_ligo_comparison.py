"""Validate UQFF vs LIGO GW150914 comparison."""
import sys
sys.path.insert(0, '.')
from CondensedPhysics import GravitationalWaveUQFFCalculator

# Create calculator
calc = GravitationalWaveUQFFCalculator(
    M_chirp=28.0,
    D_L=410.0,
    iota=0.0
)

print("=== UQFF vs LIGO GW150914 COMPARISON ===\n")

# Run comparison with phenomenological string factor
results, steps = calc.compare_to_LIGO_GW150914(
    f_TRZ=0.1,
    string_factor=0.37,  # exp(-1) ~ 0.37
    beta_m=0.01
)

print(f"Event: {results['event']} ({results['event_date']})")
print(f"Source: {results['M1_solar']:.0f} + {results['M2_solar']:.0f} M_sun → {results['M_final_solar']:.0f} M_sun")
print(f"Distance: {results['distance_Mpc']:.0f} Mpc")
print()

print("DAMPING FACTORS:")
print(f"  Aether: {results['aether_factor']:.6f}")
print(f"  SCm: {results['SCm_factor']:.6f}")
print(f"  TRZ: {results['TRZ_factor']:.4f}")
print(f"  String: {results['string_factor']:.4f}")
print(f"  Combined: {results['combined_factor']:.4f}")
print()

print("STRAIN COMPARISON:")
print(f"  Standard GR peak: {results['h_peak_GR']:.4e}")
print(f"  UQFF peak: {results['h_peak_UQFF']:.4e}")
print(f"  UQFF from observed: {results['h_UQFF_from_observed']:.4e}")
print(f"  Reduction: {results['percent_reduction']:.1f}%")
print()

print("SNR IMPACT:")
print(f"  Standard SNR: {results['SNR_observed']:.0f}")
print(f"  UQFF SNR: {results['SNR_UQFF']:.1f}")
print(f"  Still detectable: {results['detectable']}")
print()

print("OBSERVABLE SIGNATURES:")
for key, value in results['observable_signatures'].items():
    print(f"  {key}: {value}")
print()

# Verification - with phenomenological string_factor=0.37, expect ~66% reduction
# Edge case: SNR=8.0 is exactly at threshold, so detectable check may fail
expected_reduction = 60 <= results['percent_reduction'] <= 70  # ~66.7% for 0.9 × 0.37
expected_snr = results['SNR_UQFF'] < results['SNR_observed']
expected_detectable = results['SNR_UQFF'] >= 8.0  # At or above threshold

print("VERIFICATION:")
reduction_status = 'PASS' if expected_reduction else f'FAIL ({results["percent_reduction"]:.1f}%)'
print(f"  Reduction ~66%: {reduction_status}")
print(f"  SNR reduced: {'PASS' if expected_snr else 'FAIL'}")
print(f"  Still detectable: {'PASS' if expected_detectable else 'FAIL'}")
print()

if expected_reduction and expected_snr and expected_detectable:
    print("TEST PASSED!")
else:
    print("CHECK NEEDED")
