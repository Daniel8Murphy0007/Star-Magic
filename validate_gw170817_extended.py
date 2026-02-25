"""Validate GW170817 extended methods: 800Hz, tidal deformability, B-field scan."""
import sys
sys.path.insert(0, '.')
from CondensedPhysics import GravitationalWaveUQFFCalculator

# Create calculator
calc = GravitationalWaveUQFFCalculator(
    M_chirp=1.188,
    D_L=40.0,
    iota=0.0
)

print("="*70)
print("TEST 1: simulate_GW170817_800Hz() - Full frequency range")
print("="*70)

results1, _ = calc.simulate_GW170817_800Hz(
    t_duration=100.0,
    f_start=23.0,
    f_end=800.0,
    f_TRZ=0.1,
    B_NS=1e10,
    string_factor=0.37,
    Lambda_tidal=600.0
)

print(f"Frequency range: {results1['f_start']:.0f} - {results1['f_end']:.0f} Hz")
print(f"Total cycles: {results1['total_cycles']}")
print(f"Tidal Λ = {results1['Lambda_tidal']:.0f}, correction at merger: {results1['tidal_correction_at_merger']:.4f}")
print(f"Peak GR: {results1['h_peak_GR']:.4e}, Peak UQFF: {results1['h_peak_UQFF']:.4e}")
print(f"Reduction: {results1['percent_reduction']:.1f}%")
print(f"SNR: {results1['SNR_GR_final']:.1f} → {results1['SNR_UQFF_final']:.1f}")
print(f"E_GW = {results1['E_GW_solar']:.2f} M_sun c²")

# Verification
test1_pass = (
    results1['f_end'] == 800.0 and
    results1['total_cycles'] > 3000 and
    50 <= results1['percent_reduction'] <= 70
)
print(f"TEST 1: {'PASS' if test1_pass else 'FAIL'}\n")

print("="*70)
print("TEST 2: compute_tidal_deformability_constraint() - NS EOS effects")
print("="*70)

results2, _ = calc.compute_tidal_deformability_constraint(
    Lambda_range=(100, 1000),
    n_Lambda=10,
    f_TRZ=0.1,
    string_factor=0.37
)

print(f"Λ range: {results2['Lambda_range']}")
print(f"GW170817 measured: Λ = {results2['Lambda_measured']}, 90% upper: {results2['Lambda_upper_90']}")
print(f"UQFF combined factor: {results2['combined_factor']:.4f}")
print(f"Λ bias factor: {results2['Lambda_bias_factor']:.2f}")
print(f"Inferred Λ shift: {results2['Lambda_inferred_shift']:.0f}")
print(f"Mean measurability: {results2['measurability'].mean():.4f}")

# Verification
test2_pass = (
    len(results2['Lambda_values']) == 10 and
    len(results2['tidal_phase_shifts']) == 10 and
    0 < results2['Lambda_bias_factor'] < 1
)
print(f"TEST 2: {'PASS' if test2_pass else 'FAIL'}\n")

print("="*70)
print("TEST 3: scan_NS_magnetic_field_effects() - B-field parameter space")
print("="*70)

results3, _ = calc.scan_NS_magnetic_field_effects(
    B_range=(1e8, 1e14),
    n_B=12,
    f_TRZ=0.1,
    string_factor=0.37
)

print(f"B range: {results3['B_range_Gauss'][0]:.0e} - {results3['B_range_Gauss'][1]:.0e} G")
print(f"B_crit: {results3['B_crit_Gauss']:.2e} G")
print(f"Number of B values: {len(results3['B_values_Gauss'])}")
print(f"SCm factors range: {results3['SCm_factors'].min():.6f} - {results3['SCm_factors'].max():.6f}")
print(f"Detectable configurations: {results3['n_detectable']} of {len(results3['B_values_Gauss'])}")
print(f"B detection limit: {results3['B_detection_limit_Gauss']:.2e} G")

# Show a few rows
print("\nSample B-field scan:")
for i in [0, 4, 8, 11]:
    if i < len(results3['B_values_Gauss']):
        B = results3['B_values_Gauss'][i]
        SCm = results3['SCm_factors'][i]
        SNR = results3['SNR_values'][i]
        ns_type = results3['ns_types'][i]
        print(f"  B={B:.2e} G ({ns_type}): SCm={SCm:.6f}, SNR={SNR:.1f}")

# Verification - realistic NS B-fields (up to 1e15 G) never approach B_crit (4.41e17 G)
# SCm remains close to 1.0 for all realistic NS configurations
test3_pass = (
    len(results3['B_values_Gauss']) == 12 and
    results3['SCm_factors'][0] > 0.99 and  # Low B → SCm ≈ 1.0
    results3['SCm_factors'][-1] > 0.999 and  # Even magnetar B → SCm ≈ 1.0 (B << B_crit)
    results3['SCm_factors'][-1] < results3['SCm_factors'][0] and  # Higher B → slightly lower SCm
    results3['n_detectable'] > 0
)
print(f"\nTEST 3: {'PASS' if test3_pass else 'FAIL'}\n")

print("="*70)
all_pass = test1_pass and test2_pass and test3_pass
print(f"ALL TESTS: {'PASSED' if all_pass else 'SOME FAILED'}")
print("="*70)

import sys
sys.exit(0 if all_pass else 1)
