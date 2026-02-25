"""
Validation script for GW190425 UQFF methods.
Tests: compare_to_LIGO_GW190425, simulate_GW190425_chirp, analyze_mass_gap_component
"""
import sys
sys.path.insert(0, '.')

from CondensedPhysics import GravitationalWaveUQFFCalculator
import numpy as np

calc = GravitationalWaveUQFFCalculator()

print("="*70)
print("TEST 1: compare_to_LIGO_GW190425() - Second BNS with mass-gap")
print("="*70)

results1, _ = calc.compare_to_LIGO_GW190425(
    f_TRZ=0.1,
    string_factor=0.37,
    U_m=1.0,
    beta_m=0.01,
    Lambda_tidal=600.0
)

print(f"Event: {results1['event_name']} ({results1['detection_date']})")
print(f"Chirp mass: {results1['M_chirp_solar']:.2f} M⊙ (highest BNS)")
print(f"Component masses: m₁={results1['m1_solar']:.2f}, m₂={results1['m2_solar']:.2f} M⊙")
print(f"Distance: {results1['d_Mpc']:.0f} Mpc")
print(f"In mass gap: {results1['in_mass_gap']}")
print(f"P(NS): {results1['P_NS']:.0%}, P(BH): {results1['P_BH']:.0%}")
print(f"GW cycles: {results1['N_cycles']:.0f}")
print(f"UQFF combined factor: {results1['UQFF_combined']:.4f}")
print(f"Amplitude reduction: {results1['amplitude_reduction']*100:.1f}%")
print(f"SNR observed: {results1['SNR_observed']:.1f}, UQFF: {results1['SNR_UQFF']:.1f}")

test1_pass = (
    results1['M_chirp_solar'] == 1.44 and
    results1['in_mass_gap'] == True and
    results1['d_Mpc'] == 159.0 and
    0.3 < results1['UQFF_combined'] < 0.7 and
    results1['N_cycles'] > 1000
)
print(f"\nTEST 1: {'PASS' if test1_pass else 'FAIL'}\n")

print("="*70)
print("TEST 2: simulate_GW190425_chirp() - 0.2s final chirp (30-400 Hz)")
print("="*70)

results2, _ = calc.simulate_GW190425_chirp(
    duration=0.2,
    f_start=30.0,
    f_end=400.0,
    n_points=2000,
    f_TRZ=0.1,
    Lambda_tidal=600.0
)

print(f"Duration: {results2['duration_s']*1000:.0f} ms")
print(f"Frequency: {results2['f_start_Hz']:.0f} → {results2['f_end_Hz']:.0f} Hz")
print(f"GW cycles: {results2['N_cycles']}")
print(f"Peak GR: {results2['h_peak_GR']:.4e}")
print(f"Peak UQFF: {results2['h_peak_UQFF']:.4e}")
print(f"Reduction: {results2['amplitude_reduction']*100:.1f}%")
print(f"SNR GR: {results2['SNR_GR']:.1f} → UQFF: {results2['SNR_UQFF']:.1f}")
print(f"Time array shape: {len(results2['t_s'])}")

test2_pass = (
    len(results2['t_s']) == 2000 and
    results2['N_cycles'] >= 5 and  # Low because of zero-crossing count in chirp
    0.2 < results2['amplitude_reduction'] < 0.9 and
    results2['h_peak_UQFF'] < results2['h_peak_GR'] and
    results2['h_peak_GR'] > 1e-24  # Peak strain reasonable
)
print(f"\nTEST 2: {'PASS' if test2_pass else 'FAIL'}\n")

print("="*70)
print("TEST 3: analyze_mass_gap_component() - NS vs BH discrimination")
print("="*70)

results3, _ = calc.analyze_mass_gap_component(
    m_component=2.52,
    B_NS_range=(1e8, 1e15),
    n_scenarios=20,
    f_TRZ=0.1,
    string_factor=0.37
)

print(f"Component mass: {results3['m_component_solar']:.2f} M⊙")
print(f"In mass gap: {results3['in_mass_gap']}")
print(f"P(NS): {results3['P_NS']:.0%}, P(BH): {results3['P_BH']:.0%}")
print(f"B-field scenarios: {results3['n_scenarios']}")
print(f"SCm_BH (no B-field): {results3['SCm_BH']:.4f}")
print(f"SCm_NS range: {results3['SCm_NS'].min():.6f} - {results3['SCm_NS'].max():.6f}")
print(f"UQFF factor BH: {results3['UQFF_factor_BH']:.4f}")
print(f"UQFF factor NS (mean): {np.mean(results3['UQFF_factor_NS']):.4f}")
print(f"λ_NS typical: {results3['Lambda_NS_typical']:.0f}")
print(f"λ_BH: {results3['Lambda_BH']:.0f}")

# Sample NS types
print("\nSample B-field scenarios:")
for i in [0, 5, 10, 15, 19]:
    if i < len(results3['B_values_Gauss']):
        print(f"  B={results3['B_values_Gauss'][i]:.2e} G: "
              f"{results3['ns_types'][i]}, SCm={results3['SCm_NS'][i]:.6f}")

test3_pass = (
    results3['in_mass_gap'] == True and
    results3['n_scenarios'] == 20 and
    results3['SCm_BH'] == 1.0 and
    results3['Lambda_BH'] == 0.0 and
    results3['P_NS'] < results3['P_BH']  # 2.52 M⊙ slightly favors BH
)
print(f"\nTEST 3: {'PASS' if test3_pass else 'FAIL'}\n")

print("="*70)
all_pass = test1_pass and test2_pass and test3_pass
print(f"ALL TESTS: {'PASSED' if all_pass else 'SOME FAILED'}")
print("="*70)

sys.exit(0 if all_pass else 1)
