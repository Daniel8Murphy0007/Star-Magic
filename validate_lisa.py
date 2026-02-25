"""
Validation script for LISA UQFF prediction methods.
Tests: predict_LISA_SMBH_merger, simulate_LISA_EMRI, compute_LISA_detection_rates
"""
import sys
sys.path.insert(0, '.')

from CondensedPhysics import GravitationalWaveUQFFCalculator
import numpy as np

calc = GravitationalWaveUQFFCalculator()

print("="*70)
print("TEST 1: predict_LISA_SMBH_merger() - SMBH merger at z=1")
print("="*70)

results1, _ = calc.predict_LISA_SMBH_merger(
    M_total_solar=1e6,
    q=0.5,
    z=1.0,
    f_TRZ=0.1,
    string_factor=0.37,
    U_m=1.0
)

print(f"M_total: {results1['M_total_solar']:.2e} M⊙")
print(f"Chirp mass: {results1['M_chirp_solar']:.2e} M⊙")
print(f"Distance: z={results1['z']:.1f}, D_L={results1['D_L_Gpc']:.2f} Gpc")
print(f"f_ISCO (observer): {results1['f_ISCO_obs_Hz']*1000:.3f} mHz")
print(f"Time in LISA band: {results1['T_in_band_years']:.2f} years")
print(f"GW cycles: {results1['N_cycles']:.2e}")
print(f"h_peak GR: {results1['h_peak_GR']:.4e}")
print(f"h_peak UQFF: {results1['h_peak_UQFF']:.4e}")
print(f"UQFF combined factor: {results1['UQFF_combined']:.4f}")
print(f"Amplitude reduction: {results1['amplitude_reduction']*100:.1f}%")
print(f"SNR GR: {results1['SNR_GR']:.0f}, SNR UQFF: {results1['SNR_UQFF']:.1f}")
print(f"Detectable GR: {results1['detectable_GR']}, UQFF: {results1['detectable_UQFF']}")

test1_pass = (
    results1['M_total_solar'] == 1e6 and
    results1['z'] == 1.0 and
    results1['T_in_band_years'] > 0 and
    results1['UQFF_combined'] < 1.0 and
    results1['h_peak_UQFF'] < results1['h_peak_GR'] and
    results1['N_cycles'] > 1e4
)
print(f"\nTEST 1: {'PASS' if test1_pass else 'FAIL'}\n")

print("="*70)
print("TEST 2: simulate_LISA_EMRI() - Stellar BH into SMBH")
print("="*70)

results2, _ = calc.simulate_LISA_EMRI(
    M_SMBH_solar=1e6,
    m_compact_solar=10.0,
    compact_type='BH',
    z=0.5,
    T_observe_years=2.0,
    n_harmonics=5
)

print(f"SMBH: {results2['M_SMBH_solar']:.2e} M⊙")
print(f"Compact object: {results2['m_compact_solar']:.1f} M⊙ ({results2['compact_type']})")
print(f"Mass ratio: q = {results2['q']:.2e}")
print(f"Distance: z={results2['z']:.1f}, D_L={results2['D_L_Gpc']:.2f} Gpc")
print(f"f_ISCO: {results2['f_ISCO_obs_mHz']:.3f} mHz")
print(f"Observation: {results2['T_observe_years']:.1f} years")
print(f"Orbits GR: {results2['N_orbits_GR']:.2e}")
print(f"String harmonics: {results2['n_harmonics']}")
print(f"Harmonic frequencies (mHz): {results2['string_harmonics_mHz'][:3]}")
print(f"Stability factor: {results2['stability_factor']:.2f}")
print(f"h_peak UQFF: {results2['h_peak_UQFF']:.4e}")
print(f"SNR GR: {results2['SNR_GR']:.0f}, UQFF: {results2['SNR_UQFF']:.1f}")

test2_pass = (
    results2['q'] < 1e-3 and  # Extreme mass ratio
    results2['n_harmonics'] == 5 and
    len(results2['string_harmonics_mHz']) == 5 and
    results2['stability_factor'] > 1.0 and  # UQFF stabilizes
    results2['N_orbits_UQFF'] > results2['N_orbits_GR']  # More orbits
)
print(f"\nTEST 2: {'PASS' if test2_pass else 'FAIL'}\n")

print("="*70)
print("TEST 3: compute_LISA_detection_rates() - Rate predictions")
print("="*70)

results3, _ = calc.compute_LISA_detection_rates(
    f_TRZ=0.1,
    string_factor=0.37,
    U_m=1.0,
    z_max=10.0,
    n_z_bins=20
)

print(f"UQFF parameters: f_TRZ={results3['f_TRZ']}, A_TRZ={results3['A_TRZ']:.2f}")
print(f"z_max detectable GR: {results3['z_max_detectable_GR']:.1f}")
print(f"z_max detectable UQFF: {results3['z_max_detectable_UQFF']:.1f}")
print(f"Volume ratio: {results3['volume_ratio']:.2f}")

R_SMBH_GR = results3['R_SMBH_GR'][1]  # Mid estimate
R_SMBH_UQFF = results3['R_SMBH_UQFF'][1]
R_EMRI_GR = results3['R_EMRI_GR'][1]
R_EMRI_UQFF = results3['R_EMRI_UQFF'][1]

print(f"\nSMBH mergers:")
print(f"  GR: {R_SMBH_GR:.0f}/year, UQFF: {R_SMBH_UQFF:.1f}/year")
print(f"  Missing: {results3['missing_SMBH_per_year']:.0f}/year")

print(f"\nEMRIs:")
print(f"  GR: {R_EMRI_GR:.0f}/year, UQFF: {R_EMRI_UQFF:.1f}/year")
print(f"  Missing: {results3['missing_EMRI_per_year']:.0f}/year")

print(f"\nWD binaries: {results3['N_WD_GR']} (GR) → {results3['N_WD_UQFF']} (UQFF)")

print(f"\nStochastic background:")
print(f"  GR: Ω_GW = {results3['Omega_GW_GR']:.1e}")
print(f"  UQFF: Ω_GW = {results3['Omega_GW_UQFF']:.1e}")

print(f"\nTotal missing (SMBH+EMRI): {results3['total_missing_per_year']:.0f}/year")

test3_pass = (
    len(results3['z_values']) == 20 and
    results3['z_max_detectable_UQFF'] <= results3['z_max_detectable_GR'] and
    R_SMBH_UQFF < R_SMBH_GR and  # Fewer detections
    R_EMRI_UQFF < R_EMRI_GR and
    results3['N_WD_UQFF'] <= results3['N_WD_GR'] and
    results3['total_missing_per_year'] > 0
)
print(f"\nTEST 3: {'PASS' if test3_pass else 'FAIL'}\n")

print("="*70)
all_pass = test1_pass and test2_pass and test3_pass
print(f"ALL TESTS: {'PASSED' if all_pass else 'SOME FAILED'}")
print("="*70)

sys.exit(0 if all_pass else 1)
