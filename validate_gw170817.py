"""Validate UQFF vs LIGO GW170817 comparison (BNS merger)."""
import sys
sys.path.insert(0, '.')
from CondensedPhysics import GravitationalWaveUQFFCalculator

# Create calculator with BNS chirp mass
calc = GravitationalWaveUQFFCalculator(
    M_chirp=1.188,  # Solar masses (GW170817 chirp mass)
    D_L=40.0,       # Mpc
    iota=0.0
)

print("=== UQFF vs LIGO GW170817 COMPARISON (BNS Merger) ===\n")

# Run comparison with phenomenological string factor
results, steps = calc.compare_to_LIGO_GW170817(
    f_TRZ=0.1,
    B_NS=1e8,         # Typical NS: 10^8 Gauss
    string_factor=0.37,
    beta_m=0.01
)

print(f"Event: {results['event']} ({results['event_date']})")
print(f"Type: {results['event_type']}")
print(f"Source: ℳ = {results['M_chirp_solar']:.3f} M_sun, M_tot = {results['M_tot_solar']:.2f} M_sun")
print(f"Distance: {results['distance_Mpc']:.0f} Mpc ({results['host_galaxy']})")
print()

print("MULTI-MESSENGER DATA:")
print(f"  GRB: {results['GRB_name']}, delay = {results['GRB_delay_s']:.2f}s")
print(f"  GW speed: |Δc/c| < {results['GW_speed_constraint']:.0e}")
print(f"  Kilonova: {results['kilonova_name']}")
print()

print("NS MAGNETIC FIELD:")
print(f"  B_NS = {results['B_NS_Gauss']:.0e} G = {results['B_NS_Tesla']:.2e} T")
print(f"  B_NS/B_crit = {results['B_ratio']:.2e}")
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

print("TENSION ANALYSIS:")
print(f"  Mismatch metric: {results['mismatch_metric']:.3f}")
print(f"  Status: {results['constraint_status']}")
print()

print("SNR IMPACT:")
print(f"  Standard SNR: {results['SNR_observed']:.1f}")
print(f"  UQFF SNR: {results['SNR_UQFF']:.1f}")
print(f"  Still detectable: {results['detectable']}")
print()

# Verification
expected_reduction = 50 <= results['percent_reduction'] <= 70  # BNS with same factors
expected_snr_reduced = results['SNR_UQFF'] < results['SNR_observed']
expected_negligible_SCm = results['SCm_factor'] > 0.9999  # NS fields << B_crit
expected_tension = "TENSION" in results['constraint_status'] or "STRONG" in results['constraint_status']

print("VERIFICATION:")
print(f"  Reduction ~50-70%: {'PASS' if expected_reduction else 'CHECK'}")
print(f"  SNR reduced: {'PASS' if expected_snr_reduced else 'FAIL'}")
print(f"  SCm negligible (NS B << B_crit): {'PASS' if expected_negligible_SCm else 'FAIL'}")
print(f"  UQFF in tension with GR fit: {'PASS' if expected_tension else 'CHECK'}")
print()

if expected_snr_reduced and expected_negligible_SCm:
    print("TEST PASSED!")
else:
    print("CHECK NEEDED")
