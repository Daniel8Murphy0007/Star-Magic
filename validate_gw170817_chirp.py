"""Validate GW170817 BNS chirp simulation with GR vs UQFF comparison chart."""
import sys
sys.path.insert(0, '.')
from CondensedPhysics import GravitationalWaveUQFFCalculator

# Create calculator with GW170817 chirp mass
calc = GravitationalWaveUQFFCalculator(
    M_chirp=1.188,  # Solar masses
    D_L=40.0,       # Mpc
    iota=0.0
)

print("=== GW170817 BNS CHIRP SIMULATION ===\n")

# Run simulation with 0.2s 35-300 Hz chirp
results, steps = calc.simulate_GW170817_chirp(
    t_duration=0.2,
    f_start=35.0,
    f_end=300.0,
    f_TRZ=0.1,
    B_NS=1e8,
    string_factor=0.37,
    beta_m=0.01,
    n_samples=200,
    show_chart=True
)

# Print key statistics
print(f"Event: {results['event']}")
print(f"Chirp: {results['f_start']:.0f} Hz → {results['f_end']:.0f} Hz over {results['t_duration']}s")
print(f"Samples: {results['n_samples']}")
print()

print("DAMPING FACTORS:")
print(f"  Aether: {results['aether_factor']:.6f}")
print(f"  SCm: {results['SCm_factor']:.6f}")
print(f"  TRZ: {results['TRZ_factor']:.4f}")
print(f"  String: {results['string_factor']:.4f}")
print(f"  Combined: {results['combined_factor']:.4f}")
print()

print("WAVEFORM STATISTICS:")
print(f"  Peak GR: {results['h_peak_GR']:.4e}")
print(f"  Peak UQFF: {results['h_peak_UQFF']:.4e}")
print(f"  Reduction: {results['percent_reduction']:.1f}%")
print()

print("COMPARISON TO DATA:")
print(f"  GW170817 observed: h ≈ {results['h_observed']:.0e}")
print(f"  UQFF prediction: {results['h_UQFF_predicted']:.2e}")
print(f"  GR residual: {results['residual_GR']*100:.0f}%")
print(f"  UQFF mismatch: {results['mismatch_UQFF']*100:.1f}%")
print(f"  GR fits better: {results['GR_fits_better']}")
print()

# Verification
expected_reduction = 50 <= results['percent_reduction'] <= 70
expected_GR_better = results['GR_fits_better'] == True
has_arrays = len(results['h_GR']) == results['n_samples']

print("VERIFICATION:")
print(f"  Reduction ~50-70%: {'PASS' if expected_reduction else 'FAIL'}")
print(f"  GR fits better: {'PASS' if expected_GR_better else 'FAIL'}")
print(f"  Arrays populated ({results['n_samples']} samples): {'PASS' if has_arrays else 'FAIL'}")
print()

if expected_reduction and expected_GR_better and has_arrays:
    print("TEST PASSED!")
else:
    print("CHECK NEEDED")

# Print the chart from steps
print("\n" + "="*60)
print("DERIVATION EXCERPT (Chart):")
print("="*60)
# Extract chart portion from steps
chart_start = steps.find("COMPARISON CHART")
if chart_start > 0:
    chart_section = steps[chart_start:chart_start+2500]
    print(chart_section[:1800])
