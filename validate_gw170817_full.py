"""Validate GW170817 full 100s inspiral simulation."""
import sys
sys.path.insert(0, '.')
from CondensedPhysics import GravitationalWaveUQFFCalculator

# Create calculator with GW170817 chirp mass
calc = GravitationalWaveUQFFCalculator(
    M_chirp=1.188,  # Solar masses
    D_L=40.0,       # Mpc
    iota=0.0
)

print("=== GW170817 FULL 100s INSPIRAL SIMULATION ===\n")

# Run full inspiral from 23 Hz
results, steps = calc.simulate_GW170817_full_inspiral(
    t_duration=100.0,
    f_start=23.0,
    f_end=300.0,
    f_TRZ=0.1,
    B_NS=1e8,
    string_factor=0.37,
    beta_m=0.01,
    n_samples=1000,
    show_summary=True
)

print(f"Event: {results['event']}")
print(f"Simulation: {results['simulation_type']}")
print(f"Duration: {results['t_duration']:.0f}s in-band")
print(f"Frequency: {results['f_start']:.0f} Hz → {results['f_end']:.0f} Hz")
print(f"τ_chirp: {results['tau_chirp']:.0f}s")
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

print("SNR ANALYSIS:")
print(f"  SNR (GR): {results['SNR_GR_final']:.1f}")
print(f"  SNR (UQFF): {results['SNR_UQFF_final']:.1f}")
print(f"  SNR reduction: {results['SNR_reduction_percent']:.1f}%")
print(f"  Detectable (GR): {results['detectable_GR']}")
print(f"  Detectable (UQFF): {results['detectable_UQFF']}")
print()

print("PHASE ANALYSIS:")
print(f"  Total cycles: {int(results['phi_GW'][-1] / (2 * 3.14159))} cycles")
print(f"  Max phase lag: {results['max_phase_lag_rad']:.1f} rad ({results['max_phase_lag_cycles']:.1f} cycles)")
print()

print("FREQUENCY MILESTONES:")
print(f"  At 50 Hz: t = {results['t_at_50Hz']:.1f}s ({results['t_remaining_at_50Hz']:.1f}s to merger)")
print(f"  At 100 Hz: t = {results['t_at_100Hz']:.1f}s ({results['t_remaining_at_100Hz']:.1f}s to merger)")
print(f"  At 200 Hz: t = {results['t_at_200Hz']:.1f}s ({results['t_remaining_at_200Hz']:.1f}s to merger)")
print()

# Verification
expected_reduction = 50 <= results['percent_reduction'] <= 70
expected_snr_reduced = results['SNR_UQFF_final'] < results['SNR_GR_final']
has_arrays = len(results['h_GR']) == results['n_samples']
has_cumulative_snr = len(results['SNR_GR_cumulative']) == results['n_samples']
both_detectable = results['detectable_GR'] and results['detectable_UQFF']

print("VERIFICATION:")
print(f"  Reduction ~50-70%: {'PASS' if expected_reduction else 'FAIL'}")
print(f"  SNR reduced: {'PASS' if expected_snr_reduced else 'FAIL'}")
print(f"  Arrays populated ({results['n_samples']} samples): {'PASS' if has_arrays else 'FAIL'}")
print(f"  Cumulative SNR tracked: {'PASS' if has_cumulative_snr else 'FAIL'}")
print(f"  Both still detectable: {'PASS' if both_detectable else 'CHECK'}")
print()

if expected_reduction and expected_snr_reduced and has_arrays and has_cumulative_snr:
    print("TEST PASSED!")
else:
    print("CHECK NEEDED")
