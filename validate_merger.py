"""Validate UQFF BH Merger Dynamics implementation."""
import sys
sys.path.insert(0, '.')
from CondensedPhysics import BH_PHASES_MODEL

# GW150914: 36 + 29 solar masses
M_sun = 1.989e30
M1 = 36 * M_sun
M2 = 29 * M_sun

# Inspiral separation
a = 1e11  # 100 Gm

# Run merger dynamics with B_t = 10% of B_crit
results, steps = BH_PHASES_MODEL.compute_UQFF_BH_merger_dynamics(
    M1=M1, M2=M2, a=a, B_t=4.4e12, f_TRZ=0.1
)

print('=== UQFF BH MERGER DYNAMICS TEST ===')
print(f'M1: {results["M1_solar"]:.1f} M_sun')
print(f'M2: {results["M2_solar"]:.1f} M_sun')
print(f'M_tot: {results["M_tot_solar"]:.1f} M_sun')
print(f'q (mass ratio): {results["q"]:.4f}')
print()
print(f'P_GW (standard): {results["P_GW"]:.4e} W')
print(f'P_GW,UQFF: {results["P_GW_UQFF"]:.4e} W')
print(f'tau_merge: {results["tau_merge_years"]:.4e} years')
print(f'tau_UQFF: {results["tau_UQFF_years"]:.4e} years')
print(f'tau_ratio: {results["tau_ratio"]:.2f}x')
print()
print(f'E_GW (standard): {results["E_GW_standard_solar"]:.4f} M_sun')
print(f'E_GW,UQFF: {results["E_GW_UQFF_solar"]:.4f} M_sun')
print(f'Mass retention: {results["mass_retention"]:.2f}%')
print()
print(f'Aether damping: {results["aether_damping"]:.6f}')
print(f'B_factor: {results["B_factor"]:.4f}')
print(f'TRZ_factor: {results["TRZ_factor"]:.4f}')
print(f'String binding: {results["string_binding"]:.4f}')
print(f'Combined factor: {results["combined_factor"]:.6f}')
print()
print(f'Status: {results["stability_status"]}')
print()
print('TEST PASSED!' if results['tau_ratio'] > 1.1 else 'CHECK NEEDED')
