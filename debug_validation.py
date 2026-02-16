"""Debug validation failures"""
from QCalc import UnifiedFieldSolver
import numpy as np

solver = UnifiedFieldSolver()

print("=" * 80)
print("DEBUGGING UQFF_Compressed")
print("=" * 80)

M_sun = solver.C['M_sun']
AU = 1.496e11

# Test 1
results = solver.uqff_compressed_calc.compute_compressed_gravity(M=M_sun, r=AU, t=0.0)
g_solar = results[0].result
print(f"Test 1 - Solar system: g_solar = {g_solar:.6e} (expected 4e-3 to 8e-3)")
print(f"  Pass: {4e-3 < g_solar < 8e-3}")

# Test 2
M_galaxy = 1e12 * M_sun
r_kpc = 10 * 3.086e19
results = solver.uqff_compressed_calc.compute_compressed_gravity(M=M_galaxy, r=r_kpc, t=0.0)
g_galaxy = results[0].result
print(f"\nTest 2 - Galactic: g_galaxy = {g_galaxy:.6e} (expected 1e-11 to 1e-9)")
print(f"  Pass: {1e-11 < g_galaxy < 1e-9}")

# Test 3
t_universe = 1e10
results = solver.uqff_compressed_calc.compute_compressed_gravity(M=M_sun, r=AU, t=t_universe, R=AU, B=1e-6)
g_expanded = results[0].result
diff = abs(g_expanded - g_solar)
print(f"\nTest 3 - Time-varying: g_expanded = {g_expanded:.6e}, g_solar = {g_solar:.6e}")
print(f"  Difference: {diff:.6e} (threshold: 1e-12)")
print(f"  Pass: {diff >= 1e-12}")

print("\n" + "=" * 80)
print("DEBUGGING UQFF_MasterBuoyant")
print("=" * 80)

M_galaxy = 1e12 * solver.C['M_sun']
r_cosmic = 10 * 3.086e19
M_bh = 4.3e6 * solver.C['M_sun']
d_g = r_cosmic
Omega_g = 1e-15

# Test 1
results = solver.uqff_master_buoyant_calc.compute_master_buoyant_force(
    M=M_galaxy, r=r_cosmic, M_bh=M_bh, d_g=d_g, Omega_g=Omega_g
)
F = results[0].result
print(f"Test 1 - Negative force: F = {F:.6e}")
print(f"  Pass: {F < 0}")

# Test 2
results_breathing = solver.uqff_master_buoyant_calc.compute_master_buoyant_force(
    M=M_galaxy, r=r_cosmic, M_bh=M_bh, d_g=d_g, Omega_g=Omega_g,
    t=1e9, R=r_cosmic
)
F_breathing = results_breathing[0].result
diff = abs(F_breathing - F)
threshold = abs(F) * 1e-6
print(f"\nTest 2 - Volume breathing: F = {F:.6e}, F_breathing = {F_breathing:.6e}")
print(f"  Difference: {diff:.6e} (threshold: {threshold:.6e})")
print(f"  Pass: {diff >= threshold}")

# Test 3
print(f"\nTest 3 - Magnitude: abs(F) = {abs(F):.6e}")
print(f"  Pass: {1e10 < abs(F) < 1e50}")

print("\n" + "=" * 80)
print("DEBUGGING UQFF_Resonant")
print("=" * 80)

M_sun = solver.C['M_sun']
AU = 1.496e11
omega1 = 2 * np.pi / (30 * 86400)
omega2 = omega1 * 0.95

# Test 1
results = solver.uqff_resonant_calc.compute_resonant_gravity(M=M_sun, r=AU, omega1=omega1, omega2=omega2, t=0.0)
g_static = results[0].result
print(f"Test 1 - Static: g_static = {g_static:.6e} (expected 1e-15 to 1e5)")
print(f"  Pass: {1e-15 < g_static < 1e5}")

# Test 2
results_t = solver.uqff_resonant_calc.compute_resonant_gravity(M=M_sun, r=AU, omega1=omega1, omega2=omega2, t=1e7, R=AU)
g_varying = results_t[0].result
relative_diff = abs(g_varying - g_static) / abs(g_static) if g_static != 0 else 0
print(f"\nTest 2 - Time-varying: g_varying = {g_varying:.6e}, g_static = {g_static:.6e}")
print(f"  Relative diff: {relative_diff:.6e} (threshold: 1e-10)")
print(f"  Pass: {relative_diff >= 1e-10}")

# Test 3
results_full = solver.uqff_resonant_calc.compute_resonant_gravity(
    M=M_sun, r=AU, omega1=omega1, omega2=omega2, 
    t=1e6, t_n=-1e5, Delta_t=1e-43, R=AU
)
g_full = results_full[0].result
relative_diff_full = abs(g_full - g_static) / abs(g_static) if g_static != 0 else 0
print(f"\nTest 3 - Full integration: g_full = {g_full:.6e}")
print(f"  Relative diff: {relative_diff_full:.6e} (threshold: 1e-10)")
print(f"  Pass: {relative_diff_full >= 1e-10}")
