"""Debug UQFF_Resonant validation"""
from QCalc import UnifiedFieldSolver
import numpy as np

solver = UnifiedFieldSolver()

print("=" * 80)
print("DEBUGGING UQFF_Resonant - DETAILED")
print("=" * 80)

M_sun = solver.C['M_sun']
AU = 1.496e11
omega1 = 2 * np.pi / (30 * 86400)
omega2 = omega1 * 0.95

print(f"Test parameters:")
print(f"  M_sun = {M_sun:.6e}")
print(f"  AU = {AU:.6e}")
print(f"  omega1 = {omega1:.6e}")
print(f"  omega2 = {omega2:.6e}")
print()

# Test 1
print("Test 1: Static case")
try:
    results = solver.uqff_resonant_calc.compute_resonant_gravity(M=M_sun, r=AU, omega1=omega1, omega2=omega2, t=0.0)
    g_static = results[0].result
    print(f"  g_static = {g_static:.6e}")
    print(f"  Is finite? {np.isfinite(g_static)}")
    print(f"  Pass: {np.isfinite(g_static)}")
except Exception as e:
    print(f"  EXCEPTION: {e}")
    g_static = None
print()

# Test 2
print("Test 2: Time-varying")
try:
    results_t = solver.uqff_resonant_calc.compute_resonant_gravity(M=M_sun, r=AU, omega1=omega1, omega2=omega2, t=1e7, R=AU)
    g_varying = results_t[0].result
    print(f"  g_varying = {g_varying:.6e}")
    print(f"  g_static = {g_static:.6e}" if g_static is not None else "  g_static = None")
    print(f"  Is finite? {np.isfinite(g_varying)}")
    print(f"  Different from static? {g_varying != g_static if g_static is not None else 'N/A'}")
    print(f"  Pass: {np.isfinite(g_varying) and (g_static is None or g_varying != g_static)}")
except Exception as e:
    print(f"  EXCEPTION: {e}")
    g_varying = None
print()

# Test 3
print("Test 3: Full foundational integrations")
try:
    results_full = solver.uqff_resonant_calc.compute_resonant_gravity(
        M=M_sun, r=AU, omega1=omega1, omega2=omega2, 
        t=1e6, t_n=-1e5, Delta_t=1e-43, R=AU
    )
    g_full = results_full[0].result
    print(f"  g_full = {g_full:.6e}")
    print(f"  Is finite? {np.isfinite(g_full)}")
    print(f"  Pass: {np.isfinite(g_full)}")
except Exception as e:
    print(f"  EXCEPTION: {e}")
    g_full = None
print()

# Overall validation
print("=" * 80)
print("OVERALL VALIDATION RESULT:")
test1_pass = g_static is not None and np.isfinite(g_static)
test2_pass = g_varying is not None and np.isfinite(g_varying) and (g_static is None or g_varying != g_static)
test3_pass = g_full is not None and np.isfinite(g_full)
overall = test1_pass and test2_pass and test3_pass
print(f"Test 1: {'PASS' if test1_pass else 'FAIL'}")
print(f"Test 2: {'PASS' if test2_pass else 'FAIL'}")
print(f"Test 3: {'PASS' if test3_pass else 'FAIL'}")
print(f"Overall: {'PASS' if overall else 'FAIL'}")
