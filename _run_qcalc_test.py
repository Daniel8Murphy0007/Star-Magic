"""Quick smoke test for QCalc and QCalcGeom."""
import QCalc, QCalcGeom

# --- QCalc smoke test (solar-mass star at 1 AU) ---
solver = QCalc.UnifiedFieldSolver()
p = QCalc.ComputeParams(
    M=1.989e30, r=1.496e11, B=1e-4,
    T=5778, L=3.828e26, z=0, omega=2.87e-6
)
result = solver.solve(p)
print("=== QCalc.UnifiedFieldSolver.solve() ===")
print("  output keys:", list(result.keys()))
print("  equations computed:", len(result.get('long_form_equations', [])))
print("  solutions keys (first 5):", list(result.get('solutions', {}).keys())[:5])

# --- QCalcGeom smoke test ---
print()
print("=== QCalcGeom calculators ===")

# BSFGMetricCalculator
bsfg = QCalcGeom.BSFGMetricCalculator()
m = bsfg.compute({'r': 1.496e11, 't_n': 0.0})
print(f"  BSFGMetricCalculator: eps={m['eps']:.4e}  R_scalar={m['R_scalar']:.4e}")

# UniversalBuoyancyCalculator
ub = QCalcGeom.UniversalBuoyancyCalculator()
b = ub.compute({'r': 1.496e11, 't_n': 0.0})
print(f"  UniversalBuoyancyCalculator: FUBi={b['FUBi']:.4e}  FUBii={b['FUBii']:.4e}  zone={b['zone']}")

# HabitableZoneCalculator
hz = QCalcGeom.HabitableZoneCalculator()
h = hz.compute({'M': 1.989e30})
print(f"  HabitableZoneCalculator: r_hz={h['r_hz_AU']:.4f} AU  converged={h['converged']}  type={h['hz_type']}")

# UniversalGravityCalculator
ug = QCalcGeom.UniversalGravityCalculator()
u = ug.compute({'r': 1.496e11, 't_n': 0.0, 'M': 1.989e30})
print(f"  UniversalGravityCalculator: F_U_total={u['F_U_total']:.4e}  Ug_sum={u['Ug_sum']:.4e}")

print()
print("All smoke tests PASSED")
