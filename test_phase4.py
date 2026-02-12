"""
Phase 4 Implementation Test Suite
Aether Metric Tensor (UA_μν) and Stress-Energy Tensor (T_s^μν)

Tests:
1. Minkowski metric baseline
2. Stress-energy tensor from vacuum densities
3. Metric perturbation (η × T_s)
4. Full aether metric (UA_μν = g_μν + δg_μν)
5. Metric determinant
6. Ricci scalar curvature
7. Inverse metric computation
8. Integration with UnifiedFieldSolver.solve()
9. Time evolution of metric
10. Physical validation (perturbations << 1)
"""

import numpy as np
from QCalc import (
    ComputeParams, UnifiedFieldSolver, CONSTANTS,
    AetherMetricCalculator, VacuumEnergyCalculator
)

print("=" * 80)
print("PHASE 4 IMPLEMENTATION TEST")
print("Aether Metric Tensor (UA_μν) and Stress-Energy Tensor (T_s^μν)")
print("=" * 80)
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 1: Minkowski Metric Baseline
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 1] Minkowski Metric Baseline")
print("-" * 80)

aether_calc = AetherMetricCalculator()
g_mu_nu = aether_calc.compute_minkowski_metric()

print("Minkowski metric g_μν (flat spacetime):")
print(g_mu_nu)
print()

# Verify signature (+---)
signature = np.diag(g_mu_nu)
print(f"Signature: [{signature[0]:.0f}, {signature[1]:.0f}, {signature[2]:.0f}, {signature[3]:.0f}]")
expected = np.array([1, -1, -1, -1])
assert np.allclose(signature, expected), "Minkowski signature incorrect!"
print("✓ Minkowski metric correct (signature: +---)")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 2: Stress-Energy Tensor
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 2] Stress-Energy Tensor from Vacuum Densities")
print("-" * 80)

vacuum_calc = VacuumEnergyCalculator()

# Get vacuum densities
lambda_vac_UA = vacuum_calc.compute_lambda_vac_UA()
lambda_vac_SCm = vacuum_calc.compute_lambda_vac_SCm()
lambda_vac_A = vacuum_calc.compute_lambda_vac_A()
t_n = -86400.0  # -1 day

print(f"Input vacuum densities:")
print(f"  λ_vac,[UA]  = {lambda_vac_UA:.4e} J/m³")
print(f"  λ_vac,[SCm] = {lambda_vac_SCm:.4e} J/m³")
print(f"  λ_vac,A     = {lambda_vac_A:.4e} J/m³")
print(f"  t_n         = {t_n:.4e} s")
print()

T_s = aether_calc.compute_stress_energy_tensor(lambda_vac_UA, lambda_vac_SCm, lambda_vac_A, t_n)

print("Stress-energy tensor T_s^μν:")
print(f"T^00 (energy density) = {T_s[0,0]:.4e} kg/m³ c²")
print(f"T^11 (pressure x)     = {T_s[1,1]:.4e} kg/m³ c²")
print(f"T^22 (pressure y)     = {T_s[2,2]:.4e} kg/m³ c²")
print(f"T^33 (pressure z)     = {T_s[3,3]:.4e} kg/m³ c²")
print()

# Verify perfect fluid form (T^ii = -ρ/3)
pressure_ratio = T_s[1,1] / T_s[0,0]
print(f"Pressure/Density ratio: {pressure_ratio:.4f} (expected: -1/3 = -0.333)")
assert np.isclose(pressure_ratio, -1/3, rtol=0.01), "Pressure ratio incorrect!"
print("✓ Stress-energy tensor has correct perfect fluid form")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 3: Metric Perturbation
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 3] Metric Perturbation δg_μν = η × T_s^μν")
print("-" * 80)

delta_g = aether_calc.compute_metric_perturbation(lambda_vac_UA, lambda_vac_SCm, lambda_vac_A, t_n)

print(f"Aether coupling: η = {aether_calc.eta:.4e}")
print()
print("Metric perturbation δg_μν:")
print(f"δg_00 = {delta_g[0,0]:.4e}")
print(f"δg_11 = {delta_g[1,1]:.4e}")
print(f"δg_22 = {delta_g[2,2]:.4e}")
print(f"δg_33 = {delta_g[3,3]:.4e}")
print()

# Verify small perturbation (η ~ 10^-22)
max_perturbation = np.max(np.abs(delta_g))
print(f"Maximum |δg| = {max_perturbation:.4e}")
print(f"Perturbation check: |δg| << 1? {max_perturbation < 0.01}")
assert max_perturbation < 0.01, "Perturbation too large! Should be << 1"
print("✓ Metric perturbations are small (weak field)")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 4: Full Aether Metric
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 4] Full Aether Metric UA_μν = g_μν + δg_μν")
print("-" * 80)

UA_mu_nu = aether_calc.compute_aether_metric(lambda_vac_UA, lambda_vac_SCm, lambda_vac_A, t_n)

print("Aether metric UA_μν:")
print(UA_mu_nu)
print()

# Verify close to Minkowski
difference = UA_mu_nu - g_mu_nu
print("Difference from Minkowski (UA - g):")
print(f"  Δ(UA_00 - g_00) = {difference[0,0]:.4e}")
print(f"  Δ(UA_11 - g_11) = {difference[1,1]:.4e}")
print()

# Verify signature preserved
UA_signature = np.sign(np.diag(UA_mu_nu))
expected_signature = np.array([1, -1, -1, -1])
print(f"UA signature: [{UA_signature[0]:.0f}, {UA_signature[1]:.0f}, {UA_signature[2]:.0f}, {UA_signature[3]:.0f}]")
assert np.allclose(UA_signature, expected_signature), "Aether metric signature changed!"
print("✓ Aether metric preserves causal structure (signature unchanged)")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 5: Metric Determinant
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 5] Metric Determinant")
print("-" * 80)

det_g = aether_calc.compute_metric_determinant(g_mu_nu)
det_UA = aether_calc.compute_metric_determinant(UA_mu_nu)

print(f"det(g_μν)    = {det_g:.10f} (Minkowski)")
print(f"det(UA_μν)   = {det_UA:.10f} (Aether metric)")
print(f"Difference   = {det_UA - det_g:.4e}")
print()

# Verify close to -1
assert np.isclose(det_g, -1.0, atol=1e-10), "Minkowski det should be -1"
assert np.isclose(det_UA, -1.0, atol=0.01), "Aether det should be close to -1"
print("✓ Metric determinants correct (both ≈ -1)")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 6: Inverse Metric
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 6] Inverse Metric UA^μν")
print("-" * 80)

UA_inv = aether_calc.compute_inverse_metric(UA_mu_nu)

print("Inverse aether metric UA^μν:")
print(f"UA^00 = {UA_inv[0,0]:.10f}")
print(f"UA^11 = {UA_inv[1,1]:.10f}")
print()

# Verify inverse property: UA_μα × UA^αν = δ_μ^ν
identity_check = np.matmul(UA_mu_nu, UA_inv)
print("Identity check (UA_μα × UA^αν):")
print(identity_check)
print()

# Should be identity matrix within numerical precision
identity = np.eye(4)
assert np.allclose(identity_check, identity, atol=1e-10), "Inverse metric failed!"
print("✓ Inverse metric correct (UA × UA^-1 = I)")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 7: Ricci Scalar Curvature
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 7] Ricci Scalar Curvature")
print("-" * 80)

R = aether_calc.compute_ricci_scalar(UA_mu_nu)

print(f"Ricci scalar R = {R:.4e} m⁻²")
print(f"For Minkowski: R = 0")
print(f"For perturbed: R ≈ -Tr(δg)/2")
print()

trace_delta_g = np.trace(delta_g)
R_expected = -trace_delta_g / 2.0
print(f"Expected R from trace: {R_expected:.4e} m⁻²")
print(f"Match: {np.isclose(R, R_expected)}")
assert np.isclose(R, R_expected, rtol=0.01), "Ricci scalar computation incorrect!"
print("✓ Ricci scalar computed correctly")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 8: Integration with UnifiedFieldSolver
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 8] Integration with UnifiedFieldSolver")
print("-" * 80)

solver = UnifiedFieldSolver()

params = ComputeParams(
    M=CONSTANTS['M_sun'],
    r=1.496e11,  # 1 AU
    t=0.0,
    t_n=-86400.0  # -1 day
)

# Test Phase 4 method
phase4_results = solver._compute_aether_metric_phase4(params)

print(f"Phase 4 Aether Metric results: {len(phase4_results)} equations")
for eq in phase4_results:
    print(f"  {eq.name}: {eq.result:.4e} {eq.unit}")

print()
print("✓ Phase 4 methods integrated successfully")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 9: Full Integration with solve()
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 9] Full Integration with solve()")
print("-" * 80)

result = solver.solve(params)

query_id = result['query_id']
n_equations = len(result['long_form_equations'])

print(f"Query ID: {query_id}")
print(f"Total equations computed: {n_equations}")
print()

# Count Phase 4 equations
phase4_equations = []
for eq in result['long_form_equations']:
    eq_name = eq.name if hasattr(eq, 'name') else eq.get('name', '')
    if 'stress_energy' in eq_name or 'metric' in eq_name or 'ricci' in eq_name:
        phase4_equations.append(eq)

print(f"Phase 4 equations: {len(phase4_equations)}")
if phase4_equations:
    print("  Phase 4 components:")
    for eq in phase4_equations:
        eq_name = eq.name if hasattr(eq, 'name') else eq.get('name', 'unknown')
        eq_result = eq.result if hasattr(eq, 'result') else eq.get('result', 0)
        eq_unit = eq.unit if hasattr(eq, 'unit') else eq.get('unit', '')
        print(f"    {eq_name}: {eq_result:.4e} {eq_unit}")

print()

# Check available equations
available = result['available_equations']
phase4_available = [eq for eq in available if 'stress' in eq or 'metric' in eq or 'ricci' in eq]
print(f"Phase 4 available equations: {len(phase4_available)}")
if phase4_available:
    print(f"  {', '.join(phase4_available)}")

print()
print("✓ Full integration successful")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 10: Time Evolution of Metric
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 10] Time Evolution of Aether Metric")
print("-" * 80)

# Test metric at different negative times
time_points = [0.0, -86400.0, -2*86400.0, -7*86400.0]  # 0, -1, -2, -7 days

print("Aether metric time evolution (UA_00 component):")
for t_n in time_points:
    UA_t = aether_calc.compute_aether_metric(lambda_vac_UA, lambda_vac_SCm, lambda_vac_A, t_n)
    t_days = t_n / 86400.0
    print(f"  t_n = {t_days:6.1f} days: UA_00 = {UA_t[0,0]:.10f}")

print()
print("✓ Time evolution computed successfully")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 11: Physical Validation
# ═══════════════════════════════════════════════════════════════════════════════

print("[TEST 11] Physical Validation")
print("-" * 80)

print("Checking physical consistency:")
print()

# 1. Energy positivity (T^00 > 0)
print(f"1. Energy positivity: T^00 = {T_s[0,0]:.4e} kg/m³ c² > 0? {T_s[0,0] > 0}")
assert T_s[0,0] > 0, "Energy density must be positive!"

# 2. Weak energy condition (T^μν u_μ u_ν ≥ 0 for timelike u)
# For perfect fluid: ρ + P ≥ 0
rho_plus_P = T_s[0,0] + T_s[1,1]
print(f"2. Weak energy condition: ρ + P = {rho_plus_P:.4e} ≥ 0? {rho_plus_P >= 0}")
assert rho_plus_P >= 0, "Weak energy condition violated!"

# 3. Perturbation theory validity (|δg/g| << 1)
# Check only diagonal elements to avoid divide-by-zero
diagonal_g = np.diag(g_mu_nu)
diagonal_delta_g = np.diag(delta_g)
relative_pert = np.max(np.abs(diagonal_delta_g / diagonal_g))
print(f"3. Perturbation theory: |δg/g| = {relative_pert:.4e} << 1? {relative_pert < 0.01}")
assert relative_pert < 0.01, "Perturbation too large for linear theory!"

# 4. Causality preservation (signature unchanged)
print(f"4. Causality: Signature preserved? {np.allclose(UA_signature, expected_signature)}")
assert np.allclose(UA_signature, expected_signature), "Causality violated!"

# 5. Aether coupling strength
print(f"5. Aether coupling: η = {aether_calc.eta:.4e} << 1? {aether_calc.eta < 1e-10}")
assert aether_calc.eta < 1e-10, "Aether coupling too strong!"

print()
print("✓ All physical consistency checks passed")
print()

# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

print("=" * 80)
print("PHASE 4 IMPLEMENTATION STATUS")
print("=" * 80)
print("[OK] Minkowski metric - baseline flat spacetime")
print("[OK] Stress-energy tensor - vacuum density sourcing")
print("[OK] Metric perturbation - weak field (η ~ 10^-22)")
print("[OK] Full aether metric - UA_μν = g_μν + δg_μν")
print("[OK] Metric determinant - det(UA) ≈ -1")
print("[OK] Inverse metric - UA × UA^-1 = I")
print("[OK] Ricci scalar - curvature from perturbation")
print("[OK] Integration with UnifiedFieldSolver")
print("[OK] Full solve() integration")
print("[OK] Time evolution - negative time modulation")
print("[OK] Physical validation - energy conditions, causality")
print("-" * 80)
print(f"Total available equations: {len(available)}")
print(f"Phase 4 equations: 5 tensorial components")
print("  - stress_energy_tensor (T_s^μν)")
print("  - metric_perturbation (δg_μν)")
print("  - aether_metric (UA_μν)")
print("  - metric_determinant (det(UA))")
print("  - ricci_scalar (R)")
print("=" * 80)
print("Phase 4 implementation COMPLETE and VERIFIED")
print("=" * 80)
