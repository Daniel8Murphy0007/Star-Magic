#!/usr/bin/env python3
"""Test v_SCm velocity integration in SCM_VACUUM_MODEL."""

from CondensedPhysics import SCM_VACUUM_MODEL

print('=' * 70)
print('v_SCm VELOCITY INTEGRATION VALIDATION')
print('=' * 70)

# Test 1: get_v_SCm
v_SCm, explanation = SCM_VACUUM_MODEL.get_v_SCm()
print(f'\nv_SCm = {v_SCm:.0e} m/s')

# Test 2: compare to light speed
comparisons, _ = SCM_VACUUM_MODEL.compare_v_SCm_to_light_speed()
print(f'v_SCm / c = {comparisons["v_SCm_over_c"]:.4f} ≈ 1/3')
print(f'v_SCm² = {comparisons["v_SCm_squared"]:.0e} m²/s²')
print(f'Lorentz γ = {comparisons["lorentz_gamma"]:.3f}')

# Test 3: E_react calculation
E_react, _ = SCM_VACUUM_MODEL.compute_E_react_from_v_SCm(t=0)
print(f'\nE_react(t=0) = {E_react:.2f}')

# Test 4: phenomena contributions
contributions, _ = SCM_VACUUM_MODEL.compute_phenomena_v_SCm_contributions()
print(f'\n1 AU light travel: {contributions["AU_light_travel_time_s"]/60:.1f} min')
print(f'1 AU SCm travel: {contributions["AU_SCm_travel_time_s"]/60:.1f} min')
print(f'SCm/light delay ratio: {contributions["SCm_to_light_delay_ratio"]:.1f}x')

# Test 5: validate v_SCm physics
print('\n' + '=' * 70)
print('v_SCm VALIDATION TESTS')
print('=' * 70)
results = SCM_VACUUM_MODEL.validate_v_SCm_physics()
for test, status, value in results['results']:
    print(f'  {test}: {status} ({value})')
print(f'\nTotal: {results["tests_passed"]}/{results["total_tests"]} tests PASS')
print(f'Pass rate: {results["pass_rate"]:.0f}%')
