"""Smoke test for quantum-variable bundle 8 (delta_def, f_TRZ, T_s, phi_hat_j)."""
import math, importlib, inspect
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Named constants
tests.append(('F_TRZ_DEFAULT = 0.1',                  abs(m.F_TRZ_DEFAULT - 0.1) < 1e-15))
tests.append(('DELTA_DEF_AMPLITUDE = 0.01',           abs(m.DELTA_DEF_AMPLITUDE - 0.01) < 1e-15))
tests.append(('DELTA_DEF_OMEGA_PER_DAY = 0.001',      abs(m.DELTA_DEF_OMEGA_PER_DAY - 0.001) < 1e-15))
tests.append(('T_S_SUN_K = 5778',                     abs(m.T_S_SUN_K - 5778.0) < 1e-10))
tests.append(('T_S_HOT_K = 10000',                    abs(m.T_S_HOT_K - 10000.0) < 1e-10))
tests.append(('PHI_HAT_J_DEFAULT = 1.0',              abs(m.PHI_HAT_J_DEFAULT - 1.0) < 1e-15))

# F_TRZ_DEFAULT consistent with existing TRZ primitive
tests.append(('F_TRZ_DEFAULT == TRZ (alias check)',   m.F_TRZ_DEFAULT == m.TRZ))

# Defect factor at t=0 (spec eq 2)
d0 = m._defect_factor_ug1(t_days=0.0)
tests.append(('delta_def(0) = 0.0', abs(d0 - 0.0) < 1e-15))

# Defect factor at t=pi/2 / omega = 1570.7963267949 days -> exact 0.01 (spec eq 4)
t_peak = (math.pi / 2.0) / m.DELTA_DEF_OMEGA_PER_DAY
d_peak = m._defect_factor_ug1(t_days=t_peak)
tests.append(('delta_def(pi/2 omega) = 0.01 (spec eq 4)', abs(d_peak - 0.01) < 1e-12))

# Defect factor at t=pi / omega = 3141.6 days -> back to 0
t_half = math.pi / m.DELTA_DEF_OMEGA_PER_DAY
d_half = m._defect_factor_ug1(t_days=t_half)
tests.append(('delta_def(pi/omega) returns to 0.0', abs(d_half) < 1e-12))

# Period 2 pi / omega ~ 6283.2 days ~ 17.22 yr (spec)
period_days = 2.0 * math.pi / m.DELTA_DEF_OMEGA_PER_DAY
period_yr   = period_days / 365.25
tests.append(('delta_def period ~ 17.22 yr', abs(period_yr - 17.20) < 0.05))

# delta_def(period) = 0 (closure)
d_full = m._defect_factor_ug1(t_days=period_days)
tests.append(('delta_def(period) = 0', abs(d_full) < 1e-9))

# Multiplicative form (1 + delta_def) at peak should be 1.01 (spec eq 4)
tests.append(('(1 + delta_def_peak) = 1.01', abs(1.0 + d_peak - 1.01) < 1e-12))

# phi_hat_j
tests.append(('phi_hat_j(j=1) = 1.0',  abs(m._disk_unit_vector_phi_hat(j=1)  - 1.0) < 1e-15))
tests.append(('phi_hat_j(j=12) = 1.0', abs(m._disk_unit_vector_phi_hat(j=12) - 1.0) < 1e-15))

# T_s thermal scaling: B_j(5778) = 1e3 T (reference)
bj_sun = m._b_j_thermal_scaled(T_s=m.T_S_SUN_K)
tests.append(('B_j(T_s=5778) = 1e3 T (reference)', abs(bj_sun - 1.0e3) < 1e-9))

# T_s thermal scaling: B_j(10000) = 1e3 * 10000/5778 = 1730.7 T
bj_hot = m._b_j_thermal_scaled(T_s=m.T_S_HOT_K)
expected_bj_hot = 1.0e3 * (10000.0 / 5778.0)
tests.append(('B_j(T_s=10000) = 1730.7 T', abs(bj_hot - expected_bj_hot) < 1e-9))

# Composed: U_g3 at thermally scaled B_j matches spec eq 10/11 ratio
ug3_sun = m._u_g3_stellar_planetary(B_j_sum=bj_sun, t=0.0)
ug3_hot = m._u_g3_stellar_planetary(B_j_sum=bj_hot, t=0.0)
tests.append(('U_g3 at T_s=5778 = 1.8e49 (spec eq 10)',
              abs(ug3_sun - 1.8e49) / 1.8e49 < 1e-9))
# spec eq 11 says ~3.11e49 -> 1.8e49 * (10000/5778) = 3.1153e49
tests.append(('U_g3 at T_s=10000 ~ 3.11e49 (spec eq 11)',
              abs(ug3_hot - 1.8e49 * (10000.0 / 5778.0)) / (1.8e49 * 10000.0 / 5778.0) < 1e-9))
# Round-trip check vs spec quoted value (~3.11e49)
tests.append(('U_g3 at T_s=10000 within 1% of spec 3.11e49',
              abs(ug3_hot - 3.11e49) / 3.11e49 < 0.01))

# f_TRZ cross-bundle: TRZ already wired into _l96_bearden_Ui_trz param (default 0.01 there);
# spec value 0.1 obtained when caller passes F_TRZ_DEFAULT
ui_with_trz_01 = m._l96_bearden_Ui_trz(lambda_i=1.0, rho_SCm_val=7.09e-37,
                                          rho_UA_val=7.09e-36, omega_s=2.5e-6,
                                          t_n=0.0, f_TRZ=m.F_TRZ_DEFAULT)
ui_no_trz      = m._l96_bearden_Ui_trz(lambda_i=1.0, rho_SCm_val=7.09e-37,
                                          rho_UA_val=7.09e-36, omega_s=2.5e-6,
                                          t_n=0.0, f_TRZ=0.0)
# Ratio (1 + 0.1) / 1.0 = 1.1 (spec eq 7)
tests.append(('U_i(f_TRZ=0.1) / U_i(0) = 1.1 (spec eq 7)',
              abs(ui_with_trz_01 / ui_no_trz - 1.1) < 1e-10))

# Regression: bundle 7 still loads
tests.append(('Bundle 7 _step_function still present',          hasattr(m, '_step_function')))
tests.append(('Bundle 7 _b_j_from_surface_field still present', hasattr(m, '_b_j_from_surface_field')))
tests.append(('Bundle 7 T_SMUNU_DIAGONAL still present',
              abs(m.T_SMUNU_DIAGONAL - 1.123e7) < 1e-3))
tests.append(('Bundle 7 OMEGA_S_SUN still present',
              abs(m.OMEGA_S_SUN - 2.5e-6) < 1e-15))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
