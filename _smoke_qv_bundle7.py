"""Smoke test for quantum-variable bundle 7 (S, T_smunu, M_s, omega_s, B_s)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Named constants
tests.append(('T_SMUNU_DIAGONAL = 1.123e7',         abs(m.T_SMUNU_DIAGONAL - 1.123e7) < 1e-3))
tests.append(('OMEGA_S_SUN = 2.5e-6 rad/s',         abs(m.OMEGA_S_SUN - 2.5e-6) < 1e-15))
tests.append(('B_S_MIN_SUN = 1e-4 T',               abs(m.B_S_MIN_SUN - 1e-4) < 1e-15))
tests.append(('B_S_MAX_SUN = 0.4 T',                abs(m.B_S_MAX_SUN - 0.4) < 1e-15))
tests.append(('B_J_PER_B_S = 2500.0 T/T',           abs(m.B_J_PER_B_S - 2500.0) < 1e-9))

# M_s pre-existing
tests.append(('M_SUN = 1.989e30 kg', abs(m.M_SUN - 1.989e30) / 1.989e30 < 1e-10))

# Step function (spec eq 2)
tests.append(('S(r=2*R_b - R_b) = 1', abs(m._step_function(r=2*m.R_B_FIELD_BUBBLE) - 1.0) < 1e-15))
tests.append(('S(r=0.5*R_b - R_b) = 0', abs(m._step_function(r=0.5*m.R_B_FIELD_BUBBLE) - 0.0) < 1e-15))
tests.append(('S(r=R_b) = 1 (boundary)', abs(m._step_function(r=m.R_B_FIELD_BUBBLE) - 1.0) < 1e-15))
tests.append(('S(r=R_b - eps) = 0',
              abs(m._step_function(r=m.R_B_FIELD_BUBBLE * (1.0 - 1e-9)) - 0.0) < 1e-15))

# B_j(B_s) linear scaling (spec eq 16/17)
bj_max = m._b_j_from_surface_field(B_s=m.B_S_MAX_SUN)
tests.append(('B_j(B_s=0.4 T) = 1e3 T (eq 16)', abs(bj_max - 1.0e3) < 1e-9))
bj_min = m._b_j_from_surface_field(B_s=m.B_S_MIN_SUN)
tests.append(('B_j(B_s=1e-4 T) = 0.25 T (eq 17)', abs(bj_min - 0.25) < 1e-12))

# Composed: U_g3 from B_j(B_s) for both endpoints (eq 16 ~ 1.8e49, eq 17 ~ 4.5e45)
ug3_max = m._u_g3_stellar_planetary(B_j_sum=bj_max, t=0.0)
tests.append(('U_g3 at B_s=0.4 T = 1.8e49 J/m^3 (eq 16)',
              abs(ug3_max - 1.8e49) / 1.8e49 < 1e-9))
ug3_min = m._u_g3_stellar_planetary(B_j_sum=bj_min, t=0.0)
tests.append(('U_g3 at B_s=1e-4 T = 4.5e45 J/m^3 (eq 17)',
              abs(ug3_min - 4.5e45) / 4.5e45 < 1e-9))

# Linearity sanity: ratio U_g3(B_s_max) / U_g3(B_s_min) = (B_s_max/B_s_min) = 4000
tests.append(('U_g3 ratio matches B_s ratio (4000)',
              abs(ug3_max / ug3_min - (m.B_S_MAX_SUN / m.B_S_MIN_SUN)) < 1e-6))

# T_smunu wired into Aether metric (bundle 1 _aether_metric_a_munu uses 1.123e7 by default)
import inspect
sig_aether = inspect.signature(m._aether_metric_a_munu)
tests.append(('_aether_metric_a_munu default T_smunu matches T_SMUNU_DIAGONAL',
              sig_aether.parameters['T_smunu'].default == m.T_SMUNU_DIAGONAL))

# omega_s naming consistency with bundle 4 _u_g3_stellar_planetary default
sig_ug3 = inspect.signature(m._u_g3_stellar_planetary)
tests.append(('_u_g3_stellar_planetary default omega_s == OMEGA_S_SUN',
              sig_ug3.parameters['omega_s'].default == m.OMEGA_S_SUN))

# omega_s naming consistency with bundle 3 _u_g3_string_magnetism_uqff default
sig_ug3sm = inspect.signature(m._u_g3_string_magnetism_uqff)
tests.append(('_u_g3_string_magnetism_uqff default omega_s == OMEGA_S_SUN',
              sig_ug3sm.parameters['omega_s'].default == m.OMEGA_S_SUN))

# Regressions: bundle 6 still loads
tests.append(('Bundle 6 DELTA_SW_DEFAULT still present',
              abs(m.DELTA_SW_DEFAULT - 0.01) < 1e-15))
tests.append(('Bundle 6 _b_j_time_dependent still present',
              hasattr(m, '_b_j_time_dependent')))
tests.append(('Bundle 6 _solar_wind_modulation_factor still present',
              hasattr(m, '_solar_wind_modulation_factor')))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
