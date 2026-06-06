"""Smoke test for quantum-variable bundle 2 (r_j, d_g, F_U, f_feedback, Omega_g)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Constants
tests.append(('R_J_MAGNETIC = 1.496e13 (100 AU)', abs(m.R_J_MAGNETIC - 1.496e13) < 1.0))
tests.append(('D_GALACTIC_CENTER = 2.55e20 m',    abs(m.D_GALACTIC_CENTER - 2.55e20) < 1e10))
tests.append(('OMEGA_GALACTIC = 7.3e-16 rad/s',   abs(m.OMEGA_GALACTIC - 7.3e-16) < 1e-30))
tests.append(('F_FEEDBACK_DEFAULT = 0.1',         abs(m.F_FEEDBACK_DEFAULT - 0.1) < 1e-15))

# Galactic period: T = 2 pi / Omega_g -> ~8.61e15 s ~ 2.73e8 yr
T_period = m._galactic_rotation_period()
tests.append(('T = 2 pi / Omega_g ~ 8.61e15 s', abs(T_period - 8.61e15) / 8.61e15 < 0.01))
T_years = T_period / (365.25 * 86400.0)
tests.append(('T ~ 2.73e8 years',               abs(T_years - 2.73e8) / 2.73e8 < 0.02))

# F_U composer: verify literal implementation of eq 11 with controlled inputs
# Set all U_gi = 0 except check additive structure of magnetic + aether terms separately
F_zero = m._f_u_unified(t=0.0, t_n=0.0, U_g_list=[0.0, 0.0, 0.0, 0.0],
                        mu_over_r_sum=0.0, lambda_U_i_sum=0.0)
trace_only = -2.0 + 4.0 * (m.ETA_AETHER_COUPLING * 1.123e7)
tests.append(('F_U zero-input reduces to aether trace', abs(F_zero - trace_only) < 1e-6))

# F_U with zeroed inputs except magnetic term: mu/r * (1 - exp(...))
F_mag = m._f_u_unified(t=1.0, t_n=0.0, U_g_list=[0.0, 0.0, 0.0, 0.0],
                       mu_over_r_sum=2.26e10, gamma=0.00005)
expected_osc = 2.26e10 * (1.0 - math.exp(-0.00005 * 1.0))
# Plus Aether trace (Minkowski 1-1-1-1 = -2 + 4*perturbation)
trace = -2.0 + 4.0 * (m.ETA_AETHER_COUPLING * 1.123e7)
tests.append(('F_U magnetic-only branch matches', abs(F_mag - (expected_osc + trace)) < 1e-6))

# F_U with t=0 -> osc factor = 1 - exp(0) = 0 (magnetic term vanishes)
F_t0 = m._f_u_unified(t=0.0, t_n=0.0, U_g_list=[0.0, 0.0, 0.0, 0.0],
                     mu_over_r_sum=2.26e10)
tests.append(('F_U magnetic vanishes at t=0', abs(F_t0 - trace) < 1e-6))

# Regression: prior bundle still working
u_b1 = m._u_bi_solar_wind_modulated(U_gi=1.39e26, Omega_g=7.3e-16,
                                     M_bh=8.15e36, d_g=2.55e20,
                                     U_UA=1.0, t_n=0.0)
tests.append(('Regression U_bi spec ~ -1.94e27', abs(u_b1 - (-1.94e27)) / 1.94e27 < 0.01))

# Regression: existing _l96_agn_feedback_Ug4 default f_feedback=0.1 matches new constant
import inspect
sig = inspect.signature(m._l96_agn_feedback_Ug4)
tests.append(('_l96_agn_feedback_Ug4 default f_feedback == F_FEEDBACK_DEFAULT',
              sig.parameters['f_feedback'].default == m.F_FEEDBACK_DEFAULT))

# Regression: existing _l96_final_parsec_Ug4 default d_g matches new constant
sig2 = inspect.signature(m._l96_final_parsec_Ug4)
tests.append(('_l96_final_parsec_Ug4 default d_g == D_GALACTIC_CENTER',
              sig2.parameters['d_g'].default == m.D_GALACTIC_CENTER))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
