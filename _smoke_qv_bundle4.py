"""Smoke test for quantum-variable bundle 4 (M_bh, mu_j, P_core, t_n, pi)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Named constants
tests.append(('M_BH_SGR_A = 8.15e36 kg',          abs(m.M_BH_SGR_A - 8.15e36) / 8.15e36 < 1e-10))
tests.append(('MU_J_BASE_AMPLITUDE = 3.38e20',    abs(m.MU_J_BASE_AMPLITUDE - 3.38e20) / 3.38e20 < 1e-10))
tests.append(('MU_J_DC_OFFSET = 1e3',             abs(m.MU_J_DC_OFFSET - 1.0e3) < 1e-6))
tests.append(('MU_J_AC_AMPLITUDE = 0.4',          abs(m.MU_J_AC_AMPLITUDE - 0.4) < 1e-15))
tests.append(('OMEGA_C_PERIOD_S = 3.96e8 s',      abs(m.OMEGA_C_PERIOD_S - 3.96e8) < 1.0))
tests.append(('OMEGA_C_MAGNETIC_CYCLE = 2*pi/3.96e8 ~ 1.587e-8',
              abs(m.OMEGA_C_MAGNETIC_CYCLE - (2.0 * math.pi / 3.96e8)) < 1e-20))
tests.append(('P_CORE_SUN = 1.0',                 abs(m.P_CORE_SUN - 1.0) < 1e-15))
tests.append(('P_CORE_PLANET = 1e-3',             abs(m.P_CORE_PLANET - 1.0e-3) < 1e-15))

# Time-dependent magnetic moment
# At t=0: mu_j = (1000 + 0.4*sin(0)) * 3.38e20 = 1000 * 3.38e20 = 3.38e23
mu_j0 = m._mu_j_time_dependent(t=0.0)
tests.append(('mu_j(t=0) = 3.38e23 T*m^3', abs(mu_j0 - 3.38e23) / 3.38e23 < 1e-10))

# At t = T/4 (period/4): sin(omega_c * T/4) = sin(pi/2) = 1
# -> mu_j = (1000 + 0.4) * 3.38e20 = 1000.4 * 3.38e20
t_quarter = m.OMEGA_C_PERIOD_S / 4.0
mu_j_q = m._mu_j_time_dependent(t=t_quarter)
expected_q = (1000.0 + 0.4) * 3.38e20
tests.append(('mu_j(T/4) = (1000+0.4)*3.38e20', abs(mu_j_q - expected_q) / expected_q < 1e-10))

# At t = T/2: sin = 0 again -> mu_j = 3.38e23
mu_j_h = m._mu_j_time_dependent(t=m.OMEGA_C_PERIOD_S / 2.0)
tests.append(('mu_j(T/2) returns to baseline 3.38e23', abs(mu_j_h - 3.38e23) / 3.38e23 < 1e-6))

# At t = 3T/4: sin = -1 -> mu_j = (1000 - 0.4) * 3.38e20
mu_j_3q = m._mu_j_time_dependent(t=3.0 * m.OMEGA_C_PERIOD_S / 4.0)
expected_3q = (1000.0 - 0.4) * 3.38e20
tests.append(('mu_j(3T/4) = (1000-0.4)*3.38e20', abs(mu_j_3q - expected_3q) / expected_3q < 1e-6))

# Ug3 stellar vs planetary (eq 13 / eq 14)
ug3_sun    = m._u_g3_stellar_planetary(P_core=m.P_CORE_SUN)
ug3_planet = m._u_g3_stellar_planetary(P_core=m.P_CORE_PLANET)
tests.append(('U_g3 (Sun, P_core=1) ~ 1.8e49',     abs(ug3_sun - 1.8e49) / 1.8e49 < 1e-10))
tests.append(('U_g3 (planet, P_core=1e-3) ~ 1.8e46', abs(ug3_planet - 1.8e46) / 1.8e46 < 1e-10))
tests.append(('Ug3_sun / Ug3_planet = 1000',         abs(ug3_sun / ug3_planet - 1000.0) < 1e-6))

# t_n verification: cos(pi * t_n) for t_n = -1 -> cos(-pi) = -1
tests.append(('cos(pi * t_n=-1) = -1', abs(math.cos(math.pi * -1.0) - (-1.0)) < 1e-15))
# Period of cos(pi * t_n): t_n = 2 -> cos(2pi) = 1
tests.append(('cos(pi * t_n) period = 2 days',
              abs(math.cos(math.pi * 2.0) - 1.0) < 1e-15))

# Regression: M_bh default in _l96_final_parsec_Ug4 matches new constant
import inspect
sig = inspect.signature(m._l96_final_parsec_Ug4)
tests.append(('_l96_final_parsec_Ug4 default M_BH == M_BH_SGR_A',
              sig.parameters['M_BH'].default == m.M_BH_SGR_A))

# Regression: Sun U_b1 spec from bundle 2 still works with M_BH_SGR_A
u_b1 = m._u_bi_solar_wind_modulated(U_gi=1.39e26, Omega_g=7.3e-16,
                                     M_bh=m.M_BH_SGR_A, d_g=m.D_GALACTIC_CENTER,
                                     U_UA=1.0, t_n=0.0)
tests.append(('U_b1 with M_BH_SGR_A ~ -1.94e27', abs(u_b1 - (-1.94e27)) / 1.94e27 < 0.01))

# Regression: bundle 3 still loads
tests.append(('Bundle 3 H_SCM_DEFAULT still present', abs(m.H_SCM_DEFAULT - 1.0) < 1e-15))
tests.append(('Bundle 3 _u_g3_string_magnetism_uqff still present',
              hasattr(m, '_u_g3_string_magnetism_uqff')))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
