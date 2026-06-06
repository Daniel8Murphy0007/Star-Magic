"""Smoke test for quantum-variable bundle 6 (delta_sw, kappa, P_SCm, v_sw, omega_c)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Named constants
tests.append(('DELTA_SW_DEFAULT = 0.01',       abs(m.DELTA_SW_DEFAULT - 0.01) < 1e-15))
tests.append(('V_SW_DEFAULT = 5e5 m/s',        abs(m.V_SW_DEFAULT - 5.0e5) < 1.0))
tests.append(('P_SCM_SUN = 1.0',               abs(m.P_SCM_SUN - 1.0) < 1e-15))
tests.append(('P_SCM_PLANET = 1e-3',           abs(m.P_SCM_PLANET - 1.0e-3) < 1e-15))
tests.append(('B_J_DC_OFFSET_TESLA = 1e3 T',   abs(m.B_J_DC_OFFSET_TESLA - 1.0e3) < 1e-6))
tests.append(('B_J_AC_AMPLITUDE_TESLA = 0.4',  abs(m.B_J_AC_AMPLITUDE_TESLA - 0.4) < 1e-15))

# Solar wind modulation factor (spec eq 2)
sw = m._solar_wind_modulation_factor()
tests.append(('1 + delta_sw * v_sw = 5001', abs(sw - 5001.0) < 1e-9))
sw_zero = m._solar_wind_modulation_factor(delta_sw=0.0)
tests.append(('1 + 0 * v_sw = 1', abs(sw_zero - 1.0) < 1e-15))

# B_j(t=0) = 1e3 T (spec match eq 14)
bj0 = m._b_j_time_dependent(t=0.0)
tests.append(('B_j(0) = 1e3 T', abs(bj0 - 1.0e3) < 1e-10))

# B_j(t=T/4) = 1000.4 T
T = 3.96e8
bj_q = m._b_j_time_dependent(t=T/4.0)
tests.append(('B_j(T/4) = 1000.4 T', abs(bj_q - 1000.4) < 1e-10))

# B_j(t=T/2) = 1e3 (sin returns to 0)
bj_h = m._b_j_time_dependent(t=T/2.0)
tests.append(('B_j(T/2) returns to 1e3 T', abs(bj_h - 1.0e3) < 1e-6))

# B_j(t=3T/4) = 999.6 T
bj_3q = m._b_j_time_dependent(t=3.0*T/4.0)
tests.append(('B_j(3T/4) = 999.6 T', abs(bj_3q - 999.6) < 1e-6))

# Relationship: mu_j(t) / B_j(t) should equal MU_J_BASE_AMPLITUDE (3.38e20)
mu_j_q = m._mu_j_time_dependent(t=T/4.0)
ratio = mu_j_q / bj_q
tests.append(('mu_j(t)/B_j(t) = MU_J_BASE_AMPLITUDE',
              abs(ratio - m.MU_J_BASE_AMPLITUDE) / m.MU_J_BASE_AMPLITUDE < 1e-10))

# Cross-bundle: Ug2 heliosphere uses delta_sw/v_sw matching new defaults
import inspect
sig = inspect.signature(m._u_g2_heliosphere_uqff)
tests.append(('_u_g2_heliosphere_uqff default delta_sw == DELTA_SW_DEFAULT',
              sig.parameters['delta_sw'].default == m.DELTA_SW_DEFAULT))
tests.append(('_u_g2_heliosphere_uqff default v_sw == V_SW_DEFAULT',
              sig.parameters['v_sw'].default == m.V_SW_DEFAULT))

# Cross-bundle: P_SCm parameter in _l96_bearden_Um_trz can take P_SCM_SUN/PLANET
um_sun    = m._l96_bearden_Um_trz(t=1.0, t_n=0.0, mu_sum_over_r=2.26e10,
                                    gamma=0.00005, P_SCm=m.P_SCM_SUN,
                                    E_react=1.0e46, f_heaviside=0.01, f_quasi=0.01)
um_planet = m._l96_bearden_Um_trz(t=1.0, t_n=0.0, mu_sum_over_r=2.26e10,
                                    gamma=0.00005, P_SCm=m.P_SCM_PLANET,
                                    E_react=1.0e46, f_heaviside=0.01, f_quasi=0.01)
tests.append(('U_m(Sun)/U_m(planet) = 1000 (P_SCm ratio)',
              abs(um_sun / um_planet - 1000.0) < 1e-6))

# Already-captured regressions
tests.append(('kappa already captured (KAPPA_REACT_PER_DAY = 5e-4)',
              abs(m.KAPPA_REACT_PER_DAY - 5.0e-4) < 1e-15))
tests.append(('omega_c already captured (OMEGA_C_MAGNETIC_CYCLE present)',
              hasattr(m, 'OMEGA_C_MAGNETIC_CYCLE')
              and abs(m.OMEGA_C_MAGNETIC_CYCLE - 2.0*math.pi/3.96e8) < 1e-25))

# Bundle 5 still loads
tests.append(('Bundle 5 E_REACT_AMPLITUDE still present',
              abs(m.E_REACT_AMPLITUDE - 1.0e46) / 1.0e46 < 1e-10))
tests.append(('Bundle 5 _e_react_full still present', hasattr(m, '_e_react_full')))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
