"""Smoke test: Westerlund 2 super-cluster (Carina, MW) master Universal Gravity (spec 09May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_1myr = 1.0e6 * YR_S

# ---- module-level constants ----
tests.append(('M_WESTERLUND2_INIT_KG = 30000 * M_SUN',
              m.M_WESTERLUND2_INIT_KG == 3.0e4 * m.M_SUN))
tests.append(('M_WESTERLUND2_INIT_KG ~ 5.967e34 kg (spec, 0.1% tol)',
              abs(m.M_WESTERLUND2_INIT_KG - 5.967e34) / 5.967e34 < 0.001))
tests.append(('R_WESTERLUND2_M = 9.461e16 m (10 ly)',
              m.R_WESTERLUND2_M == 9.461e16))
tests.append(('B_WESTERLUND2_T = 1e-5 T (spec, 10x Tapestry)',
              m.B_WESTERLUND2_T == 1.0e-5))
tests.append(('TAU_SF_WESTERLUND2_S = 2 Myr',
              abs(m.TAU_SF_WESTERLUND2_S - 2.0e6 * YR_S) < 1e-3))
tests.append(('M_DOT0_WESTERLUND2 = 10^5/30000 ~ 3.333 (spec)',
              abs(m.M_DOT0_WESTERLUND2 - (1.0e5 / 3.0e4)) < 1e-9))
tests.append(('V_GAS_WESTERLUND2_MS = 1e5 m/s',
              m.V_GAS_WESTERLUND2_MS == 1.0e5))
tests.append(('T_WESTERLUND2_DEFAULT_S = 1 Myr',
              abs(m.T_WESTERLUND2_DEFAULT_S - t_1myr) < 1e-3))

# ---- reuse of _accretion_mass_growth_uqff for SF mass growth ----
M_t = m._accretion_mass_growth_uqff(M_initial_kg=m.M_WESTERLUND2_INIT_KG,
                                       t_s=t_1myr,
                                       M_dot_0=m.M_DOT0_WESTERLUND2,
                                       tau_acc_s=m.TAU_SF_WESTERLUND2_S)
tests.append(('M(t=1 Myr) ~ 1.803e35 kg (spec, 0.5% tol)',
              abs(M_t - 1.803e35) / 1.803e35 < 0.005))
M_dot_at_t = m.M_DOT0_WESTERLUND2 * math.exp(-0.5)
tests.append(('M_dot_at_t(1 Myr) ~ 2.021 (spec, 0.5% tol)',
              abs(M_dot_at_t - 2.021) / 2.021 < 0.005))

# ---- decomposition at t=1 Myr ----
Ug1 = m.G_NEWTON * M_t / (m.R_WESTERLUND2_M ** 2)
tests.append(('U_g1 = G*M(t)/r^2 ~ 1.344e-9 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 1.344e-9) / 1.344e-9 < 0.005))
H_factor = 1.0 + m.H0_MAGNETAR_SI * t_1myr
tests.append(('1 + H_0*t (1 Myr) ~ 1.00006893 (spec, 0.5% tol)',
              abs(H_factor - 1.00006893) / 1.00006893 < 0.005))
sc = 1.0 - m.B_WESTERLUND2_T / m.B_CRIT_MAGNETAR_T
tests.append(('1 - B/B_crit ~ 1 (B=1e-5 T, B_crit=1e11 T -> ratio=1e-16)',
              abs(sc - 1.0) < 1e-12))
grav = Ug1 * H_factor * sc
tests.append(('grav_term ~ 1.344e-9 m/s^2 (spec, 0.5% tol)',
              abs(grav - 1.344e-9) / 1.344e-9 < 0.005))
Ug_sum_trz = (Ug1 + Ug1 * sc) * 1.1
tests.append(('(Ug1+Ug4)*(1+f_TRZ) ~ 2.957e-9 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 2.957e-9) / 2.957e-9 < 0.005))

# ---- Lorentz term reuse (v=1e5, B=1e-5 -> 10x Tapestry Lorentz) ----
lor = m._lorentz_acceleration_uqff(B_T=m.B_WESTERLUND2_T, v_ms=m.V_GAS_WESTERLUND2_MS,
                                     q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Lorentz q*v*B/m_p (v=1e5, B=1e-5) * 11 * 1e-12 ~ 1.053e-3 (spec, 0.5% tol)',
              abs(lor - 1.053e-3) / 1.053e-3 < 0.005))
lor_bare = m._lorentz_acceleration_uqff(B_T=m.B_WESTERLUND2_T, v_ms=m.V_GAS_WESTERLUND2_MS,
                                          q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR,
                                          macro_scale=1.0)
tests.append(('Lorentz bare * 11 ~ 1.053e9 m/s^2 (spec pre-scaling)',
              abs(lor_bare - 1.053e9) / 1.053e9 < 0.005))

# ---- composer ----
g = m._westerlund2_g_master_uqff()
tests.append(('_westerlund2_g_master_uqff() defaults -> spec 1.053e-3 m/s^2 (0.5% tol)',
              abs(g - 1.053e-3) / 1.053e-3 < 0.005))
tests.append(('Total dominated by Lorentz term',
              (grav + Ug_sum_trz) < 0.01 * lor))

# Westerlund Lorentz is 10x Tapestry Lorentz (since B is 10x)
lor_tap = m._lorentz_acceleration_uqff(B_T=m.B_TAPESTRY_T, v_ms=m.V_GAS_TAPESTRY_MS,
                                         q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Westerlund Lorentz = 10 * Tapestry Lorentz (B 10x scaling)',
              abs(lor / lor_tap - 10.0) < 0.01))

# Setting macro_scale=0 -> total ~ grav + Ug_sum_trz
g_no_lor = m._westerlund2_g_master_uqff(macro_scale=0.0)
tests.append(('Westerlund macro_scale=0 -> total ~ grav + Ug_sum_trz (1% tol)',
              abs(g_no_lor - (grav + Ug_sum_trz)) / g_no_lor < 0.01))

# M_dot_0=0 -> M = M_initial (30000 M_sun)
g_no_sf = m._westerlund2_g_master_uqff(M_dot_0=0.0)
tests.append(('Westerlund g(M_dot_0=0) still dominated by Lorentz (~ 1.053e-3)',
              abs(g_no_sf - 1.053e-3) / 1.053e-3 < 0.01))

# Monotone in M_dot_0
g_low = m._westerlund2_g_master_uqff(M_dot_0=1.0)
g_high = m._westerlund2_g_master_uqff(M_dot_0=10.0)
tests.append(('Westerlund g monotone in M_dot_0',
              g_high > g > g_low))

# f_TRZ=0 reduces Ug_sum
g_no_trz = m._westerlund2_g_master_uqff(f_TRZ=0.0)
tests.append(('Westerlund f_TRZ=0 reduces total (sign check)',
              g_no_trz < g))

# Existing primitives untouched
tests.append(('_westerlund_g_primitive_sat untouched (returns finite saturation)',
              isinstance(m._westerlund_g_primitive_sat(), float)
                  and m._westerlund_g_primitive_sat() > 0.0))
tests.append(('westerlund_2 catalog entry untouched',
              'westerlund_2' in m.ASTRO_SYSTEMS))

# Tapestry composer untouched
g_tap = m._tapestry_g_master_uqff()
tests.append(('Tapestry composer untouched (still ~1.053e-4)',
              abs(g_tap - 1.053e-4) / 1.053e-4 < 0.005))

passed = sum(1 for _, ok in tests if ok)
total = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
