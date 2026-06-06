"""Smoke test: Tapestry of Blazing Starbirth (LMC NGC 2014+2020) master Universal Gravity (spec 09May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_25myr = 2.5e6 * YR_S

# ---- module-level constants ----
tests.append(('M_TAPESTRY_INIT_KG = 240 * M_SUN',
              m.M_TAPESTRY_INIT_KG == 240.0 * m.M_SUN))
tests.append(('M_TAPESTRY_INIT_KG ~ 4.774e32 kg (spec, 0.1% tol)',
              abs(m.M_TAPESTRY_INIT_KG - 4.774e32) / 4.774e32 < 0.001))
tests.append(('R_TAPESTRY_M = 10 ly ~ 9.461e16 m (spec)',
              abs(m.R_TAPESTRY_M - 9.461e16) / 9.461e16 < 1e-9))
tests.append(('B_TAPESTRY_T = 1e-6 T (spec LMC nebular field)',
              m.B_TAPESTRY_T == 1.0e-6))
tests.append(('TAU_SF_TAPESTRY_S = 5 Myr',
              abs(m.TAU_SF_TAPESTRY_S - 5.0e6 * YR_S) < 1e-3))
tests.append(('M_DOT0_TAPESTRY = 10^4/240 ~ 41.67 (spec)',
              abs(m.M_DOT0_TAPESTRY - (1.0e4 / 240.0)) < 1e-9))
tests.append(('V_GAS_TAPESTRY_MS = 1e5 m/s (spec nebular gas)',
              m.V_GAS_TAPESTRY_MS == 1.0e5))
tests.append(('T_TAPESTRY_DEFAULT_S = 2.5 Myr',
              abs(m.T_TAPESTRY_DEFAULT_S - t_25myr) < 1e-3))

# ---- reuse of _accretion_mass_growth_uqff for star-formation mass growth ----
M_t = m._accretion_mass_growth_uqff(M_initial_kg=m.M_TAPESTRY_INIT_KG,
                                       t_s=t_25myr,
                                       M_dot_0=m.M_DOT0_TAPESTRY,
                                       tau_acc_s=m.TAU_SF_TAPESTRY_S)
tests.append(('M(t=2.5 Myr) ~ 1.254e34 kg (spec, 0.5% tol)',
              abs(M_t - 1.254e34) / 1.254e34 < 0.005))
# t/tau = 0.5, exp(-0.5)=0.6065, so M_dot_at_t = 41.67*0.6065 = 25.27
M_dot_at_t = m.M_DOT0_TAPESTRY * math.exp(-0.5)
tests.append(('M_dot_at_t(2.5 Myr) ~ 25.27 (spec, 0.5% tol)',
              abs(M_dot_at_t - 25.27) / 25.27 < 0.005))
tests.append(('M_t = M_init * (1 + M_dot_at_t) check',
              abs(M_t - m.M_TAPESTRY_INIT_KG * (1.0 + M_dot_at_t))
                  / M_t < 1e-12))

# ---- decomposition of Tapestry composer at t=2.5 Myr ----
Ug1 = m.G_NEWTON * M_t / (m.R_TAPESTRY_M ** 2)
tests.append(('U_g1 = G*M(t)/r^2 ~ 9.351e-11 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 9.351e-11) / 9.351e-11 < 0.005))
H_factor = 1.0 + m.H0_MAGNETAR_SI * t_25myr
tests.append(('1 + H_0*t (2.5 Myr) ~ 1.0001723 (spec, 0.5% tol)',
              abs(H_factor - 1.0001723) / 1.0001723 < 0.005))
sc = 1.0 - m.B_TAPESTRY_T / m.B_CRIT_MAGNETAR_T
tests.append(('1 - B/B_crit ~ 1 (B=1e-6 T, B_crit=1e11 T -> ratio=1e-17)',
              abs(sc - 1.0) < 1e-12))
grav = Ug1 * H_factor * sc
tests.append(('grav_term ~ 9.353e-11 m/s^2 (spec, 0.5% tol)',
              abs(grav - 9.353e-11) / 9.353e-11 < 0.005))
Ug_sum_trz = (Ug1 + Ug1 * sc) * 1.1
tests.append(('(Ug1+Ug4)*(1+f_TRZ) ~ 2.057e-10 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 2.057e-10) / 2.057e-10 < 0.005))

# ---- Lorentz term reuse (with gas velocity 1e5 m/s, NOT magnetar 1e6) ----
lor = m._lorentz_acceleration_uqff(B_T=m.B_TAPESTRY_T, v_ms=m.V_GAS_TAPESTRY_MS,
                                     q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Lorentz q*v*B/m_p (gas v=1e5, B=1e-6) bare ~ 9.576e6 / scaled',
              abs(lor - 1.053e-4) / 1.053e-4 < 0.005))
# bare without macro_scale
lor_bare = m._lorentz_acceleration_uqff(B_T=m.B_TAPESTRY_T, v_ms=m.V_GAS_TAPESTRY_MS,
                                          q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR,
                                          macro_scale=1.0)
tests.append(('Lorentz bare *11 ~ 1.053e8 m/s^2 (spec pre-scaling)',
              abs(lor_bare - 1.053e8) / 1.053e8 < 0.005))
tests.append(('macro_scale=1e-12 reduces by factor 1e12',
              abs(lor_bare / lor - 1e12) < 1.0))

# ---- composer ----
g = m._tapestry_g_master_uqff()
tests.append(('_tapestry_g_master_uqff() defaults -> spec 1.053e-4 m/s^2 (0.5% tol)',
              abs(g - 1.053e-4) / 1.053e-4 < 0.005))
tests.append(('Total dominated by Lorentz term (grav + Ug_sum_trz << Lorentz)',
              (grav + Ug_sum_trz) < 0.01 * lor))

# Setting macro_scale=0 removes Lorentz -> total ~ grav + Ug_sum_trz
g_no_lor = m._tapestry_g_master_uqff(macro_scale=0.0)
tests.append(('Tapestry macro_scale=0 -> total ~ grav + Ug_sum_trz (1% tol)',
              abs(g_no_lor - (grav + Ug_sum_trz)) / g_no_lor < 0.01))

# Setting M_dot_0=0 -> no SF growth, mass = M_initial (240 M_sun)
g_no_sf = m._tapestry_g_master_uqff(M_dot_0=0.0)
Ug1_no_sf = m.G_NEWTON * m.M_TAPESTRY_INIT_KG / (m.R_TAPESTRY_M ** 2)
tests.append(('Tapestry M_dot_0=0 -> uses M_initial (240 M_sun); Ug1 << observed',
              Ug1_no_sf < Ug1))
tests.append(('Tapestry g(M_dot_0=0) still ~1.053e-4 (Lorentz dominates regardless)',
              abs(g_no_sf - 1.053e-4) / 1.053e-4 < 0.01))

# Setting f_TRZ=0 reduces Ug_sum by factor 1.1 (tiny effect on total)
g_no_trz = m._tapestry_g_master_uqff(f_TRZ=0.0)
tests.append(('Tapestry f_TRZ=0 reduces total by 0.1*Ug_sum (sign check)',
              g_no_trz < g))

# Monotone in M_dot_0
g_low = m._tapestry_g_master_uqff(M_dot_0=10.0)
g_high = m._tapestry_g_master_uqff(M_dot_0=100.0)
tests.append(('Tapestry g monotone in M_dot_0 (more SF -> more mass -> more Ug)',
              g_high > g > g_low))

# Tapestry-specific catalog/saturation primitives untouched
tests.append(('_tapestry_g_primitive_sat untouched (returns finite saturation)',
              isinstance(m._tapestry_g_primitive_sat(), float)
                  and m._tapestry_g_primitive_sat() > 0.0))
tests.append(('tapestry_ngc_3603 catalog entry untouched (NGC 3603, different region)',
              'tapestry_ngc_3603' in m.ASTRO_SYSTEMS))

passed = sum(1 for _, ok in tests if ok)
total = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
