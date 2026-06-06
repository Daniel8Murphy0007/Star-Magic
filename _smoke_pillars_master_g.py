"""Smoke test: Pillars of Creation (M16 Eagle Nebula) master Universal Gravity (spec 09May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_5e5yr = 5.0e5 * YR_S

# ---- module-level constants ----
tests.append(('M_PILLARS_INIT_KG = 10100 * M_SUN',
              m.M_PILLARS_INIT_KG == 10100.0 * m.M_SUN))
tests.append(('M_PILLARS_INIT_KG ~ 2.009e34 kg (spec, 0.1% tol)',
              abs(m.M_PILLARS_INIT_KG - 2.009e34) / 2.009e34 < 0.001))
tests.append(('R_PILLARS_M = 4.731e16 m (~5 ly)',
              m.R_PILLARS_M == 4.731e16))
tests.append(('B_PILLARS_T = 1e-6 T (spec)',
              m.B_PILLARS_T == 1.0e-6))
tests.append(('TAU_SF_PILLARS_S = 1 Myr',
              abs(m.TAU_SF_PILLARS_S - 1.0e6 * YR_S) < 1e-3))
tests.append(('M_DOT0_PILLARS = 10^4/10100 ~ 0.9901 (spec)',
              abs(m.M_DOT0_PILLARS - (1.0e4 / 10100.0)) < 1e-9))
tests.append(('E0_PILLARS = 0.1 (spec)',
              m.E0_PILLARS == 0.1))
tests.append(('TAU_ERODE_PILLARS_S = 1 Myr',
              abs(m.TAU_ERODE_PILLARS_S - 1.0e6 * YR_S) < 1e-3))
tests.append(('V_GAS_PILLARS_MS = 1e5',
              m.V_GAS_PILLARS_MS == 1.0e5))
tests.append(('T_PILLARS_DEFAULT_S = 5e5 yr',
              abs(m.T_PILLARS_DEFAULT_S - t_5e5yr) < 1e-3))

# ---- reuse: _magnetar_B_decay as generic exp decay for E(t) ----
E_t = m._magnetar_B_decay(t_5e5yr, B_0_T=m.E0_PILLARS, tau_B_s=m.TAU_ERODE_PILLARS_S)
tests.append(('E(t=5e5 yr) = 0.1*exp(-0.5) ~ 0.06065 (spec, 0.5% tol)',
              abs(E_t - 0.06065) / 0.06065 < 0.005))
erosion_factor = 1.0 - E_t
tests.append(('(1 - E(t)) ~ 0.93935 (spec, 0.5% tol)',
              abs(erosion_factor - 0.93935) / 0.93935 < 0.005))

# ---- reuse: _accretion_mass_growth_uqff for SF mass growth ----
M_t = m._accretion_mass_growth_uqff(M_initial_kg=m.M_PILLARS_INIT_KG,
                                       t_s=t_5e5yr,
                                       M_dot_0=m.M_DOT0_PILLARS,
                                       tau_acc_s=m.TAU_SF_PILLARS_S)
tests.append(('M(t=5e5 yr) ~ 3.215e34 kg (spec, 0.5% tol)',
              abs(M_t - 3.215e34) / 3.215e34 < 0.005))
M_dot_at_t = m.M_DOT0_PILLARS * math.exp(-0.5)
tests.append(('M_dot_at_t(5e5 yr) ~ 0.6005 (spec, 0.5% tol)',
              abs(M_dot_at_t - 0.6005) / 0.6005 < 0.005))

# ---- decomposition at t=5e5 yr ----
Ug1 = m.G_NEWTON * M_t / (m.R_PILLARS_M ** 2)
tests.append(('U_g1 = G*M(t)/r^2 ~ 9.588e-10 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 9.588e-10) / 9.588e-10 < 0.005))
H_factor = 1.0 + m.H0_MAGNETAR_SI * t_5e5yr
tests.append(('1 + H_0*t (5e5 yr) ~ 1.00003446 (spec, 0.5% tol)',
              abs(H_factor - 1.00003446) / 1.00003446 < 0.005))
sc = 1.0 - m.B_PILLARS_T / m.B_CRIT_MAGNETAR_T
tests.append(('1 - B/B_crit ~ 1',
              abs(sc - 1.0) < 1e-12))
grav = Ug1 * H_factor * sc * erosion_factor
tests.append(('grav_term ~ 9.005e-10 m/s^2 (erosion-suppressed, spec 0.5% tol)',
              abs(grav - 9.005e-10) / 9.005e-10 < 0.005))
Ug_sum_trz = (Ug1 + Ug1 * sc) * 1.1
tests.append(('(Ug1+Ug4)*(1+f_TRZ) ~ 2.109e-9 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 2.109e-9) / 2.109e-9 < 0.005))

# ---- Lorentz term reuse ----
lor = m._lorentz_acceleration_uqff(B_T=m.B_PILLARS_T, v_ms=m.V_GAS_PILLARS_MS,
                                     q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Lorentz q*v*B/m_p * 11 * 1e-12 ~ 1.053e-4 (spec, same as Tapestry)',
              abs(lor - 1.053e-4) / 1.053e-4 < 0.005))

# ---- composer ----
g = m._pillars_g_master_uqff()
tests.append(('_pillars_g_master_uqff() defaults -> spec 1.053e-4 m/s^2 (0.5% tol)',
              abs(g - 1.053e-4) / 1.053e-4 < 0.005))
tests.append(('Total dominated by Lorentz term',
              (grav + Ug_sum_trz) < 0.01 * lor))

# erosion_factor multiplies ONLY grav_term, not Ug_sum_trz
g_no_erosion = m._pillars_g_master_uqff(E_0=0.0)
# Without erosion, grav_term goes from 9.005e-10 to 9.588e-10 (delta ~5.83e-11)
expected_delta = Ug1 * H_factor * sc * (1.0 - 0.93935)  # what erosion was removing
tests.append(('Pillars E_0=0 -> grav_term not suppressed (delta matches expected)',
              abs((g_no_erosion - g) - expected_delta) / expected_delta < 0.05))

# E_0=1.0 -> erosion at t=0 = 1.0 -> grav_term=0; at t=5e5yr E~0.6065 -> grav suppressed by ~0.39
g_full_erode = m._pillars_g_master_uqff(E_0=1.0)
tests.append(('Pillars E_0=1.0 reduces grav_term (sign check)',
              g_full_erode < g),)

# macro_scale=0 -> total ~ grav + Ug_sum_trz
g_no_lor = m._pillars_g_master_uqff(macro_scale=0.0)
tests.append(('Pillars macro_scale=0 -> total ~ grav + Ug_sum_trz (1% tol)',
              abs(g_no_lor - (grav + Ug_sum_trz)) / g_no_lor < 0.01))

# M_dot_0=0 -> M = M_initial
g_no_sf = m._pillars_g_master_uqff(M_dot_0=0.0)
tests.append(('Pillars g(M_dot_0=0) still ~1.053e-4 (Lorentz dominates)',
              abs(g_no_sf - 1.053e-4) / 1.053e-4 < 0.01))

# Monotone in M_dot_0
g_low = m._pillars_g_master_uqff(M_dot_0=0.5)
g_high = m._pillars_g_master_uqff(M_dot_0=5.0)
tests.append(('Pillars g monotone in M_dot_0',
              g_high > g > g_low))

# f_TRZ=0 reduces Ug_sum
g_no_trz = m._pillars_g_master_uqff(f_TRZ=0.0)
tests.append(('Pillars f_TRZ=0 reduces total (sign check)',
              g_no_trz < g))

# Existing primitives untouched
tests.append(('_pillars_g_primitive_sat untouched',
              isinstance(m._pillars_g_primitive_sat(), float)
                  and m._pillars_g_primitive_sat() > 0.0))
tests.append(('_m16_g_primitive_sat untouched',
              isinstance(m._m16_g_primitive_sat(), float)
                  and m._m16_g_primitive_sat() > 0.0))
tests.append(('pillars_of_creation catalog entry untouched',
              'pillars_of_creation' in m.ASTRO_SYSTEMS))

# Other composers untouched
tests.append(('Tapestry composer untouched (~1.053e-4)',
              abs(m._tapestry_g_master_uqff() - 1.053e-4) / 1.053e-4 < 0.005))
tests.append(('Westerlund 2 composer untouched (~1.053e-3)',
              abs(m._westerlund2_g_master_uqff() - 1.053e-3) / 1.053e-3 < 0.005))

passed = sum(1 for _, ok in tests if ok)
total = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
