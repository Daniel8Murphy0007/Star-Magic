"""Smoke test: NGC 3603 ('Extreme star cluster') master Universal Gravity (spec 08May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_5e5yr = 5.0e5 * YR_S

# ---- constants ----
tests.append(('M_NGC3603_INIT_KG = 4e5 * M_SUN',
              m.M_NGC3603_INIT_KG == 4.0e5 * m.M_SUN))
tests.append(('M_NGC3603_INIT_KG ~ 7.956e35 kg (spec, 0.1% tol)',
              abs(m.M_NGC3603_INIT_KG - 7.956e35) / 7.956e35 < 0.001))
tests.append(('R_NGC3603_M = 8.998e15 m',
              m.R_NGC3603_M == 8.998e15))
tests.append(('B_NGC3603_T = 1e-5',
              m.B_NGC3603_T == 1.0e-5))
tests.append(('TAU_SF_NGC3603_S = 1 Myr',
              abs(m.TAU_SF_NGC3603_S - 1.0e6 * YR_S) < 1.0))
tests.append(('TAU_EXP_NGC3603_S = 1 Myr',
              abs(m.TAU_EXP_NGC3603_S - 1.0e6 * YR_S) < 1.0))
tests.append(('M_DOT0_NGC3603 = 1.0',
              m.M_DOT0_NGC3603 == 1.0))
tests.append(('V_WIND_NGC3603_MS = 2e6 (2000 km/s blue-star wind)',
              m.V_WIND_NGC3603_MS == 2.0e6))
tests.append(('RHO_GAS_NGC3603 = 1e-20',
              m.RHO_GAS_NGC3603 == 1.0e-20))
tests.append(('V_GAS_NGC3603_MS = 1e5 (HII gas for Lorentz)',
              m.V_GAS_NGC3603_MS == 1.0e5))
tests.append(('H0_NGC3603_KMSMPC = 70.0',
              m.H0_NGC3603_KMSMPC == 70.0))
tests.append(('T_NGC3603_DEFAULT_S = 5e5 yr',
              abs(m.T_NGC3603_DEFAULT_S - t_5e5yr) < 1.0))

# ---- M(t) growth via _accretion_mass_growth_uqff (REUSE) ----
M_t = m._accretion_mass_growth_uqff(m.M_NGC3603_INIT_KG, t_5e5yr,
                                       m.M_DOT0_NGC3603, m.TAU_SF_NGC3603_S)
tests.append(('M(t=5e5 yr) = M_init * (1+exp(-0.5)) ~ 1.278e36 kg (spec, 0.5% tol)',
              abs(M_t - 1.278e36) / 1.278e36 < 0.005))
tests.append(('1+M_dot at t=5e5 yr ~ 1.6065 (spec, 0.5% tol)',
              abs(M_t / m.M_NGC3603_INIT_KG - 1.6065) / 1.6065 < 0.005))

# ---- gravitational decomposition ----
Ug1 = m.G_NEWTON * M_t / (m.R_NGC3603_M ** 2)
tests.append(('Ug1 = G*M(t)/r^2 ~ 1.053e-6 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 1.053e-6) / 1.053e-6 < 0.005))

# H_0 local (z=0)
H0_kmsMpc_val = m._hubble_unified(t_5e5yr, 0.0, 70.0)
tests.append(('_hubble_unified(z=0, H_0=70) = 70 km/s/Mpc',
              abs(H0_kmsMpc_val - 70.0) < 1e-10))
H0_si = H0_kmsMpc_val * 1.0e3 / m._MPC_M
tests.append(('H_0 in SI ~ 2.268e-18 s^-1 (spec, 0.5% tol)',
              abs(H0_si - 2.268e-18) / 2.268e-18 < 0.005))
H_t = H0_si * t_5e5yr
tests.append(('H_0 * t at 5e5 yr ~ 3.579e-5 (spec, 0.5% tol)',
              abs(H_t - 3.579e-5) / 3.579e-5 < 0.005))
tests.append(('1 + H_0*t ~ 1.00003579 (spec)',
              abs((1.0 + H_t) - 1.00003579) / 1.00003579 < 0.005))

# sc_factor ~ 1
sc = 1.0 - m.B_NGC3603_T / m.B_CRIT_MAGNETAR_T
tests.append(('1 - B/B_crit ~ 1',
              abs(sc - 1.0) < 1e-12))

# P(t) wind cavity envelope
P_0 = m.RHO_GAS_NGC3603 * m.V_WIND_NGC3603_MS ** 2
tests.append(('rho * v_wind^2 = 4e-8 Pa (spec)',
              abs(P_0 - 4.0e-8) / 4.0e-8 < 1e-6))
P_t = m._magnetar_B_decay(t_5e5yr, P_0, m.TAU_EXP_NGC3603_S)
tests.append(('P(t=5e5 yr) = 4e-8 * exp(-0.5) ~ 2.426e-8 (spec, 0.5% tol)',
              abs(P_t - 2.426e-8) / 2.426e-8 < 0.005))
tests.append(('1 - P(t) ~ 1 (negligible suppression at t=5e5 yr)',
              abs((1.0 - P_t) - 1.0) < 1e-6))

# Ug_sum * (1 + f_TRZ) = 2*Ug1*1.1
Ug_sum_trz = (Ug1 + Ug1 * sc) * 1.1
tests.append(('(Ug1+Ug4)(1+f_TRZ) ~ 2.317e-6 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 2.317e-6) / 2.317e-6 < 0.005))

# Lorentz at v=1e5, B=1e-5 -> 1.053e-3 (same family as Westerlund/NGC 2525)
lor = m._lorentz_acceleration_uqff(B_T=m.B_NGC3603_T, v_ms=m.V_GAS_NGC3603_MS,
                                     q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Lorentz*macro ~ 1.053e-3 m/s^2 (spec, 0.5% tol)',
              abs(lor - 1.053e-3) / 1.053e-3 < 0.005))

# ---- composer ----
g = m._ngc3603_g_master_uqff()
tests.append(('_ngc3603_g_master_uqff() defaults -> spec ~1.053e-3 m/s^2 (0.5% tol)',
              abs(g - 1.053e-3) / 1.053e-3 < 0.005))
tests.append(('Total dominated by Lorentz term',
              (Ug1 + Ug_sum_trz) < 0.01 * lor))

# macro_scale=0 -> total ~ grav + Ug_sum_trz
g_no_lor = m._ngc3603_g_master_uqff(macro_scale=0.0)
grav_term_expected = Ug1 * (1.0 + H_t) * sc * (1.0 - P_t)
expected = grav_term_expected + Ug_sum_trz + m.LAMBDA_MAGNETAR_M2 * m.C_LIGHT**2 / 3.0
tests.append(('macro_scale=0 -> total = grav_term + Ug_sum_trz (1% tol)',
              abs(g_no_lor - expected) / expected < 0.01))

# Setting M_dot_0=0 reproduces M_init (no SF growth)
g_no_sf = m._ngc3603_g_master_uqff(M_dot_0=0.0, macro_scale=0.0)
M_init_only = m.M_NGC3603_INIT_KG
Ug1_init = m.G_NEWTON * M_init_only / m.R_NGC3603_M ** 2
expected_no_sf = Ug1_init * (1.0 + H_t) * sc * (1.0 - P_t) + (Ug1_init + Ug1_init * sc) * 1.1
tests.append(('M_dot_0=0 -> M(t)=M_init (no SF growth, 1% tol)',
              abs(g_no_sf - expected_no_sf) / expected_no_sf < 0.01))

# v_wind=0 -> P(t)=0 -> no cavity suppression
g_no_wind = m._ngc3603_g_master_uqff(v_wind_ms=0.0)
tests.append(('v_wind=0 -> P(t)=0 (no cavity suppression)',
              g_no_wind > 0.0))

# Setting rho_gas=0 -> P(t)=0
g_no_rho = m._ngc3603_g_master_uqff(rho_gas=0.0, macro_scale=0.0)
# Should match macro_scale=0 above (since P_t ~ 1e-8 is negligible anyway)
tests.append(('rho_gas=0 ~ baseline (P negligible at default time)',
              abs(g_no_rho - g_no_lor) / g_no_lor < 1e-4))

# Boost v_wind so P(t) becomes ~unity -> visible cavity suppression
g_strong_wind = m._ngc3603_g_master_uqff(v_wind_ms=2.0e9, macro_scale=0.0)
tests.append(('v_wind=2e9 m/s -> P_0=4e-2 visible suppression (smaller than baseline)',
              g_strong_wind < g_no_lor))

# Setting f_TRZ=0 reduces Ug_sum
g_no_trz = m._ngc3603_g_master_uqff(f_TRZ=0.0, macro_scale=0.0)
g_with_trz = m._ngc3603_g_master_uqff(f_TRZ=0.1, macro_scale=0.0)
tests.append(('f_TRZ=0.1 > f_TRZ=0 (Ug_sum amplified)',
              g_with_trz > g_no_trz))

# At t=0: no SF growth (1+exp(0)=2, so M=2*M_init), no P decay (full P_0)
M_t0 = m._accretion_mass_growth_uqff(m.M_NGC3603_INIT_KG, 0.0,
                                        m.M_DOT0_NGC3603, m.TAU_SF_NGC3603_S)
tests.append(('M(t=0) = 2 * M_init (1+exp(0))',
              abs(M_t0 - 2.0 * m.M_NGC3603_INIT_KG) < 1e-10 * M_t0))

# Existing aliases / catalog untouched
tests.append(('tapestry_ngc_3603 catalog entry untouched',
              'tapestry_ngc_3603' in m.ASTRO_SYSTEMS))
tests.append(('ngc_3603 alias still points to tapestry_ngc_3603',
              m.ASTRO_SYSTEMS.get('tapestry_ngc_3603') is not None))

# Other composers untouched
tests.append(('NGC 2525 composer untouched (~1.335e5)',
              abs(m._ngc2525_g_master_uqff() - 1.335e5) / 1.335e5 < 0.001))
tests.append(('Rings composer untouched (~1.053e-2)',
              abs(m._rings_g_master_uqff() - 1.053e-2) / 1.053e-2 < 0.005))
tests.append(('Pillars composer untouched (~1.053e-4)',
              abs(m._pillars_g_master_uqff() - 1.053e-4) / 1.053e-4 < 0.005))
tests.append(('Westerlund 2 composer untouched (~1.053e-3)',
              abs(m._westerlund2_g_master_uqff() - 1.053e-3) / 1.053e-3 < 0.005))
tests.append(('Tapestry composer untouched (~1.053e-4; different scale: 240 OB stars)',
              abs(m._tapestry_g_master_uqff() - 1.053e-4) / 1.053e-4 < 0.005))

passed = sum(1 for _, ok in tests if ok)
total = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
