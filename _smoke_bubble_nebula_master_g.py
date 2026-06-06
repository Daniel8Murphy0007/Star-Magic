"""Smoke test: Bubble Nebula NGC 7635 (BD +60 2522 Wolf-Rayet host) master Universal Gravity (spec 09May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_4Myr = 4.0e6 * YR_S

# ---- constants ----
tests.append(('M_NGC7635_STAR_KG = 45 * M_SUN',
              m.M_NGC7635_STAR_KG == 45.0 * m.M_SUN))
tests.append(('M_NGC7635_STAR_KG ~ 8.951e31 kg (spec, 0.1% tol)',
              abs(m.M_NGC7635_STAR_KG - 8.951e31) / 8.951e31 < 0.001))
tests.append(('R_NGC7635_M = 3.311e16 m (3.5 ly half-span)',
              m.R_NGC7635_M == 3.311e16))
tests.append(('B_NGC7635_T = 1e-6 (weaker than NGC 3603)',
              m.B_NGC7635_T == 1.0e-6))
tests.append(('V_WIND_NGC7635_MS = 1.789e6 (4e6 mph Wolf-Rayet)',
              m.V_WIND_NGC7635_MS == 1.789e6))
tests.append(('RHO_GAS_NGC7635 = 1e-21',
              m.RHO_GAS_NGC7635 == 1.0e-21))
tests.append(('P_0_NGC7635 = 0.1 (normalized fractional, NOT dimensional Pa)',
              m.P_0_NGC7635 == 0.1))
tests.append(('TAU_EXP_NGC7635_S = 4 Myr',
              abs(m.TAU_EXP_NGC7635_S - 4.0e6 * YR_S) < 1.0))
tests.append(('T_NGC7635_DEFAULT_S = 4 Myr',
              abs(m.T_NGC7635_DEFAULT_S - t_4Myr) < 1.0))
tests.append(('H0_NGC7635_KMSMPC = 70.0',
              m.H0_NGC7635_KMSMPC == 70.0))

# ---- gravitational decomposition ----
Ug1 = m.G_NEWTON * m.M_NGC7635_STAR_KG / (m.R_NGC7635_M ** 2)
tests.append(('G*M_star/r^2 ~ 5.449e-12 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 5.449e-12) / 5.449e-12 < 0.005))

# H_0 local (z=0)
H0_kmsMpc_val = m._hubble_unified(t_4Myr, 0.0, 70.0)
tests.append(('_hubble_unified(z=0) = 70 km/s/Mpc',
              abs(H0_kmsMpc_val - 70.0) < 1e-10))
H0_si = H0_kmsMpc_val * 1.0e3 / m._MPC_M
tests.append(('H_0 in SI ~ 2.268e-18 s^-1 (spec, 0.5% tol)',
              abs(H0_si - 2.268e-18) / 2.268e-18 < 0.005))
H_t = H0_si * t_4Myr
tests.append(('H_0 * t at 4 Myr ~ 2.863e-4 (spec, 0.5% tol)',
              abs(H_t - 2.863e-4) / 2.863e-4 < 0.005))
tests.append(('1 + H_0*t ~ 1.0002863 (spec)',
              abs((1.0 + H_t) - 1.0002863) / 1.0002863 < 0.005))

# P(t) wind cavity envelope
P_t = m._magnetar_B_decay(t_4Myr, m.P_0_NGC7635, m.TAU_EXP_NGC7635_S)
tests.append(('P(t=4 Myr) = 0.1 * exp(-1) ~ 0.03679 (spec, 0.5% tol)',
              abs(P_t - 0.03679) / 0.03679 < 0.005))
tests.append(('1 - P(t) ~ 0.96321 (spec, 0.5% tol)',
              abs((1.0 - P_t) - 0.96321) / 0.96321 < 0.005))

# Grav-term with (1+f_TRZ) factored on
grav_term = Ug1 * (1.0 + H_t) * (1.0 - P_t) * 1.1
tests.append(('grav*(1+f_TRZ) ~ 5.78e-12 m/s^2 (spec, 1% tol)',
              abs(grav_term - 5.78e-12) / 5.78e-12 < 0.01))

# ---- Lorentz with B=1e-6, v_wind=1.789e6 -> NEW dominant scaling 1.884e-3 ----
lor = m._lorentz_acceleration_uqff(B_T=m.B_NGC7635_T, v_ms=m.V_WIND_NGC7635_MS,
                                     q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Lorentz q*v*B/m_p*11*1e-12 ~ 1.884e-3 m/s^2 (spec, 0.5% tol)',
              abs(lor - 1.884e-3) / 1.884e-3 < 0.005))

# Verify Lorentz NOT 1.053e-3 (the prior family) -- this is a genuinely new scaling
lor_west = m._lorentz_acceleration_uqff(B_T=1e-5, v_ms=1e5,
                                          q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Bubble Lorentz != 1.053e-3 family (B=1e-6 vs 1e-5; v=1.789e6 vs 1e5)',
              abs(lor / lor_west - 1.789) < 0.01))  # ratio = (1e-6 * 1.789e6) / (1e-5 * 1e5) = 1.789

# ---- composer ----
g = m._bubble_nebula_g_master_uqff()
tests.append(('_bubble_nebula_g_master_uqff() defaults -> spec ~1.884e-3 m/s^2 (0.5% tol)',
              abs(g - 1.884e-3) / 1.884e-3 < 0.005))
tests.append(('Total dominated by Lorentz term (grav < 1e-7 * Lorentz)',
              grav_term < 1e-7 * lor))

# macro_scale=0 -> only grav-term remains
g_no_lor = m._bubble_nebula_g_master_uqff(macro_scale=0.0)
tests.append(('macro_scale=0 -> total = grav*(1+f_TRZ) only (1% tol)',
              abs(g_no_lor - grav_term) / grav_term < 0.01))

# f_TRZ=0 -> grav drops 10%
g_no_trz = m._bubble_nebula_g_master_uqff(f_TRZ=0.0, macro_scale=0.0)
tests.append(('f_TRZ=0 -> grav term drops to (1+H*t)(1-P)*Ug1 (no 1.1 amplification)',
              abs(g_no_trz - grav_term / 1.1) / (grav_term / 1.1) < 0.01))

# P_0=0 -> no cavity suppression, grav term increases
g_no_P = m._bubble_nebula_g_master_uqff(P_0=0.0, macro_scale=0.0)
g_with_P = m._bubble_nebula_g_master_uqff(P_0=0.1, macro_scale=0.0)
tests.append(('P_0=0 -> no suppression, grav > with-P case',
              g_no_P > g_with_P))

# At t=0: P(0)=P_0 full (max suppression), no Hubble shift
P_t0 = m._magnetar_B_decay(0.0, m.P_0_NGC7635, m.TAU_EXP_NGC7635_S)
tests.append(('P(t=0) = P_0 = 0.1 (full, no decay yet)',
              abs(P_t0 - 0.1) < 1e-10))

# t -> 10 Myr: P decays further, (1+H*t) grows
g_10Myr = m._bubble_nebula_g_master_uqff(t_s=10.0e6 * YR_S, macro_scale=0.0)
g_4Myr = m._bubble_nebula_g_master_uqff(macro_scale=0.0)
P_10Myr = m._magnetar_B_decay(10.0e6 * YR_S, m.P_0_NGC7635, m.TAU_EXP_NGC7635_S)
tests.append(('P(10 Myr) < P(4 Myr) (decay continues)',
              P_10Myr < 0.03679))
tests.append(('g(10 Myr) > g(4 Myr) (less suppression + more H_0*t)',
              g_10Myr > g_4Myr))

# v_wind=0 -> Lorentz zero (NO wind-driven Lorentz drive)
g_no_wind = m._bubble_nebula_g_master_uqff(v_wind_ms=0.0, macro_scale=m.MACROSCOPIC_SCALE_LORENTZ)
tests.append(('v_wind=0 -> Lorentz contribution -> 0',
              abs(g_no_wind - grav_term) / grav_term < 0.01))

# B=0 -> Lorentz zero
g_no_B = m._bubble_nebula_g_master_uqff(B_T=0.0)
tests.append(('B=0 -> Lorentz contribution -> 0',
              abs(g_no_B - grav_term) / grav_term < 0.01))

# Confirm composer does NOT inject extra Ug_sum or Lambda*c^2/3 leaves
# (fidelity check: spec has 4 leaves only)
expected_total = grav_term + lor
tests.append(('Composer == exactly grav*(1+f_TRZ) + Lorentz (no extra Ug_sum, no Lambda*c^2/3)',
              abs(g - expected_total) / expected_total < 1e-12))

# Existing primitives untouched
tests.append(('_bubble_g_primitive_sat (BSFG saturation triadic) untouched',
              isinstance(m._bubble_g_primitive_sat(), float)))
tests.append(('bubble_nebula catalog entry (M=44 M_sun, different mass) untouched',
              m.ASTRO_SYSTEMS['bubble_nebula']['M'] == 44.0 * m.M_SUN))
tests.append(('bubble_g alias chain untouched (bubble, bubble_nebula, ngc_7635)',
              'bubble_g' in m._LEDGER_PRIMITIVE))

# Other composers untouched
tests.append(('NGC 3603 composer untouched (~1.053e-3)',
              abs(m._ngc3603_g_master_uqff() - 1.053e-3) / 1.053e-3 < 0.005))
tests.append(('NGC 2525 composer untouched (~1.335e5)',
              abs(m._ngc2525_g_master_uqff() - 1.335e5) / 1.335e5 < 0.001))
tests.append(('Rings composer untouched (~1.053e-2)',
              abs(m._rings_g_master_uqff() - 1.053e-2) / 1.053e-2 < 0.005))
tests.append(('Pillars composer untouched (~1.053e-4)',
              abs(m._pillars_g_master_uqff() - 1.053e-4) / 1.053e-4 < 0.005))
tests.append(('Westerlund 2 composer untouched (~1.053e-3)',
              abs(m._westerlund2_g_master_uqff() - 1.053e-3) / 1.053e-3 < 0.005))
tests.append(('Tapestry composer untouched (~1.053e-4)',
              abs(m._tapestry_g_master_uqff() - 1.053e-4) / 1.053e-4 < 0.005))

passed = sum(1 for _, ok in tests if ok)
total = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
