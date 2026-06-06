"""Smoke test: Galaxy NGC 2525 (SN 2018gv host) master Universal Gravity (spec 08May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_7yr = 7.0 * YR_S
t_70Myr = 70.0e6 * YR_S

# ---- constants ----
tests.append(('M_NGC2525_STAR_KG = 1e10 * M_SUN',
              m.M_NGC2525_STAR_KG == 1.0e10 * m.M_SUN))
tests.append(('M_NGC2525_BH_KG = 2.25e7 * M_SUN',
              m.M_NGC2525_BH_KG == 2.25e7 * m.M_SUN))
tests.append(('M_NGC2525_TOTAL_KG ~ 1.993e40 kg (spec, 0.1% tol)',
              abs(m.M_NGC2525_TOTAL_KG - 1.993e40) / 1.993e40 < 0.001))
tests.append(('R_NGC2525_DISK_M = 2.836e20 m',
              m.R_NGC2525_DISK_M == 2.836e20))
tests.append(('R_NGC2525_BH_M = 1.496e11 m (1 AU)',
              m.R_NGC2525_BH_M == 1.496e11))
tests.append(('B_NGC2525_T = 1e-5',
              m.B_NGC2525_T == 1.0e-5))
tests.append(('Z_NGC2525 = 0.016',
              m.Z_NGC2525 == 0.016))
tests.append(('H0_NGC2525_KMSMPC = 70.0',
              m.H0_NGC2525_KMSMPC == 70.0))
tests.append(('V_GAS_NGC2525_MS = 1e5',
              m.V_GAS_NGC2525_MS == 1.0e5))
tests.append(('M_SN_2018GV_INIT_KG = 1.4 * M_SUN',
              m.M_SN_2018GV_INIT_KG == 1.4 * m.M_SUN))
tests.append(('TAU_SN_2018GV_S = 1 yr',
              abs(m.TAU_SN_2018GV_S - YR_S) < 1.0))
tests.append(('T_NGC2525_DEFAULT_S = 7 yr',
              abs(m.T_NGC2525_DEFAULT_S - t_7yr) < 1.0))

# ---- gravitational decomposition ----
Ug1 = m.G_NEWTON * m.M_NGC2525_TOTAL_KG / (m.R_NGC2525_DISK_M ** 2)
tests.append(('Ug1 = G*M/r^2 ~ 1.654e-11 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 1.654e-11) / 1.654e-11 < 0.005))

bh_term = m.G_NEWTON * m.M_NGC2525_BH_KG / (m.R_NGC2525_BH_M ** 2)
tests.append(('G*M_BH/r_BH^2 ~ 1.334e5 m/s^2 (spec, 0.5% tol)',
              abs(bh_term - 1.334e5) / 1.334e5 < 0.005))

# H(z=0.016, H_0=70)
H_z = m._hubble_unified(t=t_7yr, z=0.016, H_0=70.0)
tests.append(('_hubble_unified(z=0.016, H_0=70) ~ 70.51 km/s/Mpc (spec, 0.5% tol)',
              abs(H_z - 70.51) / 70.51 < 0.005))
H_si = H_z * 1.0e3 / m._MPC_M
tests.append(('H(z=0.016) in SI ~ 2.285e-18 s^-1 (spec, 0.5% tol)',
              abs(H_si - 2.285e-18) / 2.285e-18 < 0.005))

# At t=7 yr, H*t is tiny; spec uses t=70 Myr lookback for stellar evol term
H_t_7yr = H_si * t_7yr
tests.append(('H*t at 7 yr is tiny (<<1)',
              H_t_7yr < 1e-8))
H_t_70Myr = H_si * t_70Myr
tests.append(('H*t at 70 Myr ~ 5.047e-3 (spec, 0.5% tol)',
              abs(H_t_70Myr - 5.047e-3) / 5.047e-3 < 0.005))

# sc_factor ~ 1
sc = 1.0 - m.B_NGC2525_T / m.B_CRIT_MAGNETAR_T
tests.append(('1 - B/B_crit ~ 1',
              abs(sc - 1.0) < 1e-12))

# Ug_sum * (1 + f_TRZ) = 2*Ug1*1.1
Ug_sum_trz = (Ug1 + Ug1 * sc) * 1.1
tests.append(('(Ug1+Ug4)(1+f_TRZ) ~ 3.639e-11 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 3.639e-11) / 3.639e-11 < 0.005))

# Lorentz at v=1e5, B=1e-5 -> 1.053e-3 (same as Westerlund 2)
lor = m._lorentz_acceleration_uqff(B_T=m.B_NGC2525_T, v_ms=m.V_GAS_NGC2525_MS,
                                     q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Lorentz*macro ~ 1.053e-3 m/s^2 (spec, 0.5% tol)',
              abs(lor - 1.053e-3) / 1.053e-3 < 0.005))
lor_west = m._lorentz_acceleration_uqff(B_T=m.B_WESTERLUND2_T, v_ms=m.V_GAS_WESTERLUND2_MS,
                                          q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('NGC 2525 Lorentz == Westerlund 2 Lorentz (same B, same v)',
              abs(lor - lor_west) / lor_west < 1e-12))

# SN decay reuses _magnetar_B_decay
M_SN_at_7yr = m._magnetar_B_decay(t_7yr, m.M_SN_2018GV_INIT_KG, m.TAU_SN_2018GV_S)
tests.append(('M_SN(7yr) = M_SN_0 * exp(-7) ~ 2.539e27 kg (spec, 0.5% tol)',
              abs(M_SN_at_7yr - 2.539e27) / 2.539e27 < 0.005))
sn_drain = m.G_NEWTON * M_SN_at_7yr / (m.R_NGC2525_DISK_M ** 2)
tests.append(('-G*M_SN(7yr)/r^2 magnitude ~ 2.107e-24 m/s^2 (spec, 0.5% tol)',
              abs(sn_drain - 2.107e-24) / 2.107e-24 < 0.005))

# At t=0 M_SN = full Chandrasekhar mass
M_SN_t0 = m._magnetar_B_decay(0.0, m.M_SN_2018GV_INIT_KG, m.TAU_SN_2018GV_S)
tests.append(('M_SN(0) = 1.4*M_sun (no decay)',
              abs(M_SN_t0 - m.M_SN_2018GV_INIT_KG) < 1e-10))

# ---- composer ----
g = m._ngc2525_g_master_uqff()
tests.append(('_ngc2525_g_master_uqff() defaults -> spec ~1.335e5 m/s^2 (0.1% tol)',
              abs(g - 1.335e5) / 1.335e5 < 0.001))
tests.append(('Total dominated by BH term (>99.9999%)',
              bh_term / g > 0.99999))

# macro_scale=0 -> Lorentz drops out
g_no_lor = m._ngc2525_g_master_uqff(macro_scale=0.0)
tests.append(('macro_scale=0 -> total still ~ BH term (Lorentz drops)',
              abs(g_no_lor - bh_term) / bh_term < 1e-6))

# Setting M_SN_init=0 removes drain (isolate by zeroing dominant BH + Lorentz)
g_with_sn_iso = m._ngc2525_g_master_uqff(M_BH_kg=0.0, macro_scale=0.0,
                                            M_SN_init_kg=1e35)  # exaggerate to be visible
g_no_sn_iso = m._ngc2525_g_master_uqff(M_BH_kg=0.0, macro_scale=0.0,
                                          M_SN_init_kg=0.0)
tests.append(('M_SN_init=0 -> SN drain term vanishes (isolated)',
              g_no_sn_iso > g_with_sn_iso))

# Setting M_BH=0 removes the dominant term -> total drops to disk + Lorentz scale
g_no_bh = m._ngc2525_g_master_uqff(M_BH_kg=0.0)
tests.append(('M_BH=0 -> total ~ disk + Lorentz (no BH dominance)',
              g_no_bh < 1.0))

# Setting f_TRZ=0 reduces Ug_sum (tiny but consistent)
g_no_trz = m._ngc2525_g_master_uqff(f_TRZ=0.0, M_BH_kg=0.0, macro_scale=0.0)
g_with_trz = m._ngc2525_g_master_uqff(f_TRZ=0.1, M_BH_kg=0.0, macro_scale=0.0)
tests.append(('f_TRZ=0.1 > f_TRZ=0 (Ug_sum amplified)',
              g_with_trz > g_no_trz))

# At t=70 Myr (stellar evolution timescale), disk term gets H*t correction
g_70Myr = m._ngc2525_g_master_uqff(t_s=t_70Myr, M_BH_kg=0.0, macro_scale=0.0,
                                      M_SN_init_kg=0.0)
Ug1_disk_70Myr = Ug1 * (1.0 + H_t_70Myr) * sc
expected_disk = Ug1_disk_70Myr + Ug_sum_trz + m.LAMBDA_MAGNETAR_M2 * m.C_LIGHT**2 / 3.0
tests.append(('disk-only at 70 Myr matches (1+H*t) amplification (1% tol)',
              abs(g_70Myr - expected_disk) / expected_disk < 0.01))

# z=0 reduces H(z) to H_0 (sanity)
tests.append(('_hubble_unified(z=0, H_0=70) = 70',
              abs(m._hubble_unified(z=0.0, H_0=70.0) - 70.0) < 1e-10))

# Existing primitives untouched
tests.append(('_ngc_2525_g_primitive_sat (BSFG saturation triadic) untouched',
              isinstance(m._ngc_2525_g_primitive_sat(), float)))
tests.append(('ngc_2525_g (BSFG saturation) still in _LEDGER_PRIMITIVE registry',
              'ngc_2525_g' in m._LEDGER_PRIMITIVE))

# Other composers untouched
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
