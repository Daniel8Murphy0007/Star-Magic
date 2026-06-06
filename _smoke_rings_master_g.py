"""Smoke test: Rings of Relativity (GAL-CLUS-022058s Einstein ring) master Universal Gravity (spec 09May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_5gyr = 5.0e9 * YR_S

# ---- module-level constants ----
tests.append(('M_LENS_RINGS_KG = 1e14 * M_SUN',
              m.M_LENS_RINGS_KG == 1.0e14 * m.M_SUN))
tests.append(('M_LENS_RINGS_KG ~ 1.989e44 kg (spec, 0.1% tol)',
              abs(m.M_LENS_RINGS_KG - 1.989e44) / 1.989e44 < 0.001))
tests.append(('R_EINSTEIN_RINGS_M = 3.086e20 m (~10 kpc)',
              m.R_EINSTEIN_RINGS_M == 3.086e20))
tests.append(('B_RINGS_T = 1e-5 T',
              m.B_RINGS_T == 1.0e-5))
tests.append(('Z_LENS_RINGS = 0.5',
              m.Z_LENS_RINGS == 0.5))
tests.append(('Z_SOURCE_RINGS = 2.0',
              m.Z_SOURCE_RINGS == 2.0))
tests.append(('H0_RINGS_KMSMPC = 70.0 (spec, NOT 67.4)',
              m.H0_RINGS_KMSMPC == 70.0))
tests.append(('V_GAS_RINGS_MS = 1e6 m/s (10x prior nebular speeds)',
              m.V_GAS_RINGS_MS == 1.0e6))
tests.append(('T_RINGS_DEFAULT_S = 5 Gyr',
              abs(m.T_RINGS_DEFAULT_S - t_5gyr) < 1e-3))
tests.append(('_MPC_M = 3.086e22 m/Mpc',
              m._MPC_M == 3.086e22))

# ---- _D_LS_over_D_S_uqff ----
ratio = m._D_LS_over_D_S_uqff()
tests.append(('_D_LS_over_D_S_uqff(0.5, 2) = (1+0.5)/(1+2) = 0.5 (spec)',
              abs(ratio - 0.5) < 1e-12))
tests.append(('_D_LS_over_D_S_uqff(0, 0) = 1 (no redshift)',
              m._D_LS_over_D_S_uqff(0.0, 0.0) == 1.0))
tests.append(('_D_LS_over_D_S_uqff(0, 2) = 1/3',
              abs(m._D_LS_over_D_S_uqff(0.0, 2.0) - 1.0 / 3.0) < 1e-12))

# ---- _einstein_lensing_factor_uqff ----
L_t = m._einstein_lensing_factor_uqff()
GM_over_c2r = m.G_NEWTON * m.M_LENS_RINGS_KG / ((m.C_LIGHT ** 2) * m.R_EINSTEIN_RINGS_M)
tests.append(('GM/(c^2 r) ~ 4.775e-4 (spec, 0.5% tol)',
              abs(GM_over_c2r - 4.775e-4) / 4.775e-4 < 0.005))
tests.append(('_einstein_lensing_factor_uqff() = GM/(c^2 r) * 0.5 ~ 2.388e-4 (spec, 0.5% tol)',
              abs(L_t - 2.388e-4) / 2.388e-4 < 0.005))
tests.append(('L_t = GM/(c^2 r) * D_LS/D_S identity',
              abs(L_t - GM_over_c2r * 0.5) < 1e-12))

# ---- _hubble_unified at z=0.5 (REUSE check) ----
H_z = m._hubble_unified(t=t_5gyr, z=0.5, H_0=70.0)
tests.append(('_hubble_unified(z=0.5, H_0=70) ~ 91.63 km/s/Mpc (spec, 0.5% tol)',
              abs(H_z - 91.63) / 91.63 < 0.005))
H_z_si = H_z * 1.0e3 / m._MPC_M
tests.append(('H(z=0.5) in SI ~ 2.969e-18 s^-1 (spec, 0.5% tol)',
              abs(H_z_si - 2.969e-18) / 2.969e-18 < 0.005))
tests.append(('H(z=0.5)*t (5 Gyr) ~ 0.4685 (spec, 0.5% tol)',
              abs(H_z_si * t_5gyr - 0.4685) / 0.4685 < 0.005))
tests.append(('1 + H(z=0.5)*t ~ 1.4685 (spec)',
              abs((1.0 + H_z_si * t_5gyr) - 1.4685) / 1.4685 < 0.005))

# z=0 reduces H(z) to H_0 (sanity)
tests.append(('_hubble_unified(z=0, H_0=70) = 70 (no z evolution)',
              abs(m._hubble_unified(z=0.0, H_0=70.0) - 70.0) < 1e-10))

# ---- decomposition at t=5 Gyr ----
Ug1 = m.G_NEWTON * m.M_LENS_RINGS_KG / (m.R_EINSTEIN_RINGS_M ** 2)
tests.append(('U_g1 = G*M/r^2 ~ 1.394e-7 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 1.394e-7) / 1.394e-7 < 0.005))
sc = 1.0 - m.B_RINGS_T / m.B_CRIT_MAGNETAR_T
tests.append(('1 - B/B_crit ~ 1 (B=1e-5, B_crit=1e11)',
              abs(sc - 1.0) < 1e-12))
H_factor = 1.0 + H_z_si * t_5gyr
lensing_factor = 1.0 + L_t
grav = Ug1 * H_factor * sc * lensing_factor
tests.append(('grav_term ~ 2.047e-7 m/s^2 (lensing-amplified, spec 0.5% tol)',
              abs(grav - 2.047e-7) / 2.047e-7 < 0.005))
Ug_sum_trz = (Ug1 + Ug1 * sc) * 1.1
tests.append(('(Ug1+Ug4)*(1+f_TRZ) ~ 3.067e-7 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 3.067e-7) / 3.067e-7 < 0.005))

# ---- Lorentz term (v=1e6, B=1e-5 -> 10x Westerlund) ----
lor = m._lorentz_acceleration_uqff(B_T=m.B_RINGS_T, v_ms=m.V_GAS_RINGS_MS,
                                     q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Lorentz q*v*B/m_p * 11 * 1e-12 ~ 1.053e-2 (spec, 0.5% tol)',
              abs(lor - 1.053e-2) / 1.053e-2 < 0.005))
lor_west = m._lorentz_acceleration_uqff(B_T=m.B_WESTERLUND2_T, v_ms=m.V_GAS_WESTERLUND2_MS,
                                          q_C=m.EV_J, m_kg=m._M_PROTON_KG_MAGNETAR)
tests.append(('Rings Lorentz = 10 * Westerlund Lorentz (v 10x)',
              abs(lor / lor_west - 10.0) < 0.01))

# ---- composer ----
g = m._rings_g_master_uqff()
tests.append(('_rings_g_master_uqff() defaults -> spec 1.053e-2 m/s^2 (0.5% tol)',
              abs(g - 1.053e-2) / 1.053e-2 < 0.005))
tests.append(('Total dominated by Lorentz term',
              (grav + Ug_sum_trz) < 0.01 * lor))

# macro_scale=0 -> total ~ grav + Ug_sum_trz
g_no_lor = m._rings_g_master_uqff(macro_scale=0.0)
tests.append(('Rings macro_scale=0 -> total ~ grav + Ug_sum_trz (1% tol)',
              abs(g_no_lor - (grav + Ug_sum_trz)) / g_no_lor < 0.01))

# Setting z_lens = z_source -> D_LS/D_S = 1, L doubles
g_eq_z = m._rings_g_master_uqff(z_lens=1.0, z_source=1.0)
L_eq = m._einstein_lensing_factor_uqff(z_lens=1.0, z_source=1.0)
tests.append(('z_lens=z_source -> D_LS/D_S=1, L = GM/(c^2 r) full (~4.775e-4)',
              abs(L_eq - GM_over_c2r) < 1e-12))

# z_source -> inf -> D_LS/D_S -> 0, L -> 0, no lensing amplification
L_inf = m._einstein_lensing_factor_uqff(z_lens=0.5, z_source=1e10)
tests.append(('z_source -> inf -> L -> 0 (no lensing)',
              L_inf < 1e-12))

# H_0 = 67.4 (Planck) vs spec 70 -- small difference
H_planck = m._hubble_unified(z=0.5, H_0=67.4)
H_spec = m._hubble_unified(z=0.5, H_0=70.0)
tests.append(('H(z=0.5, H_0=67.4) < H(z=0.5, H_0=70) (linear in H_0)',
              abs(H_planck / H_spec - 67.4 / 70.0) < 1e-10))

# f_TRZ=0 reduces total
g_no_trz = m._rings_g_master_uqff(f_TRZ=0.0)
tests.append(('Rings f_TRZ=0 reduces total (sign check)',
              g_no_trz < g))

# Existing primitives untouched
tests.append(('_rings_of_relativity_g_primitive_sat untouched',
              isinstance(m._rings_of_relativity_g_primitive_sat(), float)
                  and m._rings_of_relativity_g_primitive_sat() > 0.0))
tests.append(('rings_of_relativity catalog entry untouched',
              'rings_of_relativity' in m.ASTRO_SYSTEMS))
tests.append(('_l95_wl_kappa_uqff (weak-lensing kappa, different physics) untouched',
              isinstance(m._l95_wl_kappa_uqff(), float)))

# Other composers untouched
tests.append(('Tapestry composer untouched (~1.053e-4)',
              abs(m._tapestry_g_master_uqff() - 1.053e-4) / 1.053e-4 < 0.005))
tests.append(('Westerlund 2 composer untouched (~1.053e-3)',
              abs(m._westerlund2_g_master_uqff() - 1.053e-3) / 1.053e-3 < 0.005))
tests.append(('Pillars composer untouched (~1.053e-4)',
              abs(m._pillars_g_master_uqff() - 1.053e-4) / 1.053e-4 < 0.005))

passed = sum(1 for _, ok in tests if ok)
total = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
