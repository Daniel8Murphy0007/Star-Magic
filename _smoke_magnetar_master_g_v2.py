"""Smoke test: magnetar evolution master Universal Gravity v2 (spec 08May2025).

v2 = v1 plus:
  - (Ug1+Ug2+Ug3+Ug4) * (1 + f_TRZ)              TRZ multiplier on Ug_sum
  - q (v x B) * (1 + rho_UA/rho_SCm) * 1e-12     Lorentz with [UA] enh + macro scaling

Spec target: g_Magnetar(t=5000 yr) ~ 4.474e12 m/s^2.
"""
import math
import uqff_pure_calculator as m

tests = []

# ---- module-level constants ----
tests.append(('V_MAGNETAR_SURFACE_MS = 1e6 m/s',
              m.V_MAGNETAR_SURFACE_MS == 1.0e6))
tests.append(('MACROSCOPIC_SCALE_LORENTZ = 1e-12',
              m.MACROSCOPIC_SCALE_LORENTZ == 1.0e-12))
tests.append(('F_TRZ_DEFAULT = 0.1 (reused)',
              m.F_TRZ_DEFAULT == 0.1))
tests.append(('EV_J reused as elementary charge magnitude (~1.602e-19)',
              abs(m.EV_J - 1.602176634e-19) < 1e-30))

# ---- _lorentz_acceleration_uqff ----
# Bare per-charge: a = q v B / m_p
# Spec: q=1.602e-19, v=1e6, B=2.865e9 T -> F=4.59e-4 N, a=2.744e23 m/s^2
B_at_5kyr = m._magnetar_B_decay(5000.0 * 365.25 * 86400.0)
tests.append(('_magnetar_B_decay(t=5000 yr) ~ 2.865e9 T (spec, 0.1% tol)',
              abs(B_at_5kyr - 2.865e9) / 2.865e9 < 0.001))

a_bare = m._lorentz_acceleration_uqff(B_T=B_at_5kyr, macro_scale=1.0,
                                       rho_UA_val=0.0, rho_SCm_val=1.0)
# a_bare = q v B / m_p * (1 + 0/1) * 1.0
expected_bare = m.EV_J * 1.0e6 * B_at_5kyr / m._M_PROTON_KG_MAGNETAR
tests.append(('_lorentz_acceleration_uqff bare (rho_UA=0, macro=1) = q v B / m_p',
              abs(a_bare - expected_bare) / expected_bare < 1e-12))

# UA-enhanced, macro-scaled, default ratio = 10 -> factor 11
a_macro = m._lorentz_acceleration_uqff(B_T=B_at_5kyr)
tests.append(('_lorentz_acceleration_uqff defaults ~ 3.018e12 m/s^2 at t=5000 yr (1% tol)',
              abs(a_macro - 3.018e12) / 3.018e12 < 0.01))
tests.append(('Lorentz UA enhancement = 11 (G-lock ratio 10)',
              abs(a_macro / (expected_bare * 1.0e-12) - 11.0) < 1e-9))

# Scales linearly in B, v, q
tests.append(('_lorentz_acceleration_uqff scales linearly in B',
              abs(m._lorentz_acceleration_uqff(B_T=2.0 * B_at_5kyr)
                  - 2.0 * a_macro) / a_macro < 1e-12))
tests.append(('_lorentz_acceleration_uqff scales linearly in v',
              abs(m._lorentz_acceleration_uqff(B_T=B_at_5kyr, v_ms=2.0e6)
                  - 2.0 * a_macro) / a_macro < 1e-12))
# Inverse in m_kg
tests.append(('_lorentz_acceleration_uqff scales 1/m_kg',
              abs(m._lorentz_acceleration_uqff(B_T=B_at_5kyr, m_kg=2.0 * m._M_PROTON_KG_MAGNETAR)
                  - 0.5 * a_macro) / a_macro < 1e-12))
# B=0 -> 0
tests.append(('_lorentz_acceleration_uqff(B=0) = 0',
              m._lorentz_acceleration_uqff(B_T=0.0) == 0.0))

# ---- _magnetar_g_master_uqff_v2 spec example (t=5000 yr) ----
g2 = m._magnetar_g_master_uqff_v2()  # defaults: t = 5000 yr
tests.append(('_magnetar_g_master_uqff_v2() defaults -> spec example 4.474e12 m/s^2 (1% tol)',
              abs(g2 - 4.474e12) / 4.474e12 < 0.01))

# Decomposition check at t=5000 yr
Ug1 = m.G_NEWTON * m.M_MAGNETAR_KG / (m.R_MAGNETAR_M ** 2)
sc_5k = 1.0 - B_at_5kyr / m.B_CRIT_MAGNETAR_T  # ~0.97135
grav = Ug1 * (1.0 + m.H0_MAGNETAR_SI * 5000.0 * 365.25 * 86400.0) * sc_5k
Ug_sum_trz = (Ug1 + Ug1 * sc_5k) * 1.1
total_dominant = grav + Ug_sum_trz + a_macro  # cosm + GW negligible
tests.append(('v2 decomposition: grav + Ug_sum*(1+TRZ) + Lorentz*macro = total (1% tol)',
              abs(g2 - total_dominant) / g2 < 0.01))
tests.append(('v2 grav_term ~ 4.512e11 m/s^2 (spec, 0.1% tol)',
              abs(grav - 4.512e11) / 4.512e11 < 0.001))
tests.append(('v2 Ug_sum*(1+TRZ) ~ 1.007e12 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 1.007e12) / 1.007e12 < 0.005))

# Setting f_TRZ=0 should reduce Ug_sum contribution by factor 1.1
g2_no_trz = m._magnetar_g_master_uqff_v2(f_TRZ=0.0)
tests.append(('v2 f_TRZ=0 reduces total by Ug_sum*(0.1) (12-decimal)',
              abs((g2 - g2_no_trz) - (Ug_sum_trz / 1.1) * 0.1) / g2 < 1e-9))

# Setting macro_scale=0 removes Lorentz term entirely
g2_no_lorentz = m._magnetar_g_master_uqff_v2(macro_scale=0.0)
tests.append(('v2 macro_scale=0 removes Lorentz term (12-decimal)',
              abs((g2 - g2_no_lorentz) - a_macro) / g2 < 1e-9))

# v2 with f_TRZ=0 and macro_scale=0 must equal v1 (at same t)
g1_at_5kyr = m._magnetar_g_master_uqff(t_s=5000.0 * 365.25 * 86400.0)
g2_stripped = m._magnetar_g_master_uqff_v2(f_TRZ=0.0, macro_scale=0.0)
tests.append(('v2(f_TRZ=0, macro=0) = v1 at same t (12-decimal)',
              abs(g2_stripped - g1_at_5kyr) / g1_at_5kyr < 1e-12))

# v2 with rho_UA=0 -> Lorentz UA factor = 1 (no enhancement)
g2_no_ua = m._magnetar_g_master_uqff_v2(rho_UA_val=0.0)
expected_no_ua_lorentz = expected_bare * 1.0 * 1.0e-12  # factor 1, not 11
tests.append(('v2 rho_UA=0 reduces Lorentz to bare * macro (1% tol)',
              abs((g2_no_ua - g2_no_lorentz) - expected_no_ua_lorentz)
                  / expected_no_ua_lorentz < 0.01))

# Caller-overridden Ug2, Ug3 must propagate through (1+f_TRZ) multiplier
g2_with_extras = m._magnetar_g_master_uqff_v2(Ug2=1.0e10, Ug3=2.0e10)
tests.append(('v2(Ug2, Ug3) extras multiplied by (1 + f_TRZ) (12-decimal)',
              abs((g2_with_extras - g2) - 3.0e10 * 1.1) < 1.0))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
