"""Smoke test: magnetar evolution master universal gravity (UQFF) leaves.

Validates the 4 new closed-form primitives + the spec composer added per
spec 'Master Universal Gravity Equation (UQFF & SM Integration)_Magnetar
Evolution_03May2025'. Spec target: g_Magnetar ~ 1.386e12 m/s^2 at t=10 kyr.
"""
import math
import uqff_pure_calculator as m

tests = []

# ---- module-level constants ----
tests.append(('B_CRIT_MAGNETAR_T = 1e11 T (=1e15 G)',
              m.B_CRIT_MAGNETAR_T == 1.0e11))
tests.append(('B0_MAGNETAR_T = 1e10 T',
              m.B0_MAGNETAR_T == 1.0e10))
tests.append(('P0_MAGNETAR_S = 5.0 s',
              m.P0_MAGNETAR_S == 5.0))
tests.append(('R_MAGNETAR_M = 2e4 m',
              m.R_MAGNETAR_M == 2.0e4))
tests.append(('M_MAGNETAR_KG = 1.4 * M_SUN',
              m.M_MAGNETAR_KG == 1.4 * m.M_SUN))
tests.append(('TAU_B_MAGNETAR_S = 4000 yr',
              abs(m.TAU_B_MAGNETAR_S - 4000.0 * 365.25 * 86400.0) < 1e-3))
tests.append(('TAU_SPIN_MAGNETAR_S = 10000 yr',
              abs(m.TAU_SPIN_MAGNETAR_S - 10000.0 * 365.25 * 86400.0) < 1e-3))
tests.append(('H0_MAGNETAR_SI = 67.4 km/s/Mpc -> 2.184e-18 s^-1 (0.1% tol)',
              abs(m.H0_MAGNETAR_SI - 2.184e-18) / 2.184e-18 < 0.001))

# ---- _magnetar_B_decay ----
tests.append(('_magnetar_B_decay(t=0) = B_0',
              m._magnetar_B_decay(0.0) == m.B0_MAGNETAR_T))
tests.append(('_magnetar_B_decay(t=tau_B) = B_0 / e',
              abs(m._magnetar_B_decay(m.TAU_B_MAGNETAR_S)
                  - m.B0_MAGNETAR_T / math.e) / (m.B0_MAGNETAR_T / math.e) < 1e-12))
tests.append(('_magnetar_B_decay(t=10 kyr) ~ 8.21e8 T (spec, 0.1% tol)',
              abs(m._magnetar_B_decay(10000.0 * 365.25 * 86400.0) - 8.21e8)
                  / 8.21e8 < 0.001))

# ---- _magnetar_spin_Omega ----
tests.append(('_magnetar_spin_Omega(t=0) = 2 pi / P_0',
              abs(m._magnetar_spin_Omega(0.0) - 2.0 * math.pi / 5.0) < 1e-12))
tests.append(('_magnetar_spin_Omega(t=tau_spin) = (2 pi / P_0) / e ~ 0.462 rad/s',
              abs(m._magnetar_spin_Omega(m.TAU_SPIN_MAGNETAR_S)
                  - (2.0 * math.pi / 5.0) / math.e) < 1e-12))

# ---- _magnetar_spin_dOmega_dt analytical derivative ----
def dOmega_numeric(t, P0=5.0, tau=10000.0 * 365.25 * 86400.0, h=1.0):
    return ((2.0 * math.pi / P0) * math.exp(-(t + h) / tau)
            - (2.0 * math.pi / P0) * math.exp(-(t - h) / tau)) / (2.0 * h)

tests.append(('_magnetar_spin_dOmega_dt(t=0) matches central-diff (1% tol)',
              abs(m._magnetar_spin_dOmega_dt(0.0) - dOmega_numeric(0.0))
                  / abs(dOmega_numeric(0.0)) < 0.01))
tests.append(('_magnetar_spin_dOmega_dt(t=tau_spin) matches central-diff (1% tol)',
              abs(m._magnetar_spin_dOmega_dt(m.TAU_SPIN_MAGNETAR_S)
                  - dOmega_numeric(m.TAU_SPIN_MAGNETAR_S))
                  / abs(dOmega_numeric(m.TAU_SPIN_MAGNETAR_S)) < 0.01))
tests.append(('_magnetar_spin_dOmega_dt(t=0) is negative (spin-down)',
              m._magnetar_spin_dOmega_dt(0.0) < 0.0))

# ---- _gw_quadrupole_spin_term ----
tests.append(('_gw_quadrupole_spin_term(dOmega=0) = 0',
              m._gw_quadrupole_spin_term(dOmega_dt=0.0) == 0.0))
# G*M^2/(c^4*r) with defaults
prefactor = m.G_NEWTON * m.M_MAGNETAR_KG * m.M_MAGNETAR_KG \
            / ((m.C_LIGHT ** 4) * m.R_MAGNETAR_M)
tests.append(('_gw_quadrupole_spin_term(dOmega=1) = G*M^2/(c^4*r) (12-decimal)',
              abs(m._gw_quadrupole_spin_term(dOmega_dt=1.0) - prefactor)
                  / prefactor < 1e-12))
tests.append(('_gw_quadrupole_spin_term scales as (dOmega/dt)^2',
              abs(m._gw_quadrupole_spin_term(dOmega_dt=2.0)
                  - 4.0 * m._gw_quadrupole_spin_term(dOmega_dt=1.0)) < 1e-30))

# ---- _magnetar_g_master_uqff composer (spec target = 1.386e12 m/s^2) ----
g = m._magnetar_g_master_uqff()
tests.append(('_magnetar_g_master_uqff() defaults -> spec example 1.386e12 m/s^2 (0.5% tol)',
              abs(g - 1.386e12) / 1.386e12 < 0.005))
# Decomposition: U_g1 (G*M/r^2) ~ 4.645e11
Ug1 = m.G_NEWTON * m.M_MAGNETAR_KG / (m.R_MAGNETAR_M ** 2)
tests.append(('Spec U_g1 = G*M/r^2 ~ 4.645e11 m/s^2 (0.1% tol)',
              abs(Ug1 - 4.645e11) / 4.645e11 < 0.001))
# At t=0: B=B_0, 1-B/Bcrit = 0.9, total = Ug1*0.9 + Ug1 + 0 + 0 + Ug1*0.9 + neg = Ug1*2.8
g_t0 = m._magnetar_g_master_uqff(t_s=0.0)
tests.append(('_magnetar_g_master_uqff(t=0) ~ 2.8 * U_g1 (0.1% tol)',
              abs(g_t0 - 2.8 * Ug1) / (2.8 * Ug1) < 0.001))
# Composer is monotone in t over [0, 5*tau_B] (B decays -> sc_factor -> 1)
g_late = m._magnetar_g_master_uqff(t_s=5.0 * m.TAU_B_MAGNETAR_S)
tests.append(('_magnetar_g_master_uqff(t=5*tau_B) > _magnetar_g_master_uqff(t=0) (B decays)',
              g_late > g_t0))
# Caller can override U_g2, U_g3 (additive)
g_with_extras = m._magnetar_g_master_uqff(Ug2=1.0e10, Ug3=2.0e10)
tests.append(('_magnetar_g_master_uqff(Ug2,Ug3) is additive (12-decimal)',
              abs((g_with_extras - g) - 3.0e10) < 1.0))
# Lambda*c^2/3 contribution magnitude
cosm = m.LAMBDA_MAGNETAR_M2 * (m.C_LIGHT ** 2) / 3.0
tests.append(('Lambda*c^2/3 ~ 3.3e-36 m/s^2 (spec, 1% tol)',
              abs(cosm - 3.3e-36) / 3.3e-36 < 0.01))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
