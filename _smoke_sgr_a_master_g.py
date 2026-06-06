"""Smoke test: Sgr A* (SMBH) evolution master Universal Gravity (spec 09May2025)."""
import math
import uqff_pure_calculator as m

tests = []
YR_S = 365.25 * 86400.0
t_45gyr = 4.5e9 * YR_S

# ---- module-level constants ----
tests.append(('M_SGRA_KG = 4.3e6 * M_SUN',
              m.M_SGRA_KG == 4.3e6 * m.M_SUN))
tests.append(('M_SGRA_KG ~ 8.552e36 kg (spec, 0.1% tol)',
              abs(m.M_SGRA_KG - 8.552e36) / 8.552e36 < 0.001))
tests.append(('R_SGRA_M = 1.27e10 m (spec)',
              m.R_SGRA_M == 1.27e10))
tests.append(('B0_SGRA_T = 1.0 T (=10^4 G)',
              m.B0_SGRA_T == 1.0))
tests.append(('TAU_B_SGRA_S = 1 Myr',
              abs(m.TAU_B_SGRA_S - 1.0e6 * YR_S) < 1e-3))
tests.append(('TAU_ACC_SGRA_S = 9 Gyr',
              abs(m.TAU_ACC_SGRA_S - 9.0e9 * YR_S) < 1e-3))
tests.append(('M_DOT0_SGRA = 0.01 (1% amplitude)',
              m.M_DOT0_SGRA == 0.01))
tests.append(('OMEGA_PREFACTOR_SGRA = 0.3 (eta = 0.3)',
              m.OMEGA_PREFACTOR_SGRA == 0.3))
tests.append(('SGRA_PRECESSION_DEG = 30.0',
              m.SGRA_PRECESSION_DEG == 30.0))

# ---- _accretion_mass_growth_uqff ----
tests.append(('_accretion_mass_growth_uqff(t=0) = M_initial * (1 + M_dot_0)',
              abs(m._accretion_mass_growth_uqff(t_s=0.0)
                  - m.M_SGRA_KG * (1.0 + 0.01)) / (m.M_SGRA_KG * 1.01) < 1e-12))
M_45 = m._accretion_mass_growth_uqff(t_s=t_45gyr)
tests.append(('M(t=4.5 Gyr) ~ 8.604e36 kg (spec, 0.1% tol)',
              abs(M_45 - 8.604e36) / 8.604e36 < 0.001))
tests.append(('M(t->inf) -> M_initial (no accretion residual)',
              abs(m._accretion_mass_growth_uqff(t_s=1e25) - m.M_SGRA_KG)
                  / m.M_SGRA_KG < 1e-12))

# ---- _sgra_spin_Omega_uqff ----
Omega_0 = 0.3 * m.C_LIGHT / m.R_SGRA_M
tests.append(('_sgra_spin_Omega_uqff(t=0) = 0.3 c / r ~ 7.087e-3 rad/s (spec)',
              abs(m._sgra_spin_Omega_uqff(t_s=0.0) - Omega_0) < 1e-15))
tests.append(('Omega(t=0) ~ 7.087e-3 rad/s (spec literal, 0.1% tol)',
              abs(m._sgra_spin_Omega_uqff(t_s=0.0) - 7.087e-3) / 7.087e-3 < 0.001))
tests.append(('Omega(t=4.5 Gyr) ~ 4.299e-3 rad/s (spec, 0.5% tol)',
              abs(m._sgra_spin_Omega_uqff(t_s=t_45gyr) - 4.299e-3) / 4.299e-3 < 0.005))

# ---- _sgra_spin_dOmega_dt_uqff (analytical) ----
def dO_numeric(t, r=m.R_SGRA_M, eta=0.3, tau=9.0e9 * YR_S, h=1e10):
    return ((eta * m.C_LIGHT / r) * math.exp(-(t + h) / tau)
            - (eta * m.C_LIGHT / r) * math.exp(-(t - h) / tau)) / (2.0 * h)

tests.append(('dOmega/dt(t=0) matches central-diff (1% tol)',
              abs(m._sgra_spin_dOmega_dt_uqff(t_s=0.0) - dO_numeric(0.0))
                  / abs(dO_numeric(0.0)) < 0.01))
tests.append(('dOmega/dt(t=0) is negative (spin-down)',
              m._sgra_spin_dOmega_dt_uqff(t_s=0.0) < 0.0))
tests.append(('dOmega/dt(t=4.5 Gyr) ~ -1.512e-20 rad/s^2 (spec, 1% tol)',
              abs(m._sgra_spin_dOmega_dt_uqff(t_s=t_45gyr) - (-1.512e-20))
                  / 1.512e-20 < 0.01))

# ---- _dm_perturbation_precession_uqff ----
dm = m._dm_perturbation_precession_uqff(
        M_visible_kg=M_45, M_DM_kg=0.1 * M_45,
        delta_rho_over_rho=1e-5,
        M_central_kg=M_45, r_m=m.R_SGRA_M,
        theta_deg=30.0)
tests.append(('_dm_perturbation_precession_uqff(theta=30, Sgr A* defaults) ~ 4.076e33 kg/m (spec, 5% tol)',
              abs(dm - 4.076e33) / 4.076e33 < 0.05))
# sin(0) -> term reduces to (M_v+M_DM) * delta_rho/rho
dm0 = m._dm_perturbation_precession_uqff(
        M_visible_kg=M_45, M_DM_kg=0.1 * M_45,
        delta_rho_over_rho=1e-5,
        M_central_kg=M_45, r_m=m.R_SGRA_M,
        theta_deg=0.0)
tests.append(('_dm_perturbation_precession_uqff(theta=0) = (M_v+M_DM) * delta_rho/rho',
              abs(dm0 - (1.1 * M_45) * 1e-5) / abs((1.1 * M_45) * 1e-5) < 1e-12))
# sin(90) -> maximum precession
dm90 = m._dm_perturbation_precession_uqff(
        M_visible_kg=M_45, M_DM_kg=0.1 * M_45,
        delta_rho_over_rho=1e-5,
        M_central_kg=M_45, r_m=m.R_SGRA_M,
        theta_deg=90.0)
tests.append(('_dm_perturbation_precession_uqff(theta=90) ~ 2 * theta=30 term (since sin(90)/sin(30)=2)',
              abs((dm90 - 1.1 * M_45 * 1e-5) / (dm - 1.1 * M_45 * 1e-5) - 2.0) < 0.01))

# ---- _sgr_a_g_master_uqff (spec example) ----
g = m._sgr_a_g_master_uqff()
tests.append(('_sgr_a_g_master_uqff() defaults -> spec example 1.250e7 m/s^2 (1% tol)',
              abs(g - 1.250e7) / 1.250e7 < 0.01))

# Decomposition at t=4.5 Gyr
Ug1 = m.G_NEWTON * M_45 / (m.R_SGRA_M ** 2)
tests.append(('U_g1 = G*M(t)/r^2 ~ 3.561e6 m/s^2 (spec, 0.5% tol)',
              abs(Ug1 - 3.561e6) / 3.561e6 < 0.005))
B_45 = m._magnetar_B_decay(t_45gyr, B_0_T=1.0, tau_B_s=1.0e6 * YR_S)
tests.append(('B(t=4.5 Gyr) ~ 0 T (decayed; spec)',
              B_45 < 1e-100))
sc = 1.0 - B_45 / m.B_CRIT_MAGNETAR_T
H_factor = 1.0 + m.H0_MAGNETAR_SI * t_45gyr
tests.append(('1 + H_0 * t (4.5 Gyr) ~ 1.3101 (spec, 0.5% tol)',
              abs(H_factor - 1.3101) / 1.3101 < 0.005))
grav = Ug1 * H_factor * sc
tests.append(('grav_term ~ 4.665e6 m/s^2 (spec, 0.5% tol)',
              abs(grav - 4.665e6) / 4.665e6 < 0.005))
Ug_sum_trz = (Ug1 + Ug1 * sc) * 1.1
tests.append(('(Ug1+Ug4)*(1+f_TRZ) ~ 7.834e6 m/s^2 (spec, 0.5% tol)',
              abs(Ug_sum_trz - 7.834e6) / 7.834e6 < 0.005))

# Setting f_TRZ=0 reduces Ug_sum by factor 1.1
g_no_trz = m._sgr_a_g_master_uqff(f_TRZ=0.0)
tests.append(('Sgr A* f_TRZ=0 reduces total by (Ug_sum)*0.1',
              abs((g - g_no_trz) - (Ug1 + Ug1 * sc) * 0.1) / g < 1e-9))

# M_dot_0=0 (no accretion) -> M(t) = M_initial; total uses M_SGRA_KG instead of M_45
g_no_acc = m._sgr_a_g_master_uqff(M_dot_0=0.0)
Ug1_no_acc = m.G_NEWTON * m.M_SGRA_KG / (m.R_SGRA_M ** 2)
grav_no_acc = Ug1_no_acc * H_factor * sc
Ug_sum_trz_no_acc = (Ug1_no_acc + Ug1_no_acc * sc) * 1.1
total_no_acc_expected = grav_no_acc + Ug_sum_trz_no_acc  # cosm + GW negligible
tests.append(('Sgr A* M_dot_0=0 -> uses M_initial (1% tol)',
              abs(g_no_acc - total_no_acc_expected) / g_no_acc < 0.01))

# Composer should be monotone in M_dot_0 (more accretion -> more mass -> more g)
g_low = m._sgr_a_g_master_uqff(M_dot_0=0.001)
g_high = m._sgr_a_g_master_uqff(M_dot_0=0.1)
tests.append(('Sgr A* g increases with M_dot_0 (more accretion -> more g)',
              g_high > g > g_low > g_no_acc))

# Schwarzschild radius primitive matches spec r
r_s_from_primitive = m._l29_r_schwarzschild(m.M_SGRA_KG)
tests.append(('_l29_r_schwarzschild(M_SGRA) ~ 1.27e10 m (spec, 0.5% tol)',
              abs(r_s_from_primitive - 1.27e10) / 1.27e10 < 0.005))

passed = sum(1 for _, ok in tests if ok)
total = len(tests)
for name, ok in tests:
    if not ok:
        print('FAIL:', name)
print('PASS %d/%d' % (passed, total))
