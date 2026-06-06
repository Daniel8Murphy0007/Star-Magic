"""Smoke test for Horsehead Nebula Barnard 33 erosion evolution master Universal Gravity.

Spec: 'Master Universal Gravity Equation_Infrared view of the Horsehead Nebula Evolution_09May2025'
Composer: _horsehead_g_master_uqff (4-leaf clean form, P_rad NOT negligible per fidelity directive)
ONE NEW primitive validated: _radiation_pressure_acceleration_uqff
REUSE win: _merger_progress_saturating_uqff (Antennae) -> E(t) erosion (exact 1-exp shape)
Spec example total: 1.097e-3 m/s^2 at t=1 Myr (Lorentz 96%, P_rad 4%, grav <<).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _horsehead_g_master_uqff,
    _radiation_pressure_acceleration_uqff,
    _merger_progress_saturating_uqff,       # reused for E(t)
    _hubble_unified,
    _lorentz_acceleration_uqff,
    _horsehead_g_primitive_sat,             # existing BSFG triadic -- regression check
    _L_SUN_W_HORSEHEAD,
    M_HORSEHEAD_KG,
    R_HORSEHEAD_M,
    L_SIGMA_ORIONIS_W,
    RHO_HORSEHEAD_KGM3,
    _M_H_HORSEHEAD_KG,
    E_0_HORSEHEAD,
    TAU_ERODE_HORSEHEAD_S,
    Z_HORSEHEAD,
    H0_HORSEHEAD_KMSMPC,
    B_HORSEHEAD_T,
    V_HORSEHEAD_MS,
    T_HORSEHEAD_DEFAULT_S,
    G_NEWTON,
    M_SUN,
    EV_J,
    C_LIGHT,
    _MPC_M,
    _YEAR_S_MAGNETAR,
    _F_TRZ_DEFAULT_MAGNETAR,
    _M_PROTON_KG_MAGNETAR,
    MACROSCOPIC_SCALE_LORENTZ,
)

passed = 0
failed = 0

def check(cond, msg):
    global passed, failed
    if cond:
        passed += 1
        print(f"  PASS: {msg}")
    else:
        failed += 1
        print(f"  FAIL: {msg}")

def close(a, b, tol):
    return abs(a - b) <= tol * max(abs(b), 1e-300)

print("=" * 78)
print("Horsehead Nebula Barnard 33 master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(_L_SUN_W_HORSEHEAD, 3.826e26, 1e-12),                "_L_SUN_W_HORSEHEAD = 3.826e26 W (spec value)")
check(close(M_HORSEHEAD_KG, 120.0 * M_SUN, 1e-12),               "M_HORSEHEAD_KG = 120 M_sun = 2.387e32 kg")
check(close(R_HORSEHEAD_M, 1.182e16, 1e-12),                     "R_HORSEHEAD_M = 1.182e16 m (half 2.5 ly span)")
check(close(L_SIGMA_ORIONIS_W, 1.0e5 * 3.826e26, 1e-12),         "L_SIGMA_ORIONIS_W = 1e5 L_sun = 3.826e31 W")
check(close(RHO_HORSEHEAD_KGM3, 1.0e-21, 1e-12),                 "RHO_HORSEHEAD_KGM3 = 1e-21 kg/m^3")
check(close(_M_H_HORSEHEAD_KG, 1.67e-27, 1e-12),                 "_M_H_HORSEHEAD_KG = 1.67e-27 (spec value, NOT m_p precision)")
check(close(E_0_HORSEHEAD, 0.2, 1e-12),                          "E_0_HORSEHEAD = 0.2 (20% erosion amplitude)")
check(close(TAU_ERODE_HORSEHEAD_S, 5.0e6 * _YEAR_S_MAGNETAR, 1e-12),
      "TAU_ERODE_HORSEHEAD_S = 5 Myr (1.578e14 s)")
check(close(Z_HORSEHEAD, 0.0003, 1e-12),                         "Z_HORSEHEAD = 0.0003 (1500 ly)")
check(close(H0_HORSEHEAD_KMSMPC, 70.0, 1e-12),                   "H0_HORSEHEAD_KMSMPC = 70")
check(close(B_HORSEHEAD_T, 1.0e-5, 1e-12),                       "B_HORSEHEAD_T = 1e-5 T (nebula field)")
check(close(V_HORSEHEAD_MS, 1.0e5, 1e-12),                       "V_HORSEHEAD_MS = 1e5 m/s (gas velocity)")
check(close(T_HORSEHEAD_DEFAULT_S, 1.0e6 * _YEAR_S_MAGNETAR, 1e-12),
      "T_HORSEHEAD_DEFAULT_S = 1 Myr (3.156e13 s)")

# --- 2. Radiation pressure NEW primitive ---
P_rad = _radiation_pressure_acceleration_uqff()
check(close(P_rad, 4.347e-5, 0.005),                             f"P_rad = {P_rad:.4e} ~ 4.347e-5 (spec example, 0.5% tol)")
# Decompose
rad_pa = L_SIGMA_ORIONIS_W / (4.0 * math.pi * R_HORSEHEAD_M ** 2 * C_LIGHT)
check(close(rad_pa, 7.26e-11, 0.01),                             f"L/(4 pi r^2 c) = {rad_pa:.4e} ~ 7.26e-11 Pa (1% tol)")
n_density = RHO_HORSEHEAD_KGM3 / _M_H_HORSEHEAD_KG
check(close(n_density, 5.988e5, 0.005),                          f"rho/m_H = {n_density:.4e} ~ 5.988e5 (0.5% tol)")
check(close(P_rad, rad_pa * n_density, 1e-12),                   "P_rad = rad_pa * n_density exactly")
# Boundary cases for new primitive
check(close(_radiation_pressure_acceleration_uqff(L_star_W=0.0), 0.0, 1e-12),
      "P_rad(L=0) = 0")
check(close(_radiation_pressure_acceleration_uqff(rho_kgm3=0.0), 0.0, 1e-12),
      "P_rad(rho=0) = 0")
# Inverse r^2 scaling
P_2r = _radiation_pressure_acceleration_uqff(r_m=2.0 * R_HORSEHEAD_M)
check(close(P_2r, P_rad / 4.0, 1e-9),                            "P_rad(2r) = P_rad(r)/4 (inverse r^2)")

# --- 3. Saturating erosion E(t) -- REUSED from Antennae primitive ---
E_t = _merger_progress_saturating_uqff(T_HORSEHEAD_DEFAULT_S, E_0_HORSEHEAD, TAU_ERODE_HORSEHEAD_S)
expected_E = 0.2 * (1.0 - math.exp(-0.2))
check(close(E_t, expected_E, 1e-6),                              f"E(1 Myr) = {E_t:.5f} = 0.2*(1-e^-0.2) = {expected_E:.5f}")
check(close(E_t, 0.03626, 0.005),                                "E(1 Myr) ~ 0.03626 (0.5% tol)")
erosion_factor = 1.0 - E_t
check(close(erosion_factor, 0.96374, 0.005),                     f"(1-E) = {erosion_factor:.5f} ~ 0.96374")
# Boundary
check(close(_merger_progress_saturating_uqff(0.0, E_0_HORSEHEAD, TAU_ERODE_HORSEHEAD_S), 0.0, 1e-12),
      "E(t=0) = 0 (no erosion yet)")

# --- 4. Hubble at z=0.0003 ---
H_z_kmsMpc = _hubble_unified(T_HORSEHEAD_DEFAULT_S, Z_HORSEHEAD, H0_HORSEHEAD_KMSMPC)
check(close(H_z_kmsMpc, 70.0095, 0.005),                         f"H(z=0.0003) = {H_z_kmsMpc:.4f} km/s/Mpc ~ 70.0095 (0.5% tol)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_HORSEHEAD_DEFAULT_S
check(close(H_t, 7.161e-5, 0.01),                                f"H(z)*t = {H_t:.4e} ~ 7.161e-5 (1% tol)")
check(close(1.0 + H_t, 1.00007161, 1e-7),                        "(1+H t) ~ 1.00007161")

# --- 5. Gravitational term ---
Ug1 = G_NEWTON * M_HORSEHEAD_KG / (R_HORSEHEAD_M ** 2)
check(close(Ug1, 1.141e-10, 0.005),                              f"Ug1 = G*M/r^2 = {Ug1:.4e} ~ 1.141e-10 m/s^2 (0.5% tol)")
grav_full = Ug1 * (1.0 + H_t) * erosion_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_full, 1.21e-10, 0.01),                          f"grav*(1+H t)*(1-E)*(1+f_TRZ) = {grav_full:.4e} ~ 1.21e-10 (1% tol, CARRIED FULL PRECISION per fidelity directive)")

# --- 6. Lorentz term (1.053e-3 family -- same as NGC 3603/Westerlund) ---
lorentz = _lorentz_acceleration_uqff(B_T=B_HORSEHEAD_T, v_ms=V_HORSEHEAD_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 1.053e-3, 0.005),                           f"Lorentz q*v*B/m_p*11*1e-12 = {lorentz:.4e} ~ 1.053e-3 (1.053e-3 family, 0.5% tol)")

# --- 7. Composer total = 1.097e-3 m/s^2 at t=1 Myr ---
g_total = _horsehead_g_master_uqff()
check(close(g_total, 1.097e-3, 0.005),                           f"_horsehead_g_master_uqff() = {g_total:.4e} m/s^2 ~ 1.097e-3 (spec example, 0.5% tol)")

# --- 8. Fidelity check: composer == grav + P_rad + Lorentz exactly (3-leaf additive sum) ---
manual = grav_full + P_rad + lorentz
check(close(g_total, manual, 1e-9),                              "Composer == grav*(1+H t)*(1-E)*(1+f_TRZ) + P_rad + Lorentz exactly (no extra leaves)")

# --- 9. P_rad contribution fraction (per fidelity directive: 4% NOT negligible) ---
P_rad_frac = P_rad / g_total
check(0.035 < P_rad_frac < 0.045,                                f"P_rad contributes {P_rad_frac*100:.2f}% (~4%, NOT negligible per directive)")
lorentz_frac = lorentz / g_total
check(0.94 < lorentz_frac < 0.98,                                f"Lorentz contributes {lorentz_frac*100:.2f}% (~96% dominant)")
grav_frac = grav_full / g_total
check(grav_frac < 1e-5,                                          f"grav contributes {grav_frac:.2e} (carried per directive, even though << 1)")

# --- 10. Boundary cases ---
g_no_macro = _horsehead_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_full + P_rad, 1e-9),                "macro_scale=0 -> Lorentz=0; grav + P_rad remain")

g_no_L = _horsehead_g_master_uqff(L_star_W=0.0)
check(close(g_no_L, grav_full + lorentz, 1e-9),                  "L_star=0 -> P_rad=0; grav + Lorentz remain")

g_no_rho = _horsehead_g_master_uqff(rho_kgm3=0.0)
check(close(g_no_rho, grav_full + lorentz, 1e-9),                "rho=0 -> P_rad=0; grav + Lorentz remain")

g_no_erode = _horsehead_g_master_uqff(E_0=0.0)
grav_no_erode = Ug1 * (1.0 + H_t) * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_erode, grav_no_erode + P_rad + lorentz, 1e-9),  "E_0=0 -> erosion_factor=1, no erosion suppression")

g_no_trz = _horsehead_g_master_uqff(f_TRZ=0.0)
grav_no_trz = Ug1 * (1.0 + H_t) * erosion_factor * 1.0
check(close(g_no_trz, grav_no_trz + P_rad + lorentz, 1e-9),      "f_TRZ=0 -> grav drops factor (1+f_TRZ)")

g_no_B = _horsehead_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_full + P_rad, 1e-9),                    "B=0 -> Lorentz=0; grav + P_rad remain")

# --- 11. Time evolution ---
g_5Myr = _horsehead_g_master_uqff(t_s=5.0e6 * _YEAR_S_MAGNETAR)
E_5Myr = _merger_progress_saturating_uqff(5.0e6 * _YEAR_S_MAGNETAR, 0.2, TAU_ERODE_HORSEHEAD_S)
check(E_5Myr > E_t,                                              "E(5 Myr) > E(1 Myr) (erosion progresses)")
expected_E_5 = 0.2 * (1.0 - math.exp(-1.0))
check(close(E_5Myr, expected_E_5, 1e-6),                         f"E(5 Myr = tau_erode) = 0.2*(1-e^-1) = {expected_E_5:.5f}")
# At late times, E saturates to E_0=0.2 -> (1-E)=0.8
E_inf = _merger_progress_saturating_uqff(1e20, 0.2, TAU_ERODE_HORSEHEAD_S)
check(close(E_inf, 0.2, 1e-6),                                   "E(t->inf) = E_0 = 0.2 (saturates)")

# --- 12. Existing primitives UNTOUCHED ---
sat_val = _horsehead_g_primitive_sat()
check(sat_val is not None and sat_val > 0,                       f"_horsehead_g_primitive_sat() = {sat_val:.4e} unchanged (BSFG triadic, different observable)")

# --- 13. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
    _westerlund2_g_master_uqff,
)
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),          "Regression: Antennae composer = 1.053e-1 (NEW Lorentz family)")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),     "Regression: Bubble Nebula composer = 1.884e-3")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),           "Regression: NGC 3603 composer = 1.053e-3")
check(close(_westerlund2_g_master_uqff(), 1.053e-3, 0.01),       "Regression: Westerlund 2 composer = 1.053e-3 family")

print("=" * 78)
print(f"Horsehead smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
