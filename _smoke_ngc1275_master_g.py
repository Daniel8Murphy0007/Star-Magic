"""Smoke test for NGC 1275 Perseus A 'Magnetic Monster' AGN master Universal Gravity.

Spec: 'Master Universal Gravity Equation_Magnetic Monster NGC 1275 Evolution_09May2025'
Composer: _ngc1275_g_master_uqff (3-leaf additive: grav*(...)+a_fil+Lorentz)
ONE NEW primitive validated: _magnetic_tension_acceleration_uqff
REUSE WIN: _merger_progress_saturating_uqff -> F_BH(t) feedback (3rd use of saturating shape)
Spec example total: 3.160e-5 m/s^2 at t=50 Myr (Lorentz 99.99% dominant, NEW family).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _ngc1275_g_master_uqff,
    _magnetic_tension_acceleration_uqff,
    _merger_progress_saturating_uqff,         # reused for F_BH(t)
    _hubble_unified,
    _lorentz_acceleration_uqff,
    _ngc_1275_g_primitive_sat,                # existing BSFG triadic -- regression check
    _MU_0_NGC1275,
    M_NGC1275_STELLAR_KG,
    M_NGC1275_SMBH_KG,
    M_NGC1275_TOTAL_KG,
    R_NGC1275_M,
    F_0_NGC1275,
    TAU_BH_NGC1275_S,
    Z_NGC1275,
    H0_NGC1275_KMSMPC,
    B_NGC1275_T,
    V_NGC1275_FIL_M3,
    M_NGC1275_FIL_KG,
    V_HVS_NGC1275_MS,
    T_NGC1275_DEFAULT_S,
    G_NEWTON,
    M_SUN,
    EV_J,
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
print("NGC 1275 Perseus A 'Magnetic Monster' master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(_MU_0_NGC1275, 4.0 * math.pi * 1.0e-7, 1e-12),       "_MU_0_NGC1275 = 4 pi 1e-7 (SI exact)")
check(close(M_NGC1275_STELLAR_KG, 1.0e12 * M_SUN, 1e-12),        "M_NGC1275_STELLAR_KG = 10^12 M_sun")
check(close(M_NGC1275_SMBH_KG, 8.0e8 * M_SUN, 1e-12),            "M_NGC1275_SMBH_KG = 8e8 M_sun (800 M SMBH)")
check(close(M_NGC1275_TOTAL_KG, (1.0e12 + 8.0e8) * M_SUN, 1e-12),"M_NGC1275_TOTAL_KG = (10^12 + 8e8) M_sun = 1.991e42 kg")
check(close(R_NGC1275_M, 9.46e20, 1e-12),                        "R_NGC1275_M = 9.46e20 m (half 200 kly span)")
check(close(F_0_NGC1275, 0.1, 1e-12),                            "F_0_NGC1275 = 0.1 (10% SMBH feedback)")
check(close(TAU_BH_NGC1275_S, 100.0e6 * _YEAR_S_MAGNETAR, 1e-12),"TAU_BH_NGC1275_S = 100 Myr (3.156e15 s)")
check(close(Z_NGC1275, 0.0176, 1e-12),                           "Z_NGC1275 = 0.0176 (237 Mly)")
check(close(H0_NGC1275_KMSMPC, 70.0, 1e-12),                     "H0_NGC1275_KMSMPC = 70")
check(close(B_NGC1275_T, 1.0e-8, 1e-12),                         "B_NGC1275_T = 1e-8 T (weak ICM field)")
check(close(V_NGC1275_FIL_M3, 1.42e50, 1e-12),                   "V_NGC1275_FIL_M3 = 1.42e50 m^3 (spec LITERAL declared)")
check(close(M_NGC1275_FIL_KG, 1.0e6 * M_SUN, 1e-12),             "M_NGC1275_FIL_KG = 10^6 M_sun per filament")
check(close(V_HVS_NGC1275_MS, 3.0e6, 1e-12),                     "V_HVS_NGC1275_MS = 3e6 m/s (3000 km/s HVS)")
check(close(T_NGC1275_DEFAULT_S, 50.0e6 * _YEAR_S_MAGNETAR, 1e-12),"T_NGC1275_DEFAULT_S = 50 Myr (1.578e15 s)")

# --- 2. Magnetic tension NEW primitive ---
a_fil = _magnetic_tension_acceleration_uqff()
check(close(a_fil, 2.840e-9, 0.005),                             f"a_fil = {a_fil:.4e} ~ 2.840e-9 m/s^2 (spec example, 0.5% tol)")
# Decompose
P_mag = (B_NGC1275_T ** 2) / (2.0 * _MU_0_NGC1275)
check(close(P_mag, 3.978e-11, 0.005),                            f"B^2/(2 mu_0) = {P_mag:.4e} ~ 3.978e-11 Pa (0.5% tol)")
Force = P_mag * V_NGC1275_FIL_M3
check(close(Force, 5.649e39, 0.005),                             f"F = P*V_fil = {Force:.4e} ~ 5.649e39 N (0.5% tol)")
accel_raw = Force / M_NGC1275_FIL_KG
check(close(accel_raw, 2.840e3, 0.005),                          f"F/M_fil = {accel_raw:.4e} ~ 2.840e3 m/s^2 (0.5% tol)")
check(close(a_fil, accel_raw * 1e-12, 1e-12),                    "a_fil = F/M * 1e-12 exactly (macro scaling)")
# Boundary cases
check(close(_magnetic_tension_acceleration_uqff(B_T=0.0), 0.0, 1e-12),
      "a_fil(B=0) = 0")
check(close(_magnetic_tension_acceleration_uqff(V_m3=0.0), 0.0, 1e-12),
      "a_fil(V=0) = 0")
check(close(_magnetic_tension_acceleration_uqff(macro_scale=0.0), 0.0, 1e-12),
      "a_fil(macro=0) = 0")
# B^2 scaling
a_2B = _magnetic_tension_acceleration_uqff(B_T=2.0 * B_NGC1275_T)
check(close(a_2B, 4.0 * a_fil, 1e-9),                            "a_fil(2B) = 4 * a_fil (B^2 scaling)")
# Linear V scaling
a_2V = _magnetic_tension_acceleration_uqff(V_m3=2.0 * V_NGC1275_FIL_M3)
check(close(a_2V, 2.0 * a_fil, 1e-9),                            "a_fil(2V) = 2 * a_fil (linear V)")
# Inverse M scaling
a_2M = _magnetic_tension_acceleration_uqff(M_kg=2.0 * M_NGC1275_FIL_KG)
check(close(a_2M, 0.5 * a_fil, 1e-9),                            "a_fil(2M) = a_fil/2 (inverse M)")

# --- 3. Saturating feedback F_BH(t) -- REUSED ---
F_BH_t = _merger_progress_saturating_uqff(T_NGC1275_DEFAULT_S, F_0_NGC1275, TAU_BH_NGC1275_S)
expected_F = 0.1 * (1.0 - math.exp(-0.5))
check(close(F_BH_t, expected_F, 1e-6),                           f"F_BH(50 Myr) = {F_BH_t:.5f} = 0.1*(1-e^-0.5) = {expected_F:.5f}")
check(close(F_BH_t, 0.03935, 0.005),                             "F_BH(50 Myr) ~ 0.03935 (0.5% tol)")
feedback_factor = 1.0 - F_BH_t
check(close(feedback_factor, 0.96065, 0.005),                    f"(1-F_BH) = {feedback_factor:.5f} ~ 0.96065")
check(close(_merger_progress_saturating_uqff(0.0, F_0_NGC1275, TAU_BH_NGC1275_S), 0.0, 1e-12),
      "F_BH(t=0) = 0 (no feedback yet)")

# --- 4. Hubble at z=0.0176 ---
H_z_kmsMpc = _hubble_unified(T_NGC1275_DEFAULT_S, Z_NGC1275, H0_NGC1275_KMSMPC)
check(close(H_z_kmsMpc, 70.56, 0.01),                            f"H(z=0.0176) = {H_z_kmsMpc:.4f} km/s/Mpc ~ 70.56 (1% tol)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_NGC1275_DEFAULT_S
check(close(H_t, 3.609e-3, 0.01),                                f"H(z)*t = {H_t:.4e} ~ 3.609e-3 (1% tol)")
check(close(1.0 + H_t, 1.003609, 1e-5),                          "(1+H t) ~ 1.003609")

# --- 5. Gravitational term ---
Ug1 = G_NEWTON * M_NGC1275_TOTAL_KG / (R_NGC1275_M ** 2)
check(close(Ug1, 1.485e-10, 0.01),                               f"Ug1 = G*M/r^2 = {Ug1:.4e} ~ 1.485e-10 m/s^2 (1% tol)")
grav_full = Ug1 * (1.0 + H_t) * feedback_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_full, 1.576e-10, 0.01),                         f"grav*(1+H t)*(1-F)*(1+f_TRZ) = {grav_full:.4e} ~ 1.576e-10 (CARRIED FULL per directive)")

# --- 6. Lorentz term (NEW 3.160e-5 family) ---
lorentz = _lorentz_acceleration_uqff(B_T=B_NGC1275_T, v_ms=V_HVS_NGC1275_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 3.160e-5, 0.005),                           f"Lorentz q*v*B/m_p*11*1e-12 = {lorentz:.4e} ~ 3.160e-5 (NEW family, 0.5% tol)")
# Verify NEW family (distinct from prior 4 families)
check(lorentz < 1e-4,                                            "Lorentz < 1e-4 (NEW family, distinct from 1.053e-3, 1.884e-3, 1.053e-1)")
check(lorentz > 1e-5,                                            "Lorentz > 1e-5 (NEW 3.160e-5 family)")

# --- 7. Composer total = 3.160e-5 m/s^2 at t=50 Myr ---
g_total = _ngc1275_g_master_uqff()
check(close(g_total, 3.160e-5, 0.005),                           f"_ngc1275_g_master_uqff() = {g_total:.4e} m/s^2 ~ 3.160e-5 (spec example, 0.5% tol)")

# --- 8. Fidelity check: composer == grav + a_fil + Lorentz exactly ---
manual = grav_full + a_fil + lorentz
check(close(g_total, manual, 1e-9),                              "Composer == grav*(1+H t)*(1-F)*(1+f_TRZ) + a_fil + Lorentz exactly (no extra leaves)")

# --- 9. Contribution fractions per fidelity directive ---
a_fil_frac = a_fil / g_total
check(5e-5 < a_fil_frac < 1.5e-4,                                f"a_fil contributes {a_fil_frac:.2e} (~9e-5, carried per directive)")
lorentz_frac = lorentz / g_total
check(0.999 < lorentz_frac < 1.001,                              f"Lorentz contributes {lorentz_frac*100:.4f}% (~99.99% dominant)")
grav_frac = grav_full / g_total
check(grav_frac < 1e-4,                                          f"grav contributes {grav_frac:.2e} (carried per directive, even though << 1)")

# --- 10. Boundary cases ---
g_no_macro = _ngc1275_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_full, 1e-9),                        "macro_scale=0 -> a_fil=0 AND Lorentz=0; only grav remains")

g_no_B = _ngc1275_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_full, 1e-9),                            "B=0 -> a_fil=0 AND Lorentz=0; only grav remains")

g_no_V = _ngc1275_g_master_uqff(V_fil_m3=0.0)
check(close(g_no_V, grav_full + lorentz, 1e-9),                  "V_fil=0 -> a_fil=0; grav + Lorentz remain")

g_no_F = _ngc1275_g_master_uqff(F_0=0.0)
grav_no_F = Ug1 * (1.0 + H_t) * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_F, grav_no_F + a_fil + lorentz, 1e-9),          "F_0=0 -> feedback_factor=1, no SMBH feedback suppression")

g_no_trz = _ngc1275_g_master_uqff(f_TRZ=0.0)
grav_no_trz = Ug1 * (1.0 + H_t) * feedback_factor * 1.0
check(close(g_no_trz, grav_no_trz + a_fil + lorentz, 1e-9),      "f_TRZ=0 -> grav drops factor (1+f_TRZ)")

g_no_v = _ngc1275_g_master_uqff(v_hvs_ms=0.0)
check(close(g_no_v, grav_full + a_fil, 1e-9),                    "v=0 -> Lorentz=0; grav + a_fil remain")

# --- 11. Time evolution ---
g_100Myr = _ngc1275_g_master_uqff(t_s=100.0e6 * _YEAR_S_MAGNETAR)
F_100 = _merger_progress_saturating_uqff(100.0e6 * _YEAR_S_MAGNETAR, 0.1, TAU_BH_NGC1275_S)
check(F_100 > F_BH_t,                                            "F_BH(100 Myr) > F_BH(50 Myr) (feedback progresses)")
expected_F_100 = 0.1 * (1.0 - math.exp(-1.0))
check(close(F_100, expected_F_100, 1e-6),                        f"F_BH(100 Myr = tau_BH) = 0.1*(1-e^-1) = {expected_F_100:.5f}")
F_inf = _merger_progress_saturating_uqff(1e20, 0.1, TAU_BH_NGC1275_S)
check(close(F_inf, 0.1, 1e-6),                                   "F_BH(t->inf) = F_0 = 0.1 (saturates)")

# --- 12. Existing primitives UNTOUCHED ---
sat_val = _ngc_1275_g_primitive_sat()
check(sat_val is not None and sat_val > 0,                       f"_ngc_1275_g_primitive_sat() = {sat_val:.4e} unchanged (BSFG^2 triadic, different observable)")

# --- 13. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _horsehead_g_master_uqff,
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
    _westerlund2_g_master_uqff,
)
check(close(_horsehead_g_master_uqff(), 1.097e-3, 0.01),         "Regression: Horsehead composer = 1.097e-3")
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),          "Regression: Antennae composer = 1.053e-1")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),     "Regression: Bubble Nebula composer = 1.884e-3")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),           "Regression: NGC 3603 composer = 1.053e-3")
check(close(_westerlund2_g_master_uqff(), 1.053e-3, 0.01),       "Regression: Westerlund 2 composer = 1.053e-3 family")

print("=" * 78)
print(f"NGC 1275 smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
