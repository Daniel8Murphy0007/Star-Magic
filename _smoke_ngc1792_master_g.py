"""Smoke test for NGC 1792 'The Stellar Forge' starburst spiral master Universal Gravity.

Spec: 'Master Universal Gravity Equation for NGC 1792 Evolution_09May2025'
Composer: _ngc1792_g_master_uqff (5-leaf clean form: grav*(1+H z t)*(1+M_sf)*(1-F_sn)*(1+f_TRZ) + Lorentz)
ZERO NEW primitives -- 2ND CONSECUTIVE PURE REUSE WIN.
NEW LORENTZ FAMILY: 1.053e-2 (v*B=10, distinct from 1, 100, 0.03 families).
5TH CONSECUTIVE USE of _merger_progress_saturating_uqff.
Spec example total: 1.053e-2 m/s^2 at t=100 Myr (Lorentz dominates ~99.99999990%).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _ngc1792_g_master_uqff,
    _merger_progress_saturating_uqff,         # reused for F_sn(t) (5TH USE!)
    _sfr_linear_mass_growth_uqff,             # reused for M_sf growth fraction
    _hubble_unified,
    _lorentz_acceleration_uqff,
    M_NGC1792_KG,
    M_NGC1792_MSUN,
    R_NGC1792_M,
    SFR_NGC1792_MSUN_YR,
    F_SN_0_NGC1792,
    TAU_SN_NGC1792_S,
    Z_NGC1792,
    H0_NGC1792_KMSMPC,
    B_NGC1792_T,
    V_NGC1792_MS,
    T_NGC1792_DEFAULT_S,
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
print("NGC 1792 'The Stellar Forge' master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(M_NGC1792_KG, 1.0e10 * M_SUN, 1e-12),                  "M_NGC1792_KG = 10^10 M_sun = 1.989e40 kg")
check(close(M_NGC1792_MSUN, 1.0e10, 1e-12),                        "M_NGC1792_MSUN = 1e10 (for dimensionless ratio)")
check(close(R_NGC1792_M, 3.78e20, 1e-12),                          "R_NGC1792_M = 3.78e20 m (half 80 kly diameter)")
check(close(SFR_NGC1792_MSUN_YR, 10.0, 1e-12),                     "SFR_NGC1792_MSUN_YR = 10 M_sun/yr (10x Milky Way)")
check(close(F_SN_0_NGC1792, 0.05, 1e-12),                          "F_SN_0_NGC1792 = 0.05 (5% supernova/wind feedback)")
check(close(TAU_SN_NGC1792_S, 100.0e6 * _YEAR_S_MAGNETAR, 1e-12),  "TAU_SN_NGC1792_S = 100 Myr starburst duration")
check(close(Z_NGC1792, 0.0095, 1e-12),                             "Z_NGC1792 = 0.0095 (from 42 Mly)")
check(close(H0_NGC1792_KMSMPC, 70.0, 1e-12),                       "H0_NGC1792_KMSMPC = 70")
check(close(B_NGC1792_T, 1.0e-5, 1e-12),                           "B_NGC1792_T = 1e-5 T (galactic ISM field)")
check(close(V_NGC1792_MS, 1.0e6, 1e-12),                           "V_NGC1792_MS = 1e6 m/s (supernova wind)")
check(close(T_NGC1792_DEFAULT_S, 100.0e6 * _YEAR_S_MAGNETAR, 1e-12),
                                                                    "T_NGC1792_DEFAULT_S = 100 Myr (3.156e15 s)")

# --- 2. M_sf dimensionless starburst growth (REUSE of linear SFR shape) ---
t_yr = T_NGC1792_DEFAULT_S / _YEAR_S_MAGNETAR
M_sf_frac = SFR_NGC1792_MSUN_YR * t_yr / M_NGC1792_MSUN
check(close(M_sf_frac, 0.1, 1e-9),                                 f"M_sf(100 Myr) = SFR*t/M_0 = {M_sf_frac:.4f} = 0.1")
check(close(1.0 + M_sf_frac, 1.1, 1e-9),                           f"(1+M_sf) = {1.0+M_sf_frac:.4f} = 1.1 (growth factor)")
# Mathematical equivalence: _sfr_linear_mass_growth_uqff yields same M(t)
SFR_kgs = SFR_NGC1792_MSUN_YR * M_SUN / _YEAR_S_MAGNETAR
M_t_linear = _sfr_linear_mass_growth_uqff(M_NGC1792_KG, T_NGC1792_DEFAULT_S, SFR_kgs)
M_sf_from_linear = (M_t_linear - M_NGC1792_KG) / M_NGC1792_KG
check(close(M_sf_from_linear, M_sf_frac, 1e-9),                    "M_sf dimensionless ratio derives from _sfr_linear_mass_growth_uqff exactly")

# --- 3. F_sn(t) saturating REUSE (5TH consecutive use of saturating primitive!) ---
F_sn_t = _merger_progress_saturating_uqff(T_NGC1792_DEFAULT_S, F_SN_0_NGC1792, TAU_SN_NGC1792_S)
expected_F_sn = 0.05 * (1.0 - math.exp(-1.0))
check(close(F_sn_t, expected_F_sn, 1e-9),                          f"F_sn(100 Myr) = {F_sn_t:.6f} = 0.05*(1-e^-1) = {expected_F_sn:.6f}")
check(close(F_sn_t, 0.031606, 1e-4),                               f"F_sn(100 Myr) ~ 0.03161 (t/tau=1)")
feedback_factor = 1.0 - F_sn_t
check(close(feedback_factor, 0.968394, 1e-4),                      f"(1-F_sn) = {feedback_factor:.6f} ~ 0.96840")
check(close(_merger_progress_saturating_uqff(0.0, F_SN_0_NGC1792, TAU_SN_NGC1792_S), 0.0, 1e-12),
                                                                    "F_sn(t=0) = 0")
# Saturation behavior
F_sn_inf = _merger_progress_saturating_uqff(1.0e12 * _YEAR_S_MAGNETAR, F_SN_0_NGC1792, TAU_SN_NGC1792_S)
check(close(F_sn_inf, F_SN_0_NGC1792, 1e-9),                       f"F_sn(t->inf) = {F_sn_inf:.6f} = F_0 = 0.05 (saturated)")

# --- 4. Hubble at z=0.0095 (nearby) ---
H_z_kmsMpc = _hubble_unified(T_NGC1792_DEFAULT_S, Z_NGC1792, H0_NGC1792_KMSMPC)
expected_H = 70.0 * math.sqrt(0.3 * (1.0095 ** 3) + 0.7)
check(close(H_z_kmsMpc, expected_H, 1e-9),                         f"H(z=0.0095) = {H_z_kmsMpc:.4f} km/s/Mpc = 70*sqrt(0.3*1.0287+0.7) = {expected_H:.4f}")
check(close(H_z_kmsMpc, 70.30, 0.002),                             f"H(z=0.0095) ~ 70.30 km/s/Mpc (spec 70.301)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_NGC1792_DEFAULT_S
check(close(H_t, 7.189e-3, 0.01),                                  f"H(z)*t(100 Myr) = {H_t:.4e} ~ 7.189e-3")
expansion_factor = 1.0 + H_t
check(close(expansion_factor, 1.00719, 0.001),                     f"(1+H t) = {expansion_factor:.6f} ~ 1.00719")

# --- 5. Gravitational term ---
Ug1 = G_NEWTON * M_NGC1792_KG / (R_NGC1792_M ** 2)
check(close(Ug1, 9.293e-12, 0.005),                                f"Ug1 = G*M/r^2 = {Ug1:.4e} ~ 9.293e-12 m/s^2 (0.5% tol)")
grav_full = Ug1 * expansion_factor * (1.0 + M_sf_frac) * feedback_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_full, 1.097e-11, 0.01),                           f"grav*(1+H t)*(1+M_sf)*(1-F_sn)*(1+f_TRZ) = {grav_full:.4e} ~ 1.097e-11 (CARRIED FULL per directive)")

# --- 6. Lorentz term (NEW 1.053e-2 FAMILY -- v*B=10) ---
lorentz = _lorentz_acceleration_uqff(B_T=B_NGC1792_T, v_ms=V_NGC1792_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 1.053e-2, 0.005),                             f"Lorentz q*v*B/m_p*11*1e-12 = {lorentz:.4e} ~ 1.053e-2 (NEW 1.053e-2 FAMILY, 0.5% tol)")
# Verify v*B product = 10 (distinct from prior families)
check(close(B_NGC1792_T * V_NGC1792_MS, 10.0, 1e-12),              "B*v = 10.0 (NEW family -- distinct from 1.0 HUDF/NGC3603, 100 Antennae, 0.03 NGC 1275)")

# --- 7. Composer total = 1.053e-2 m/s^2 at t=100 Myr ---
g_total = _ngc1792_g_master_uqff()
check(close(g_total, 1.053e-2, 0.005),                             f"_ngc1792_g_master_uqff() = {g_total:.4e} m/s^2 ~ 1.053e-2 (spec example, 0.5% tol)")

# --- 8. Fidelity check: composer == grav + Lorentz exactly ---
manual = grav_full + lorentz
check(close(g_total, manual, 1e-9),                                "Composer == grav*(1+H t)*(1+M_sf)*(1-F_sn)*(1+f_TRZ) + Lorentz exactly")

# --- 9. Contribution fractions per fidelity directive ---
grav_frac = grav_full / g_total
check(grav_frac < 1e-8,                                            f"grav contributes {grav_frac:.2e} (carried per directive, even though << 1)")
lorentz_frac = lorentz / g_total
check(0.99999 < lorentz_frac < 1.00001,                            f"Lorentz contributes {lorentz_frac*100:.7f}% (~99.999999% dominant)")

# --- 10. Boundary cases ---
g_no_macro = _ngc1792_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_full, 1e-9),                          "macro_scale=0 -> Lorentz=0; only grav remains")

g_no_SFR = _ngc1792_g_master_uqff(SFR_msun_yr=0.0)
grav_no_growth = Ug1 * expansion_factor * 1.0 * feedback_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_SFR, grav_no_growth + lorentz, 1e-9),             "SFR=0 -> M_sf=0, (1+M_sf)=1, no starburst growth")

g_no_fb = _ngc1792_g_master_uqff(F_sn_0=0.0)
grav_no_fb = Ug1 * expansion_factor * (1.0 + M_sf_frac) * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_fb, grav_no_fb + lorentz, 1e-9),                  "F_sn_0=0 -> feedback_factor=1, no supernova suppression")

g_no_trz = _ngc1792_g_master_uqff(f_TRZ=0.0)
grav_no_trz = Ug1 * expansion_factor * (1.0 + M_sf_frac) * feedback_factor * 1.0
check(close(g_no_trz, grav_no_trz + lorentz, 1e-9),                "f_TRZ=0 -> grav drops factor (1+f_TRZ)")

g_no_B = _ngc1792_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_full, 1e-9),                              "B=0 -> Lorentz=0; only grav remains")

g_no_v = _ngc1792_g_master_uqff(v_wind_ms=0.0)
check(close(g_no_v, grav_full, 1e-9),                              "v=0 -> Lorentz=0; only grav remains")

# --- 11. Time evolution ---
# At t=0: all evolution factors unity except (1+f_TRZ)
g_t0 = _ngc1792_g_master_uqff(t_s=0.0)
grav_t0 = Ug1 * 1.0 * 1.0 * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_t0, grav_t0 + lorentz, 1e-9),                        "t=0 -> all evolution factors unity (no H*t, no M_sf, no F_sn)")
# At t=500 Myr (late starburst): M_sf grows to 0.5, F_sn approaches F_0
g_500Myr = _ngc1792_g_master_uqff(t_s=500.0e6 * _YEAR_S_MAGNETAR)
M_sf_500 = 10.0 * 500.0e6 / 1.0e10
check(close(M_sf_500, 0.5, 1e-9),                                  f"M_sf(500 Myr) = {M_sf_500:.4f} = 0.5 (linear in t)")
F_sn_500 = 0.05 * (1.0 - math.exp(-5.0))
check(F_sn_500 > 0.04,                                             f"F_sn(500 Myr) = {F_sn_500:.6f} approaches F_0=0.05 (e^-5=0.0067)")

# --- 12. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _hudf_g_master_uqff,
    _ngc1275_g_master_uqff,
    _horsehead_g_master_uqff,
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
)
check(close(_hudf_g_master_uqff(), 1.053e-3, 0.01),                "Regression: HUDF composer = 1.053e-3 (v*B=1 family)")
check(close(_ngc1275_g_master_uqff(), 3.160e-5, 0.01),             "Regression: NGC 1275 composer = 3.160e-5 (v*B=0.03 family)")
check(close(_horsehead_g_master_uqff(), 1.097e-3, 0.01),           "Regression: Horsehead composer = 1.097e-3 (v*B=1 family)")
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),            "Regression: Antennae composer = 1.053e-1 (v*B=100 family)")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),       "Regression: Bubble Nebula composer = 1.884e-3")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),             "Regression: NGC 3603 composer = 1.053e-3 (v*B=1 family)")

# --- 13. Lorentz family separation: NGC 1792 1.053e-2 is exactly 10x HUDF 1.053e-3 ---
ratio_to_hudf = lorentz / 1.053e-3
check(close(ratio_to_hudf, 10.0, 0.005),                           f"NGC 1792 Lorentz / HUDF Lorentz = {ratio_to_hudf:.4f} = 10 (B factor 10x)")
# And exactly 10x smaller than Antennae 1.053e-1
ratio_to_antennae = 1.053e-1 / lorentz
check(close(ratio_to_antennae, 10.0, 0.005),                       f"Antennae Lorentz / NGC 1792 Lorentz = {ratio_to_antennae:.4f} = 10 (B factor 10x)")

print("=" * 78)
print(f"NGC 1792 smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
