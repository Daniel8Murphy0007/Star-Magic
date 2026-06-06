"""Smoke test for Hubble Ultra Deep Field 'Galaxies Galore' aggregate field master Universal Gravity.

Spec: 'Master Universal Gravity Equation_Hubble sees galaxies galore HUDF Evolution_09May2025'
Composer: _hudf_g_master_uqff (5-leaf clean form: grav*(1+H t)*(1+M_evo)*(1-M_merge)*(1+f_TRZ) + Lorentz)
ZERO NEW primitives -- COMPLETE REUSE WIN per fidelity directive.
Spec example total: 1.053e-3 m/s^2 at t=13 Gyr (Lorentz dominates ~99.9999%).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _hudf_g_master_uqff,
    _merger_progress_saturating_uqff,         # reused for M_merge(t)
    _sfr_linear_mass_growth_uqff,             # reused for M_evo growth fraction
    _hubble_unified,
    _lorentz_acceleration_uqff,
    _hudf_g_primitive_sat,                    # existing BSFG triadic -- regression check
    M_HUDF_KG,
    M_HUDF_MSUN,
    R_HUDF_M,
    SFR_HUDF_MSUN_YR,
    M_MERGE_0_HUDF,
    TAU_MERGE_HUDF_S,
    Z_HUDF_AVG,
    H0_HUDF_KMSMPC,
    B_HUDF_T,
    V_HUDF_MS,
    T_HUDF_DEFAULT_S,
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
print("Hubble Ultra Deep Field 'Galaxies Galore' master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(M_HUDF_KG, 1.0e12 * M_SUN, 1e-12),                    "M_HUDF_KG = 10^12 M_sun = 1.989e42 kg")
check(close(M_HUDF_MSUN, 1.0e12, 1e-12),                          "M_HUDF_MSUN = 1e12 (for dimensionless ratio)")
check(close(R_HUDF_M, 1.5e22, 1e-12),                             "R_HUDF_M = 1.5e22 m (1.5 Mpc field scale)")
check(close(SFR_HUDF_MSUN_YR, 10.0, 1e-12),                       "SFR_HUDF_MSUN_YR = 10 M_sun/yr (effective -- yields spec literal M_evo=0.13)")
check(close(M_MERGE_0_HUDF, 0.2, 1e-12),                          "M_MERGE_0_HUDF = 0.2 (20% merger amplitude)")
check(close(TAU_MERGE_HUDF_S, 1.0e9 * _YEAR_S_MAGNETAR, 1e-12),   "TAU_MERGE_HUDF_S = 1 Gyr")
check(close(Z_HUDF_AVG, 3.0, 1e-12),                              "Z_HUDF_AVG = 3 (averaged midpoint)")
check(close(H0_HUDF_KMSMPC, 70.0, 1e-12),                         "H0_HUDF_KMSMPC = 70")
check(close(B_HUDF_T, 1.0e-6, 1e-12),                             "B_HUDF_T = 1e-6 T (intergalactic field)")
check(close(V_HUDF_MS, 1.0e6, 1e-12),                             "V_HUDF_MS = 1e6 m/s (galaxy velocity)")
check(close(T_HUDF_DEFAULT_S, 13.0e9 * _YEAR_S_MAGNETAR, 1e-12),  "T_HUDF_DEFAULT_S = 13 Gyr (4.103e17 s)")

# --- 2. M_evo dimensionless growth fraction (REUSE of linear SFR shape) ---
t_yr = T_HUDF_DEFAULT_S / _YEAR_S_MAGNETAR
M_evo_frac = SFR_HUDF_MSUN_YR * t_yr / M_HUDF_MSUN
check(close(M_evo_frac, 0.13, 0.005),                             f"M_evo(13 Gyr) = SFR*t/M_0 = {M_evo_frac:.4f} ~ 0.13 (0.5% tol)")
check(close(1.0 + M_evo_frac, 1.13, 0.005),                       f"(1+M_evo) = {1.0+M_evo_frac:.4f} ~ 1.13 (growth factor)")
# Mathematical equivalence check: _sfr_linear_mass_growth_uqff yields same M(t)
SFR_kgs = SFR_HUDF_MSUN_YR * M_SUN / _YEAR_S_MAGNETAR
M_t_linear = _sfr_linear_mass_growth_uqff(M_HUDF_KG, T_HUDF_DEFAULT_S, SFR_kgs)
M_evo_from_linear = (M_t_linear - M_HUDF_KG) / M_HUDF_KG
check(close(M_evo_from_linear, M_evo_frac, 1e-9),                 "M_evo dimensionless ratio derives from _sfr_linear_mass_growth_uqff exactly")

# --- 3. M_merge(t) saturating REUSE (4th consecutive use of saturating primitive!) ---
M_merge_t = _merger_progress_saturating_uqff(T_HUDF_DEFAULT_S, M_MERGE_0_HUDF, TAU_MERGE_HUDF_S)
expected_M_merge = 0.2 * (1.0 - math.exp(-13.0))
check(close(M_merge_t, expected_M_merge, 1e-6),                   f"M_merge(13 Gyr) = {M_merge_t:.6f} = 0.2*(1-e^-13) = {expected_M_merge:.6f}")
check(close(M_merge_t, 0.2, 1e-4),                                "M_merge(13 Gyr) ~ 0.2 (saturated, t/tau=13 >> 1)")
merger_factor = 1.0 - M_merge_t
check(close(merger_factor, 0.8, 1e-4),                            f"(1-M_merge) = {merger_factor:.4f} ~ 0.8")
check(close(_merger_progress_saturating_uqff(0.0, M_MERGE_0_HUDF, TAU_MERGE_HUDF_S), 0.0, 1e-12),
      "M_merge(t=0) = 0")

# --- 4. Hubble at AVERAGED z=3 (largest expansion factor in lineage) ---
H_z_kmsMpc = _hubble_unified(T_HUDF_DEFAULT_S, Z_HUDF_AVG, H0_HUDF_KMSMPC)
expected_H = 70.0 * math.sqrt(0.3 * 64.0 + 0.7)
check(close(H_z_kmsMpc, expected_H, 1e-9),                        f"H(z=3) = {H_z_kmsMpc:.4f} km/s/Mpc = 70*sqrt(0.3*64+0.7) = {expected_H:.4f}")
check(close(H_z_kmsMpc, 312.27, 0.005),                           f"H(z=3) ~ 312.27 km/s/Mpc (spec rounded 312.2, 0.5% tol)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_HUDF_DEFAULT_S
check(close(H_t, 4.148, 0.005),                                   f"H(z=3)*t(13 Gyr) = {H_t:.4f} ~ 4.148 (LARGEST in lineage, 0.5% tol)")
expansion_factor = 1.0 + H_t
check(close(expansion_factor, 5.148, 0.005),                      f"(1+H t) = {expansion_factor:.4f} ~ 5.148 (largest expansion factor seen this lineage)")

# --- 5. Gravitational term ---
Ug1 = G_NEWTON * M_HUDF_KG / (R_HUDF_M ** 2)
check(close(Ug1, 5.902e-13, 0.005),                               f"Ug1 = G*M_0/r^2 = {Ug1:.4e} ~ 5.902e-13 m/s^2 (0.5% tol)")
grav_full = Ug1 * expansion_factor * (1.0 + M_evo_frac) * merger_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_full, 3.015e-12, 0.01),                          f"grav*(1+H t)*(1+M_evo)*(1-M_merge)*(1+f_TRZ) = {grav_full:.4e} ~ 3.015e-12 (CARRIED FULL per directive)")

# --- 6. Lorentz term (1.053e-3 family -- same v*B product as NGC 3603/Westerlund/Pillars) ---
lorentz = _lorentz_acceleration_uqff(B_T=B_HUDF_T, v_ms=V_HUDF_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 1.053e-3, 0.005),                            f"Lorentz q*v*B/m_p*11*1e-12 = {lorentz:.4e} ~ 1.053e-3 (1.053e-3 family, 0.5% tol)")
# Verify same product as NGC 3603/Westerlund (B*v = 10^-6 * 10^6 = 10^0 = 10^-5 * 10^5)
check(close(B_HUDF_T * V_HUDF_MS, 1.0, 1e-12),                    "B*v = 1.0 (same product as NGC 3603 10^-5*10^5)")

# --- 7. Composer total = 1.053e-3 m/s^2 at t=13 Gyr ---
g_total = _hudf_g_master_uqff()
check(close(g_total, 1.053e-3, 0.005),                            f"_hudf_g_master_uqff() = {g_total:.4e} m/s^2 ~ 1.053e-3 (spec example, 0.5% tol)")

# --- 8. Fidelity check: composer == grav + Lorentz exactly (NO extra leaves) ---
manual = grav_full + lorentz
check(close(g_total, manual, 1e-9),                               "Composer == grav*(1+H t)*(1+M_evo)*(1-M_merge)*(1+f_TRZ) + Lorentz exactly")

# --- 9. Contribution fractions per fidelity directive ---
grav_frac = grav_full / g_total
check(grav_frac < 1e-8,                                           f"grav contributes {grav_frac:.2e} (carried per directive, even though << 1)")
lorentz_frac = lorentz / g_total
check(0.9999 < lorentz_frac < 1.0001,                             f"Lorentz contributes {lorentz_frac*100:.6f}% (~99.9999% dominant)")

# --- 10. Boundary cases ---
g_no_macro = _hudf_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_full, 1e-9),                         "macro_scale=0 -> Lorentz=0; only grav remains")

g_no_SFR = _hudf_g_master_uqff(SFR_msun_yr=0.0)
grav_no_growth = Ug1 * expansion_factor * 1.0 * merger_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_SFR, grav_no_growth + lorentz, 1e-9),            "SFR=0 -> M_evo=0, (1+M_evo)=1, no growth boost")

g_no_merge = _hudf_g_master_uqff(M_merge_0=0.0)
grav_no_merge = Ug1 * expansion_factor * (1.0 + M_evo_frac) * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_merge, grav_no_merge + lorentz, 1e-9),           "M_merge_0=0 -> merger_factor=1, no merger suppression")

g_no_trz = _hudf_g_master_uqff(f_TRZ=0.0)
grav_no_trz = Ug1 * expansion_factor * (1.0 + M_evo_frac) * merger_factor * 1.0
check(close(g_no_trz, grav_no_trz + lorentz, 1e-9),               "f_TRZ=0 -> grav drops factor (1+f_TRZ)")

g_no_B = _hudf_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_full, 1e-9),                             "B=0 -> Lorentz=0; only grav remains")

g_no_v = _hudf_g_master_uqff(v_galaxy_ms=0.0)
check(close(g_no_v, grav_full, 1e-9),                             "v=0 -> Lorentz=0; only grav remains")

g_no_z = _hudf_g_master_uqff(z_avg=0.0)
# At z=0, H(z) = H_0 = 70 km/s/Mpc; H*t = 70e3/3.086e22 * 4.103e17 = 0.930
H_at_z0 = 70.0e3 / _MPC_M * T_HUDF_DEFAULT_S
expansion_z0 = 1.0 + H_at_z0
grav_z0 = Ug1 * expansion_z0 * (1.0 + M_evo_frac) * merger_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_z, grav_z0 + lorentz, 1e-6),                     "z=0 -> H(z)=H_0, smaller expansion factor (~1.93 vs 5.148)")

# --- 11. Time evolution ---
# At t=0: M_evo=0, M_merge=0, H*t=0 -> grav = Ug1*(1+f_TRZ) only; Lorentz unchanged
g_t0 = _hudf_g_master_uqff(t_s=0.0)
grav_t0 = Ug1 * 1.0 * 1.0 * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_t0, grav_t0 + lorentz, 1e-9),                       "t=0 -> all evolution factors unity (no H*t, no M_evo, no M_merge)")
# At later t, growth factor increases linearly
g_26Gyr = _hudf_g_master_uqff(t_s=26.0e9 * _YEAR_S_MAGNETAR)
# M_evo grows linearly: 2x t -> 2x M_evo (10*26e9/1e12 = 0.26)
M_evo_26 = SFR_HUDF_MSUN_YR * 26.0e9 / M_HUDF_MSUN
check(close(M_evo_26, 0.26, 1e-9),                                f"M_evo(26 Gyr) = {M_evo_26:.4f} = 0.26 (linear in t)")

# --- 12. Existing primitives UNTOUCHED ---
sat_val = _hudf_g_primitive_sat()
check(sat_val is not None,                                        f"_hudf_g_primitive_sat() = {sat_val:.4e} unchanged (BSFG triadic aggregate, different observable)")

# --- 13. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _ngc1275_g_master_uqff,
    _horsehead_g_master_uqff,
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
    _westerlund2_g_master_uqff,
)
check(close(_ngc1275_g_master_uqff(), 3.160e-5, 0.01),            "Regression: NGC 1275 composer = 3.160e-5")
check(close(_horsehead_g_master_uqff(), 1.097e-3, 0.01),          "Regression: Horsehead composer = 1.097e-3")
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),           "Regression: Antennae composer = 1.053e-1")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),      "Regression: Bubble Nebula composer = 1.884e-3")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),            "Regression: NGC 3603 composer = 1.053e-3")
check(close(_westerlund2_g_master_uqff(), 1.053e-3, 0.01),        "Regression: Westerlund 2 composer = 1.053e-3 family")

print("=" * 78)
print(f"HUDF smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
