"""Smoke test for Crab Nebula (M1, NGC 1952) pulsar-wind-driven SNR master Universal Gravity.

Spec: 'Master Universal Gravity Equation (UQFF & SM Integration)_Crab Nebula Evolution_03May2025'
      + DeepSearch overlay 09May2025 04:00 EDT
Composer: _crab_g_master_uqff (3 ADDITIVE leaves: grav_chain + a_wind + M_mag)
ONE NEW primitive: _relativistic_wind_acceleration_uqff (pulsar wind pressure with (1+v/c) boost).
PURE-REUSE STREAK ENDS at 4 (HUDF/NGC1792/Sombrero+1/Saturn/M16-Eagle -> Crab+1).
1st SNR composer; 1st NEUTRON STAR central engine; 1st ELECTRON mass in Lorentz;
1st UA-DISABLED Lorentz configuration; 1st relativistic correction (1+v_shock/c).
Spec example total: 1.481e6 m/s^2 at t=971 yr (pulsar wind DOMINATES 99.99999%).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _crab_g_master_uqff,
    _relativistic_wind_acceleration_uqff,   # NEW primitive
    _lorentz_acceleration_uqff,             # REUSE with electron mass + UA disabled
    _hubble_unified,
    M_CRAB_KG,
    R_CRAB_M,
    P_PULSAR_CRAB_W,
    V_SHOCK_CRAB_MS,
    RHO_NEBULA_CRAB_KGM3,
    B_CRAB_T,
    V_GAS_CRAB_MS,
    Z_CRAB,
    H0_CRAB_KMSMPC,
    T_CRAB_DEFAULT_S,
    _M_ELECTRON_KG_CRAB,
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
print("Crab Nebula (M1, NGC 1952) pulsar-wind-driven SNR master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(M_CRAB_KG, 4.6 * 1.989e30, 1e-12),                     "M_CRAB_KG = 4.6 M_sun = 9.149e30 kg (1.4 pulsar + 3.2 ejecta)")
check(close(R_CRAB_M, 5.2e16, 1e-12),                              "R_CRAB_M = 5.2e16 m (5.5 ly radius)")
check(close(P_PULSAR_CRAB_W, 5.0e31, 1e-12),                       "P_PULSAR_CRAB = 5e31 W (pulsar luminosity)")
check(close(V_SHOCK_CRAB_MS, 1.5e6, 1e-12),                        "V_SHOCK_CRAB = 1.5e6 m/s (1500 km/s expansion)")
check(close(RHO_NEBULA_CRAB_KGM3, 1.0e-21, 1e-12),                 "RHO_NEBULA_CRAB = 1e-21 kg/m^3 (filament density)")
check(close(B_CRAB_T, 1.0e-8, 1e-12),                              "B_CRAB = 1e-8 T (avg nebula B field)")
check(close(V_GAS_CRAB_MS, 1.5e6, 1e-12),                          "V_GAS_CRAB = 1.5e6 m/s (shock velocity for M_mag)")
check(close(Z_CRAB, 0.0015, 1e-12),                                "Z_CRAB = 0.0015 (from 6500 ly)")
check(close(H0_CRAB_KMSMPC, 70.0, 1e-12),                          "H_0 = 70 km/s/Mpc")
check(close(T_CRAB_DEFAULT_S, 971.0 * _YEAR_S_MAGNETAR, 1e-12),    "T_CRAB default = 971 yr = 3.064e10 s (since SN 1054 AD)")
check(close(_M_ELECTRON_KG_CRAB, 9.11e-31, 1e-12),                 "m_e = 9.11e-31 kg (1st composer with electron mass in Lorentz)")

# --- 2. Gravitational base ---
Ug1 = G_NEWTON * M_CRAB_KG / (R_CRAB_M ** 2)
check(close(Ug1, 2.258e-13, 0.005),                                f"G*M/r^2 = {Ug1:.4e} ~ 2.258e-13 (spec)")

# --- 3. Hubble at z=0.0015 (same redshift family as M16-Eagle) ---
H_z_kmsMpc = _hubble_unified(T_CRAB_DEFAULT_S, Z_CRAB, H0_CRAB_KMSMPC)
check(close(H_z_kmsMpc, 70.0473, 0.001),                           f"H(z=0.0015) = {H_z_kmsMpc:.4f} ~ 70.047 km/s/Mpc (spec, same as M16-Eagle)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_CRAB_DEFAULT_S
check(close(H_t, 6.952e-8, 0.01),                                  f"H(z=0.0015)*t(971 yr) = {H_t:.4e} ~ 6.952e-8 (negligible but carried)")
expansion_factor = 1.0 + H_t
check(close(expansion_factor, 1.0000001, 1e-6),                    f"(1+H t) = {expansion_factor:.10f} ~ 1.0000001 (essentially 1)")

# --- 4. Gravitational chain ---
grav_chain = Ug1 * expansion_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_chain, 2.484e-13, 0.005),                         f"grav_chain = {grav_chain:.4e} ~ 2.484e-13 (CARRIED FULL precision per directive)")

# --- 5. NEW primitive: relativistic pulsar wind acceleration ---
# Component check: F_wind_flux = P/(4 pi r^2)
F_flux = P_PULSAR_CRAB_W / (4.0 * math.pi * R_CRAB_M * R_CRAB_M)
check(close(F_flux, 1.471e-3, 0.005),                              f"F_wind_flux = P/(4 pi r^2) = {F_flux:.4e} ~ 1.471e-3 W/m^2 (spec)")
# Relativistic boost
rel_boost = 1.0 + V_SHOCK_CRAB_MS / C_LIGHT
check(close(rel_boost, 1.005, 0.001),                              f"(1+v_shock/c) = {rel_boost:.6f} ~ 1.005 (1st relativistic correction in calculator)")
# Wind pressure
F_wind = F_flux * rel_boost
check(close(F_wind, 1.479e-3, 0.01),                               f"F_wind = F_flux*(1+v/c) = {F_wind:.4e} ~ 1.479e-3 (spec)")
# Pre-macro acceleration
a_pre = F_wind / RHO_NEBULA_CRAB_KGM3
check(close(a_pre, 1.479e18, 0.01),                                f"a_pre_macro = F_wind/rho = {a_pre:.4e} ~ 1.479e18 (spec, before macro scale)")
# Full primitive call
a_wind = _relativistic_wind_acceleration_uqff(P_pulsar_W=P_PULSAR_CRAB_W,
                                                r_m=R_CRAB_M,
                                                v_shock_ms=V_SHOCK_CRAB_MS,
                                                rho_kgm3=RHO_NEBULA_CRAB_KGM3,
                                                c_ms=C_LIGHT,
                                                macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(a_wind, 1.479e6, 0.01),                                f"a_wind primitive = {a_wind:.4e} ~ 1.479e6 m/s^2 (DOMINANT 99.99999%, spec)")

# --- 6. NEW primitive boundary cases ---
check(close(_relativistic_wind_acceleration_uqff(0, R_CRAB_M, V_SHOCK_CRAB_MS, RHO_NEBULA_CRAB_KGM3), 0.0, 1e-12),
                                                                    "Primitive: P=0 -> a_wind=0")
check(close(_relativistic_wind_acceleration_uqff(P_PULSAR_CRAB_W, 0, V_SHOCK_CRAB_MS, RHO_NEBULA_CRAB_KGM3), 0.0, 1e-12),
                                                                    "Primitive: r=0 -> a_wind=0 (safe handling)")
check(close(_relativistic_wind_acceleration_uqff(P_PULSAR_CRAB_W, R_CRAB_M, V_SHOCK_CRAB_MS, 0.0), 0.0, 1e-12),
                                                                    "Primitive: rho=0 -> a_wind=0 (safe handling)")
check(close(_relativistic_wind_acceleration_uqff(P_PULSAR_CRAB_W, R_CRAB_M, 0.0, RHO_NEBULA_CRAB_KGM3), 
            P_PULSAR_CRAB_W/(4*math.pi*R_CRAB_M*R_CRAB_M)/RHO_NEBULA_CRAB_KGM3*MACROSCOPIC_SCALE_LORENTZ, 1e-9),
                                                                    "Primitive: v_shock=0 -> (1+v/c)=1 (no rel boost)")

# --- 7. M_mag: REUSE Lorentz primitive with electron mass + UA disabled ---
M_mag = _lorentz_acceleration_uqff(B_T=B_CRAB_T, v_ms=V_GAS_CRAB_MS,
                                     q_C=EV_J, m_kg=_M_ELECTRON_KG_CRAB,
                                     rho_UA_val=0.0,             # disable UA
                                     rho_SCm_val=1.0,            # any non-zero
                                     macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(M_mag, 2.638e-3, 0.005),                               f"M_mag (Lorentz reuse, electron mass, UA disabled) = {M_mag:.4e} ~ 2.638e-3 (spec)")
check(close(B_CRAB_T * V_GAS_CRAB_MS, 0.015, 1e-9),                "B*v = 0.015 (NEW 8th Lorentz family with electron mass)")

# --- 8. Lorentz UA disable verification: check 11x absence ---
# If UA were active (rho_UA=7.09e-36, rho_SCm=7.09e-37 -> 11x), M_mag would be 2.902e-2
M_mag_with_UA = _lorentz_acceleration_uqff(B_T=B_CRAB_T, v_ms=V_GAS_CRAB_MS,
                                             q_C=EV_J, m_kg=_M_ELECTRON_KG_CRAB,
                                             rho_UA_val=None,        # default: UA active
                                             rho_SCm_val=None,
                                             macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(M_mag_with_UA, M_mag * 11.0, 0.01),                    f"With UA: {M_mag_with_UA:.4e} = 11x = {M_mag*11.0:.4e} (UA enhancement confirmed; we use disabled)")

# --- 9. Composer total ~ 1.481e6 m/s^2 at t=971 yr ---
g_total = _crab_g_master_uqff()
check(close(g_total, 1.481e6, 0.01),                               f"_crab_g_master_uqff() = {g_total:.4e} m/s^2 ~ 1.481e6 (spec)")

# --- 10. Fidelity: composer == 3 leaves summed exactly ---
manual = grav_chain + a_wind + M_mag
check(close(g_total, manual, 1e-9),                                "Composer == grav_chain + a_wind + M_mag exactly (3 additive leaves)")

# --- 11. Contribution fractions per fidelity directive (a_wind dwarfs everything) ---
a_wind_frac = a_wind / g_total
check(a_wind_frac > 0.9999999,                                     f"a_wind contributes {a_wind_frac*100:.7f}% (DOMINANT >99.99999%)")
grav_frac = grav_chain / g_total
check(grav_frac < 1e-18,                                           f"grav_chain contributes {grav_frac:.3e} (carried full precision per directive)")
M_mag_frac = M_mag / g_total
check(M_mag_frac < 1e-8,                                           f"M_mag contributes {M_mag_frac:.3e} (carried full precision per directive)")

# --- 12. Boundary cases ---
g_no_wind = _crab_g_master_uqff(P_pulsar_W=0.0)
check(close(g_no_wind, grav_chain + M_mag, 1e-9),                  "P_pulsar=0 -> a_wind=0; only grav_chain + M_mag remain")

g_no_mag = _crab_g_master_uqff(B_T=0.0)
check(close(g_no_mag, grav_chain + a_wind, 1e-9),                  "B=0 -> M_mag=0; only grav_chain + a_wind remain")

g_no_macro_wind = _crab_g_master_uqff(macro_scale_wind=0.0)
check(close(g_no_macro_wind, grav_chain + M_mag, 1e-9),            "macro_scale_wind=0 -> a_wind=0")

g_no_macro_mag = _crab_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro_mag, grav_chain + a_wind, 1e-9),            "macro_scale (mag)=0 -> M_mag=0")

g_no_trz = _crab_g_master_uqff(f_TRZ=0.0)
grav_no_trz = Ug1 * expansion_factor * 1.0
check(close(g_no_trz, grav_no_trz + a_wind + M_mag, 1e-9),         "f_TRZ=0 -> grav_chain drops (1+f_TRZ) factor")

g_no_rho = _crab_g_master_uqff(rho_nebula_kgm3=0.0)
check(close(g_no_rho, grav_chain + M_mag, 1e-9),                   "rho_nebula=0 -> a_wind=0 (primitive safe-handles)")

g_t0 = _crab_g_master_uqff(t_s=0.0)
expected_g_t0_grav = Ug1 * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_t0, expected_g_t0_grav + a_wind + M_mag, 1e-9),      "t=0 -> H*t=0; grav_chain = Ug1*1.1; a_wind & M_mag unaffected")

# --- 13. Relativistic boost amplitude verification ---
# At v_shock=1.5e6, c=3e8: (1+v/c) = 1.005, boosting wind by 0.5%
g_no_rel = _crab_g_master_uqff(v_shock_ms=0.0)
# When v_shock=0, a_wind_no_rel = P/(4 pi r^2) / rho * macro (no 1.005 boost)
a_wind_no_rel = (P_PULSAR_CRAB_W / (4 * math.pi * R_CRAB_M * R_CRAB_M)) / RHO_NEBULA_CRAB_KGM3 * MACROSCOPIC_SCALE_LORENTZ
check(close(g_no_rel, grav_chain + a_wind_no_rel + M_mag, 1e-9),   "v_shock=0 -> (1+v/c)=1; a_wind drops relativistic boost")
ratio_rel = a_wind / a_wind_no_rel
check(close(ratio_rel, 1.005, 1e-4),                               f"Relativistic boost ratio = {ratio_rel:.6f} ~ 1.005 (0.5% enhancement)")

# --- 14. Pulsar wind energetics sanity check ---
# E_pulsar over 971 yr = P_pulsar * t = 5e31 * 3.064e10 = 1.532e42 J
E_total = P_PULSAR_CRAB_W * T_CRAB_DEFAULT_S
check(close(E_total, 1.532e42, 0.01),                              f"E_pulsar(971 yr) = P*t = {E_total:.4e} J ~ 1.53e42 (total energy injected by pulsar)")

# --- 15. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _m16_eagle_g_master_uqff,
    _saturn_g_master_uqff,
    _sombrero_g_master_uqff,
    _ngc1792_g_master_uqff,
    _hudf_g_master_uqff,
    _ngc1275_g_master_uqff,
    _horsehead_g_master_uqff,
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
    _pillars_g_master_uqff,
)
check(close(_m16_eagle_g_master_uqff(), 1.053e-3, 0.01),           "Regression: M16-Eagle composer = 1.053e-3 (v*B=1 family)")
check(close(_saturn_g_master_uqff(), 10.44, 0.01),                 "Regression: Saturn composer = 10.44 m/s^2 (g_planet DOMINANT)")
check(close(_sombrero_g_master_uqff(), 5.351e-1, 0.01),            "Regression: Sombrero composer = 5.351e-1 (v*B=2 family)")
check(close(_ngc1792_g_master_uqff(), 1.053e-2, 0.01),             "Regression: NGC 1792 composer = 1.053e-2 (v*B=10 family)")
check(close(_hudf_g_master_uqff(), 1.053e-3, 0.01),                "Regression: HUDF composer = 1.053e-3 (v*B=1 family)")
check(close(_ngc1275_g_master_uqff(), 3.160e-5, 0.01),             "Regression: NGC 1275 composer = 3.160e-5 (v*B=0.03 family)")
check(close(_horsehead_g_master_uqff(), 1.097e-3, 0.01),           "Regression: Horsehead composer = 1.097e-3 (v*B=1 family)")
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),            "Regression: Antennae composer = 1.053e-1 (v*B=100 family)")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),       "Regression: Bubble Nebula composer = 1.884e-3 (v*B=1.789)")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),             "Regression: NGC 3603 composer = 1.053e-3 (v*B=1 family)")
check(close(_pillars_g_master_uqff(), 1.053e-4, 0.01),             "Regression: Pillars composer = 1.053e-4 (sub-structure, DISTINCT from M16-Eagle)")

# --- 16. NEW Lorentz family separation: Crab v*B=0.015 with m_e -- 8th distinct ---
# Saturn Lorentz used m_p; Crab uses m_e (1836x lighter -> 1836x larger acceleration per unit v*B)
saturn_lorentz_check = _lorentz_acceleration_uqff(B_T=1e-7, v_ms=500.0,
                                                    q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                                    rho_UA_val=None, rho_SCm_val=None,
                                                    macro_scale=MACROSCOPIC_SCALE_LORENTZ)
# Saturn ratio: m_p/m_e = 1836; v*B ratio Saturn/Crab = 5e-5/0.015 = 3.333e-3
# Saturn = 5.27e-8 (proton, UA on=11x), Crab = 2.638e-3 (electron, UA off=1x)
# Crab/Saturn = 2.638e-3 / 5.27e-8 = 50057 -> v*B ratio (300) * m_p/m_e (1836) / 11 (UA disable) = 300*1836/11 = 50073
ratio = M_mag / saturn_lorentz_check
expected_ratio = (0.015 / 5e-5) * (_M_PROTON_KG_MAGNETAR / _M_ELECTRON_KG_CRAB) / 11.0
check(close(ratio, expected_ratio, 0.01),                          f"Crab/Saturn Lorentz ratio = {ratio:.4e} = (v*B ratio 300)*(m_p/m_e 1836)/(UA 11) = {expected_ratio:.4e}")

print("=" * 78)
print(f"Crab Nebula smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
