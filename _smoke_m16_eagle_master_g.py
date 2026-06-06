"""Smoke test for M16 Eagle Nebula entire 70 ly region master Universal Gravity.

Spec: 'Master Universal Gravity Equation for M16 Evolution_09May2025' (03:40 EDT)
Composer: _m16_eagle_g_master_uqff (2 ADDITIVE leaves: grav_chain + Lorentz)
ZERO NEW primitives -- 4TH CONSECUTIVE PURE REUSE WIN.
SATURATING PRIMITIVE RETURNS after Sombrero+Saturn 2-composer non-use streak.
6th v*B=1 Lorentz family composer (joins HUDF/NGC3603/Westerlund/Horsehead).
DISTINCT from _pillars_g_master_uqff (sub-pillars vs whole-region scale).
Spec example total: 1.053e-3 m/s^2 at t=5 Myr (Lorentz dominates 99.9999%).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _m16_eagle_g_master_uqff,
    _merger_progress_saturating_uqff,       # REUSE (returns after 2-composer streak)
    _hubble_unified,
    _lorentz_acceleration_uqff,
    M_M16_EAGLE_KG,
    R_M16_EAGLE_M,
    SFR_M16_MSUN_PER_YR,
    M_SF_FACTOR_M16,
    E_0_M16_EAGLE,
    TAU_ERODE_M16_EAGLE_S,
    Z_M16_EAGLE,
    H0_M16_EAGLE_KMSMPC,
    B_M16_EAGLE_T,
    V_GAS_M16_EAGLE_MS,
    T_M16_EAGLE_DEFAULT_S,
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
print("M16 Eagle Nebula entire 70 ly region master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(M_M16_EAGLE_KG, 1200.0 * 1.989e30, 1e-12),             "M_M16_EAGLE_KG = 1200 M_sun = 2.387e33 kg (spec)")
check(close(R_M16_EAGLE_M, 3.31e17, 1e-12),                        "R_M16_EAGLE_M = 3.31e17 m (35 ly half-span)")
check(close(SFR_M16_MSUN_PER_YR, 1.0, 1e-12),                      "SFR_M16 = 1 M_sun/yr per pillar (spec)")
check(close(M_SF_FACTOR_M16, 4.472, 1e-12),                        "M_SF_FACTOR_M16 = 4.472 (spec LITERAL declared)")
check(close(E_0_M16_EAGLE, 0.3, 1e-12),                            "E_0 = 0.3 (30% peak erosion amplitude)")
check(close(TAU_ERODE_M16_EAGLE_S, 3.0e6 * _YEAR_S_MAGNETAR, 1e-12), "TAU_ERODE_M16 = 3 Myr = 9.468e13 s (spec)")
check(close(Z_M16_EAGLE, 0.0015, 1e-12),                           "Z_M16 = 0.0015 (from 6500 ly distance)")
check(close(H0_M16_EAGLE_KMSMPC, 70.0, 1e-12),                     "H_0 = 70 km/s/Mpc")
check(close(B_M16_EAGLE_T, 1.0e-5, 1e-12),                         "B_M16 = 10^-5 T (nebular B field)")
check(close(V_GAS_M16_EAGLE_MS, 1.0e5, 1e-12),                     "V_GAS_M16 = 10^5 m/s")
check(close(T_M16_EAGLE_DEFAULT_S, 5.0e6 * _YEAR_S_MAGNETAR, 1e-12),"T_M16 default = 5 Myr = 1.578e14 s (young stars age)")

# --- 2. Gravitational base ---
Ug1 = G_NEWTON * M_M16_EAGLE_KG / (R_M16_EAGLE_M ** 2)
check(close(Ug1, 1.454e-12, 0.01),                                 f"G*M/r^2 = {Ug1:.4e} ~ 1.454e-12 (spec)")

# --- 3. Hubble at z=0.0015 (NEW redshift -- smallest non-zero z so far) ---
H_z_kmsMpc = _hubble_unified(T_M16_EAGLE_DEFAULT_S, Z_M16_EAGLE, H0_M16_EAGLE_KMSMPC)
check(close(H_z_kmsMpc, 70.0473, 0.001),                           f"H(z=0.0015) = {H_z_kmsMpc:.4f} ~ 70.047 km/s/Mpc (spec)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_M16_EAGLE_DEFAULT_S
check(close(H_t, 3.581e-4, 0.01),                                  f"H(z=0.0015)*t(5 Myr) = {H_t:.4e} ~ 3.581e-4 (spec)")
expansion_factor = 1.0 + H_t
check(close(expansion_factor, 1.0003581, 1e-6),                    f"(1+H t) = {expansion_factor:.7f} ~ 1.0003581")

# --- 4. M_sf factor (spec literal -- arithmetic chain has unit-cancel inconsistency) ---
check(close(M_SF_FACTOR_M16, 4.472, 1e-9),                         "(1+M_sf) = 4.472 (spec LITERAL declared per fidelity hierarchy)")
# Spec procedure: (SFR_Msun_per_yr * t_yr / M_0_Msun) / M_0_Msun + 1
t_yr = T_M16_EAGLE_DEFAULT_S / _YEAR_S_MAGNETAR
M_0_Msun = M_M16_EAGLE_KG / M_SUN
spec_M_sf_step1 = (SFR_M16_MSUN_PER_YR * t_yr) / M_0_Msun           # = 4167 (spec's M_sf intermediate)
spec_M_sf_factor = 1.0 + spec_M_sf_step1 / M_0_Msun                 # = 4.472 (spec's final factor)
check(close(spec_M_sf_step1, 4166.67, 0.001),                      f"Spec M_sf step1 = SFR*t/M_0_Msun = {spec_M_sf_step1:.2f} ~ 4167")
check(close(spec_M_sf_factor, 4.472, 0.001),                       f"Spec procedure (1 + step1/M_0_Msun) = {spec_M_sf_factor:.4f} ~ 4.472 (reproduces spec chain)")

# --- 5. Erosion (saturating REUSE -- returns after Sombrero+Saturn 2-composer non-use streak) ---
E_rad = _merger_progress_saturating_uqff(T_M16_EAGLE_DEFAULT_S, E_0_M16_EAGLE, TAU_ERODE_M16_EAGLE_S)
# Expected: E_0 * (1 - exp(-t/tau)) = 0.3 * (1 - exp(-5/3))
expected_E = 0.3 * (1.0 - math.exp(-5.0/3.0))
check(close(E_rad, expected_E, 1e-9),                              f"E_rad(5 Myr) = {E_rad:.4f} = 0.3*(1-exp(-5/3)) (saturating primitive REUSE)")
check(close(E_rad, 0.2433, 0.005),                                 f"E_rad(5 Myr) = {E_rad:.4f} ~ 0.2433 (spec)")
erosion_factor = 1.0 - E_rad
check(close(erosion_factor, 0.7567, 0.005),                        f"(1-E_rad) = {erosion_factor:.4f} ~ 0.7567 (spec)")

# --- 6. Gravitational chain ---
grav_chain = Ug1 * expansion_factor * M_SF_FACTOR_M16 * erosion_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_chain, 5.413e-12, 0.01),                          f"grav_chain = {grav_chain:.4e} ~ 5.413e-12 (CARRIED FULL precision per directive)")

# --- 7. Lorentz term (v*B=1 family 6th composer REUSE) ---
lorentz = _lorentz_acceleration_uqff(B_T=B_M16_EAGLE_T, v_ms=V_GAS_M16_EAGLE_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 1.053e-3, 0.005),                             f"Lorentz q*v*B/m_p*11*1e-12 = {lorentz:.4e} ~ 1.053e-3 (v*B=1 family DOMINANT)")
check(close(B_M16_EAGLE_T * V_GAS_M16_EAGLE_MS, 1.0, 1e-9),        "B*v = 1.0 (v*B=1 family -- HUDF/NGC3603/Westerlund/Horsehead)")

# --- 8. Composer total ~ 1.053e-3 m/s^2 at t=5 Myr ---
g_total = _m16_eagle_g_master_uqff()
check(close(g_total, 1.053e-3, 0.005),                             f"_m16_eagle_g_master_uqff() = {g_total:.4e} m/s^2 ~ 1.053e-3 (spec)")

# --- 9. Fidelity: composer == 2 leaves summed exactly ---
manual = grav_chain + lorentz
check(close(g_total, manual, 1e-9),                                "Composer == grav_chain + Lorentz exactly (2 additive leaves)")

# --- 10. Contribution fractions per fidelity directive (Lorentz dwarfs grav_chain) ---
lorentz_frac = lorentz / g_total
check(lorentz_frac > 0.999999,                                     f"Lorentz contributes {lorentz_frac*100:.7f}% (DOMINANT >99.9999%)")
grav_frac = grav_chain / g_total
check(grav_frac < 1e-5,                                            f"grav_chain contributes {grav_frac:.3e} (carried full precision per directive)")

# --- 11. Boundary cases ---
g_no_macro = _m16_eagle_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_chain, 1e-9),                         "macro_scale=0 -> Lorentz=0; only grav_chain remains")

g_no_B = _m16_eagle_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_chain, 1e-9),                             "B=0 -> Lorentz=0; only grav_chain remains")

g_no_trz = _m16_eagle_g_master_uqff(f_TRZ=0.0)
grav_chain_no_trz = Ug1 * expansion_factor * M_SF_FACTOR_M16 * erosion_factor * 1.0
check(close(g_no_trz, grav_chain_no_trz + lorentz, 1e-9),          "f_TRZ=0 -> grav_chain drops (1+f_TRZ) factor")

g_t0 = _m16_eagle_g_master_uqff(t_s=0.0)
# At t=0: H*t=0 so (1+H t)=1; E_rad(0)=0 so (1-E)=1; grav = Ug1 * 1 * 4.472 * 1 * 1.1
expected_g_t0_grav = Ug1 * 1.0 * M_SF_FACTOR_M16 * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_t0, expected_g_t0_grav + lorentz, 1e-9),             "t=0 -> H*t=0 & E_rad=0; grav_chain = Ug1 * 4.472 * 1.1")

g_no_E = _m16_eagle_g_master_uqff(E_0=0.0)
grav_no_E = Ug1 * expansion_factor * M_SF_FACTOR_M16 * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_E, grav_no_E + lorentz, 1e-9),                    "E_0=0 -> no erosion; (1-E)=1")

g_no_sf = _m16_eagle_g_master_uqff(m_sf_factor=1.0)
grav_no_sf = Ug1 * expansion_factor * 1.0 * erosion_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_sf, grav_no_sf + lorentz, 1e-9),                  "m_sf_factor=1 -> no star formation growth")

# --- 12. Erosion saturation behavior (saturating primitive REUSE validation) ---
E_at_tau = _merger_progress_saturating_uqff(TAU_ERODE_M16_EAGLE_S, E_0_M16_EAGLE, TAU_ERODE_M16_EAGLE_S)
# At t=tau: E = E_0*(1-e^-1) = E_0*0.6321
check(close(E_at_tau, 0.3 * (1.0 - math.exp(-1.0)), 1e-9),         f"E_rad(t=tau) = E_0*(1-e^-1) = {E_at_tau:.4f} (saturating shape verified)")

E_long = _merger_progress_saturating_uqff(100*TAU_ERODE_M16_EAGLE_S, E_0_M16_EAGLE, TAU_ERODE_M16_EAGLE_S)
check(close(E_long, E_0_M16_EAGLE, 1e-12),                         f"E_rad(t->inf) -> E_0 = {E_long:.4f} (asymptote)")

# --- 13. DISTINCT from Pillars composer ---
from uqff_pure_calculator import _pillars_g_master_uqff, R_PILLARS_M, M_PILLARS_INIT_KG, B_PILLARS_T
check(R_PILLARS_M != R_M16_EAGLE_M,                                f"M16-Eagle r={R_M16_EAGLE_M:.3e} != Pillars r={R_PILLARS_M:.3e} (distinct scales)")
check(M_PILLARS_INIT_KG != M_M16_EAGLE_KG,                         f"M16-Eagle M={M_M16_EAGLE_KG:.3e} != Pillars M={M_PILLARS_INIT_KG:.3e} (distinct mass)")
check(B_PILLARS_T != B_M16_EAGLE_T,                                f"M16-Eagle B={B_M16_EAGLE_T:.3e} != Pillars B={B_PILLARS_T:.3e} (distinct B fields: 1e-5 vs 1e-6)")
g_pillars = _pillars_g_master_uqff()
g_m16 = g_total
check(g_pillars != g_m16,                                          f"Pillars composer ({g_pillars:.4e}) != M16-Eagle composer ({g_m16:.4e}) -- DISTINCT astrophysical objects")

# --- 14. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _saturn_g_master_uqff,
    _sombrero_g_master_uqff,
    _ngc1792_g_master_uqff,
    _hudf_g_master_uqff,
    _ngc1275_g_master_uqff,
    _horsehead_g_master_uqff,
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
)
check(close(_saturn_g_master_uqff(), 10.44, 0.01),                 "Regression: Saturn composer = 10.44 m/s^2 (g_planet DOMINANT)")
check(close(_sombrero_g_master_uqff(), 5.351e-1, 0.01),            "Regression: Sombrero composer = 5.351e-1 (v*B=2 family)")
check(close(_ngc1792_g_master_uqff(), 1.053e-2, 0.01),             "Regression: NGC 1792 composer = 1.053e-2 (v*B=10 family)")
check(close(_hudf_g_master_uqff(), 1.053e-3, 0.01),                "Regression: HUDF composer = 1.053e-3 (v*B=1 family)")
check(close(_ngc1275_g_master_uqff(), 3.160e-5, 0.01),             "Regression: NGC 1275 composer = 3.160e-5 (v*B=0.03 family)")
check(close(_horsehead_g_master_uqff(), 1.097e-3, 0.01),           "Regression: Horsehead composer = 1.097e-3 (v*B=1 family)")
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),            "Regression: Antennae composer = 1.053e-1 (v*B=100 family)")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),       "Regression: Bubble Nebula composer = 1.884e-3 (v*B=1.789)")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),             "Regression: NGC 3603 composer = 1.053e-3 (v*B=1 family)")

# --- 15. v*B=1 family 6th composer arithmetic verification ---
# All v*B=1 composers should produce IDENTICAL Lorentz term (1.053e-3) since
# primitive doesn't depend on system, only on B*v product, q, m, and macro_scale.
lorentz_hudf_check = _lorentz_acceleration_uqff(B_T=1.0e-6, v_ms=1.0e6,
                                                  q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                                  rho_UA_val=None, rho_SCm_val=None,
                                                  macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, lorentz_hudf_check, 1e-9),                    f"M16-Eagle Lorentz == HUDF Lorentz (both v*B=1 family) -- {lorentz:.4e} = {lorentz_hudf_check:.4e}")

print("=" * 78)
print(f"M16-Eagle smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
