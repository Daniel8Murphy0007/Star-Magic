"""Smoke test for Saturn gas-giant planetary master Universal Gravity.

Spec: 'Master Universal Gravity Equation for Saturn Evolution_09May2025'
Composer: _saturn_g_master_uqff (5 ADDITIVE leaves: grav_sun_term + g_planet + T_ring + a_wind + Lorentz)
ZERO NEW primitives -- 3RD CONSECUTIVE PURE REUSE WIN.
NEW LORENTZ FAMILY: 5.268e-7 (v*B=5e-5, 7th distinct family, SMALLEST in lineage).
FIRST PLANETARY composer (all 6 prior were galactic/nebular).
Saturn surface gravity (10.44 m/s^2) dominates 99.9994pct.
Spec example total: 10.44 m/s^2 at t=4.5 Gyr (Solar System age).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _saturn_g_master_uqff,
    _dynamical_friction_acceleration_uqff,     # REUSE Sombrero primitive (2nd use!)
    _hubble_unified,
    _lorentz_acceleration_uqff,
    M_SUN_SATURN_KG,
    R_ORBIT_SATURN_M,
    M_SATURN_KG,
    R_SATURN_M,
    M_RING_SATURN_KG,
    R_RING_SATURN_M,
    RHO_ATM_SATURN,
    V_WIND_SATURN_MS,
    B_SATURN_T,
    Z_SATURN,
    H0_SATURN_KMSMPC,
    T_SATURN_DEFAULT_S,
    G_NEWTON,
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
print("Saturn gas-giant planetary master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(M_SUN_SATURN_KG, 1.989e30, 1e-12),                     "M_SUN_SATURN_KG = 1.989e30 kg (spec)")
check(close(R_ORBIT_SATURN_M, 1.43e12, 1e-12),                     "R_ORBIT_SATURN_M = 1.43e12 m (9.58 AU)")
check(close(M_SATURN_KG, 5.683e26, 1e-12),                         "M_SATURN_KG = 5.683e26 kg")
check(close(R_SATURN_M, 6.0268e7, 1e-12),                          "R_SATURN_M = 6.0268e7 m (60268 km equatorial)")
check(close(M_RING_SATURN_KG, 1.5e19, 1e-12),                      "M_RING_SATURN_KG = 1.5e19 kg (water ice)")
check(close(R_RING_SATURN_M, 7.0e7, 1e-12),                        "R_RING_SATURN_M = 7e7 m (~70000 km)")
check(close(RHO_ATM_SATURN, 2.0e-4, 1e-12),                        "RHO_ATM_SATURN = 2e-4 kg/m^3")
check(close(V_WIND_SATURN_MS, 500.0, 1e-12),                       "V_WIND_SATURN_MS = 500 m/s (~1800 km/h)")
check(close(B_SATURN_T, 1.0e-7, 1e-12),                            "B_SATURN_T = 1e-7 T (cloud-top field)")
check(close(Z_SATURN, 0.0, 1e-12),                                 "Z_SATURN = 0 (Solar System scale)")
check(close(H0_SATURN_KMSMPC, 70.0, 1e-12),                        "H0_SATURN_KMSMPC = 70")
check(close(T_SATURN_DEFAULT_S, 4.5e9 * _YEAR_S_MAGNETAR, 1e-12),  "T_SATURN_DEFAULT_S = 4.5 Gyr (1.420e17 s)")

# --- 2. Sun's gravitational pull at Saturn's orbit ---
g_sun_orbit = G_NEWTON * M_SUN_SATURN_KG / (R_ORBIT_SATURN_M ** 2)
check(close(g_sun_orbit, 6.494e-5, 0.005),                         f"g_sun_orbit = G*M_Sun/r_orbit^2 = {g_sun_orbit:.4e} ~ 6.494e-5 (spec)")

# --- 3. Hubble at z=0 ---
H_z_kmsMpc = _hubble_unified(T_SATURN_DEFAULT_S, Z_SATURN, H0_SATURN_KMSMPC)
check(close(H_z_kmsMpc, 70.0, 1e-9),                               f"H(z=0) = {H_z_kmsMpc:.4f} = H_0 = 70 (sqrt(0.3+0.7)=1)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_SATURN_DEFAULT_S
check(close(H_t, 0.3221, 0.01),                                    f"H(z=0)*t(4.5 Gyr) = {H_t:.4f} ~ 0.3221 (spec)")
expansion_factor = 1.0 + H_t
check(close(expansion_factor, 1.3221, 0.01),                       f"(1+H t) = {expansion_factor:.4f} ~ 1.3221")

grav_sun_term = g_sun_orbit * expansion_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_sun_term, 9.443e-5, 0.005),                       f"g_sun_orbit*(1+H t)*(1+f_TRZ) = {grav_sun_term:.4e} ~ 9.443e-5 (CARRIED FULL per directive)")

# --- 4. Saturn surface gravity (DOMINANT leaf) ---
g_planet = G_NEWTON * M_SATURN_KG / (R_SATURN_M ** 2)
check(close(g_planet, 10.44, 0.005),                               f"g_planet = G*M/r^2 = {g_planet:.4f} ~ 10.44 m/s^2 (DOMINANT, spec)")

# --- 5. Ring tidal gravitational leaf ---
T_ring = G_NEWTON * M_RING_SATURN_KG / (R_RING_SATURN_M ** 2)
check(close(T_ring, 2.043e-7, 0.005),                              f"T_ring = G*M_ring/r_ring^2 = {T_ring:.4e} ~ 2.043e-7 (spec, CARRIED per directive)")

# --- 6. Wind dynamical friction (REUSE Sombrero primitive -- DEGENERATE case rho_drag=rho_medium) ---
a_wind = _dynamical_friction_acceleration_uqff(rho_drag_kgm3=RHO_ATM_SATURN,
                                                 v_orbit_ms=V_WIND_SATURN_MS,
                                                 rho_medium_kgm3=RHO_ATM_SATURN,
                                                 macro_scale=1.0e-12)
# DEGENERATE: rho_atm*v^2/rho_atm = v^2; then *1e-12
expected_a_wind = (V_WIND_SATURN_MS ** 2) * 1.0e-12
check(close(a_wind, expected_a_wind, 1e-12),                       f"a_wind = v_wind^2*1e-12 = {a_wind:.4e} (rho_drag=rho_medium degenerate cancels)")
check(close(a_wind, 2.5e-7, 1e-9),                                 f"a_wind = {a_wind:.4e} ~ 2.5e-7 m/s^2 (spec literal)")

# --- 7. Lorentz term (NEW 7TH FAMILY 5.27e-8 -- SMALLEST in lineage) ---
# FIDELITY NOTE: spec declares Lorentz=5.268e-7 from chain '1.602e-19*500*1e-7=8.01e-23'
# but 1.602e-19*500*1e-7 = 8.01e-24 (spec arithmetic error 10x in intermediate q*v*B).
# Our primitive evaluates spec's formula q(v x B)*(1+rho)*1e-12 CORRECTLY to 5.27e-8.
# Spec's FINAL composer total 10.44 is preserved either way since Lorentz<<10.44.
lorentz = _lorentz_acceleration_uqff(B_T=B_SATURN_T, v_ms=V_WIND_SATURN_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 5.27e-8, 0.005),                              f"Lorentz q*v*B/m_p*11*1e-12 = {lorentz:.4e} ~ 5.27e-8 (NEW 7TH FAMILY SMALLEST -- spec arith error declared 5.268e-7, correct=5.27e-8)")
# Linear scaling check: Saturn / HUDF should be exactly v*B ratio (5e-5/1 = 5e-5)
# HUDF Lorentz = 1.0537e-3, so Saturn = 1.0537e-3 * 5e-5 = 5.27e-8 (verified)
from uqff_pure_calculator import _hudf_g_master_uqff as _hudf_check
_ = _hudf_check
hudf_lorentz_check = _lorentz_acceleration_uqff(B_T=1.0e-6, v_ms=1.0e6,
                                                  q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                                  rho_UA_val=None, rho_SCm_val=None,
                                                  macro_scale=MACROSCOPIC_SCALE_LORENTZ)
scaling_ratio = lorentz / hudf_lorentz_check
expected_scaling = (B_SATURN_T * V_WIND_SATURN_MS) / (1.0e-6 * 1.0e6)
check(close(scaling_ratio, expected_scaling, 1e-9),               f"Saturn/HUDF Lorentz ratio = {scaling_ratio:.4e} = v*B ratio {expected_scaling:.4e} (linear scaling verified)")
check(close(B_SATURN_T * V_WIND_SATURN_MS, 5.0e-5, 1e-12),         "B*v = 5e-5 (NEW SMALLEST family -- distinct from 1, 100, 10, 0.03, 1.789, 2)")

# --- 8. Composer total ~ 10.44 m/s^2 at t=4.5 Gyr ---
g_total = _saturn_g_master_uqff()
check(close(g_total, 10.44, 0.005),                                f"_saturn_g_master_uqff() = {g_total:.4f} m/s^2 ~ 10.44 (spec, 0.5% tol)")

# --- 9. Fidelity check: composer == 5 leaves summed exactly ---
manual = grav_sun_term + g_planet + T_ring + a_wind + lorentz
check(close(g_total, manual, 1e-9),                                "Composer == grav_sun_term + g_planet + T_ring + a_wind + Lorentz exactly (5 additive leaves)")

# --- 10. Contribution fractions per fidelity directive (Saturn dwarfs everything) ---
planet_frac = g_planet / g_total
check(planet_frac > 0.9999,                                        f"g_planet contributes {planet_frac*100:.5f}% (DOMINANT >99.9994%)")
sun_frac = grav_sun_term / g_total
check(sun_frac < 1e-4,                                             f"grav_sun_term contributes {sun_frac:.3e} (carried per directive)")
ring_frac = T_ring / g_total
check(ring_frac < 1e-7,                                            f"T_ring contributes {ring_frac:.3e} (carried per directive)")
wind_frac = a_wind / g_total
check(wind_frac < 1e-7,                                            f"a_wind contributes {wind_frac:.3e} (carried per directive)")
lorentz_frac = lorentz / g_total
check(lorentz_frac < 1e-7,                                         f"Lorentz contributes {lorentz_frac:.3e} (SMALLEST in lineage, carried per directive)")

# --- 11. Boundary cases ---
g_no_macro = _saturn_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_sun_term + g_planet + T_ring + a_wind, 1e-9),
                                                                    "macro_scale=0 -> Lorentz=0; 4 other leaves remain")

g_no_macro_wind = _saturn_g_master_uqff(macro_scale_wind=0.0)
check(close(g_no_macro_wind, grav_sun_term + g_planet + T_ring + lorentz, 1e-9),
                                                                    "macro_scale_wind=0 -> a_wind=0; 4 other leaves remain")

g_no_ring = _saturn_g_master_uqff(M_ring_kg=0.0)
check(close(g_no_ring, grav_sun_term + g_planet + a_wind + lorentz, 1e-9),
                                                                    "M_ring=0 -> T_ring=0; 4 other leaves remain")

g_no_atm = _saturn_g_master_uqff(rho_atm=0.0)
# When rho_atm=0, both rho_drag and rho_medium become 0; primitive returns 0 (rho_medium<=0 branch)
check(close(g_no_atm, grav_sun_term + g_planet + T_ring + lorentz, 1e-9),
                                                                    "rho_atm=0 -> a_wind=0 (primitive handles rho_medium=0 safely)")

g_no_B = _saturn_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_sun_term + g_planet + T_ring + a_wind, 1e-9),
                                                                    "B=0 -> Lorentz=0; 4 other leaves remain")

g_no_trz = _saturn_g_master_uqff(f_TRZ=0.0)
grav_sun_no_trz = g_sun_orbit * expansion_factor * 1.0
check(close(g_no_trz, grav_sun_no_trz + g_planet + T_ring + a_wind + lorentz, 1e-9),
                                                                    "f_TRZ=0 -> grav_sun_term drops (1+f_TRZ) factor")

# --- 12. Time evolution: only grav_sun_term has H*t dependence; planet/ring/wind/Lorentz constant ---
g_t0 = _saturn_g_master_uqff(t_s=0.0)
grav_sun_t0 = g_sun_orbit * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_t0, grav_sun_t0 + g_planet + T_ring + a_wind + lorentz, 1e-9),
                                                                    "t=0 -> H*t=0; only grav_sun_term changes")
# At present (t=0 from now), grav_sun_term should still be tiny
check(close(grav_sun_t0, 7.143e-5, 0.01),                          f"grav_sun_term(t=0) = {grav_sun_t0:.4e} = g_sun*(1+f_TRZ) only")

# --- 13. Wind dynamical friction degeneracy: 2x rho_atm should still give SAME a_wind ---
a_wind_2rho = _dynamical_friction_acceleration_uqff(rho_drag_kgm3=2*RHO_ATM_SATURN,
                                                      v_orbit_ms=V_WIND_SATURN_MS,
                                                      rho_medium_kgm3=2*RHO_ATM_SATURN)
check(close(a_wind_2rho, a_wind, 1e-9),                            "Wind degenerate: 2x rho_atm (both drag and medium) -> SAME a_wind (cancels)")

# --- 14. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _sombrero_g_master_uqff,
    _ngc1792_g_master_uqff,
    _hudf_g_master_uqff,
    _ngc1275_g_master_uqff,
    _horsehead_g_master_uqff,
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
)
check(close(_sombrero_g_master_uqff(), 5.351e-1, 0.01),            "Regression: Sombrero composer = 5.351e-1 (v*B=2 family)")
check(close(_ngc1792_g_master_uqff(), 1.053e-2, 0.01),             "Regression: NGC 1792 composer = 1.053e-2 (v*B=10 family)")
check(close(_hudf_g_master_uqff(), 1.053e-3, 0.01),                "Regression: HUDF composer = 1.053e-3 (v*B=1 family)")
check(close(_ngc1275_g_master_uqff(), 3.160e-5, 0.01),             "Regression: NGC 1275 composer = 3.160e-5 (v*B=0.03 family)")
check(close(_horsehead_g_master_uqff(), 1.097e-3, 0.01),           "Regression: Horsehead composer = 1.097e-3 (v*B=1 family)")
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),            "Regression: Antennae composer = 1.053e-1 (v*B=100 family)")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),       "Regression: Bubble Nebula composer = 1.884e-3 (v*B=1.789)")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),             "Regression: NGC 3603 composer = 1.053e-3 (v*B=1 family)")

# --- 15. Lorentz family separation: Saturn 5.27e-8 -- SMALLEST in lineage ---
# Saturn / HUDF = 5e-5 (v*B ratio 5e-5)
ratio_to_hudf = lorentz / 1.053e-3
check(close(ratio_to_hudf, 5.0e-5, 0.01),                          f"Saturn Lorentz / HUDF Lorentz = {ratio_to_hudf:.4e} = 5e-5 (v*B ratio 5e-5)")
# NGC 1275 / Saturn = 600 (v*B ratio 0.03/5e-5 = 600)
ratio_ngc1275 = 3.160e-5 / lorentz
check(close(ratio_ngc1275, 600.0, 0.01),                           f"NGC 1275 Lorentz / Saturn Lorentz = {ratio_ngc1275:.4f} ~ 600 (v*B ratio 0.03/5e-5 = 600x exact)")

print("=" * 78)
print(f"Saturn smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
