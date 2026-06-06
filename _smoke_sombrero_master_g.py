"""Smoke test for Sombrero Galaxy M104 (NGC 4594) spiral+SMBH+dust-lane master Universal Gravity.

Spec: 'Master Universal Gravity Equation for Sombrero Galaxy Evolution_09May2025'
Composer: _sombrero_g_master_uqff (4 ADDITIVE leaves: grav_galaxy + g_BH + a_dust + Lorentz)
1 NEW primitive (_dynamical_friction_acceleration_uqff) + multiple reuses.
NEW LORENTZ FAMILY: 2.107e-3 (v*B=2, 6th distinct family, exactly 2x of v*B=1).
FIRST composer where Lorentz is NOT dominant (a_dust 75% + g_BH 25%).
Spec example total: 5.351e-1 m/s^2 at t=10 Gyr.
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _sombrero_g_master_uqff,
    _dynamical_friction_acceleration_uqff,    # NEW primitive
    _hubble_unified,
    _lorentz_acceleration_uqff,
    M_SOMBRERO_KG,
    M_BH_SOMBRERO_KG,
    M_TOTAL_SOMBRERO_KG,
    R_SOMBRERO_M,
    R_BH_SOMBRERO_M,
    RHO_DUST_SOMBRERO,
    RHO_ISM_SOMBRERO,
    V_ORBIT_SOMBRERO_MS,
    Z_SOMBRERO,
    H0_SOMBRERO_KMSMPC,
    B_SOMBRERO_T,
    T_SOMBRERO_DEFAULT_S,
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
print("Sombrero Galaxy M104 (NGC 4594) master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(M_SOMBRERO_KG, 1.0e11 * M_SUN, 1e-12),                 "M_SOMBRERO_KG = 10^11 M_sun = 1.989e41 kg")
check(close(M_BH_SOMBRERO_KG, 1.0e9 * M_SUN, 1e-12),               "M_BH_SOMBRERO_KG = 10^9 M_sun = 1.989e39 kg")
check(close(M_TOTAL_SOMBRERO_KG, 1.01e11 * M_SUN, 1e-9),           "M_TOTAL_SOMBRERO_KG = 1.01e11 M_sun = 2.009e41 kg")
check(close(R_SOMBRERO_M, 2.36e20, 1e-12),                         "R_SOMBRERO_M = 2.36e20 m (half 50 kly)")
check(close(R_BH_SOMBRERO_M, 1.0e15, 1e-12),                       "R_BH_SOMBRERO_M = 1e15 m (~0.1 pc SMBH influence)")
check(close(RHO_DUST_SOMBRERO, 1.0e-20, 1e-12),                    "RHO_DUST_SOMBRERO = 1e-20 kg/m^3")
check(close(RHO_ISM_SOMBRERO, 1.0e-21, 1e-12),                     "RHO_ISM_SOMBRERO = 1e-21 kg/m^3 (denominator)")
check(close(V_ORBIT_SOMBRERO_MS, 2.0e5, 1e-12),                    "V_ORBIT_SOMBRERO_MS = 2e5 m/s (orbital)")
check(close(Z_SOMBRERO, 0.0063, 1e-12),                            "Z_SOMBRERO = 0.0063 (28 Mly distance)")
check(close(H0_SOMBRERO_KMSMPC, 70.0, 1e-12),                      "H0_SOMBRERO_KMSMPC = 70")
check(close(B_SOMBRERO_T, 1.0e-5, 1e-12),                          "B_SOMBRERO_T = 1e-5 T (galactic ISM field)")
check(close(T_SOMBRERO_DEFAULT_S, 10.0e9 * _YEAR_S_MAGNETAR, 1e-12),
                                                                    "T_SOMBRERO_DEFAULT_S = 10 Gyr (3.156e17 s)")

# --- 2. NEW PRIMITIVE: _dynamical_friction_acceleration_uqff (spec literal) ---
a_dust = _dynamical_friction_acceleration_uqff(rho_drag_kgm3=RHO_DUST_SOMBRERO,
                                                 v_orbit_ms=V_ORBIT_SOMBRERO_MS,
                                                 rho_medium_kgm3=RHO_ISM_SOMBRERO,
                                                 macro_scale=1.0e-12)
# Spec: 10^-20 * (2e5)^2 / 10^-21 * 1e-12 = 10^-20 * 4e10 / 10^-21 * 1e-12
#     = 4e-10 / 1e-21 * 1e-12 = 4e11 * 1e-12 = 4e-1
expected_a_dust = (1.0e-20 * (2.0e5 ** 2) / 1.0e-21) * 1.0e-12
check(close(a_dust, expected_a_dust, 1e-9),                        f"a_dust = {a_dust:.4e} = (rho*v^2/rho_med)*1e-12 exact")
check(close(a_dust, 4.0e-1, 1e-9),                                 f"a_dust = {a_dust:.4e} = 0.4 m/s^2 (spec literal)")
# Boundary: rho_medium=0 -> 0
check(_dynamical_friction_acceleration_uqff(1e-20, 2e5, 0.0) == 0.0, "Dynamical friction: rho_medium=0 -> 0")
# Boundary: rho_drag=0 -> 0
check(_dynamical_friction_acceleration_uqff(0.0, 2e5, 1e-21) == 0.0, "Dynamical friction: rho_drag=0 -> 0")
# Boundary: v=0 -> 0
check(_dynamical_friction_acceleration_uqff(1e-20, 0.0, 1e-21) == 0.0, "Dynamical friction: v=0 -> 0")
# Boundary: macro_scale=0 -> 0
check(_dynamical_friction_acceleration_uqff(1e-20, 2e5, 1e-21, macro_scale=0.0) == 0.0,
                                                                    "Dynamical friction: macro_scale=0 -> 0")
# Scaling: 2x v -> 4x a (v^2 dependence)
a_2v = _dynamical_friction_acceleration_uqff(1e-20, 4e5, 1e-21)
check(close(a_2v / a_dust, 4.0, 1e-9),                             "Dynamical friction: 2x v -> 4x a (v^2 scaling)")
# Scaling: 2x rho_drag -> 2x a
a_2rho = _dynamical_friction_acceleration_uqff(2e-20, 2e5, 1e-21)
check(close(a_2rho / a_dust, 2.0, 1e-9),                           "Dynamical friction: 2x rho_drag -> 2x a (linear)")
# Scaling: 2x rho_medium -> 0.5x a
a_2med = _dynamical_friction_acceleration_uqff(1e-20, 2e5, 2e-21)
check(close(a_2med / a_dust, 0.5, 1e-9),                           "Dynamical friction: 2x rho_medium -> 0.5x a (inverse)")

# --- 3. SMBH separate gravitational leaf (inline -- SAME shape as Ug1, no new primitive) ---
g_BH = G_NEWTON * M_BH_SOMBRERO_KG / (R_BH_SOMBRERO_M ** 2)
expected_g_BH = 6.6743e-11 * 1.989e39 / (1.0e15 ** 2)
check(close(g_BH, expected_g_BH, 1e-6),                            f"g_BH = G*M_BH/r_BH^2 = {g_BH:.4e}")
check(close(g_BH, 1.327e-1, 0.005),                                f"g_BH ~ 1.327e-1 m/s^2 (spec, 0.5% tol)")

# --- 4. Hubble at z=0.0063 ---
H_z_kmsMpc = _hubble_unified(T_SOMBRERO_DEFAULT_S, Z_SOMBRERO, H0_SOMBRERO_KMSMPC)
expected_H = 70.0 * math.sqrt(0.3 * (1.0063 ** 3) + 0.7)
check(close(H_z_kmsMpc, expected_H, 1e-9),                         f"H(z=0.0063) = {H_z_kmsMpc:.4f} km/s/Mpc = 70*sqrt(0.3*1.0190+0.7)")
check(close(H_z_kmsMpc, 70.196, 0.002),                            f"H(z=0.0063) ~ 70.196 km/s/Mpc (spec)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
H_t = H_si * T_SOMBRERO_DEFAULT_S
check(close(H_t, 0.7178, 0.01),                                    f"H(z)*t(10 Gyr) = {H_t:.4f} ~ 0.7178")
expansion_factor = 1.0 + H_t
check(close(expansion_factor, 1.7178, 0.005),                      f"(1+H t) = {expansion_factor:.4f} ~ 1.7178")

# --- 5. Galaxy gravitational leaf (with H expansion and TRZ) ---
Ug1_galaxy = G_NEWTON * M_TOTAL_SOMBRERO_KG / (R_SOMBRERO_M ** 2)
check(close(Ug1_galaxy, 2.408e-10, 0.01),                          f"Ug1_galaxy = G*M_total/r^2 = {Ug1_galaxy:.4e} ~ 2.408e-10 m/s^2")
grav_galaxy = Ug1_galaxy * expansion_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_galaxy, 4.552e-10, 0.01),                         f"grav_galaxy*(1+H t)*(1+f_TRZ) = {grav_galaxy:.4e} ~ 4.552e-10 (CARRIED FULL per directive)")

# --- 6. Lorentz term (NEW 6TH FAMILY 2.107e-3, v*B=2) ---
lorentz = _lorentz_acceleration_uqff(B_T=B_SOMBRERO_T, v_ms=V_ORBIT_SOMBRERO_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 2.107e-3, 0.005),                             f"Lorentz q*v*B/m_p*11*1e-12 = {lorentz:.4e} ~ 2.107e-3 (NEW 6TH FAMILY, 0.5% tol)")
# Verify v*B=2 (distinct family)
check(close(B_SOMBRERO_T * V_ORBIT_SOMBRERO_MS, 2.0, 1e-12),       "B*v = 2.0 (NEW family -- distinct from 1, 100, 10, 0.03, 1.789)")

# --- 7. Composer total = 5.351e-1 m/s^2 at t=10 Gyr ---
g_total = _sombrero_g_master_uqff()
check(close(g_total, 5.351e-1, 0.005),                             f"_sombrero_g_master_uqff() = {g_total:.4e} m/s^2 ~ 5.351e-1 (spec, 0.5% tol)")

# --- 8. Fidelity check: composer == grav_galaxy + g_BH + a_dust + Lorentz exactly ---
manual = grav_galaxy + g_BH + a_dust + lorentz
check(close(g_total, manual, 1e-9),                                "Composer == grav_galaxy + g_BH + a_dust + Lorentz exactly (4 additive leaves)")

# --- 9. Contribution fractions per fidelity directive (FIRST composer where Lorentz NOT dominant) ---
dust_frac = a_dust / g_total
check(0.7 < dust_frac < 0.78,                                      f"a_dust contributes {dust_frac*100:.2f}% (MOST DOMINANT ~75%)")
bh_frac = g_BH / g_total
check(0.22 < bh_frac < 0.28,                                       f"g_BH contributes {bh_frac*100:.2f}% (~25%)")
lorentz_frac = lorentz / g_total
check(0.003 < lorentz_frac < 0.005,                                f"Lorentz contributes {lorentz_frac*100:.4f}% (only ~0.4% -- BREAKS Lorentz dominance pattern!)")
grav_galaxy_frac = grav_galaxy / g_total
check(grav_galaxy_frac < 1e-8,                                     f"grav_galaxy contributes {grav_galaxy_frac:.2e} (carried per directive)")

# --- 10. Boundary cases ---
g_no_macro = _sombrero_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_galaxy + g_BH + a_dust, 1e-9),        "macro_scale=0 -> Lorentz=0; grav + g_BH + a_dust remain")

g_no_macro_dust = _sombrero_g_master_uqff(macro_scale_dust=0.0)
check(close(g_no_macro_dust, grav_galaxy + g_BH + lorentz, 1e-9),  "macro_scale_dust=0 -> a_dust=0; only grav + g_BH + Lorentz")

g_no_BH = _sombrero_g_master_uqff(M_BH_kg=0.0)
check(close(g_no_BH, grav_galaxy + a_dust + lorentz, 1e-9),        "M_BH=0 -> g_BH=0; only grav + a_dust + Lorentz")

g_no_dust = _sombrero_g_master_uqff(rho_dust=0.0)
check(close(g_no_dust, grav_galaxy + g_BH + lorentz, 1e-9),        "rho_dust=0 -> a_dust=0; only grav + g_BH + Lorentz")

g_no_B = _sombrero_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_galaxy + g_BH + a_dust, 1e-9),            "B=0 -> Lorentz=0; grav + g_BH + a_dust remain")

g_no_trz = _sombrero_g_master_uqff(f_TRZ=0.0)
grav_galaxy_no_trz = Ug1_galaxy * expansion_factor * 1.0
check(close(g_no_trz, grav_galaxy_no_trz + g_BH + a_dust + lorentz, 1e-9),
                                                                    "f_TRZ=0 -> grav_galaxy drops (1+f_TRZ) factor")

# --- 11. Time evolution: SMBH and dust independent of t ---
g_t0 = _sombrero_g_master_uqff(t_s=0.0)
grav_galaxy_t0 = Ug1_galaxy * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_t0, grav_galaxy_t0 + g_BH + a_dust + lorentz, 1e-9), "t=0 -> H*t=0; g_BH/a_dust/Lorentz unchanged (no t-dependence)")

# --- 12. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _ngc1792_g_master_uqff,
    _hudf_g_master_uqff,
    _ngc1275_g_master_uqff,
    _horsehead_g_master_uqff,
    _antennae_g_master_uqff,
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
)
check(close(_ngc1792_g_master_uqff(), 1.053e-2, 0.01),             "Regression: NGC 1792 composer = 1.053e-2 (v*B=10 family)")
check(close(_hudf_g_master_uqff(), 1.053e-3, 0.01),                "Regression: HUDF composer = 1.053e-3 (v*B=1 family)")
check(close(_ngc1275_g_master_uqff(), 3.160e-5, 0.01),             "Regression: NGC 1275 composer = 3.160e-5 (v*B=0.03 family)")
check(close(_horsehead_g_master_uqff(), 1.097e-3, 0.01),           "Regression: Horsehead composer = 1.097e-3 (v*B=1 family)")
check(close(_antennae_g_master_uqff(), 1.053e-1, 0.01),            "Regression: Antennae composer = 1.053e-1 (v*B=100 family)")
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),       "Regression: Bubble Nebula composer = 1.884e-3 (v*B=1.789)")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),             "Regression: NGC 3603 composer = 1.053e-3 (v*B=1 family)")

# --- 13. Lorentz family separation: Sombrero 2.107e-3 is exactly 2x HUDF 1.053e-3 ---
ratio_to_hudf = lorentz / 1.053e-3
check(close(ratio_to_hudf, 2.0, 0.005),                            f"Sombrero Lorentz / HUDF Lorentz = {ratio_to_hudf:.4f} = 2 (v*B ratio 2x)")
# And exactly 5x smaller than NGC 1792 1.053e-2
ratio_to_ngc1792 = 1.053e-2 / lorentz
check(close(ratio_to_ngc1792, 5.0, 0.005),                         f"NGC 1792 Lorentz / Sombrero Lorentz = {ratio_to_ngc1792:.4f} = 5 (v*B ratio 5x)")

print("=" * 78)
print(f"Sombrero smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
