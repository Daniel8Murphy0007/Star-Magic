"""Smoke test for Antennae Galaxies NGC 4038/4039 merger evolution master Universal Gravity.

Spec: 'Master Universal Gravity Equation_Antennae Galaxies reloaded Evolution_09May2025'
Composer: _antennae_g_master_uqff (4-leaf clean form per fidelity directive)
Two NEW primitives validated: _sfr_linear_mass_growth_uqff + _merger_progress_saturating_uqff
Spec example total: 1.053e-1 m/s^2 at t=300 Myr (Lorentz dominates by ~10^8).
"""
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uqff_pure_calculator import (
    _antennae_g_master_uqff,
    _sfr_linear_mass_growth_uqff,
    _merger_progress_saturating_uqff,
    _hubble_unified,
    _lorentz_acceleration_uqff,
    _antennae_g_primitive_sat,         # existing BSFG triadic -- regression check
    M_ANTENNAE_INIT_KG,
    R_ANTENNAE_M,
    B_ANTENNAE_T,
    SFR_ANTENNAE_KGS,
    TAU_MERGE_ANTENNAE_S,
    M_COLL_0_ANTENNAE,
    Z_ANTENNAE,
    H0_ANTENNAE_KMSMPC,
    V_GAS_ANTENNAE_MS,
    T_ANTENNAE_DEFAULT_S,
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
print("Antennae Galaxies NGC 4038/4039 master g -- smoke")
print("=" * 78)

# --- 1. Constants exact ---
check(close(M_ANTENNAE_INIT_KG, 2.0e11 * M_SUN, 1e-12),       "M_ANTENNAE_INIT_KG = 2e11 M_sun = 3.978e41 kg")
check(close(R_ANTENNAE_M, 2.838e20, 1e-12),                    "R_ANTENNAE_M = 2.838e20 m (30 kly core sep)")
check(close(B_ANTENNAE_T, 1.0e-4, 1e-12),                      "B_ANTENNAE_T = 1e-4 T (enhanced starburst)")
check(close(SFR_ANTENNAE_KGS, 20.0 * M_SUN / _YEAR_S_MAGNETAR, 1e-12),
      "SFR_ANTENNAE_KGS = 20 M_sun/yr in kg/s")
check(close(TAU_MERGE_ANTENNAE_S, 400.0e6 * _YEAR_S_MAGNETAR, 1e-12),
      "TAU_MERGE_ANTENNAE_S = 400 Myr in seconds")
check(close(M_COLL_0_ANTENNAE, 0.5, 1e-12),                    "M_COLL_0_ANTENNAE = 0.5")
check(close(Z_ANTENNAE, 0.0105, 1e-12),                        "Z_ANTENNAE = 0.0105 (45 Mly)")
check(close(H0_ANTENNAE_KMSMPC, 70.0, 1e-12),                  "H0_ANTENNAE_KMSMPC = 70 km/s/Mpc")
check(close(V_GAS_ANTENNAE_MS, 1.0e6, 1e-12),                  "V_GAS_ANTENNAE_MS = 1e6 m/s (starburst outflow)")
check(close(T_ANTENNAE_DEFAULT_S, 300.0e6 * _YEAR_S_MAGNETAR, 1e-12),
      "T_ANTENNAE_DEFAULT_S = 300 Myr in seconds")

# --- 2. Linear SFR mass growth primitive ---
M_t_300Myr = _sfr_linear_mass_growth_uqff()
check(close(M_t_300Myr, 4.097e41, 0.005),                      f"M(300 Myr) = {M_t_300Myr:.4e} kg ~ 4.097e41 kg (0.5% tol)")
M_merge_frac = (M_t_300Myr - M_ANTENNAE_INIT_KG) / M_ANTENNAE_INIT_KG
check(close(M_merge_frac, 0.0300, 0.01),                       f"SFR*t/M_init = {M_merge_frac:.4f} ~ 0.0300 (1% tol)")
check(close(_sfr_linear_mass_growth_uqff(t_s=0.0), M_ANTENNAE_INIT_KG, 1e-12),
      "M(t=0) = M_init exactly")
check(_sfr_linear_mass_growth_uqff(t_s=1e17) > M_ANTENNAE_INIT_KG,
      "M(t>>0) > M_init (linear growth)")
# Linearity check
M_at_t = _sfr_linear_mass_growth_uqff(t_s=1e15)
M_at_2t = _sfr_linear_mass_growth_uqff(t_s=2e15)
delta1 = M_at_t - M_ANTENNAE_INIT_KG
delta2 = M_at_2t - M_ANTENNAE_INIT_KG
check(close(delta2, 2.0 * delta1, 1e-9),                       "Linear: delta(2t) = 2*delta(t) exactly")

# --- 3. Saturating merger progress primitive ---
M_coll_300Myr = _merger_progress_saturating_uqff()
expected_M_coll = 0.5 * (1.0 - math.exp(-0.75))
check(close(M_coll_300Myr, expected_M_coll, 1e-6),             f"M_coll(300 Myr) = {M_coll_300Myr:.4f} = 0.5*(1-e^-0.75) = {expected_M_coll:.4f}")
check(close(M_coll_300Myr, 0.2638, 0.005),                     "M_coll(300 Myr) ~ 0.2638 (0.5% tol)")
check(close(_merger_progress_saturating_uqff(t_s=0.0), 0.0, 1e-12),
      "M_coll(t=0) = 0 exactly")
check(_merger_progress_saturating_uqff(t_s=1e20) > 0.4999,
      "M_coll(t>>tau) -> M_0 = 0.5 (saturates)")
M_coll_inf = _merger_progress_saturating_uqff(t_s=1e20)
check(close(M_coll_inf, 0.5, 1e-6),                            f"M_coll(t->inf) = {M_coll_inf:.6f} ~ 0.5")

# --- 4. Hubble at z=0.0105 ---
H_z_kmsMpc = _hubble_unified(T_ANTENNAE_DEFAULT_S, Z_ANTENNAE, H0_ANTENNAE_KMSMPC)
check(close(H_z_kmsMpc, 70.34, 0.005),                         f"H(z=0.0105) = {H_z_kmsMpc:.4f} km/s/Mpc ~ 70.34 (0.5% tol)")
H_si = H_z_kmsMpc * 1.0e3 / _MPC_M
check(close(H_si, 2.279e-18, 0.005),                           f"H(z) SI = {H_si:.4e} s^-1 ~ 2.279e-18")
H_t = H_si * T_ANTENNAE_DEFAULT_S
check(close(H_t, 2.158e-2, 0.005),                             f"H(z)*t = {H_t:.4e} ~ 2.158e-2 (1+H*t = 1.0216)")

# --- 5. Gravitational term ---
Ug1 = G_NEWTON * M_t_300Myr / (R_ANTENNAE_M * R_ANTENNAE_M)
check(close(Ug1, 3.395e-10, 0.005),                            f"Ug1 = G*M(t)/r^2 = {Ug1:.4e} ~ 3.395e-10 m/s^2 (0.5% tol)")
merger_factor = 1.0 - M_coll_300Myr
check(close(merger_factor, 0.7362, 0.005),                     f"(1-M_coll) = {merger_factor:.4f} ~ 0.7362")
grav_full = Ug1 * (1.0 + H_t) * merger_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(grav_full, 2.808e-10, 0.01),                       f"grav*(1+H t)*(1-M_coll)*(1+f_TRZ) = {grav_full:.4e} ~ 2.81e-10 (1% tol)")

# --- 6. Lorentz term (NEW 1.053e-1 family) ---
lorentz = _lorentz_acceleration_uqff(B_T=B_ANTENNAE_T, v_ms=V_GAS_ANTENNAE_MS,
                                       q_C=EV_J, m_kg=_M_PROTON_KG_MAGNETAR,
                                       rho_UA_val=None, rho_SCm_val=None,
                                       macro_scale=MACROSCOPIC_SCALE_LORENTZ)
check(close(lorentz, 1.053e-1, 0.005),                         f"Lorentz q*v*B/m_p*(1+rho_UA/rho_SCm)*1e-12 = {lorentz:.4e} ~ 1.053e-1 (NEW family, 0.5% tol)")
# Verify it's 100x stronger than NGC 3603/Westerlund family (1.053e-3)
check(lorentz > 0.01,                                          "Lorentz > 1e-2 (100x stronger than 1.053e-3 family)")

# --- 7. Composer total = 1.053e-1 m/s^2 at t=300 Myr ---
g_total = _antennae_g_master_uqff()
check(close(g_total, 1.053e-1, 0.005),                         f"_antennae_g_master_uqff() = {g_total:.4e} m/s^2 ~ 1.053e-1 (spec example, 0.5% tol)")

# --- 8. Fidelity check: composer == grav + Lorentz exactly (NO extra leaves) ---
manual = grav_full + lorentz
check(close(g_total, manual, 1e-12),                           "Composer == grav*(1+H t)*(1-M_coll)*(1+f_TRZ) + Lorentz exactly (no extra leaves)")

# --- 9. Lorentz dominates by ~10^8 ---
ratio = grav_full / lorentz
check(ratio < 1e-7,                                            f"Grav/Lorentz = {ratio:.3e} < 1e-7 (Lorentz dominates by ~10^8)")

# --- 10. Boundary cases ---
g_no_macro = _antennae_g_master_uqff(macro_scale=0.0)
check(close(g_no_macro, grav_full, 1e-9),                      "macro_scale=0 -> Lorentz=0, only grav remains")

g_no_trz = _antennae_g_master_uqff(f_TRZ=0.0)
grav_no_trz = Ug1 * (1.0 + H_t) * merger_factor * 1.0
check(close(g_no_trz, grav_no_trz + lorentz, 1e-6),            "f_TRZ=0 -> grav drops factor (1+f_TRZ)")

g_no_merge = _antennae_g_master_uqff(M_coll_0=0.0)
grav_no_merge = Ug1 * (1.0 + H_t) * 1.0 * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_merge, grav_no_merge + lorentz, 1e-6),        "M_coll_0=0 -> merger_factor=1, no merger suppression")

g_no_sfr = _antennae_g_master_uqff(SFR_kgs=0.0)
Ug1_no_sfr = G_NEWTON * M_ANTENNAE_INIT_KG / (R_ANTENNAE_M * R_ANTENNAE_M)
grav_no_sfr = Ug1_no_sfr * (1.0 + H_t) * merger_factor * (1.0 + _F_TRZ_DEFAULT_MAGNETAR)
check(close(g_no_sfr, grav_no_sfr + lorentz, 1e-6),            "SFR=0 -> M(t)=M_init (static mass)")

g_no_B = _antennae_g_master_uqff(B_T=0.0)
check(close(g_no_B, grav_full, 1e-9),                          "B=0 -> Lorentz=0, only grav remains")

g_no_v = _antennae_g_master_uqff(v_gas_ms=0.0)
check(close(g_no_v, grav_full, 1e-9),                          "v=0 -> Lorentz=0, only grav remains")

# --- 11. Time evolution monotonicity ---
g_1Gyr = _antennae_g_master_uqff(t_s=1.0e9 * _YEAR_S_MAGNETAR)
g_100Myr = _antennae_g_master_uqff(t_s=100.0e6 * _YEAR_S_MAGNETAR)
# Lorentz dominates so g_total is approx constant; check grav-only path
M_1G = _sfr_linear_mass_growth_uqff(t_s=1.0e9 * _YEAR_S_MAGNETAR)
M_100M = _sfr_linear_mass_growth_uqff(t_s=100.0e6 * _YEAR_S_MAGNETAR)
check(M_1G > M_100M,                                            "M(1 Gyr) > M(100 Myr) (mass grows linearly)")
M_coll_1G = _merger_progress_saturating_uqff(t_s=1.0e9 * _YEAR_S_MAGNETAR)
M_coll_100M = _merger_progress_saturating_uqff(t_s=100.0e6 * _YEAR_S_MAGNETAR)
check(M_coll_1G > M_coll_100M,                                  "M_coll(1 Gyr) > M_coll(100 Myr) (merger progresses)")
check(M_coll_1G < 0.5,                                          "M_coll(1 Gyr) < 0.5 (saturates below M_0)")

# --- 12. Existing primitives UNTOUCHED ---
# _antennae_g_primitive_sat is the BSFG saturation triadic (different observable)
sat_val = _antennae_g_primitive_sat()
check(sat_val is not None and sat_val > 0,                      f"_antennae_g_primitive_sat() = {sat_val:.4e} unchanged (BSFG triadic, different observable)")

# --- 13. Regression: prior composers still importable + sane ---
from uqff_pure_calculator import (
    _bubble_nebula_g_master_uqff,
    _ngc3603_g_master_uqff,
    _ngc2525_g_master_uqff,
    _rings_g_master_uqff,
    _pillars_g_master_uqff,
    _westerlund2_g_master_uqff,
    _tapestry_g_master_uqff,
    _sgr_a_g_master_uqff,
    _magnetar_g_master_uqff,
    _magnetar_g_master_uqff_v2,
)
check(close(_bubble_nebula_g_master_uqff(), 1.884e-3, 0.01),    "Regression: Bubble Nebula composer = 1.884e-3")
check(close(_ngc3603_g_master_uqff(), 1.053e-3, 0.01),          "Regression: NGC 3603 composer = 1.053e-3")
check(close(_westerlund2_g_master_uqff(), 1.053e-3, 0.01),      "Regression: Westerlund 2 composer = 1.053e-3 family")

print("=" * 78)
print(f"Antennae smoke: {passed} passed, {failed} failed")
print("=" * 78)
sys.exit(0 if failed == 0 else 1)
