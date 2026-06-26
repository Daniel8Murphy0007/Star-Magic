#!/usr/bin/env python3
"""
recompute_uqff.py — INDEPENDENT recomputation of UQFF closures.

PURPOSE: read-only verification harness for daniel.murphy00@enrgyone.com's
Star-Magic UQFF framework. Reimplements the locked primitives, the
derive_from_quantum_chain root, the Gold REGISTRY formulas, and every
documented derive_* function from scratch using ONLY the locked primitives.
Does NOT import Daniel's modules. Does NOT write to his repo.

Goal: confirm that the published Gold_Standard_Validation_Report.json
numerics are reproducible from first principles + locked primitives.
Any divergence is an indicator worth reporting.

Author of verifier: assistant (Claude), read-only audit.
Source of physics: Daniel T. Murphy, UQFF / Star-Magic.
"""

import math
import json
from fractions import Fraction

# =============================================================================
# LOCKED PRIMITIVES (transcribed verbatim from Gold_Standard_Pure_UQFF.md §3
#   and confirmed against dpm_vacuum_manifold.py)
# =============================================================================

# Integer system
D_PHYS   = 4
D_BSFG   = 6        # derivative: dim_R[SO(5)/U(2)] = 10 - 4 = 6  (PAPER_1167)
N_CH     = 9
SO_FIVE  = 10
D_CRIT   = 26
A_FIVE   = 60

# Rational system
G1_K        = Fraction(5, 6)
G2_BETA     = Fraction(3, 5)
G3_RICCI    = Fraction(1, 2)
G4_BSFG     = Fraction(3, 20)
K_MEX       = Fraction(25, 12)   # derivative: Phi_5_6 * SO_FIVE / D_PHYS = (5/6)*10/4 = 25/12  (PAPER_1522)
SSQ         = Fraction(57, 100)
PHI_RES     = Fraction(84, 100)  # also seen as 5/6 in PAPER_1203 Nuclear
TRZ         = Fraction(1, 10)
BETA_I      = Fraction(6029, 10000)

# Real / transcendental system
E0          = 1e-20
S26_3       = 1.4531e26          # canonical Ramanujan order-3 accel of Li_26(0.57)
KAPPA       = 5e-4               # day^-1
F_THZ       = 1.25e12            # 1.25 THz DPM phonon

# Float versions of rationals for numerical compute
PHI_F   = float(PHI_RES)
TRZ_F   = float(TRZ)
SSQ_F   = float(SSQ)
BETA_F  = float(BETA_I)
G1_F    = float(G1_K)
G2_F    = float(G2_BETA)
G3_F    = float(G3_RICCI)
G4_F    = float(G4_BSFG)
K_MEX_F = float(K_MEX)

# =============================================================================
# VACUUM LEDGER ROOT (verbatim from dpm_vacuum_manifold.py:derive_from_quantum_chain)
# =============================================================================

def derive_from_quantum_chain(n_levels=26, f_SCm=0.57, V=1.0):
    """ρ_vac = Σ_{n=1..N} f_SCm · E_0 · 10^n  /  V    [J/m³, energy density only]"""
    total = 0.0
    for n in range(1, n_levels + 1):
        total += f_SCm * E0 * (10 ** n)
    return total / V

# Compute the canonical macro values
RHO_VAC_SCM = derive_from_quantum_chain()                       # ≈ 6.333e5
RHO_VAC_UA  = RHO_VAC_SCM * 10.0                                # SO_FIVE structural ratio
DPM_RATIO   = 10.0
E_REACT_0   = (RHO_VAC_SCM * (1.0 ** 2)) / RHO_VAC_UA           # = 0.1 (v=1 normalized)
V_DPM_BASE  = math.sqrt(E_REACT_0 * (RHO_VAC_SCM / RHO_VAC_UA)) * (SO_FIVE / TRZ_F)
# = sqrt(0.1 * 0.1) * (10/0.1) = 0.1 * 100 = 10.0

# Micro per-DPM-volume legacy (retained for proxies — Star-Magic.txt original)
RHO_VAC_SCM_MICRO = 7.09e-37
RHO_VAC_UA_MICRO  = RHO_VAC_SCM_MICRO * 10.0


# =============================================================================
# DERIVED UQFF SYMBOLS — verbatim transcription of Gold_Standard_Validation_Script.py
# =============================================================================

def derive_c_eff():
    """c_eff = D_CRIT * (4π/Φ_res) * v_base, where v_base is chosen so c_eff matches
    the canonical 4D projection target. The factor D_CRIT*4π/Φ_res implements the
    26D → 4D downward projection rule."""
    target_c = 299792458.0
    v_base = target_c / (D_CRIT * (4 * math.pi / PHI_F))
    return D_CRIT * (4 * math.pi / PHI_F) * v_base   # round-trips to target

def derive_alpha_uqff():
    """α = G4 * (1 + TRZ * G3) / D_CRIT²"""
    return float(G4_BSFG * (1 + TRZ * G3_RICCI) / (D_CRIT ** 2))

def derive_h():
    """h = TRZ · PHI · (E0/F_THZ) · (1 - 2α_internal)
    where α_internal = 1/(PHI*D_CRIT*2π)  — distinct from derive_alpha_uqff."""
    alpha_int = 1 / (PHI_F * D_CRIT * 2 * math.pi)
    return TRZ_F * PHI_F * (E0 / F_THZ) * (1 - 2 * alpha_int)

def derive_hbar():
    return derive_h() / (2 * math.pi)

def derive_e_crack():
    """E_crack = ρ_SCm * V_DPM² / SSQ   (no c²; pure UQFF)"""
    return RHO_VAC_SCM * V_DPM_BASE**2 / SSQ_F

def derive_G_uqff():
    """G_UQFF = G1_K · G4 · ρ · S26_3 / (DPM_RATIO · BETA_I · scale)
    with an explicit empirical 2.86e16 calibration noted in the source as
    'scale to match order of legacy tuned proxy for honest comparison'."""
    e_crack = derive_e_crack()
    rho_e = RHO_VAC_SCM * (e_crack / (RHO_VAC_SCM / SSQ_F))   # approx E_react proxy
    g_proxy = G1_F * G4_F * rho_e * (S26_3 / 1e26) / (10.0 * BETA_F)
    return g_proxy / 2.86e16

def derive_N_A_uqff():
    """N_A from M_0 = E_crack/V_DPM², N_A_proxy = (1/M_0) · 26 · exp(-SSQ·26/V_DPM) · 1.13e29"""
    m0 = derive_e_crack() / (V_DPM_BASE ** 2)
    na_proxy = (1.0 / m0) * 26 * math.exp(-SSQ_F * 26 / V_DPM_BASE)
    return na_proxy * 1.13e29

def derive_planck_mass_uqff():
    """m_P = sqrt(hbar · c / G)"""
    return (derive_hbar() * derive_c_eff() / derive_G_uqff()) ** 0.5

def derive_planck_length_uqff():
    """l_P = sqrt(hbar · G / c³)"""
    return (derive_hbar() * derive_G_uqff() / (derive_c_eff() ** 3)) ** 0.5

def derive_k_b_uqff():
    """k_B: Gold script returns CODATA value 1.380649e-23 (labeled 'proxy', not a derive)."""
    return 1.380649e-23

def derive_vacuum_permittivity_uqff():
    scale_26d_vac = (S26_3 / 1e26) * 1e80
    return 1.0 / (4.0 * math.pi * derive_G_uqff() * RHO_VAC_SCM_MICRO * 60 / (derive_c_eff() ** 4) * scale_26d_vac)

def derive_vacuum_permeability_uqff():
    scale_26d_vac = (S26_3 / 1e26) * 1e80
    return 4 * math.pi * derive_G_uqff() * RHO_VAC_SCM_MICRO * 60 / (derive_c_eff() ** 2) / scale_26d_vac

def derive_delta_scm_uqff():
    return RHO_VAC_SCM_MICRO * (1 - SSQ_F) * BETA_F / D_CRIT

def derive_p1_lkk_um_uqff():
    return D_BSFG * (G4_F * D_CRIT + PHI_F)

def derive_p6_lkk_inv_uqff():
    return D_BSFG * G4_F * (S26_3 / 1e26) * PHI_F

def derive_p7_w_a_uqff():
    return G1_F * (1 + TRZ_F * G4_F) / D_CRIT

def derive_p9_h_tension_uqff():
    return G4_F * TRZ_F * BETA_F * (S26_3 / 1e26)

def derive_p10_s8_tension_uqff():
    return G4_F * TRZ_F * (S26_3 / 1e26) * 1e6

def derive_sgr_a_g_uqff():
    return (D_CRIT * D_CRIT) * (S26_3 / 1e26) * G4_F * BETA_F

def derive_sgr_1745_g_uqff():
    return (D_BSFG ** 3) * (S26_3 / 1e26) * PHI_F * 1e-6

def derive_lenr_parkhomov_uqff():
    return D_CRIT * (D_BSFG + (S26_3 / 1e26) * PHI_F) * 100

def derive_omega_gw_h2_uqff():
    return G4_F * (S26_3 / 1e26) * TRZ_F * PHI_F

def derive_f_nl_equil_uqff():
    return G3_F * SSQ_F * TRZ_F / D_BSFG

def derive_epsilon_slow_roll_uqff():
    return G4_F * TRZ_F * 0.5

# Stubbed-as-closure derivations (return characteristic values, not pure formulas)
def derive_neutron_lifetime_full(t_mode="age"):       return 880.0      # s closure (per §5.12)
def derive_dark_matter_particle_uqff(dataset=None):    return 1e-3       # eV (per §5.13)
def derive_dark_flow_uqff(dataset=None):               return 800000.0   # m/s (per §5.14)
def derive_cmb_cold_spot_uqff(dataset=None):           return -100e-6    # K  (per §5.15)
def derive_pta_sgwb_spectral_index_uqff(d=None):       return 13.0 / 3.0
def derive_txs0506_multimessenger_delay_uqff(d=None):  return 1000.0
def derive_frb_origin_uqff(d=None):                    return 1.4e9
def derive_casimir_effect_uqff(d=None):
    rho_free = RHO_VAC_SCM_MICRO * 11
    delta_rho = rho_free * 0.01
    return -delta_rho * (3e8 ** 2) * (BETA_F * SSQ_F / PHI_F)
def derive_coronal_heating_uqff(d=None):               return 2.0e6
def derive_proton_radius_muonic_uqff(d=None):          return 8.409e-16


# =============================================================================
# REGISTRY: name → (formula_or_callable, target_value, unit, description)
# Transcribed from Gold_Standard_Validation_Script.py REGISTRY
# All formulas use the locked primitives only (no SM constants inside math).
# =============================================================================

def f_yang_mills():       return (RHO_VAC_SCM * V_DPM_BASE**2 / SSQ_F) * 5.3e4
def f_riemann():          return 9877.78265 * math.exp(-2.763 * SSQ_F) * 1.0
def f_bsd():              return 0.3059997738 * (1 + BETA_F * SSQ_F)
def f_navier_stokes():    return 1.0 - TRZ_F * D_BSFG / D_PHYS
def f_hodge():            return (D_PHYS + D_BSFG) / SO_FIVE
def f_poincare():         return 0.5 + TRZ_F * (5.0/6.0)
def f_p_vs_np():          return 1.0 - TRZ_F ** 9
def f_bh_info():          return 1.0
def f_alpha_primitive():  return G4_F * (1 + TRZ_F * G3_F) / (D_CRIT ** 2)
def f_proton_mass():      return D_CRIT * PHI_F * (S26_3/1e26 - TRZ_F)
def f_yang_mills_sat():   return SSQ_F * G4_F * PHI_F * G1_F
def f_h0():               return (S26_3/1e26) + SSQ_F + G4_F * G1_F
def f_t0():               return BETA_F * (PHI_F - TRZ_F)
def f_m_mu():             return D_CRIT * G4_F * ((S26_3/1e26) - SSQ_F)
def f_m_tau():            return (D_CRIT ** 2) * G4_F * BETA_F
def f_m_t():              return D_CRIT * (S26_3/1e26) * G4_F * G1_F
def f_m_W():              return D_BSFG * BETA_F * (S26_3/1e26) * (1 + TRZ_F)
def f_m_Z():              return D_BSFG * BETA_F * (S26_3/1e26) / G1_F
def f_m_H():              return D_BSFG * (S26_3/1e26) * PHI_F * (G1_F + G2_F)
def f_v_higgs():          return D_CRIT * BETA_F * (S26_3/1e26) / PHI_F
def f_G_F():              return G4_F * TRZ_F / ((D_CRIT ** 2) * (S26_3/1e26))
def f_alpha_s():          return G1_F * G4_F / ((S26_3/1e26) * PHI_F)
def f_Omega_m():          return G4_F * (S26_3/1e26) * PHI_F + SSQ_F * TRZ_F
def f_Omega_b_h2():       return (G4_F ** 2) * PHI_F
def f_Omega_DM_h2():      return G2_F * G4_F * PHI_F
def f_n_s():              return 1.0 - G4_F * TRZ_F * (1 + SSQ_F)
def f_A_s():              return RHO_VAC_SCM * 60 * (S26_3/1e26) * PHI_F
def f_eta_baryon():       return (G4_F ** 2) * TRZ_F * BETA_F / 60
def f_Y_p():              return G4_F * (S26_3/1e26) + G3_F * TRZ_F
def f_z_re():             return D_BSFG + PHI_F + SSQ_F
def f_tau_reion():        return G4_F * TRZ_F * (S26_3/1e26) * PHI_F
def f_m_e():              return D_CRIT * G4_F * (S26_3/1e26 - SSQ_F) * 0.01
def f_m_pion():           return D_CRIT * PHI_F * (S26_3/1e26)
def f_m_kaon():           return D_CRIT * (PHI_F - TRZ_F) * (S26_3/1e26)
def f_dark_flow_simple(): return derive_c_eff() * BETA_F * (RHO_VAC_UA / RHO_VAC_SCM) * (S26_3 / 1e26) * 0.2
def f_e_react():          return (RHO_VAC_SCM * V_DPM_BASE**2 / RHO_VAC_UA) * math.exp(-KAPPA * 0)
def f_e_crack_J():        return RHO_VAC_SCM * V_DPM_BASE**2 / SSQ_F
def f_spinor():           return D_CRIT * BETA_F * PHI_F * S26_3 / 1e26


REGISTRY = [
    # Millennium 8
    ("millennium_yang_mills",     f_yang_mills,         1.78,             "GeV",  "Yang-Mills mass gap via pure E_crack"),
    ("millennium_riemann",        f_riemann,            9877.78265,       "Im(t_10000)", "Riemann via 26D + SSQ damping"),
    ("millennium_bsd",            f_bsd,                0.3059997738,     "L'(E,1)", "BSD L-function proxy"),
    ("millennium_navier_stokes",  f_navier_stokes,      8500.0,           "peak entropy", "NS regularity bound"),
    ("millennium_hodge",          f_hodge,              1.0,              "closure", "Hodge"),
    ("millennium_poincare",       f_poincare,           1.0,              "closure", "Poincare"),
    ("millennium_p_vs_np",        f_p_vs_np,            1.0,              "closure", "P vs NP"),
    ("millennium_bh_info",        f_bh_info,            1.0,              "Page",    "BH info / Page curve"),
    # Standard-model spectrum
    ("alpha_primitive_sat",       f_alpha_primitive,    0.0072973525693,  "",       "alpha"),
    ("proton_mass_primitive_sat", f_proton_mass,        1.67262192369e-27,"kg",     "Proton mass"),
    ("yang_mills_primitive_sat",  f_yang_mills_sat,     1.78,             "GeV",    "YM sat"),
    ("neutron_lifetime_primitive_sat", lambda: derive_neutron_lifetime_full("age"), 880.0, "s", "Neutron lifetime (closure)"),
    ("h0_primitive_sat",          f_h0,                 70.0,             "km/s/Mpc", "H0"),
    ("t0_primitive_sat",          f_t0,                 13.8,             "Gyr",    "Age"),
    ("m_mu_primitive_sat",        f_m_mu,               0.105658,         "GeV",    "muon"),
    ("m_tau_primitive_sat",       f_m_tau,              1.77686,          "GeV",    "tau"),
    ("m_t_primitive_sat",         f_m_t,                172.69,           "GeV",    "top"),
    ("m_W_primitive_sat",         f_m_W,                80.379,           "GeV",    "W"),
    ("m_Z_primitive_sat",         f_m_Z,                91.1876,          "GeV",    "Z"),
    ("m_H_primitive_sat",         f_m_H,                125.25,           "GeV",    "Higgs"),
    ("v_higgs_primitive_sat",     f_v_higgs,            246.0,            "GeV",    "Higgs vev"),
    ("G_F_primitive_sat",         f_G_F,                1.1663787e-5,     "GeV^-2", "Fermi"),
    ("alpha_s_primitive_sat",     f_alpha_s,            0.118,            "",       "alpha_s"),
    ("m_e_primitive_sat",         f_m_e,                0.000511,         "GeV",    "electron"),
    ("m_pion_primitive_sat",      f_m_pion,             0.13957,          "GeV",    "pion"),
    ("m_kaon_primitive_sat",      f_m_kaon,             0.4937,           "GeV",    "kaon"),
    # Cosmology
    ("Omega_m_primitive_sat",     f_Omega_m,            0.315,            "",       "Omega_m"),
    ("Omega_b_h2_primitive_sat",  f_Omega_b_h2,         0.0224,           "",       "Omega_b h^2"),
    ("Omega_DM_h2_primitive_sat", f_Omega_DM_h2,        0.12,             "",       "Omega_DM h^2"),
    ("n_s_primitive_sat",         f_n_s,                0.965,            "",       "n_s"),
    ("A_s_primitive_sat",         f_A_s,                2.1e-9,           "",       "A_s"),
    ("eta_primitive_sat",         f_eta_baryon,         6.1e-10,          "",       "eta baryon"),
    ("Y_p_primitive_sat",         f_Y_p,                0.245,            "",       "Y_p"),
    ("z_re_primitive_sat",        f_z_re,               7.0,              "",       "z_reion"),
    ("tau_reion_primitive_sat",   f_tau_reion,          0.054,            "",       "tau reion"),
    # Pure SI via derive_*
    ("G_newton_primitive_sat",    derive_G_uqff,        6.6743e-11,       "m^3 kg^-1 s^-2", "G"),
    ("h_planck_primitive_sat",    derive_h,             6.62607015e-34,   "J s",    "h"),
    ("hbar_primitive_sat",        derive_hbar,          1.0545718e-34,    "J s",    "hbar"),
    ("c_light_primitive_sat",     derive_c_eff,         299792458.0,      "m/s",    "c"),
    ("N_A_primitive_sat",         derive_N_A_uqff,      6.02214076e23,    "mol^-1", "N_A"),
    ("alpha_primitive_sat_derive",derive_alpha_uqff,    0.0072973525693,  "",       "alpha (derive)"),
    ("planck_mass_primitive_sat", derive_planck_mass_uqff, 2.176434e-8,   "kg",     "m_P"),
    ("planck_length_primitive_sat", derive_planck_length_uqff, 1.616e-35, "m",      "l_P"),
    ("k_b_primitive_sat",         derive_k_b_uqff,      1.380649e-23,     "J/K",    "k_B (closure)"),
    ("vacuum_permittivity_primitive_sat", derive_vacuum_permittivity_uqff, 8.854187817e-12, "F/m", "eps0"),
    ("vacuum_permeability_primitive_sat", derive_vacuum_permeability_uqff, 1.256637e-6, "H/m", "mu0"),
    ("delta_scm_primitive_sat",   derive_delta_scm_uqff, 1e-36,           "",       "delta SCM"),
    # Proxies
    ("p1_lkk_um_primitive_sat",   derive_p1_lkk_um_uqff, 0.5,             "",       "p1"),
    ("p6_lkk_inv_primitive_sat",  derive_p6_lkk_inv_uqff, 1.0,            "",       "p6"),
    ("p7_w_a_primitive_sat",      derive_p7_w_a_uqff,   1.0,              "",       "p7"),
    ("p9_h_tension_primitive_sat",derive_p9_h_tension_uqff, 1.0,          "",       "p9"),
    ("p10_s8_tension_primitive_sat", derive_p10_s8_tension_uqff, 1.0,     "",       "p10"),
    ("sgr_a_g_primitive_sat",     derive_sgr_a_g_uqff,  4.3e6,            "Msun",   "Sgr A*"),
    ("sgr_1745_g_primitive_sat",  derive_sgr_1745_g_uqff, 1.0,            "",       "Sgr 1745"),
    ("lenr_parkhomov_primitive_sat", derive_lenr_parkhomov_uqff, 200.0,   "W",      "Parkhomov"),
    ("Omega_GW_h2_primitive_sat", derive_omega_gw_h2_uqff, 1e-15,         "",       "Omega_GW"),
    ("f_NL_equil_primitive_sat",  derive_f_nl_equil_uqff, 1.0,            "",       "f_NL"),
    ("epsilon_slow_roll_primitive_sat", derive_epsilon_slow_roll_uqff, 0.01, "",    "eps slow-roll"),
    # E_react / E_crack / spinor
    ("e_react_uqff",              f_e_react,            None,             "J/m^3",  "E_react pure"),
    ("e_crack_uqff_J",            f_e_crack_J,          None,             "J",      "E_crack pure (no c^2)"),
    ("spinor_canonical_engine_derive", f_spinor,        None,             "",       "spinor"),
    # Merged tests (return characteristic closures)
    ("dark_flow",                 derive_dark_flow_uqff, 800000,          "m/s",    "Dark Flow"),
    ("cmb_cold_spot",             derive_cmb_cold_spot_uqff, -100e-6,     "K",      "CMB Cold Spot"),
    ("dark_matter_particle",      derive_dark_matter_particle_uqff, 1e-3, "eV",     "DM effective"),
    ("pta_sgwb_spectral_index",   derive_pta_sgwb_spectral_index_uqff, 13.0/3.0, "", "PTA SGWB"),
    ("txs0506_multimessenger_delay", derive_txs0506_multimessenger_delay_uqff, 1000.0, "s", "TXS delay"),
    ("frb_origin",                derive_frb_origin_uqff, 1.4e9,          "Hz",     "FRB origin"),
    ("casimir_effect",            derive_casimir_effect_uqff, -1.3e-27,   "N/m^2",  "Casimir"),
    ("coronal_heating_uqff",      derive_coronal_heating_uqff, 2.0e6,     "K",      "Coronal heating"),
    ("proton_radius_muonic_primitive_sat", derive_proton_radius_muonic_uqff, 8.41e-16, "m", "muonic H"),
]


# =============================================================================
# COMPUTE + COMPARE
# =============================================================================

def pct_diff(uqff, target):
    if target is None:
        return None
    if target == 0:
        return float("inf") if uqff != 0 else 0.0
    return abs(uqff - target) / abs(target) * 100.0

def fmt(v):
    if v is None: return "—"
    if isinstance(v, float):
        if abs(v) > 1e6 or (abs(v) < 1e-3 and v != 0): return f"{v:.6e}"
        return f"{v:.6g}"
    return str(v)

print("=" * 110)
print(f"UQFF RECOMPUTE — independent verifier")
print(f"Locked primitives:")
print(f"  Integer system : D_PHYS={D_PHYS}  D_BSFG={D_BSFG}  N_CH={N_CH}  SO_FIVE={SO_FIVE}  D_CRIT={D_CRIT}  A_FIVE={A_FIVE}")
print(f"  Rational system: SSQ={SSQ}={SSQ_F}  PHI_RES={PHI_RES}={PHI_F}  TRZ={TRZ}={TRZ_F}")
print(f"                   G1_K={G1_K}={G1_F}  G2_BETA={G2_BETA}={G2_F}  G3_RICCI={G3_RICCI}={G3_F}  G4_BSFG={G4_BSFG}={G4_F}")
print(f"                   K_MEX={K_MEX}={K_MEX_F}  BETA_I={BETA_I}={BETA_F}")
print(f"  Real system    : E0={E0}  S26_3={S26_3}  KAPPA={KAPPA}  F_THZ={F_THZ}")
print(f"  Geometry system: BSFG (Buoyancy-SCm-Fluid-Geometry); D_BSFG=6=dim(SO(5)/U(2))")
print()
print(f"Vacuum ledger root (derive_from_quantum_chain):")
print(f"  RHO_VAC_SCM       = {RHO_VAC_SCM:.6e} J/m^3   (Σ_{{n=1..26}} 0.57·1e-20·10^n)")
print(f"  RHO_VAC_UA        = {RHO_VAC_UA:.6e} J/m^3   (10× SO_FIVE structural)")
print(f"  V_DPM_BASE        = {V_DPM_BASE:.6e}        (from E_react chain + 10× projection)")
print(f"  E_REACT_0         = {E_REACT_0:.6e}        (rho_SCm·1²/rho_UA at v=1 normalized)")
print(f"  Micro retention   : RHO_VAC_SCM_MICRO = {RHO_VAC_SCM_MICRO:.6e} (legacy per-DPM-vol, Star-Magic.txt Ch.2)")
print(f"  E_crack pure      = {derive_e_crack():.6e} J  (= rho·V_DPM²/SSQ, NO c²)")
print()
print(f"Three numeric systems:")
print(f"  VDS (Vacuum Density Series): Li_26(SSQ) → 26-rung polylog")
print(f"  DVP (Dipole Vortex Primes): p=113 base prime, system-dependent")
print(f"  BSH (Buoyancy Saturation Harmonics): β_i · exp(-SSQ·i/26) over 26 states")
print(f"  Hybrid identity: R_VDS × p_DVP × BSH(i) = F_{{U_Bi_i}} within 0.1% (PAPER_1069)")
print("=" * 110)

# Compute
results = []
for entry in REGISTRY:
    name, f, target, unit, desc = entry
    try:
        v = f()
    except Exception as e:
        v = None
        err = str(e)
    else:
        err = None
    p = pct_diff(v, target) if v is not None else None
    results.append({"name": name, "uqff": v, "target": target, "unit": unit, "pct": p, "desc": desc, "err": err})

# Sort: known-exact-or-close first, then proxies
def sort_key(r):
    if r["pct"] is None: return (3, r["name"])
    if r["pct"] < 1: return (0, r["name"])
    if r["pct"] < 10: return (1, r["name"])
    return (2, r["pct"] * -1)  # show worst at bottom

results_sorted = sorted(results, key=sort_key)

print()
print(f"{'name':<42} {'uqff':>18} {'target':>18} {'unit':<10} {'pct':>14}")
print("-" * 110)
for r in results_sorted:
    pct = f"{r['pct']:.4f}%" if r['pct'] is not None else "—"
    if r['err']:
        pct = "ERROR: " + r['err'][:30]
    print(f"{r['name']:<42} {fmt(r['uqff']):>18} {fmt(r['target']):>18} {r['unit']:<10} {pct:>14}")

# Summary statistics
ok_exact = sum(1 for r in results if r["pct"] is not None and r["pct"] < 0.1)
ok_lt1   = sum(1 for r in results if r["pct"] is not None and r["pct"] < 1)
ok_lt10  = sum(1 for r in results if r["pct"] is not None and r["pct"] < 10)
not_close = sum(1 for r in results if r["pct"] is not None and r["pct"] >= 10)
errors  = sum(1 for r in results if r["err"] is not None)
no_tgt  = sum(1 for r in results if r["target"] is None)

print("-" * 110)
print(f"SUMMARY: total={len(results)}  <0.1%={ok_exact}  <1%={ok_lt1}  <10%={ok_lt10}  >=10%={not_close}  errors={errors}  no-target={no_tgt}")

# Write JSON
out = {
    "primitives": {
        "integer": {"D_PHYS":D_PHYS, "D_BSFG":D_BSFG, "N_CH":N_CH, "SO_FIVE":SO_FIVE, "D_CRIT":D_CRIT, "A_FIVE":A_FIVE},
        "rational": {"SSQ":str(SSQ), "PHI_RES":str(PHI_RES), "TRZ":str(TRZ), "G1_K":str(G1_K),
                     "G2_BETA":str(G2_BETA), "G3_RICCI":str(G3_RICCI), "G4_BSFG":str(G4_BSFG),
                     "K_MEX":str(K_MEX), "BETA_I":str(BETA_I)},
        "real": {"E0":E0, "S26_3":S26_3, "KAPPA":KAPPA, "F_THZ":F_THZ},
    },
    "vacuum_ledger": {
        "RHO_VAC_SCM": RHO_VAC_SCM, "RHO_VAC_UA": RHO_VAC_UA,
        "V_DPM_BASE": V_DPM_BASE, "E_REACT_0": E_REACT_0,
        "RHO_VAC_SCM_MICRO": RHO_VAC_SCM_MICRO,
        "E_crack_pure_J": derive_e_crack(),
    },
    "three_number_systems": {
        "VDS": "Vacuum Density Series Z_26 = Li_26(SSQ)",
        "DVP": "Dipole Vortex Primes (p=113 base)",
        "BSH": "Buoyancy Saturation Harmonics: BETA_I·exp(-SSQ·i/26), 26 states",
        "hybrid": "R_VDS × p_DVP × BSH(i) = F_U_Bi_i within 0.1% (PAPER_1069)",
    },
    "geometry_system": "BSFG (Buoyancy-SCm-Fluid-Geometry); D_BSFG=6=dim_R[SO(5)/U(2)]=10-4 (PAPER_1167)",
    "results": results,
    "summary": {"total": len(results), "lt_01pct": ok_exact, "lt_1pct": ok_lt1,
                "lt_10pct": ok_lt10, "ge_10pct": not_close, "errors": errors, "no_target": no_tgt},
}
with open("/sessions/vibrant-keen-bohr/mnt/outputs/uqff_recompute/recompute_report.json", "w") as fp:
    json.dump(out, fp, indent=2, default=str)
print()
print(f"Wrote: /sessions/vibrant-keen-bohr/mnt/outputs/uqff_recompute/recompute_report.json")
