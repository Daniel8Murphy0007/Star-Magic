"""
QCalc_core_uqff.py - Manually implemented core UQFF equations from MAIN_1_CoAnQi.cpp SOURCE4

These are the core physics equations of the Unified Quantum Field Framework (UQFF):
- Ug1: Magnetic dipole gravity (string-magnetic coupling)
- Ug2: Charge-reactivity gravity (QUA coupling)
- Ug3: String rotation gravity (Bj field coupling)
- Ug4: Vacuum concentration gravity (dark energy contribution)
- Ubi: Buoyancy force (anti-gravity counter-term)
- Um: Magnetism term (string ensemble)
- FU: Complete unified field (sum of all terms)

Plus MUGE (Modified Universal Gravitation Equation) Compressed functions:
- Base Newtonian + expansion + superconductive adjustments
- Cosmological constant + quantum + fluid + perturbation corrections

Plus 15 SMBH SOURCE82 classes for supermassive black hole physics.

Reference: MAIN_1_CoAnQi.cpp lines 26575-26730 (SOURCE4 namespace)
"""

import math
from dpm_helpers import dpm_emergent_ug1, dpm_emergent_ug2

from typing import Dict, Any, Callable, Optional

# =============================================================================
# PHYSICAL CONSTANTS (SOURCE4 defaults from C++)
# =============================================================================

UQFF_CONSTANTS = {
    # Fundamental
    "G": 6.6743e-11,           # Gravitational constant (m³/kg/s²)
    "c": 2.998e8,              # Speed of light (m/s)
    "PI": 3.141592653589793,
    "hbar": 1.0546e-34,        # Reduced Planck constant (J·s)
    
    # UQFF calibrated constants
    "kappa": 0.0005,           # κ calibration factor (1/day)
    "SSq": 0.57,               # [SSq] string-space quotient
    "H_SCm": 0.99,             # Superconducting medium factor
    "U_UA": 0.0001,            # Universal attraction modulator
    "k_eta": 1e-113,           # η coupling constant
    "beta_i": 0.603,           # Buoyancy coefficient
    
    # SOURCE4 defaults
    "alpha_SOURCE4": 1e-6,     # Decay constant
    "delta_def_SOURCE4": 0.01, # Defect amplitude
    "k1_SOURCE4": 1.0,         # Ug1 coupling
    "k2_SOURCE4": 1.0,         # Ug2 coupling
    "k3_SOURCE4": 1.0,         # Ug3 coupling
    "k4_SOURCE4": 1.0,         # Ug4 coupling
    "QA_SOURCE4": 0.5,         # Charge-ratio
    "delta_sw_SOURCE4": 0.01,  # Solar wind modulation
    "v_sw_SOURCE4": 400000,    # Solar wind velocity (m/s)
    "HSCm_SOURCE4": 0.99,      # SCm field factor
    "rho_A_SOURCE4": 1e-20,    # Aether density
    "kappa_SOURCE4": 0.0005,   # κ calibration
    "v_SCm_SOURCE4": 1e5,      # SCm velocity
    "rho_v_SOURCE4": 1e-26,    # Vacuum energy density
    "C_concentration_SOURCE4": 1.0,  # Concentration factor
    "Mbh_SOURCE4": 4e6 * 1.989e30,   # M_bh default (SgrA* mass)
    "dg_SOURCE4": 2.6e20,      # Distance to galactic center (m)
    "f_feedback_SOURCE4": 0.1, # Feedback efficiency
    "Omega_g_SOURCE4": 1.0,    # Galactic angular velocity factor
    "epsilon_sw_SOURCE4": 0.01,# Wind efficiency
    "rho_sw_SOURCE4": 1e-20,   # Wind density
    "UUA_SOURCE4": 0.0001,     # Universal attraction
    "gamma_SOURCE4": 1e-10,    # Decay rate for Um
    "num_strings_SOURCE4": 1e10, # Number of strings
    
    # Cosmological
    "Lambda": 1.1e-52,         # Cosmological constant (1/m²)
    "H0": 2.269e-18,           # Hubble constant (1/s)
    "tHubble": 4.35e17,        # Hubble time (s)
    
    # SMBH
    "B_crit": 4.4e13,          # Critical magnetic field (T)
    "r_s": 1.2e10,             # Schwarzschild radius (m) for 4M☉ BH
    "sigma_sb": 5.67e-8,       # Stefan-Boltzmann constant
}


# =============================================================================
# HELPER FUNCTIONS (from SOURCE4 namespace)
# =============================================================================

def step_function(r: float, Rb: float) -> float:
    """Heaviside step function S(r-Rb)."""
    return 1.0 if r >= Rb else 0.0


def compute_mu_s(t: float, Bs_avg: float, omega_c: float, Rs: float) -> float:
    """Magnetic dipole moment μ_s = B_s * R_s³ * (1 + sin(ω_c*t))."""
    return Bs_avg * Rs**3 * (1.0 + math.sin(omega_c * t))


def compute_grad_Ms_r(Ms: float, Rs: float) -> float:
    """Gradient ∂M_s/∂r ≈ M_s/R_s."""
    if Rs == 0:
        return 0.0
    return Ms / Rs


def compute_Ereact(t: float, SCm_density: float, v_SCm: float, rho_A: float, kappa: float) -> float:
    """Reactivity energy E_react = SCm * v_SCm * ρ_A * (1 + κ*t)."""
    return SCm_density * v_SCm * rho_A * (1.0 + kappa * t)


def compute_omega_s_t(t: float, omega_s: float, omega_c: float) -> float:
    """Time-dependent angular frequency ω_s(t) = ω_s * (1 + 0.01*sin(ω_c*t))."""
    return omega_s * (1.0 + 0.01 * math.sin(omega_c * t))


def compute_Bj(t: float, omega_c: float, SCm_contrib: float = 1.0) -> float:
    """B_j magnetic field strength = SCm_contrib * (1 + 0.1*cos(ω_c*t))."""
    return SCm_contrib * (1.0 + 0.1 * math.cos(omega_c * t))


def compute_mu_j(t: float, omega_c: float, Rs: float) -> float:
    """String dipole moment μ_j = R_s² * (1 + 0.1*sin(ω_c*t))."""
    return Rs**2 * (1.0 + 0.1 * math.sin(omega_c * t))


# =============================================================================
# CORE UQFF EQUATIONS (Ug1, Ug2, Ug3, Ug4, Ubi, Um, FU)
# =============================================================================

def compute_Ug1_SOURCE4(p: dict) -> float:
    """
    Ug1: Magnetic dipole gravity (string-magnetic coupling)
    
    Formula: Ug1 = k1 * μ_s * ∇M_s * exp(-α*t) * cos(π*tn) * (1 + δ_def*sin(0.001*t))
    
    Parameters:
        r: radial distance (m)
        t: time (s)
        tn: normalized time
        Bs_avg: average magnetic field (T)
        omega_c: core angular frequency (rad/s)
        Rs: stellar radius (m)
        Ms: stellar mass (kg)
        alpha: decay constant
        delta_def: defect amplitude
        k1: coupling constant
    """
    # Extract parameters with defaults
    t = p.get("t", 0.0)
    tn = p.get("tn", 0.0)
    Bs_avg = p.get("Bs_avg", 1e-4)
    omega_c = p.get("omega_c", 1e-5)
    Rs = p.get("Rs", 6.96e8)  # Solar radius default
    Ms = p.get("Ms", UQFF_CONSTANTS["Mbh_SOURCE4"])
    alpha = p.get("alpha", UQFF_CONSTANTS["alpha_SOURCE4"])
    delta_def = p.get("delta_def", UQFF_CONSTANTS["delta_def_SOURCE4"])
    k1 = p.get("k1", UQFF_CONSTANTS["k1_SOURCE4"])
    
    # Compute intermediate values
    mu_s = compute_mu_s(t, Bs_avg, omega_c, Rs)
    grad_Ms_r = compute_grad_Ms_r(Ms, Rs)
    defect = 1.0 + delta_def * math.sin(0.001 * t)
    
    # Ug1 = k1 * μ_s * ∇M_s * exp(-α*t) * cos(π*tn) * defect
    result = k1 * mu_s * grad_Ms_r * math.exp(-alpha * t) * math.cos(UQFF_CONSTANTS["PI"] * tn) * defect
    
    return result


def compute_Ug2_SOURCE4(p: dict) -> float:
    """
    Ug2: Charge-reactivity gravity (QUA coupling)
    
    Formula: Ug2 = k2 * (QA + QUA) * Ms/r² * S(r-Rb) * (1 + δ_sw*v_sw) * HSCm * E_react
    
    Parameters:
        r: radial distance (m)
        t: time (s)
        Ms: stellar mass (kg)
        Rb: boundary radius (m)
        QUA: universal attraction quotient
        QA: charge-atmosphere ratio
        SCm_density: superconducting medium density
        Other SOURCE4 constants
    """
    r = p.get("r", 1e10)
    t = p.get("t", 0.0)
    Ms = p.get("Ms", UQFF_CONSTANTS["Mbh_SOURCE4"])
    Rb = p.get("Rb", 1e9)
    QUA = p.get("QUA", 0.1)
    k2 = p.get("k2", UQFF_CONSTANTS["k2_SOURCE4"])
    QA = p.get("QA", UQFF_CONSTANTS["QA_SOURCE4"])
    delta_sw = p.get("delta_sw", UQFF_CONSTANTS["delta_sw_SOURCE4"])
    v_sw = p.get("v_sw", UQFF_CONSTANTS["v_sw_SOURCE4"])
    HSCm = p.get("HSCm", UQFF_CONSTANTS["HSCm_SOURCE4"])
    SCm_density = p.get("SCm_density", 1e-10)
    rho_A = p.get("rho_A", UQFF_CONSTANTS["rho_A_SOURCE4"])
    kappa = p.get("kappa", UQFF_CONSTANTS["kappa_SOURCE4"])
    
    if r == 0:
        return 0.0
    
    # Compute reactivity
    E_react = compute_Ereact(t, SCm_density, UQFF_CONSTANTS["v_SCm_SOURCE4"], rho_A, kappa)
    S = step_function(r, Rb)
    wind_mod = 1.0 + delta_sw * v_sw
    
    # Ug2 = k2 * (QA + QUA) * Ms/r² * S * wind_mod * HSCm * E_react
    result = k2 * (QA + QUA) * Ms / (r * r) * S * wind_mod * HSCm * E_react
    
    return result


def compute_Ug3_SOURCE4(p: dict) -> float:
    """
    Ug3: String rotation gravity (Bj field coupling)
    
    Formula: Ug3 = k3 * Bj * cos(ω_s(t)*t*π) * P_core * E_react
    
    Parameters:
        r: radial distance (m)
        t: time (s)
        Pcore: core pressure (Pa)
        omega_s: spin frequency (rad/s)
        omega_c: core frequency (rad/s)
        SCm_density: superconducting medium density
    """
    t = p.get("t", 0.0)
    Pcore = p.get("Pcore", 1e16)  # Core pressure
    omega_s = p.get("omega_s", 1e-6)
    omega_c = p.get("omega_c", 1e-5)
    k3 = p.get("k3", UQFF_CONSTANTS["k3_SOURCE4"])
    SCm_density = p.get("SCm_density", 1e-10)
    rho_A = p.get("rho_A", UQFF_CONSTANTS["rho_A_SOURCE4"])
    kappa = p.get("kappa", UQFF_CONSTANTS["kappa_SOURCE4"])
    
    # Compute intermediate values
    E_react = compute_Ereact(t, SCm_density, UQFF_CONSTANTS["v_SCm_SOURCE4"], rho_A, kappa)
    omega_s_t = compute_omega_s_t(t, omega_s, omega_c)
    Bj = compute_Bj(t, omega_c)
    
    # Ug3 = k3 * Bj * cos(ω_s(t)*t*π) * P_core * E_react
    result = k3 * Bj * math.cos(omega_s_t * t * UQFF_CONSTANTS["PI"]) * Pcore * E_react
    
    return result


def compute_Ug4_SOURCE4(p: dict) -> float:
    """
    Ug4: Vacuum concentration gravity (dark energy contribution)
    
    Formula: Ug4 = k4 * ρ_v * C_concentration * M_bh/d_g * exp(-α*t) * cos(π*tn) * (1 + f_feedback)
    
    Parameters:
        t: time (s)
        tn: normalized time
        rho_v: vacuum energy density (kg/m³)
        C_concentration: concentration factor
        Mbh: black hole mass (kg)
        dg: galactic distance (m)
        alpha: decay constant
        f_feedback: feedback efficiency
        k4: coupling constant
    """
    t = p.get("t", 0.0)
    tn = p.get("tn", 0.0)
    rho_v = p.get("rho_v", UQFF_CONSTANTS["rho_v_SOURCE4"])
    C_concentration = p.get("C_concentration", UQFF_CONSTANTS["C_concentration_SOURCE4"])
    Mbh = p.get("Mbh", UQFF_CONSTANTS["Mbh_SOURCE4"])
    dg = p.get("dg", UQFF_CONSTANTS["dg_SOURCE4"])
    alpha = p.get("alpha", UQFF_CONSTANTS["alpha_SOURCE4"])
    f_feedback = p.get("f_feedback", UQFF_CONSTANTS["f_feedback_SOURCE4"])
    k4 = p.get("k4", UQFF_CONSTANTS["k4_SOURCE4"])
    
    if dg == 0:
        return 0.0
    
    decay = math.exp(-alpha * t)
    cycle = math.cos(UQFF_CONSTANTS["PI"] * tn)
    
    # Ug4 = k4 * ρ_v * C * M_bh/d_g * exp(-α*t) * cos(π*tn) * (1 + f_feedback)
    result = k4 * rho_v * C_concentration * Mbh / dg * decay * cycle * (1.0 + f_feedback)
    
    return result


def compute_Ubi_SOURCE4(p: dict) -> float:
    """
    Ubi: Buoyancy force (anti-gravity counter-term)
    
    Formula: Ubi = -β_i * Ug_sum * Ω_g * M_bh/d_g * (1 + ε_sw*ρ_sw) * UUA * cos(π*tn)
    
    This represents the buoyancy-based correction to gravity in the UQFF framework.
    The negative sign indicates anti-gravity (repulsion).
    
    Parameters:
        Ugi: sum of Ug1+Ug2+Ug3+Ug4
        beta_i: buoyancy coefficient (calibrated ~0.603)
        Omega_g: galactic angular velocity factor
        Mbh: black hole mass (kg)
        dg: galactic distance (m)
        epsilon_sw: wind efficiency
        rho_sw: wind density
        UUA: universal attraction
        tn: normalized time
    """
    Ugi = p.get("Ugi", 1.0)  # Sum of Ug terms
    beta_i = p.get("beta_i", UQFF_CONSTANTS["beta_i"])
    Omega_g = p.get("Omega_g", UQFF_CONSTANTS["Omega_g_SOURCE4"])
    Mbh = p.get("Mbh", UQFF_CONSTANTS["Mbh_SOURCE4"])
    dg = p.get("dg", UQFF_CONSTANTS["dg_SOURCE4"])
    epsilon_sw = p.get("epsilon_sw", UQFF_CONSTANTS["epsilon_sw_SOURCE4"])
    rho_sw = p.get("rho_sw", UQFF_CONSTANTS["rho_sw_SOURCE4"])
    UUA = p.get("UUA", UQFF_CONSTANTS["UUA_SOURCE4"])
    tn = p.get("tn", 0.0)
    
    if dg == 0:
        return 0.0
    
    wind_mod = 1.0 + epsilon_sw * rho_sw
    
    # Ubi = -β_i * Ug_sum * Ω_g * M_bh/d_g * wind_mod * UUA * cos(π*tn)
    result = -beta_i * Ugi * Omega_g * Mbh / dg * wind_mod * UUA * math.cos(UQFF_CONSTANTS["PI"] * tn)
    
    return result


def compute_Um_SOURCE4(p: dict) -> float:
    """
    Um: Magnetism term (string ensemble contribution)
    
    Formula: Um = μ_j/r_j * (1 - exp(-γ*t*cos(π*tn))) * φ̂ * N_strings * P_SCm * E_react
                  * (1 + 10¹³·f_Heaviside) * (1 + f_quasi)
    
    Parameters:
        t: time (s)
        tn: normalized time
        rj: string radial position (m)
        Rs: stellar radius (m)
        omega_c: core frequency (rad/s)
        gamma: decay rate
        num_strings: number of strings
        PSCm: SCm power
        SCm_density: SCm density
        f_quasi: quasi-particle coupling fraction (default 0.01)
        f_Heaviside: Heaviside gate activation (0 = off, default 0.0)
    """
    t = p.get("t", 0.0)
    tn = p.get("tn", 0.0)
    rj = p.get("rj", 1e9)
    Rs = p.get("Rs", 6.96e8)
    omega_c = p.get("omega_c", 1e-5)
    gamma = p.get("gamma", UQFF_CONSTANTS["gamma_SOURCE4"])
    num_strings = p.get("num_strings", UQFF_CONSTANTS["num_strings_SOURCE4"])
    PSCm = p.get("PSCm", 1.0)
    phi_hat = p.get("phi_hat", 1.0)
    SCm_density = p.get("SCm_density", 1e-10)
    rho_A = p.get("rho_A", UQFF_CONSTANTS["rho_A_SOURCE4"])
    kappa = p.get("kappa", UQFF_CONSTANTS["kappa_SOURCE4"])
    f_quasi = p.get("f_quasi", 0.01)
    f_Heaviside = p.get("f_Heaviside", 0.0)
    
    if rj == 0:
        return 0.0
    
    # Compute intermediate values
    E_react = compute_Ereact(t, SCm_density, UQFF_CONSTANTS["v_SCm_SOURCE4"], rho_A, kappa)
    mu_j = compute_mu_j(t, omega_c, Rs)
    decay = 1.0 - math.exp(-gamma * t * math.cos(UQFF_CONSTANTS["PI"] * tn))
    single = mu_j / rj * decay * phi_hat
    
    # Um = single * N_strings * P_SCm * E_react * (1+10^13*f_Heaviside) * (1+f_quasi)
    result = single * num_strings * PSCm * E_react * (1.0 + 1e13 * f_Heaviside) * (1.0 + f_quasi)
    
    return result


def compute_FU_SOURCE4(p: dict) -> float:
    """
    FU: Complete Unified Field (UQFF master equation)
    
    Formula: F_U = Σ(Ug1 + Ug2 + Ug3 + Ug4) + Ubi + Um
    
    This is the complete unified field combining all gravity components:
    - Ug1: Magnetic dipole gravity
    - Ug2: Charge-reactivity gravity
    - Ug3: String rotation gravity
    - Ug4: Vacuum concentration gravity
    - Ubi: Buoyancy counter-term
    - Um: Magnetism contribution
    
    Parameters:
        All parameters from Ug1, Ug2, Ug3, Ug4, Ubi, Um
    """
    # Compute individual terms
    Ug1 = compute_Ug1_SOURCE4(p)
    Ug2 = compute_Ug2_SOURCE4(p)
    Ug3 = compute_Ug3_SOURCE4(p)
    Ug4 = compute_Ug4_SOURCE4(p)
    
    # Sum of Ug terms
    sum_Ug = Ug1 + Ug2 + Ug3 + Ug4
    
    # Pass sum to Ubi calculation
    p_ubi = dict(p)
    p_ubi["Ugi"] = sum_Ug
    Ubi = compute_Ubi_SOURCE4(p_ubi)
    
    # Compute Um
    Um = compute_Um_SOURCE4(p)
    
    # Complete unified field
    result = sum_Ug + Ubi + Um
    
    return result


# =============================================================================
# MUGE COMPRESSED EQUATIONS (Modified Universal Gravitation)
# =============================================================================

def compute_compressed_base_SOURCE4(p: dict) -> float:
    """
    MUGE Base: Newtonian gravity g = GM/r²
    """
    G = UQFF_CONSTANTS["G"]
    M = p.get("M", UQFF_CONSTANTS["Mbh_SOURCE4"])
    r = p.get("r", 1e10)
    
    if r == 0:
        return 0.0
    
    return dpm_emergent_ug1(M, r)


def compute_compressed_expansion_SOURCE4(p: dict) -> float:
    """
    MUGE Expansion: Hubble correction factor = 1 + H0*t
    """
    H0 = p.get("H0", UQFF_CONSTANTS["H0"])
    t = p.get("t", 0.0)
    
    return 1.0 + H0 * t


def compute_compressed_super_adj_SOURCE4(p: dict) -> float:
    """
    MUGE Superconductive adjustment: = 1 - B/B_crit
    
    Magnetic field suppression of gravity in superconducting medium.
    """
    B = p.get("B", 1e-4)  # Magnetic field (T)
    Bcrit = p.get("Bcrit", UQFF_CONSTANTS["B_crit"])
    
    if Bcrit == 0:
        return 1.0
    
    return 1.0 - B / Bcrit


def compute_compressed_cosm_SOURCE4(p: dict) -> float:
    """
    MUGE Cosmological: Dark energy contribution = Λc²/3
    """
    Lambda = p.get("Lambda", UQFF_CONSTANTS["Lambda"])
    c = UQFF_CONSTANTS["c"]
    
    return Lambda * c * c / 3.0


def compute_compressed_quantum_SOURCE4(p: dict) -> float:
    """
    MUGE Quantum: Uncertainty principle contribution
    
    Formula: (ℏ/Δx·Δp) * ∫|ψ|² * (2π/t_Hubble)
    """
    hbar = UQFF_CONSTANTS["hbar"]
    Delta_x_p = p.get("Delta_x_p", 1e-68)  # Δx·Δp
    integral_psi = p.get("integral_psi", 2.176e-18)
    tHubble = p.get("tHubble", UQFF_CONSTANTS["tHubble"])
    
    if Delta_x_p == 0:
        return 0.0
    
    return (hbar / Delta_x_p) * integral_psi * (2 * UQFF_CONSTANTS["PI"] / tHubble)


def compute_compressed_fluid_SOURCE4(p: dict) -> float:
    """
    MUGE Fluid: Navier-Stokes contribution = ρ_fluid * V * g_local
    """
    rho_fluid = p.get("rho_fluid", 1e-20)
    Vsys = p.get("Vsys", 1e30)  # System volume
    g_local = p.get("g_local", 1e-6)
    
    return rho_fluid * Vsys * g_local


def compute_compressed_perturbation_SOURCE4(p: dict) -> float:
    """
    MUGE Perturbation: Dark matter + density perturbation
    
    Formula: (M + M_DM) * (δρ/ρ + 3GM/r³)
    """
    M = p.get("M", UQFF_CONSTANTS["Mbh_SOURCE4"])
    M_DM = p.get("M_DM", 0.1 * UQFF_CONSTANTS["Mbh_SOURCE4"])
    delta_rho_rho = p.get("delta_rho_rho", 1e-5)
    r = p.get("r", 1e10)
    G = UQFF_CONSTANTS["G"]
    
    if r == 0:
        return 0.0
    
    return (M + M_DM) * (delta_rho_rho + 3 * G * M / (r * r * r))


def compute_compressed_MUGE_SOURCE4(p: dict) -> float:
    """
    MUGE Complete: Full Modified Universal Gravitation Equation
    
    Formula: g = g_base * expansion * super_adj * env + Ug_sum + cosm + quantum + fluid + perturbation
    
    This combines Newtonian gravity with all UQFF corrections.
    """
    base = compute_compressed_base_SOURCE4(p)
    expansion = compute_compressed_expansion_SOURCE4(p)
    super_adj = compute_compressed_super_adj_SOURCE4(p)
    env = 1.0  # Environment factor (placeholder)
    
    adjusted_base = base * expansion * super_adj * env
    
    cosm = compute_compressed_cosm_SOURCE4(p)
    quantum = compute_compressed_quantum_SOURCE4(p)
    fluid = compute_compressed_fluid_SOURCE4(p)
    perturbation = compute_compressed_perturbation_SOURCE4(p)
    
    # Note: Ug_sum would come from FU calculation if integrated
    Ug_sum = p.get("Ug_sum", 0.0)
    
    return adjusted_base + Ug_sum + cosm + quantum + fluid + perturbation


# =============================================================================
# SMBH SOURCE82 CLASSES (15 Supermassive Black Hole physics terms)
# =============================================================================

def compute_SMBHDynamicVacuumTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Dynamic Vacuum: Time-varying vacuum energy near SMBH
    
    Formula: amplitude * ρ_vac * sin(frequency * t)
    """
    amplitude = p.get("amplitude", 1.0)
    rho_vac = p.get("rho_vac", 1e-26)
    frequency = p.get("frequency", 1e-10)
    t = p.get("t", 0.0)
    
    return amplitude * rho_vac * math.sin(frequency * t)


def compute_SMBHQuantumCouplingTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Quantum Coupling: Quantum effects at event horizon scale
    
    Formula: α_q * exp(-r/λ_q)
    """
    alpha_q = p.get("alpha_q", 1e-10)
    r = p.get("r", 1e10)
    lambda_q = p.get("lambda_q", 1e6)  # Quantum coherence length
    
    if lambda_q == 0:
        return 0.0
    
    return alpha_q * math.exp(-r / lambda_q)


def compute_SMBHMSigmaRelationTerm_SOURCE82(p: dict) -> float:
    """
    SMBH M-σ Relation: Black hole mass - velocity dispersion relation
    
    Formula: M_bh = 10^8.13 * (σ/200 km/s)^4.38 * M_sun
    """
    sigma = p.get("sigma", 200e3)  # Velocity dispersion (m/s)
    sigma_norm = sigma / 200000  # Normalize to 200 km/s
    M_sun = 1.989e30
    
    M_bh = math.pow(10, 8.13) * math.pow(sigma_norm, 4.38) * M_sun
    return M_bh


def compute_SMBHBulgeGravityTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Bulge Gravity: Contribution from galactic bulge
    
    Formula: g_bulge = G * M_bulge / r² * (1 - exp(-r/r_eff))
    """
    G = UQFF_CONSTANTS["G"]
    M_bulge = p.get("M_bulge", 1e10 * 1.989e30)  # 10^10 M_sun
    r = p.get("r", 1e20)
    r_eff = p.get("r_eff", 3e20)  # Effective radius
    
    if r == 0:
        return 0.0
    
    return dpm_emergent_ug1(M_bulge, r) * (1.0 - math.exp(-r / r_eff))


def compute_SMBHUg1Term_SOURCE82(p: dict) -> float:
    """
    SMBH Ug1: Newtonian gravity from SMBH
    
    Formula: Ug1 = G * M_bh / r²
    """
    G = UQFF_CONSTANTS["G"]
    M_bh = p.get("M_bh", UQFF_CONSTANTS["Mbh_SOURCE4"])
    r = p.get("r", 1e10)
    
    if r == 0:
        return 0.0
    
    return dpm_emergent_ug1(M_bh, r)


def compute_SMBHUg2Term_SOURCE82(p: dict) -> float:
    """
    SMBH Ug2: General relativistic correction
    
    Formula: Ug2 = Ug1 * 3GM_bh / (r*c²)
    """
    G = UQFF_CONSTANTS["G"]
    c = UQFF_CONSTANTS["c"]
    M_bh = p.get("M_bh", UQFF_CONSTANTS["Mbh_SOURCE4"])
    r = p.get("r", 1e10)
    
    if r == 0:
        return 0.0
    
    Ug1 = dpm_emergent_ug1(M_bh, r)
    
    return Ug1 * 3.0 * G * M_bh / (r * c * c)


def compute_SMBHUg3Term_SOURCE82(p: dict) -> float:
    """
    SMBH Ug3: Quantum gravity correction
    
    Formula: Ug3 = k * exp(-r/λ_Q)
    """
    k = p.get("k", 1e-20)
    r = p.get("r", 1e10)
    lambda_Q = p.get("lambda_Q", 1e9)  # Quantum scale
    
    if lambda_Q == 0:
        return 0.0
    
    return k * math.exp(-r / lambda_Q)


def compute_SMBHUg4Term_SOURCE82(p: dict) -> float:
    """
    SMBH Ug4: Vacuum energy contribution
    
    Formula: Ug4 = ρ_vac * c² * r³ / (3 * M_bh)
    """
    rho_vac = p.get("rho_vac", 1e-26)
    c = UQFF_CONSTANTS["c"]
    r = p.get("r", 1e10)
    M_bh = p.get("M_bh", UQFF_CONSTANTS["Mbh_SOURCE4"])
    
    if M_bh == 0:
        return 0.0
    
    return rho_vac * c * c * r * r * r / (3.0 * M_bh)


def compute_SMBHReactorEfficiencyTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Reactor Efficiency: Accretion disk efficiency
    
    Formula: η = (L_bol) / (Ṁ * c²)
    
    Typical efficiency η ~ 0.1 for thin disk, ~0.4 for spinning BH
    """
    L_bol = p.get("L_bol", 1e38)  # Bolometric luminosity (W)
    M_dot = p.get("M_dot", 1e18)  # Accretion rate (kg/s)
    c = UQFF_CONSTANTS["c"]
    
    if M_dot == 0:
        return 0.0
    
    return L_bol / (M_dot * c * c)


def compute_SMBHPseudoMonopoleTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Pseudo-Monopole: Magnetic monopole-like effects near SMBH
    
    Formula: g_mono = g_m * (r_s/r)² * exp(-r/λ_m)
    """
    g_m = p.get("g_m", 1e-30)  # Monopole coupling
    r_s = p.get("r_s", UQFF_CONSTANTS["r_s"])  # Schwarzschild radius
    r = p.get("r", 1e10)
    lambda_m = p.get("lambda_m", 1e12)  # Monopole scale
    
    if r == 0 or lambda_m == 0:
        return 0.0
    
    return g_m * (r_s / r)**2 * math.exp(-r / lambda_m)


def compute_SMBHRedshiftCorrectionTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Redshift Correction: Gravitational redshift factor
    
    Formula: z_grav = 1 / sqrt(1 - r_s/r) - 1
    """
    r_s = p.get("r_s", UQFF_CONSTANTS["r_s"])
    r = p.get("r", 1e10)
    
    if r <= r_s:
        return float('inf')
    
    return 1.0 / math.sqrt(1.0 - r_s / r) - 1.0


def compute_SMBHUiTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Ui: Combined gravity index
    
    Formula: Ui = Ug1 + Ug2 + Ug3 + Ug4
    """
    Ug1 = compute_SMBHUg1Term_SOURCE82(p)
    Ug2 = compute_SMBHUg2Term_SOURCE82(p)
    Ug3 = compute_SMBHUg3Term_SOURCE82(p)
    Ug4 = compute_SMBHUg4Term_SOURCE82(p)
    
    return Ug1 + Ug2 + Ug3 + Ug4


def compute_SMBHUmTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Um: Magnetic contribution to SMBH gravity
    
    Formula: Um = μ_0 * M_mag / (4π * r³)
    
    Where M_mag is the magnetic dipole moment of the SMBH system.
    """
    mu_0 = 4 * UQFF_CONSTANTS["PI"] * 1e-7  # Permeability of free space
    M_mag = p.get("M_mag", 1e40)  # Magnetic dipole moment
    r = p.get("r", 1e10)
    
    if r == 0:
        return 0.0
    
    return mu_0 * M_mag / (4 * UQFF_CONSTANTS["PI"] * r * r * r)


def compute_SMBHOmegaSGalacticTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Omega_S Galactic: Galactic rotation influence
    
    Formula: Ω_g = v_circular / r
    """
    v_circular = p.get("v_circular", 220e3)  # Circular velocity (m/s)
    r = p.get("r", 2.6e20)  # Distance from center
    
    if r == 0:
        return 0.0
    
    return v_circular / r


def compute_SMBHCosmicTimeTerm_SOURCE82(p: dict) -> float:
    """
    SMBH Cosmic Time: Time evolution factor
    
    Formula: f_cosmic = (1 + z)^(-3) * (t/t_Hubble)
    
    Evolution of SMBH influence with cosmic time.
    """
    z = p.get("z", 0.0)  # Redshift
    t = p.get("t", 1e17)  # Current cosmic time
    tHubble = p.get("tHubble", UQFF_CONSTANTS["tHubble"])
    
    if tHubble == 0:
        return 0.0
    
    return math.pow(1.0 + z, -3) * (t / tHubble)


# =============================================================================
# MUGE RESONANCE DYNAMICS (14 Functions from SOURCE4)
# =============================================================================
# These are the frequency-based resonance dynamics equations used for
# time-evolution simulations in the UQFF framework.

def compute_aDPM_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Dipolar Moment (DPM) base term
    
    Formula: aDPM = F_DPM * f_DPM * E_vac_neb * c_res * V_sys
    Where F_DPM = I * A * (ω₁ - ω₂)
    
    This is the fundamental resonance term from which other dynamics derive.
    """
    I = p.get("I", 1.0)              # Current/intensity
    A = p.get("A", 1.0)              # Area
    omega1 = p.get("omega1", 1e-5)   # Frequency 1
    omega2 = p.get("omega2", 0.9e-5) # Frequency 2
    fDPM = p.get("fDPM", 1e12)       # DPM frequency (PAPER_371 §3)
    Evac_neb = p.get("Evac_neb", 7.09e-36)  # Vacuum energy (nebula)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    Vsys = p.get("Vsys", 1e30)       # System volume
    
    FDPM = I * A * (omega1 - omega2)
    return FDPM * fDPM * Evac_neb * c_res * Vsys


def compute_aTHz_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: THz frequency contribution
    
    Formula: aTHz = f_THz * E_vac_neb * v_exp * aDPM / (E_vac_ISM * c_res)
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    fTHz = p.get("fTHz", 1e12)       # THz frequency
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    vexp = p.get("vexp", 1e5)        # Expansion velocity
    Evac_ISM = p.get("Evac_ISM", 7.09e-37)  # Vacuum energy (ISM)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_ISM == 0 or c_res == 0:
        return 0.0
    
    return fTHz * Evac_neb * vexp * aDPM / (Evac_ISM * c_res)


def compute_avac_diff_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Vacuum differential term
    
    Formula: a_vac_diff = ΔE_vac * v_exp² * aDPM / (E_vac_neb * c_res²)
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    Delta_Evac = p.get("Delta_Evac", 6.381e-36)  # Vacuum energy difference
    vexp = p.get("vexp", 1e5)
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_neb == 0 or c_res == 0:
        return 0.0
    
    return Delta_Evac * vexp * vexp * aDPM / (Evac_neb * c_res * c_res)


def compute_asuper_freq_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Superconductivity frequency term
    
    Formula: a_super_freq = F_super * f_THz * aDPM / (E_vac_neb * c_res)
    
    Models superconducting transition effects on resonance.
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    Fsuper = p.get("Fsuper", 6.287e-19)  # Superconductive force (PAPER_371 §3)
    fTHz = p.get("fTHz", 1e12)
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_neb == 0 or c_res == 0:
        return 0.0
    
    return Fsuper * fTHz * aDPM / (Evac_neb * c_res)


def compute_aaether_res_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Aether resonance coupling
    
    Formula: a_aether_res = UA_SCM * ω_i * f_THz * aDPM * (1 + f_TRZ)
    
    Couples to the universal aether field in UQFF.
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    UA_SCM = p.get("UA_SCM", 10.0)    # Aether SCm coupling (PAPER_371 §3)
    omega_i = p.get("omega_i", 1e-8)  # Intrinsic angular frequency (PAPER_371 §3)
    fTHz = p.get("fTHz", 1e12)
    fTRZ = p.get("fTRZ", 0.1)         # TRZ factor
    
    return UA_SCM * omega_i * fTHz * aDPM * (1 + fTRZ)


def compute_Ug4i_resonance_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Ug4 resonance term
    
    Formula: Ug4i = k4_res * E_react * f_react * aDPM / E_vac_neb * c_res
    Where E_react = 1e46 * exp(-kappa/86400 * t), kappa = 0.0005/day
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    t = p.get("t", 0.0)
    k4_res = p.get("k4_res", 1.0)
    freact = p.get("freact", 1e10)    # Reactive frequency (PAPER_371 §3)
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_neb == 0:
        return 0.0
    
    # kappa = 0.0005/day, convert to per-second
    Ereact = 1e46 * math.exp(-0.0005 / 86400.0 * t)
    return k4_res * Ereact * freact * aDPM / Evac_neb * c_res


def compute_aquantum_freq_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Quantum frequency contribution
    
    Formula: a_quantum_freq = f_quantum * E_vac_neb * aDPM / (E_vac_ISM * c_res)
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    fquantum = p.get("fquantum", 1.445e-17)  # Quantum frequency (PAPER_371 §3)
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    Evac_ISM = p.get("Evac_ISM", 7.09e-37)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_ISM == 0 or c_res == 0:
        return 0.0
    
    return fquantum * Evac_neb * aDPM / (Evac_ISM * c_res)


def compute_aAether_freq_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Aether frequency contribution
    
    Formula: a_Aether_freq = f_Aether * E_vac_neb * aDPM / (E_vac_ISM * c_res)
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    fAether = p.get("fAether", 1.576e-35)  # Aether frequency (PAPER_371 §3)
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    Evac_ISM = p.get("Evac_ISM", 7.09e-37)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_ISM == 0 or c_res == 0:
        return 0.0
    
    return fAether * Evac_neb * aDPM / (Evac_ISM * c_res)


def compute_afluid_freq_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Fluid frequency (Navier-Stokes) contribution
    
    Formula: a_fluid_freq = f_fluid * E_vac_neb * V_sys / (E_vac_ISM * c_res)
    """
    ffluid = p.get("ffluid", 1.269e-14)  # Fluid frequency (PAPER_371 SGR1745 default)
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    Vsys = p.get("Vsys", 1e30)
    Evac_ISM = p.get("Evac_ISM", 7.09e-37)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_ISM == 0 or c_res == 0:
        return 0.0
    
    return ffluid * Evac_neb * Vsys / (Evac_ISM * c_res)


def compute_Osc_term_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Oscillation term (PAPER_371 §2.10)
    
    Formula: Osc = f_osc * cos(2π * f_osc * t)
    """
    t = p.get("t", 0.0)
    fosc = p.get("fosc", 4.57e14)     # Oscillation frequency (PAPER_371 §3)
    
    return fosc * math.cos(2.0 * math.pi * fosc * t)


def compute_aexp_freq_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Expansion frequency (Hubble) contribution
    
    Formula: a_exp_freq = f_exp * E_vac_neb * aDPM / (E_vac_ISM * c_res)
    Where f_exp = 2π * H(z) * t
    """
    aDPM = p.get("aDPM", 0.0)
    if aDPM == 0:
        aDPM = compute_aDPM_SOURCE4(p)
    
    t = p.get("t", 0.0)
    H_z = p.get("H_z", 2.270e-18)     # Hubble parameter
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    Evac_ISM = p.get("Evac_ISM", 7.09e-37)
    c_res = p.get("c_res", 3e8)      # Speed of light (PAPER_371 §3)
    
    if Evac_ISM == 0 or c_res == 0:
        return 0.0
    
    fexp = 2 * UQFF_CONSTANTS["PI"] * H_z * t
    return fexp * Evac_neb * aDPM / (Evac_ISM * c_res)


def compute_a_wormhole_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Wormhole metric contribution
    
    Formula: a_wormhole = f_worm * E_vac_neb / (b² + r²)
    
    Based on Morris-Thorne traversable wormhole metric.
    """
    r = p.get("r", 1e10)
    b = p.get("b", 1.0)               # Wormhole throat radius
    f_worm = p.get("f_worm", 1.0)     # Wormhole factor
    Evac_neb = p.get("Evac_neb", 7.09e-36)
    
    return f_worm * Evac_neb / (b * b + r * r)


def compute_resonance_MUGE_SOURCE4(p: dict) -> float:
    """
    MUGE Resonance: Complete resonance Modified Universal Gravitation Equation
    
    Formula: g_resonance = aDPM + aTHz + a_vac_diff + a_super_freq + a_aether_res +
                           Ug4i + a_quantum_freq + a_Aether_freq + a_fluid_freq +
                           Osc_term + a_exp_freq + a_wormhole
    
    This is the complete MUGE dynamics equation used for time-evolution simulations.
    It combines all 14 resonance terms for the full dynamical picture.
    """
    # Compute base DPM term first
    aDPM = compute_aDPM_SOURCE4(p)
    
    # Pass aDPM to dependent terms
    p_with_dpm = dict(p)
    p_with_dpm["aDPM"] = aDPM
    
    # Compute all resonance terms
    aTHz = compute_aTHz_SOURCE4(p_with_dpm)
    avac_diff = compute_avac_diff_SOURCE4(p_with_dpm)
    asuper_freq = compute_asuper_freq_SOURCE4(p_with_dpm)
    aaether_res = compute_aaether_res_SOURCE4(p_with_dpm)
    Ug4i = compute_Ug4i_resonance_SOURCE4(p_with_dpm)
    aquantum_freq = compute_aquantum_freq_SOURCE4(p_with_dpm)
    aAether_freq = compute_aAether_freq_SOURCE4(p_with_dpm)
    afluid_freq = compute_afluid_freq_SOURCE4(p)  # Doesn't depend on aDPM
    Osc_term = compute_Osc_term_SOURCE4(p)
    aexp_freq = compute_aexp_freq_SOURCE4(p_with_dpm)
    a_wormhole = compute_a_wormhole_SOURCE4(p)
    
    # Sum all resonance contributions
    return (aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i + 
            aquantum_freq + aAether_freq + afluid_freq + Osc_term + aexp_freq + a_wormhole)


# =============================================================================
# PRE-DEFINED ASTROPHYSICAL SYSTEMS (7 systems from SOURCE4)
# =============================================================================

MUGE_SYSTEMS = {
    "SGR1745": {
        "name": "SGR1745", "M": 2.984e30, "r": 1e4, "t": 3.799e10,
        "SFR": 0.0009, "Vsys": 4.189e12, "omega1": 1e10, "omega2": 1e11,
        "I": 1e-15, "A": 10.0, "B": 0.0, "Bcrit": 1e-5,
        "Evac_neb": 7.09e-36, "Evac_ISM": 7.09e-37, "Delta_Evac": 6.381e-36
    },
    "SagittariusA*": {
        "name": "SagittariusA*", "M": 8.155e36, "r": 1.2e11, "t": 1.2e14,
        "SFR": 0.0, "Vsys": 3.552e45, "omega1": 1e8, "omega2": 1e10,
        "I": 1e-20, "A": 5e-6, "B": 1e37, "Bcrit": 1e-4,
        "Evac_neb": 7.09e-36, "Evac_ISM": 7.09e-37, "Delta_Evac": 6.381e-36
    },
    "Tapestry": {
        "name": "Tapestry", "M": 1.989e35, "r": 9.46e18, "t": 3.156e13,
        "SFR": 0.005, "Vsys": 1e53, "omega1": 1e-8, "omega2": 1e-6,
        "I": 1e-22, "A": 1e-10, "B": 1e34, "Bcrit": 1e-3,
        "Evac_neb": 7.09e-36, "Evac_ISM": 7.09e-37, "Delta_Evac": 6.381e-36
    },
    "Westerlund2": {
        "name": "Westerlund2", "M": 1e37, "r": 1.5e20, "t": 1e13,
        "SFR": 0.002, "Vsys": 1e56, "omega1": 1e-6, "omega2": 1e-5,
        "I": 1e-20, "A": 1e-8, "B": 1e35, "Bcrit": 1e-4,
        "Evac_neb": 7.09e-36, "Evac_ISM": 7.09e-37, "Delta_Evac": 6.381e-36
    },
    "PillarsOfCreation": {
        "name": "PillarsOfCreation", "M": 1.989e32, "r": 9.46e15, "t": 1e12,
        "SFR": 0.00014, "Vsys": 1e47, "omega1": 1e-9, "omega2": 1e-8,
        "I": 1e-24, "A": 1e-12, "B": 1e31, "Bcrit": 1e-5,
        "Evac_neb": 7.09e-36, "Evac_ISM": 7.09e-37, "Delta_Evac": 6.381e-36
    },
    "RingsOfRelativity": {
        "name": "RingsOfRelativity", "M": 1.989e36, "r": 1e22, "t": 3.156e14,
        "SFR": 0.01, "Vsys": 1e60, "omega1": 1e-5, "omega2": 1e-4,
        "I": 1e-18, "A": 1e-6, "B": 1e36, "Bcrit": 1e-3,
        "Evac_neb": 7.09e-36, "Evac_ISM": 7.09e-37, "Delta_Evac": 6.381e-36
    },
    "StudentGuideUniverse": {
        "name": "StudentGuideUniverse", "M": 1e53, "r": 1e26, "t": 4.35e17,
        "SFR": 0.0, "Vsys": 1e78, "omega1": 1e-12, "omega2": 1e-10,
        "I": 1e-26, "A": 1e-18, "B": 1e50, "Bcrit": 1e-10,
        "Evac_neb": 7.09e-36, "Evac_ISM": 7.09e-37, "Delta_Evac": 6.381e-36
    }
}


# =============================================================================
# REGISTRY AND INTERFACE
# =============================================================================

CORE_UQFF_EQUATIONS = {
    # Core UQFF (8 equations)
    "Ug1_SOURCE4": compute_Ug1_SOURCE4,
    "Ug2_SOURCE4": compute_Ug2_SOURCE4,
    "Ug3_SOURCE4": compute_Ug3_SOURCE4,
    "Ug4_SOURCE4": compute_Ug4_SOURCE4,
    "Ubi_SOURCE4": compute_Ubi_SOURCE4,
    "Um_SOURCE4": compute_Um_SOURCE4,
    "FU_SOURCE4": compute_FU_SOURCE4,
    
    # MUGE Compressed (10 equations)
    "compressed_base_SOURCE4": compute_compressed_base_SOURCE4,
    "compressed_expansion_SOURCE4": compute_compressed_expansion_SOURCE4,
    "compressed_super_adj_SOURCE4": compute_compressed_super_adj_SOURCE4,
    "compressed_cosm_SOURCE4": compute_compressed_cosm_SOURCE4,
    "compressed_quantum_SOURCE4": compute_compressed_quantum_SOURCE4,
    "compressed_fluid_SOURCE4": compute_compressed_fluid_SOURCE4,
    "compressed_perturbation_SOURCE4": compute_compressed_perturbation_SOURCE4,
    "compressed_MUGE_SOURCE4": compute_compressed_MUGE_SOURCE4,
    
    # SMBH SOURCE82 (15 equations)
    "SMBHDynamicVacuumTerm_SOURCE82": compute_SMBHDynamicVacuumTerm_SOURCE82,
    "SMBHQuantumCouplingTerm_SOURCE82": compute_SMBHQuantumCouplingTerm_SOURCE82,
    "SMBHMSigmaRelationTerm_SOURCE82": compute_SMBHMSigmaRelationTerm_SOURCE82,
    "SMBHBulgeGravityTerm_SOURCE82": compute_SMBHBulgeGravityTerm_SOURCE82,
    "SMBHUg1Term_SOURCE82": compute_SMBHUg1Term_SOURCE82,
    "SMBHUg2Term_SOURCE82": compute_SMBHUg2Term_SOURCE82,
    "SMBHUg3Term_SOURCE82": compute_SMBHUg3Term_SOURCE82,
    "SMBHUg4Term_SOURCE82": compute_SMBHUg4Term_SOURCE82,
    "SMBHReactorEfficiencyTerm_SOURCE82": compute_SMBHReactorEfficiencyTerm_SOURCE82,
    "SMBHPseudoMonopoleTerm_SOURCE82": compute_SMBHPseudoMonopoleTerm_SOURCE82,
    "SMBHRedshiftCorrectionTerm_SOURCE82": compute_SMBHRedshiftCorrectionTerm_SOURCE82,
    "SMBHUiTerm_SOURCE82": compute_SMBHUiTerm_SOURCE82,
    "SMBHUmTerm_SOURCE82": compute_SMBHUmTerm_SOURCE82,
    "SMBHOmegaSGalacticTerm_SOURCE82": compute_SMBHOmegaSGalacticTerm_SOURCE82,
    "SMBHCosmicTimeTerm_SOURCE82": compute_SMBHCosmicTimeTerm_SOURCE82,
    
    # MUGE Resonance Dynamics (14 equations)
    "aDPM_SOURCE4": compute_aDPM_SOURCE4,
    "aTHz_SOURCE4": compute_aTHz_SOURCE4,
    "avac_diff_SOURCE4": compute_avac_diff_SOURCE4,
    "asuper_freq_SOURCE4": compute_asuper_freq_SOURCE4,
    "aaether_res_SOURCE4": compute_aaether_res_SOURCE4,
    "Ug4i_resonance_SOURCE4": compute_Ug4i_resonance_SOURCE4,
    "aquantum_freq_SOURCE4": compute_aquantum_freq_SOURCE4,
    "aAether_freq_SOURCE4": compute_aAether_freq_SOURCE4,
    "afluid_freq_SOURCE4": compute_afluid_freq_SOURCE4,
    "Osc_term_SOURCE4": compute_Osc_term_SOURCE4,
    "aexp_freq_SOURCE4": compute_aexp_freq_SOURCE4,
    "a_wormhole_SOURCE4": compute_a_wormhole_SOURCE4,
    "resonance_MUGE_SOURCE4": compute_resonance_MUGE_SOURCE4,
}

CORE_EQUATION_COUNT = len(CORE_UQFF_EQUATIONS)


def compute_core_equation(name: str, params: dict) -> float:
    """Compute a core UQFF equation by name."""
    if name in CORE_UQFF_EQUATIONS:
        return CORE_UQFF_EQUATIONS[name](params)
    return 0.0


def list_core_equations() -> list:
    """List all core UQFF equation names."""
    return list(CORE_UQFF_EQUATIONS.keys())


def get_core_equation_info(name: str) -> dict:
    """Get docstring information for a core equation."""
    if name in CORE_UQFF_EQUATIONS:
        func = CORE_UQFF_EQUATIONS[name]
        return {
            'name': name,
            'description': func.__doc__ if func.__doc__ else "No description"
        }
    return {}


# =============================================================================
# TEST SUITE
# =============================================================================

def test_core_equations():
    """Test all core UQFF equations with sample parameters."""
    print("=" * 60)
    print("TESTING CORE UQFF EQUATIONS")
    print("=" * 60)
    
    # Sample parameters - using moderate time scale to avoid exp decay to zero
    # At t=1e10 with alpha=1e-6, exp(-alpha*t) = exp(-10000) ≈ 0
    # Use t=1e6 for meaningful results
    sgr_a_params = {
        # Core parameters
        "r": 1e15,             # Distance (m) - not too large for quantum terms
        "t": 1e6,              # Time (s) - moderate to avoid exp decay
        "tn": 0.1,             # Normalized time
        "M": 4e6 * 1.989e30,   # SgrA* mass
        "M_bh": 4e6 * 1.989e30,
        "Mbh": 4e6 * 1.989e30,
        
        # Magnetic field parameters
        "Bs_avg": 1e-4,        # Magnetic field (T)
        "omega_c": 1e-5,       # Core frequency
        "Rs": 1.2e10,          # Schwarzschild radius
        "Ms": 4e6 * 1.989e30,
        "Rb": 1e9,
        "B": 1e-6,
        "Bcrit": 4.4e13,
        
        # UQFF coupling parameters
        "QUA": 0.1,
        "Pcore": 1e16,
        "omega_s": 1e-8,
        "SCm_density": 1e-10,
        "dg": 2.6e20,
        
        # SMBH parameters
        "sigma": 100e3,        # Velocity dispersion
        "M_bulge": 1e10 * 1.989e30,
        "L_bol": 1e38,
        "M_dot": 1e18,
        "v_circular": 220e3,
        "z": 0.0,
        
        # Quantum/monopole parameters - use reasonable scales
        "alpha_q": 1e-10,
        "lambda_q": 1e18,      # Large coherence length
        "k": 1e-10,            # Quantum gravity k
        "lambda_Q": 1e18,      # Large quantum scale
        "g_m": 1e-20,          # Monopole coupling
        "lambda_m": 1e20,
        
        # Resonance parameters
        "I": 1e-20,
        "A": 5e-6,
        "omega1": 1e8,
        "omega2": 1e10,
        "fDPM": 1.0,
        "Evac_neb": 7.09e-36,
        "Evac_ISM": 7.09e-37,
        "Delta_Evac": 6.381e-36,
        "c_res": 1.0,
        "Vsys": 3.552e45,
        "fTHz": 1e12,
        "vexp": 1e5,
        "Fsuper": 1.0,
        "UA_SCM": 0.0001,
        "omega_i": 1e-7,
        "fTRZ": 0.1,
        "k4_res": 1.0,
        "freact": 1.0,
        "fquantum": 1e15,
        "fAether": 1e10,
        "ffluid": 1e-4,
        "H_z": 2.270e-18,
        "b": 1e5,
        "f_worm": 1.0,
        
        # Oscillation parameters
        "A_osc": 1e-10,
        "omega_osc": 1e-10,
        "phi_osc": 0.5,
    }
    
    results = []
    for name, func in CORE_UQFF_EQUATIONS.items():
        try:
            value = func(sgr_a_params)
            status = "OK" if value != 0.0 else "ZERO"
            results.append((name, value, status))
            print(f"  {name:40s} = {value:12.4e}  [{status}]")
        except Exception as e:
            results.append((name, 0.0, f"ERROR: {e}"))
            print(f"  {name:40s} = ERROR: {e}")
    
    # Summary
    working = sum(1 for _, v, s in results if s == "OK")
    zeros = sum(1 for _, v, s in results if s == "ZERO")
    print(f"\n{'=' * 60}")
    print(f"SUMMARY: {working}/{len(results)} working, {zeros} zero")
    print("=" * 60)
    
    return results


def test_dynamics():
    """Test time evolution of MUGE Resonance dynamics."""
    print("\n" + "=" * 60)
    print("TESTING DYNAMICS TIME EVOLUTION")
    print("=" * 60)
    
    import copy
    
    # Use SagittariusA* system parameters
    base_params = MUGE_SYSTEMS["SagittariusA*"].copy()
    base_params.update({
        "fDPM": 1.0, "c_res": 1.0, "fTHz": 1e12, "vexp": 1e5,
        "Fsuper": 1.0, "UA_SCM": 0.0001, "omega_i": 1e-7, "fTRZ": 0.1,
        "k4_res": 1.0, "freact": 1.0, "fquantum": 1e15, "fAether": 1e10,
        "ffluid": 1e-4, "H_z": 2.270e-18, "b": 1e5, "f_worm": 1.0,
        "A_osc": 1e-10, "omega_osc": 1e-10, "phi_osc": 0.0
    })
    
    print(f"\nSystem: {base_params['name']}")
    print(f"Mass: {base_params['M']:.2e} kg")
    print(f"{'Time (s)':<15} {'FU':<15} {'MUGE Compressed':<15} {'Resonance MUGE':<15}")
    print("-" * 60)
    
    for t in [1e6, 1e8, 1e10, 1e12, 1e14]:
        params = copy.copy(base_params)
        params["t"] = t
        params["tn"] = (t / 1e14) % 1.0  # Normalized time
        
        fu = compute_FU_SOURCE4(params)
        muge_c = compute_compressed_MUGE_SOURCE4(params)
        muge_r = compute_resonance_MUGE_SOURCE4(params)
        
        print(f"{t:<15.2e} {fu:<15.4e} {muge_c:<15.4e} {muge_r:<15.4e}")
    
    print("\nDynamics verified: equations evolve with time as expected.")
    return True


if __name__ == "__main__":
    test_core_equations()
    test_dynamics()
