from dpm_helpers import dpm_emergent_ug1, dpm_emergent_ug2
"""
Standard Astrophysics Equations Module (100+ Equations)
══════════════════════════════════════════════════════════════════════════════════

Extracted from Grok Conversation b4469997f5324be48bc0697cdeaf21f9
Complete formula library for astrophysical calculations

Categories:
  01-10: Cosmology & Expansion
  11-20: Stellar Structure & Evolution
  21-30: Supernovae & Remnants
  31-40: Gravitational Waves
  41-50: Black Holes & Compact Objects
  51-60: Galaxy Dynamics
  61-70: Magnetohydrodynamics
  71-80: Particle Acceleration & Radiation
  81-90: Nuclear & Quantum Physics
  91-100: Miscellaneous Astrophysics

Author: Star-Magic Framework
Date: March 1, 2026
"""

import numpy as np
from typing import Dict, Optional, Tuple, Callable
from dataclasses import dataclass


# ═══════════════════════════════════════════════════════════════════════════════
# PHYSICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

C = {
    'G': 6.674e-11,           # Gravitational constant (m³/(kg·s²))
    'c': 2.998e8,             # Speed of light (m/s)
    'hbar': 1.055e-34,        # Reduced Planck constant (J·s)
    'k_B': 1.381e-23,         # Boltzmann constant (J/K)
    'sigma_sb': 5.670e-8,     # Stefan-Boltzmann constant (W/(m²·K⁴))
    'm_p': 1.673e-27,         # Proton mass (kg)
    'm_e': 9.109e-31,         # Electron mass (kg)
    'e': 1.602e-19,           # Elementary charge (C)
    'epsilon_0': 8.854e-12,   # Vacuum permittivity (F/m)
    'mu_0': 4 * np.pi * 1e-7, # Vacuum permeability (H/m)
    'M_sun': 1.989e30,        # Solar mass (kg)
    'R_sun': 6.96e8,          # Solar radius (m)
    'L_sun': 3.828e26,        # Solar luminosity (W)
    'AU': 1.496e11,           # Astronomical unit (m)
    'pc': 3.086e16,           # Parsec (m)
    'ly': 9.461e15,           # Light year (m)
    'yr': 3.156e7,            # Year (s)
    'H_0': 67.4e3/3.086e22,   # Hubble constant (s⁻¹) ~67.4 km/s/Mpc
    'rho_c': 9.47e-27,        # Critical density (kg/m³)
    'Lambda': 1.1e-52,        # Cosmological constant (m⁻²)
}


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 01-10: COSMOLOGY & EXPANSION
# ═══════════════════════════════════════════════════════════════════════════════

def eq_01_friedmann_1(H: float, rho: float, k: float = 0, a: float = 1) -> Dict:
    """
    Eq 01: First Friedmann Equation
    H² = (8πG/3)ρ - kc²/a² + Λc²/3
    """
    G, c, Lambda = C['G'], C['c'], C['Lambda']
    H_sq = (8 * np.pi * G / 3) * rho - k * c**2 / a**2 + Lambda * c**2 / 3
    return {
        'H_squared': H_sq,
        'H': np.sqrt(abs(H_sq)),
        'equation': "H² = (8πG/3)ρ - kc²/a² + Λc²/3",
        'description': "Expansion rate from total energy density"
    }


def eq_02_friedmann_2(rho: float, P: float) -> Dict:
    """
    Eq 02: Second Friedmann Equation (Acceleration)
    ä/a = -(4πG/3)(ρ + 3P/c²) + Λc²/3
    """
    G, c, Lambda = C['G'], C['c'], C['Lambda']
    a_ddot_over_a = -(4 * np.pi * G / 3) * (rho + 3 * P / c**2) + Lambda * c**2 / 3
    return {
        'a_ddot_over_a': a_ddot_over_a,
        'equation': "ä/a = -(4πG/3)(ρ + 3P/c²) + Λc²/3",
        'description': "Cosmic acceleration from matter+pressure+dark energy"
    }


def eq_03_hubble_law(v: float = None, d: float = None, H: float = None) -> Dict:
    """
    Eq 03: Hubble's Law
    v = H₀ × d
    """
    H_0 = C['H_0']
    if H is None:
        H = H_0
    if v is None and d is not None:
        v = H * d
    elif d is None and v is not None:
        d = v / H
    return {
        'v': v,
        'd': d,
        'H': H,
        'equation': "v = H₀ × d",
        'description': "Recession velocity proportional to distance"
    }


def eq_04_cosmological_redshift(z: float, a: float = None) -> Dict:
    """
    Eq 04: Cosmological Redshift
    1 + z = 1/a = λ_obs/λ_emit
    """
    if a is None:
        a = 1 / (1 + z)
    return {
        'z': z,
        'a': a,
        'equation': "1 + z = 1/a",
        'description': "Redshift from scale factor"
    }


def eq_05_luminosity_distance(z: float, Omega_m: float = 0.3, Omega_L: float = 0.7) -> Dict:
    """
    Eq 05: Luminosity Distance (flat universe)
    d_L = (c/H₀)(1+z)∫(dz'/E(z'))
    E(z) = √(Ω_m(1+z)³ + Ω_Λ)
    """
    c, H_0 = C['c'], C['H_0']
    from scipy.integrate import quad
    E = lambda zp: np.sqrt(Omega_m * (1 + zp)**3 + Omega_L)
    integral, _ = quad(lambda zp: 1 / E(zp), 0, z)
    d_L = (c / H_0) * (1 + z) * integral
    return {
        'd_L': d_L,
        'd_L_Mpc': d_L / (3.086e22),
        'equation': "d_L = (c/H₀)(1+z)∫dz'/E(z')",
        'description': "Relationship between flux and intrinsic luminosity"
    }


def eq_06_age_of_universe(H_0: float = None, Omega_m: float = 0.3) -> Dict:
    """
    Eq 06: Age of the Universe (approximation)
    t_0 ≈ (2/3)/H₀ × f(Ω_m, Ω_Λ)
    """
    if H_0 is None:
        H_0 = C['H_0']
    t_H = 1 / H_0  # Hubble time
    t_0 = (2/3) * t_H * 1.0  # Simplified for flat ΛCDM
    return {
        't_0': t_0,
        't_0_Gyr': t_0 / (3.156e16),
        'equation': "t_0 ≈ (2/3)/H₀ × f(Ω)",
        'description': "Current age ~13.8 Gyr"
    }


def eq_07_inflation_slow_roll(V: float, dV_dphi: float, M_pl: float = None) -> Dict:
    """
    Eq 07: Slow-Roll Parameters
    ε = (M_pl²/2)(V'/V)², η = M_pl²(V''/V)
    """
    if M_pl is None:
        M_pl = np.sqrt(C['hbar'] * C['c'] / C['G'])
    epsilon = (M_pl**2 / 2) * (dV_dphi / V)**2
    return {
        'epsilon': epsilon,
        'inflating': epsilon < 1,
        'equation': "ε = (M_pl²/2)(V'/V)²",
        'description': "Inflation requires ε < 1"
    }


def eq_08_matter_radiation_equality(Omega_m: float = 0.3, Omega_r: float = 9e-5) -> Dict:
    """
    Eq 08: Matter-Radiation Equality
    z_eq = Ω_m/Ω_r - 1 ≈ 3400
    """
    z_eq = Omega_m / Omega_r - 1
    return {
        'z_eq': z_eq,
        'equation': "z_eq = Ω_m/Ω_r - 1",
        'description': "Transition from radiation to matter domination"
    }


def eq_09_cmb_temperature(z: float, T_0: float = 2.725) -> Dict:
    """
    Eq 09: CMB Temperature Evolution
    T(z) = T₀(1 + z)
    """
    T_z = T_0 * (1 + z)
    return {
        'T_z': T_z,
        'equation': "T(z) = T₀(1 + z)",
        'description': "CMB photon temperature at redshift z"
    }


def eq_10_baryonic_acoustic_oscillation(z: float = 1100) -> Dict:
    """
    Eq 10: BAO Sound Horizon
    r_s = ∫c_s dt/a ≈ 147 Mpc (comoving)
    """
    r_s = 147 * 3.086e22  # Mpc to m
    theta_s = r_s / ((1 + z) * C['c'] / C['H_0'])  # angular size
    return {
        'r_s': r_s,
        'r_s_Mpc': 147,
        'equation': "r_s = ∫c_s dt/a",
        'description': "Standard ruler for cosmology"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 11-20: STELLAR STRUCTURE & EVOLUTION
# ═══════════════════════════════════════════════════════════════════════════════

def eq_11_hydrostatic_equilibrium(M_r: float, rho: float, r: float) -> Dict:
    """
    Eq 11: Hydrostatic Equilibrium
    dP/dr = -GM_r ρ/r²
    """
    G = C['G']
    dP_dr = -dpm_emergent_ug1(M_r, r) * rho
    return {
        'dP_dr': dP_dr,
        'equation': "dP/dr = -GM_r ρ/r²",
        'description': "Pressure gradient balances gravity"
    }


def eq_12_mass_continuity(rho: float, r: float) -> Dict:
    """
    Eq 12: Mass Continuity
    dM_r/dr = 4πr²ρ
    """
    dM_dr = 4 * np.pi * r**2 * rho
    return {
        'dM_dr': dM_dr,
        'equation': "dM_r/dr = 4πr²ρ",
        'description': "Mass enclosed in shells"
    }


def eq_13_radiative_transport(L_r: float, kappa: float, rho: float, r: float, T: float) -> Dict:
    """
    Eq 13: Radiative Temperature Gradient
    dT/dr = -3κρL_r/(16πacr²T³)
    """
    a = 4 * C['sigma_sb'] / C['c']  # radiation constant
    c = C['c']
    dT_dr = -3 * kappa * rho * L_r / (16 * np.pi * a * c * r**2 * T**3)
    return {
        'dT_dr': dT_dr,
        'equation': "dT/dr = -3κρL_r/(16πacr²T³)",
        'description': "Radiative diffusion"
    }


def eq_14_luminosity_generation(epsilon: float, rho: float, r: float) -> Dict:
    """
    Eq 14: Luminosity Generation
    dL_r/dr = 4πr²ρε
    """
    dL_dr = 4 * np.pi * r**2 * rho * epsilon
    return {
        'dL_dr': dL_dr,
        'equation': "dL_r/dr = 4πr²ρε",
        'description': "Energy generation rate"
    }


def eq_15_eddington_luminosity(M: float) -> Dict:
    """
    Eq 15: Eddington Luminosity
    L_Edd = 4πGMm_p c/σ_T = 1.26×10³¹(M/M_☉) W
    """
    G, c, m_p = C['G'], C['c'], C['m_p']
    sigma_T = 6.652e-29  # Thomson cross-section (m²)
    L_Edd = 4 * np.pi * G * M * m_p * c / sigma_T
    return {
        'L_Edd': L_Edd,
        'L_Edd_solar': L_Edd / C['L_sun'],
        'equation': "L_Edd = 4πGMm_p c/σ_T",
        'description': "Maximum luminosity before radiation pressure ejects matter"
    }


def eq_16_schwarzschild_criterion(nabla_rad: float, nabla_ad: float) -> Dict:
    """
    Eq 16: Convection Criterion
    Convection if: ∇_rad > ∇_ad
    """
    convective = nabla_rad > nabla_ad
    return {
        'convective': convective,
        'nabla_rad': nabla_rad,
        'nabla_ad': nabla_ad,
        'equation': "∇_rad > ∇_ad → convective",
        'description': "Schwarzschild criterion for convection"
    }


def eq_17_main_sequence_lifetime(M: float) -> Dict:
    """
    Eq 17: Main Sequence Lifetime
    τ_MS ≈ 10¹⁰ yr × (M/M_☉)^(-2.5)
    """
    M_sun, yr = C['M_sun'], C['yr']
    tau_MS = 1e10 * yr * (M / M_sun)**(-2.5)
    return {
        'tau_MS': tau_MS,
        'tau_MS_yr': tau_MS / yr,
        'equation': "τ_MS ≈ 10¹⁰(M/M_☉)^(-2.5) yr",
        'description': "Time on main sequence"
    }


def eq_18_mass_luminosity_relation(M: float) -> Dict:
    """
    Eq 18: Mass-Luminosity Relation
    L ∝ M^3.5 (for M ~ 0.5-10 M_☉)
    """
    M_sun, L_sun = C['M_sun'], C['L_sun']
    L = L_sun * (M / M_sun)**3.5
    return {
        'L': L,
        'L_solar': L / L_sun,
        'equation': "L ∝ M^3.5",
        'description': "Higher mass = much higher luminosity"
    }


def eq_19_chandrasekhar_limit() -> Dict:
    """
    Eq 19: Chandrasekhar Mass Limit
    M_Ch = (ℏc/G)^(3/2) / (μ_e m_p)² ≈ 1.4 M_☉
    """
    G, c, hbar, m_p = C['G'], C['c'], C['hbar'], C['m_p']
    mu_e = 2  # mean molecular weight per electron
    M_Ch = (hbar * c / G)**(3/2) / (mu_e * m_p)**2
    return {
        'M_Ch': M_Ch,
        'M_Ch_solar': M_Ch / C['M_sun'],
        'equation': "M_Ch = (ℏc/G)^(3/2)/(μ_e m_p)²",
        'description': "Maximum white dwarf mass"
    }


def eq_20_virial_theorem(K: float = None, W: float = None) -> Dict:
    """
    Eq 20: Virial Theorem
    2K + W = 0 → K = -W/2
    """
    if K is None and W is not None:
        K = -W / 2
    elif W is None and K is not None:
        W = -2 * K
    return {
        'K': K,
        'W': W,
        'equation': "2K + W = 0",
        'description': "Kinetic energy = -(1/2) potential energy in equilibrium"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 21-30: SUPERNOVAE & REMNANTS
# ═══════════════════════════════════════════════════════════════════════════════

def eq_21_sedov_taylor(E: float, rho_0: float, t: float) -> Dict:
    """
    Eq 21: Sedov-Taylor Blast Wave
    R(t) = (Et²/ρ₀)^(1/5) × ξ₀
    """
    xi_0 = 1.15  # dimensionless constant for γ=5/3
    R = xi_0 * (E * t**2 / rho_0)**(1/5)
    v = (2/5) * R / t  # velocity
    return {
        'R': R,
        'v': v,
        'equation': "R(t) = ξ₀(Et²/ρ₀)^(1/5)",
        'description': "SNR expansion in energy-conserving phase"
    }


def eq_22_supernova_energy() -> Dict:
    """
    Eq 22: Supernova Energy Budget
    E_total ≈ 3×10⁴⁶ J (mostly neutrinos)
    E_kinetic ≈ 10⁴⁴ J
    """
    return {
        'E_total': 3e46,
        'E_kinetic': 1e44,
        'E_light': 1e43,
        'equation': "E_total ≈ 3×10⁴⁶ J",
        'description': "Core-collapse SN energy breakdown"
    }


def eq_23_snr_transition_time(E: float, n: float) -> Dict:
    """
    Eq 23: Free Expansion → Sedov Transition
    t_trans ≈ M_ej/(ρ₀ v³) ≈ 200 yr
    """
    yr = C['yr']
    m_p = C['m_p']
    # Estimate for 10 M_☉ ejecta, 10⁴ km/s, n=1 cm⁻³
    M_ej = 10 * C['M_sun']
    v = 1e7  # m/s
    rho_0 = n * m_p
    t_trans = (M_ej / rho_0)**(1/3) / v
    return {
        't_trans': t_trans,
        't_trans_yr': t_trans / yr,
        'equation': "t_trans ≈ M_ej^(1/3)/(ρ₀^(1/3) v)",
        'description': "Transition when swept mass ~ ejecta mass"
    }


def eq_24_pulsar_spindown(P: float, P_dot: float) -> Dict:
    """
    Eq 24: Pulsar Spin-down Age
    τ = P/(2Ṗ)
    """
    tau = P / (2 * P_dot)
    return {
        'tau': tau,
        'tau_yr': tau / C['yr'],
        'equation': "τ = P/(2Ṗ)",
        'description': "Characteristic age assuming n=3"
    }


def eq_25_pulsar_magnetic_field(P: float, P_dot: float) -> Dict:
    """
    Eq 25: Pulsar Surface Magnetic Field
    B = 3.2×10¹⁹ √(P·Ṗ) T
    """
    B = 3.2e19 * np.sqrt(P * P_dot)
    return {
        'B': B,
        'B_G': B * 1e4,  # Gauss
        'equation': "B = 3.2×10¹⁹ √(P·Ṗ) T",
        'description': "Dipole B from spin-down"
    }


def eq_26_neutron_star_radius() -> Dict:
    """
    Eq 26: Neutron Star Radius (TOV)
    R ≈ 10-12 km for M ~ 1.4 M_☉
    """
    return {
        'R_typical': 11e3,  # m
        'M_typical': 1.4 * C['M_sun'],
        'equation': "TOV → R ≈ 10-12 km",
        'description': "From solving Tolman-Oppenheimer-Volkoff equation"
    }


def eq_27_tov_equation(P: float, rho: float, M_r: float, r: float) -> Dict:
    """
    Eq 27: Tolman-Oppenheimer-Volkoff Equation
    dP/dr = -Gρ(M_r + 4πr³P/c²)/[r(r - 2GM_r/c²)]
    """
    G, c = C['G'], C['c']
    numerator = -G * rho * (M_r + 4 * np.pi * r**3 * P / c**2)
    denominator = r * (r - 2 * dpm_emergent_ug1(M_r, c))
    dP_dr = numerator / denominator
    return {
        'dP_dr': dP_dr,
        'equation': "dP/dr = -Gρ(M_r + 4πr³P/c²)/[r(r - 2GM_r/c²)]",
        'description': "Relativistic hydrostatic equilibrium"
    }


def eq_28_glitch_model(delta_nu: float, nu: float) -> Dict:
    """
    Eq 28: Pulsar Glitch
    Δν/ν = I_s ΔΩ_s / (I_c Ω_c)
    """
    fractional = delta_nu / nu
    return {
        'Delta_nu_over_nu': fractional,
        'equation': "Δν/ν = I_s ΔΩ_s/(I_c Ω_c)",
        'description': "Angular momentum transfer from crust superfluid"
    }


def eq_29_magnetar_field() -> Dict:
    """
    Eq 29: Magnetar Critical Field
    B_QED = m_e²c³/(eℏ) ≈ 4.4×10⁹ T
    """
    m_e, c, e, hbar = C['m_e'], C['c'], C['e'], C['hbar']
    B_QED = m_e**2 * c**3 / (e * hbar)
    return {
        'B_QED': B_QED,
        'B_QED_G': B_QED * 1e4,
        'equation': "B_QED = m_e²c³/(eℏ)",
        'description': "Field where Landau levels ~ m_e c²"
    }


def eq_30_braking_index(nu: float, nu_dot: float, nu_ddot: float) -> Dict:
    """
    Eq 30: Pulsar Braking Index
    n = ν ν̈/ν̇²
    """
    n = nu * nu_ddot / nu_dot**2
    return {
        'n': n,
        'equation': "n = ν ν̈/ν̇²",
        'description': "n=3 for pure dipole radiation"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 31-40: GRAVITATIONAL WAVES
# ═══════════════════════════════════════════════════════════════════════════════

def eq_31_gw_chirp_mass(m1: float, m2: float) -> Dict:
    """
    Eq 31: Chirp Mass
    M_c = (m₁m₂)^(3/5)/(m₁+m₂)^(1/5)
    """
    M_c = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    return {
        'M_c': M_c,
        'M_c_solar': M_c / C['M_sun'],
        'equation': "M_c = (m₁m₂)^(3/5)/(m₁+m₂)^(1/5)",
        'description': "Mass combination determining GW amplitude"
    }


def eq_32_gw_strain(M_c: float, f: float, d: float) -> Dict:
    """
    Eq 32: GW Strain Amplitude
    h = (4/d)(GM_c/c²)^(5/3)(πf/c)^(2/3)
    """
    G, c = C['G'], C['c']
    h = (4 / d) * (dpm_emergent_ug1(M_c, c))**(5/3) * (np.pi * f / c)**(2/3)
    return {
        'h': h,
        'equation': "h = (4/d)(GM_c/c²)^(5/3)(πf/c)^(2/3)",
        'description': "Dimensionless strain at detector"
    }


def eq_33_gw_frequency_evolution(M_c: float, f: float) -> Dict:
    """
    Eq 33: GW Frequency Evolution (Inspiral)
    ḟ = (96/5)π^(8/3)(GM_c/c³)^(5/3)f^(11/3)
    """
    G, c = C['G'], C['c']
    f_dot = (96/5) * np.pi**(8/3) * (G * M_c / c**3)**(5/3) * f**(11/3)
    return {
        'f_dot': f_dot,
        'equation': "ḟ = (96/5)π^(8/3)(GM_c/c³)^(5/3)f^(11/3)",
        'description': "Frequency increases as orbit shrinks"
    }


def eq_34_gw_merger_frequency(M_total: float) -> Dict:
    """
    Eq 34: Merger Frequency (ISCO)
    f_merger ≈ c³/(6^(3/2)πGM)
    """
    G, c = C['G'], C['c']
    f_merger = c**3 / (6**(3/2) * np.pi * G * M_total)
    return {
        'f_merger': f_merger,
        'equation': "f_merger ≈ c³/(6^(3/2)πGM)",
        'description': "GW frequency at ISCO"
    }


def eq_35_gw_ringdown_frequency(M: float, a: float = 0) -> Dict:
    """
    Eq 35: Ringdown QNM Frequency
    f_QNM ≈ c³/(2πGM) × [1 - 0.63(1-a)^(3/10)]
    """
    G, c = C['G'], C['c']
    f_QNM = c**3 / (2 * np.pi * G * M) * (1 - 0.63 * (1 - a)**(3/10))
    return {
        'f_QNM': f_QNM,
        'spin': a,
        'equation': "f_QNM ≈ c³/(2πGM) × [1 - 0.63(1-a)^(3/10)]",
        'description': "Quasi-normal mode frequency"
    }


def eq_36_gw_luminosity(M_c: float, f: float) -> Dict:
    """
    Eq 36: GW Luminosity
    L_GW = (32/5)(c⁵/G)(πGM_c f/c³)^(10/3)
    """
    G, c = C['G'], C['c']
    L_GW = (32/5) * (c**5 / G) * (np.pi * G * M_c * f / c**3)**(10/3)
    return {
        'L_GW': L_GW,
        'equation': "L_GW = (32/5)(c⁵/G)(πGM_c f/c³)^(10/3)",
        'description': "Power radiated in gravitational waves"
    }


def eq_37_gw_coalescence_time(M_c: float, f: float) -> Dict:
    """
    Eq 37: Time to Coalescence
    t_c = (5/256)(c⁵/G)(πGM_c/c³)^(-8/3)f^(-8/3)
    """
    G, c = C['G'], C['c']
    t_c = (5/256) * (c**5 / G) * (np.pi * G * M_c / c**3)**(-8/3) * f**(-8/3)
    return {
        't_c': t_c,
        't_c_yr': t_c / C['yr'],
        'equation': "t_c = (5/256)(c⁵/G)(πGM_c/c³)^(-8/3)f^(-8/3)",
        'description': "Time from current frequency to merger"
    }


def eq_38_gw_energy_radiated(M_total: float, M_final: float) -> Dict:
    """
    Eq 38: Energy Radiated in GWs
    E_GW = (M_total - M_final)c²
    """
    c = C['c']
    E_GW = (M_total - M_final) * c**2
    return {
        'E_GW': E_GW,
        'Delta_M': M_total - M_final,
        'efficiency': (M_total - M_final) / M_total,
        'equation': "E_GW = (M_total - M_final)c²",
        'description': "Mass loss converted to GWs"
    }


def eq_39_gw_kick_velocity(eta: float, delta: float) -> Dict:
    """
    Eq 39: BH Merger Kick Velocity
    v_kick ≈ A η²(1 - q)/(1 + q) [1 + B sin(...)], A ~ 5000 km/s
    """
    A = 5000e3  # m/s
    v_kick = A * eta**2 * delta
    return {
        'v_kick': v_kick,
        'v_kick_kmps': v_kick / 1000,
        'equation': "v_kick ≈ A η²(1-q)/(1+q)...",
        'description': "Recoil from asymmetric GW emission"
    }


def eq_40_gw_antenna_pattern(theta: float, phi: float, psi: float) -> Dict:
    """
    Eq 40: Detector Antenna Pattern
    F+ = (1/2)(1 + cos²θ)cos2φ cos2ψ - cosθ sin2φ sin2ψ
    """
    F_plus = (1/2) * (1 + np.cos(theta)**2) * np.cos(2*phi) * np.cos(2*psi) \
             - np.cos(theta) * np.sin(2*phi) * np.sin(2*psi)
    F_cross = (1/2) * (1 + np.cos(theta)**2) * np.cos(2*phi) * np.sin(2*psi) \
              + np.cos(theta) * np.sin(2*phi) * np.cos(2*psi)
    return {
        'F_plus': F_plus,
        'F_cross': F_cross,
        'equation': "F+ = ...(see full expression)",
        'description': "Detector response to polarizations"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 41-50: BLACK HOLES & COMPACT OBJECTS
# ═══════════════════════════════════════════════════════════════════════════════

def eq_41_schwarzschild_radius(M: float) -> Dict:
    """
    Eq 41: Schwarzschild Radius
    r_s = 2GM/c²
    """
    G, c = C['G'], C['c']
    r_s = 2 * dpm_emergent_ug1(M, c)
    return {
        'r_s': r_s,
        'r_s_km': r_s / 1000,
        'equation': "r_s = 2GM/c²",
        'description': "Event horizon radius"
    }


def eq_42_hawking_temperature(M: float) -> Dict:
    """
    Eq 42: Hawking Temperature
    T_H = ℏc³/(8πGM k_B)
    """
    G, c, hbar, k_B = C['G'], C['c'], C['hbar'], C['k_B']
    T_H = hbar * c**3 / (8 * np.pi * G * M * k_B)
    return {
        'T_H': T_H,
        'equation': "T_H = ℏc³/(8πGM k_B)",
        'description': "Black hole thermal radiation"
    }


def eq_43_hawking_luminosity(M: float) -> Dict:
    """
    Eq 43: Hawking Luminosity
    L_H = ℏc⁶/(15360 π G² M²)
    """
    G, c, hbar = C['G'], C['c'], C['hbar']
    L_H = hbar * c**6 / (15360 * np.pi * G**2 * M**2)
    return {
        'L_H': L_H,
        'equation': "L_H = ℏc⁶/(15360π G² M²)",
        'description': "Power radiated via Hawking radiation"
    }


def eq_44_bh_evaporation_time(M: float) -> Dict:
    """
    Eq 44: Black Hole Evaporation Time
    t_evap = 5120 π G² M³ / (ℏc⁴)
    """
    G, c, hbar = C['G'], C['c'], C['hbar']
    t_evap = 5120 * np.pi * G**2 * M**3 / (hbar * c**4)
    return {
        't_evap': t_evap,
        't_evap_yr': t_evap / C['yr'],
        'equation': "t_evap = 5120π G² M³/(ℏc⁴)",
        'description': "Lifetime due to Hawking radiation"
    }


def eq_45_isco_radius(M: float, a: float = 0) -> Dict:
    """
    Eq 45: ISCO Radius
    r_ISCO = 6GM/c² (Schwarzschild), r_ISCO = GM/c² (extreme Kerr prograde)
    """
    G, c = C['G'], C['c']
    r_s = 2 * dpm_emergent_ug1(M, c)
    # Approximation for general spin
    if a == 0:
        r_ISCO = 3 * r_s
    elif a == 1:
        r_ISCO = 0.5 * r_s
    else:
        r_ISCO = r_s * (3 + 2 - a)  # simplified
    return {
        'r_ISCO': r_ISCO,
        'spin': a,
        'equation': "r_ISCO = 6GM/c² (a=0) to GM/c² (a=1)",
        'description': "Innermost stable circular orbit"
    }


def eq_46_bondi_accretion(M: float, rho_inf: float, c_s_inf: float) -> Dict:
    """
    Eq 46: Bondi Accretion Rate
    Ṁ = 4πρ_∞(GM)²/c_s³
    """
    G = C['G']
    M_dot = 4 * np.pi * rho_inf * (G * M)**2 / c_s_inf**3
    return {
        'M_dot': M_dot,
        'M_dot_solar_yr': M_dot / C['M_sun'] * C['yr'],
        'equation': "Ṁ = 4πρ_∞(GM)²/c_s³",
        'description': "Spherical accretion rate"
    }


def eq_47_eddington_accretion(M: float, eta: float = 0.1) -> Dict:
    """
    Eq 47: Eddington Accretion Rate
    Ṁ_Edd = L_Edd/(ηc²) = 4πGMm_p/(ηcσ_T)
    """
    G, c, m_p = C['G'], C['c'], C['m_p']
    sigma_T = 6.652e-29
    M_dot_Edd = 4 * np.pi * G * M * m_p / (eta * c * sigma_T)
    return {
        'M_dot_Edd': M_dot_Edd,
        'M_dot_Edd_solar_yr': M_dot_Edd / C['M_sun'] * C['yr'],
        'equation': "Ṁ_Edd = L_Edd/(ηc²)",
        'description': "Maximum sustainable accretion rate"
    }


def eq_48_blandford_znajek(M: float, a: float, B: float) -> Dict:
    """
    Eq 48: Blandford-Znajek Jet Power
    L_BZ ≈ (1/32) ω_H² B² r_H⁴ c
    ω_H = a c/(2 r_H), r_H = (1 + √(1-a²)) r_s/2
    """
    G, c = C['G'], C['c']
    r_s = 2 * dpm_emergent_ug1(M, c)
    r_H = (1 + np.sqrt(1 - a**2)) * r_s / 2
    omega_H = a * c / (2 * r_H)
    L_BZ = (1/32) * omega_H**2 * B**2 * r_H**4 * c
    return {
        'L_BZ': L_BZ,
        'omega_H': omega_H,
        'r_H': r_H,
        'equation': "L_BZ ≈ (1/32) ω_H² B² r_H⁴ c",
        'description': "Electromagnetic power extraction from rotating BH"
    }


def eq_49_kerr_metric_ergosphere(M: float, a: float, theta: float) -> Dict:
    """
    Eq 49: Ergosphere Radius (Kerr)
    r_erg(θ) = r_s/2 + √((r_s/2)² - a²cos²θ)
    """
    G, c = C['G'], C['c']
    r_s = 2 * dpm_emergent_ug1(M, c)
    r_g = r_s / 2
    r_erg = r_g + np.sqrt(r_g**2 - a**2 * np.cos(theta)**2)
    return {
        'r_erg': r_erg,
        'theta': theta,
        'equation': "r_erg(θ) = r_g + √(r_g² - a²cos²θ)",
        'description': "Boundary where stationary observers impossible"
    }


def eq_50_penrose_process_efficiency(a: float) -> Dict:
    """
    Eq 50: Penrose Process Efficiency
    η_max = (1 - √((1 + √(1-a²))/2))
    """
    eta_max = 1 - np.sqrt((1 + np.sqrt(1 - a**2)) / 2)
    return {
        'eta_max': eta_max,
        'spin': a,
        'equation': "η_max = 1 - √((1 + √(1-a²))/2)",
        'description': "Maximum energy extraction from rotating BH"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 51-60: GALAXY DYNAMICS
# ═══════════════════════════════════════════════════════════════════════════════

def eq_51_rotation_curve(M: float, r: float) -> Dict:
    """
    Eq 51: Circular Rotation Velocity
    v_c = √(GM/r)
    """
    G = C['G']
    v_c = np.sqrt(G * M / r)
    return {
        'v_c': v_c,
        'v_c_kmps': v_c / 1000,
        'equation': "v_c = √(GM/r)",
        'description': "Keplerian rotation"
    }


def eq_52_nfw_profile(rho_s: float, r_s: float, r: float) -> Dict:
    """
    Eq 52: NFW Dark Matter Profile
    ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
    """
    x = r / r_s
    rho = rho_s / (x * (1 + x)**2)
    return {
        'rho': rho,
        'equation': "ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]",
        'description': "Universal dark matter halo profile"
    }


def eq_53_tully_fisher(L: float) -> Dict:
    """
    Eq 53: Tully-Fisher Relation
    L ∝ v_max⁴
    """
    # L = A × v^4, for I-band A ≈ 3e28 for v in km/s, L in L_☉
    v_max = (L / C['L_sun'] / 3e28)**(1/4) * 1000  # m/s
    return {
        'v_max': v_max,
        'v_max_kmps': v_max / 1000,
        'equation': "L ∝ v_max⁴",
        'description': "Spiral galaxy luminosity-velocity relation"
    }


def eq_54_faber_jackson(L: float) -> Dict:
    """
    Eq 54: Faber-Jackson Relation
    L ∝ σ⁴
    """
    # L = B × σ^4 for ellipticals
    sigma = (L / C['L_sun'] / 1e28)**(1/4) * 1000  # m/s
    return {
        'sigma': sigma,
        'sigma_kmps': sigma / 1000,
        'equation': "L ∝ σ⁴",
        'description': "Elliptical galaxy luminosity-dispersion relation"
    }


def eq_55_dynamical_mass(sigma: float, R_e: float) -> Dict:
    """
    Eq 55: Dynamical Mass
    M_dyn = k σ² R_e / G
    """
    G = C['G']
    k = 5  # geometric factor
    M_dyn = k * sigma**2 * R_e / G
    return {
        'M_dyn': M_dyn,
        'M_dyn_solar': M_dyn / C['M_sun'],
        'equation': "M_dyn = k σ² R_e / G",
        'description': "Total mass from velocity dispersion"
    }


def eq_56_jeans_mass(T: float, mu: float, rho: float) -> Dict:
    """
    Eq 56: Jeans Mass
    M_J = (5kT/Gμm_p)^(3/2) × (3/(4πρ))^(1/2)
    """
    G, k_B, m_p = C['G'], C['k_B'], C['m_p']
    M_J = (5 * k_B * T / (G * mu * m_p))**(3/2) * (3 / (4 * np.pi * rho))**(1/2)
    return {
        'M_J': M_J,
        'M_J_solar': M_J / C['M_sun'],
        'equation': "M_J = (5kT/Gμm_p)^(3/2) (3/(4πρ))^(1/2)",
        'description': "Minimum mass for gravitational collapse"
    }


def eq_57_toomre_Q(sigma_r: float, kappa: float, Sigma: float) -> Dict:
    """
    Eq 57: Toomre Stability Parameter
    Q = κ σ_r / (3.36 G Σ)
    Q > 1 stable, Q < 1 unstable to axisymmetric perturbations
    """
    G = C['G']
    Q = kappa * sigma_r / (3.36 * G * Sigma)
    return {
        'Q': Q,
        'stable': Q > 1,
        'equation': "Q = κ σ_r / (3.36 G Σ)",
        'description': "Disk gravitational stability"
    }


def eq_58_dynamical_friction(M: float, v: float, rho: float, ln_Lambda: float = 10) -> Dict:
    """
    Eq 58: Chandrasekhar Dynamical Friction
    F_df = -4π G² M² ρ ln(Λ) / v²
    """
    G = C['G']
    F_df = 4 * np.pi * G**2 * M**2 * rho * ln_Lambda / v**2
    return {
        'F_df': F_df,
        'a_df': F_df / M,
        'equation': "F_df = 4π G² M² ρ ln(Λ) / v²",
        'description': "Drag from gravitational wake"
    }


def eq_59_oort_constants(R: float, v_c: float, dv_dR: float) -> Dict:
    """
    Eq 59: Oort Constants
    A = (1/2)(v_c/R - dv_c/dR)
    B = -(1/2)(v_c/R + dv_c/dR)
    """
    A = (1/2) * (v_c / R - dv_dR)
    B = -(1/2) * (v_c / R + dv_dR)
    return {
        'A': A,
        'B': B,
        'Omega': A - B,
        'equation': "A = (1/2)(v/R - dv/dR)",
        'description': "Local rotation curve parameters"
    }


def eq_60_relaxation_time(N: float, R: float, v: float) -> Dict:
    """
    Eq 60: Two-Body Relaxation Time
    t_relax = N/(8 ln N) × R / v
    """
    t_relax = N / (8 * np.log(N)) * R / v
    return {
        't_relax': t_relax,
        't_relax_yr': t_relax / C['yr'],
        'equation': "t_relax = N/(8 ln N) × R/v",
        'description': "Time for stellar encounters to redistribute energy"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 61-70: MAGNETOHYDRODYNAMICS
# ═══════════════════════════════════════════════════════════════════════════════

def eq_61_alfven_velocity(B: float, rho: float) -> Dict:
    """
    Eq 61: Alfvén Velocity
    v_A = B / √(μ₀ρ)
    """
    mu_0 = C['mu_0']
    v_A = B / np.sqrt(mu_0 * rho)
    return {
        'v_A': v_A,
        'v_A_kmps': v_A / 1000,
        'equation': "v_A = B / √(μ₀ρ)",
        'description': "MHD wave speed along field lines"
    }


def eq_62_plasma_beta(P: float, B: float) -> Dict:
    """
    Eq 62: Plasma Beta
    β = 2μ₀P/B² = P_thermal/P_magnetic
    """
    mu_0 = C['mu_0']
    P_mag = B**2 / (2 * mu_0)
    beta = P / P_mag
    return {
        'beta': beta,
        'P_magnetic': P_mag,
        'equation': "β = 2μ₀P/B²",
        'description': "Ratio of thermal to magnetic pressure"
    }


def eq_63_magnetic_reynolds(v: float, L: float, eta: float) -> Dict:
    """
    Eq 63: Magnetic Reynolds Number
    R_m = vL/η  (η = magnetic diffusivity)
    """
    R_m = v * L / eta
    return {
        'R_m': R_m,
        'frozen_in': R_m >> 1,
        'equation': "R_m = vL/η",
        'description': "R_m >> 1: flux frozen in plasma"
    }


def eq_64_sweet_parker_reconnection(v_A: float, R_m: float) -> Dict:
    """
    Eq 64: Sweet-Parker Reconnection Rate
    v_in/v_A = 1/√R_m
    """
    v_in = v_A / np.sqrt(R_m)
    return {
        'v_in': v_in,
        'v_in_over_vA': 1 / np.sqrt(R_m),
        'equation': "v_in/v_A = 1/√R_m",
        'description': "Inflow speed in steady reconnection"
    }


def eq_65_synchrotron_frequency(B: float, gamma: float) -> Dict:
    """
    Eq 65: Synchrotron Characteristic Frequency
    ν_s = (3/4π)(eB/m_e c)γ²
    """
    e, m_e, c = C['e'], C['m_e'], C['c']
    nu_s = (3/(4 * np.pi)) * (e * B / (m_e * c)) * gamma**2
    return {
        'nu_s': nu_s,
        'equation': "ν_s = (3/4π)(eB/m_e c)γ²",
        'description': "Peak emission frequency for relativistic electron"
    }


def eq_66_synchrotron_cooling(B: float, gamma: float) -> Dict:
    """
    Eq 66: Synchrotron Cooling Time
    t_cool = 6πm_e c/(σ_T B² γ)
    """
    m_e, c = C['m_e'], C['c']
    sigma_T = 6.652e-29
    t_cool = 6 * np.pi * m_e * c / (sigma_T * B**2 * gamma)
    return {
        't_cool': t_cool,
        't_cool_yr': t_cool / C['yr'],
        'equation': "t_cool = 6πm_e c/(σ_T B² γ)",
        'description': "Energy loss time for synchrotron radiation"
    }


def eq_67_faraday_rotation(B_parallel: float, n_e: float, L: float, lambda_: float) -> Dict:
    """
    Eq 67: Faraday Rotation Measure
    RM = 0.81 ∫n_e B_∥ dl  [rad/m²]
    Δψ = RM × λ²
    """
    RM = 0.81 * n_e * B_parallel * L  # simplified
    Delta_psi = RM * lambda_**2
    return {
        'RM': RM,
        'Delta_psi': Delta_psi,
        'equation': "RM = 0.81 ∫n_e B_∥ dl",
        'description': "Polarization rotation through magnetized plasma"
    }


def eq_68_parker_instability(H: float, v_A: float) -> Dict:
    """
    Eq 68: Parker Instability Growth Rate
    γ_PI ≈ v_A / H
    """
    gamma_PI = v_A / H
    return {
        'gamma_PI': gamma_PI,
        't_PI': 1 / gamma_PI,
        'equation': "γ_PI ≈ v_A / H",
        'description': "Magnetic buoyancy instability in disk"
    }


def eq_69_ideal_mhd_induction(B: float, v: float, L: float) -> Dict:
    """
    Eq 69: Ideal MHD Induction Equation
    ∂B/∂t = ∇×(v×B) → B grows at rate v/L
    """
    growth_rate = v / L
    e_folding_time = L / v
    return {
        'growth_rate': growth_rate,
        'e_folding': e_folding_time,
        'equation': "∂B/∂t = ∇×(v×B)",
        'description': "Field line stretching amplifies B"
    }


def eq_70_magnetic_pressure_gradient(B: float, L: float) -> Dict:
    """
    Eq 70: Magnetic Pressure Gradient Force
    f_B = -∇(B²/2μ₀) + (B·∇)B/μ₀
    """
    mu_0 = C['mu_0']
    P_B = B**2 / (2 * mu_0)
    f_B = P_B / L  # gradient estimate
    return {
        'P_B': P_B,
        'f_B': f_B,
        'equation': "f_B = -∇(B²/2μ₀) + (B·∇)B/μ₀",
        'description': "Magnetic pressure + tension forces"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 71-80: PARTICLE ACCELERATION & RADIATION
# ═══════════════════════════════════════════════════════════════════════════════

def eq_71_dsa_momentum_gain(u_1: float, u_2: float, p: float) -> Dict:
    """
    Eq 71: Diffusive Shock Acceleration (DSA)
    dp/dt = (4/3)(u_1 - u_2)/r × p  where r = u_1/u_2
    """
    r = u_1 / u_2
    dp_dt = (4/3) * (u_1 - u_2) / r * p
    return {
        'dp_dt': dp_dt,
        'compression': r,
        'equation': "dp/dt = (4/3)(u_1 - u_2)p/r",
        'description': "First-order Fermi acceleration"
    }


def eq_72_dsa_spectral_index(r: float) -> Dict:
    """
    Eq 72: DSA Power Law Index
    q = 3r/(r-1) → f(p) ∝ p^(-q), N(E) ∝ E^(-s), s = (r+2)/(2(r-1))
    """
    q = 3 * r / (r - 1)
    s = (r + 2) / (2 * (r - 1))
    return {
        'q': q,
        's': s,
        'equation': "q = 3r/(r-1)",
        'description': "q=4, s=2 for strong shock r=4"
    }


def eq_73_cosmic_ray_energy_density() -> Dict:
    """
    Eq 73: Galactic Cosmic Ray Energy Density
    u_CR ≈ 1 eV/cm³ ≈ 1.6×10⁻¹³ J/m³
    """
    return {
        'u_CR': 1.6e-13,
        'u_CR_eV_cm3': 1.0,
        'equation': "u_CR ≈ 1 eV/cm³",
        'description': "Local cosmic ray energy density"
    }


def eq_74_hillas_criterion(E: float, Z: int, B: float, R: float) -> Dict:
    """
    Eq 74: Hillas Criterion for Maximum CR Energy
    E_max = ZeBR (gyroradius ≤ source size)
    """
    e = C['e']
    E_max = Z * e * B * R
    return {
        'E_max': E_max,
        'E_max_eV': E_max / (1.6e-19),
        'equation': "E_max = ZeBR",
        'description': "Confinement limit for acceleration"
    }


def eq_75_inverse_compton(E_photon: float, gamma: float) -> Dict:
    """
    Eq 75: Inverse Compton Scattered Photon Energy
    E_IC ≈ (4/3)γ²E_photon (Thomson limit)
    """
    E_IC = (4/3) * gamma**2 * E_photon
    return {
        'E_IC': E_IC,
        'boost_factor': (4/3) * gamma**2,
        'equation': "E_IC ≈ (4/3)γ²E_photon",
        'description': "Upscattered photon energy"
    }


def eq_76_bremsstrahlung_spectrum(T: float) -> Dict:
    """
    Eq 76: Thermal Bremsstrahlung Spectrum
    ε_ff ∝ n_e² T^(1/2) exp(-hν/kT)
    """
    k_B = C['k_B']
    E_cutoff = k_B * T
    return {
        'E_cutoff': E_cutoff,
        'E_cutoff_keV': E_cutoff / (1.6e-16),
        'equation': "ε_ff ∝ n_e² T^(1/2) exp(-hν/kT)",
        'description': "Free-free emission spectrum"
    }


def eq_77_photon_photon_pair_production(E1: float, E2: float, theta: float) -> Dict:
    """
    Eq 77: Photon-Photon Pair Production Threshold
    E₁E₂(1 - cosθ) ≥ 2(m_e c²)²
    """
    m_e, c = C['m_e'], C['c']
    threshold = 2 * (m_e * c**2)**2
    product = E1 * E2 * (1 - np.cos(theta))
    return {
        'product': product,
        'threshold': threshold,
        'pairs_possible': product >= threshold,
        'equation': "E₁E₂(1-cosθ) ≥ 2(m_e c²)²",
        'description': "γγ → e⁺e⁻ threshold"
    }


def eq_78_compton_wavelength() -> Dict:
    """
    Eq 78: Compton Wavelength
    λ_C = h/(m_e c) = 2.43×10⁻¹² m
    """
    hbar, m_e, c = C['hbar'], C['m_e'], C['c']
    lambda_C = 2 * np.pi * hbar / (m_e * c)
    return {
        'lambda_C': lambda_C,
        'equation': "λ_C = h/(m_e c)",
        'description': "Quantum scale for electron"
    }


def eq_79_gyrofrequency(B: float, m: float, q: float) -> Dict:
    """
    Eq 79: Cyclotron/Gyrofrequency
    ω_c = qB/m
    """
    omega_c = q * B / m
    nu_c = omega_c / (2 * np.pi)
    return {
        'omega_c': omega_c,
        'nu_c': nu_c,
        'equation': "ω_c = qB/m",
        'description': "Particle gyration frequency"
    }


def eq_80_gyroradius(v_perp: float, B: float, m: float, q: float) -> Dict:
    """
    Eq 80: Larmor/Gyroradius
    r_L = mv_⊥/(qB)
    """
    r_L = m * v_perp / (q * B)
    return {
        'r_L': r_L,
        'equation': "r_L = mv_⊥/(qB)",
        'description': "Particle gyration radius"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 81-90: NUCLEAR & QUANTUM PHYSICS
# ═══════════════════════════════════════════════════════════════════════════════

def eq_81_pp_chain_rate(T: float) -> Dict:
    """
    Eq 81: p-p Chain Reaction Rate
    ε_pp ∝ ρ X² T⁴ (for T ~ 10⁷ K)
    """
    epsilon_pp_prefactor = 1.08e-12  # approximate
    return {
        'T_exponent': 4,
        'equation': "ε_pp ∝ ρ X² T⁴",
        'description': "Hydrogen burning in Sun"
    }


def eq_82_cno_cycle_rate(T: float) -> Dict:
    """
    Eq 82: CNO Cycle Reaction Rate
    ε_CNO ∝ ρ X X_CNO T^16 (for T ~ 2×10⁷ K)
    """
    return {
        'T_exponent': 16,
        'equation': "ε_CNO ∝ ρ X X_CNO T^16",
        'description': "Hydrogen burning via CNO catalysts"
    }


def eq_83_triple_alpha_rate(T: float) -> Dict:
    """
    Eq 83: Triple-Alpha Reaction Rate
    ε_3α ∝ ρ² Y³ T^40 (highly temperature sensitive)
    """
    return {
        'T_exponent': 40,
        'equation': "ε_3α ∝ ρ² Y³ T^40",
        'description': "Helium burning to carbon"
    }


def eq_84_nuclear_binding_energy(A: int, Z: int) -> Dict:
    """
    Eq 84: Bethe-Weizsäcker Mass Formula
    B(A,Z) = a_V A - a_S A^(2/3) - a_C Z²/A^(1/3) - a_A(A-2Z)²/A + δ
    """
    a_V, a_S, a_C, a_A = 15.8, 18.3, 0.714, 23.2  # MeV
    # Pairing term
    if A % 2 == 1:
        delta = 0
    elif Z % 2 == 0:
        delta = 12 / A**(1/2)
    else:
        delta = -12 / A**(1/2)
    
    B = (a_V * A - a_S * A**(2/3) - a_C * Z**2 / A**(1/3) 
         - a_A * (A - 2*Z)**2 / A + delta)
    return {
        'B': B,
        'B_per_nucleon': B / A,
        'equation': "B(A,Z) = a_V A - a_S A^(2/3) - ...",
        'description': "Nuclear binding energy (MeV)"
    }


def eq_85_gamow_peak(T: float, Z1: int, Z2: int, A: float) -> Dict:
    """
    Eq 85: Gamow Peak Energy
    E_0 = (b k_B T / 2)^(2/3) where b = √(2μ) π e² Z₁Z₂/ℏ
    """
    k_B = C['k_B']
    E_0 = k_B * T * (Z1 * Z2)**(2/3)  # simplified
    return {
        'E_0': E_0,
        'E_0_keV': E_0 / (1.6e-16),
        'equation': "E_0 = (b k_B T/2)^(2/3)",
        'description': "Most probable tunneling energy"
    }


def eq_86_neutrino_luminosity_collapse() -> Dict:
    """
    Eq 86: Core-Collapse Neutrino Luminosity
    L_ν ≈ 3×10⁴⁶ W for ~10 seconds
    """
    return {
        'L_nu': 3e46,
        't_nu': 10,
        'E_total_nu': 3e53,
        'equation': "L_ν ≈ 3×10⁴⁶ W",
        'description': "Neutrino energy dominates SN output"
    }


def eq_87_electron_degeneracy_pressure(rho: float, mu_e: float = 2) -> Dict:
    """
    Eq 87: Electron Degeneracy Pressure
    P = K ρ^(5/3) (non-relativistic)
    """
    m_e, hbar = C['m_e'], C['hbar']
    K = (hbar**2 / (15 * np.pi**2 * m_e)) * (3 * np.pi**2 / (mu_e * C['m_p']))**(5/3)
    P = K * rho**(5/3)
    return {
        'P': P,
        'K': K,
        'equation': "P = K ρ^(5/3)",
        'description': "White dwarf support (NR limit)"
    }


def eq_88_neutron_degeneracy_pressure(rho: float) -> Dict:
    """
    Eq 88: Neutron Degeneracy Pressure
    P_n = K_n ρ^(5/3) (non-relativistic)
    """
    m_n, hbar = C['m_p'], C['hbar']  # m_n ≈ m_p
    K_n = (hbar**2 / (15 * np.pi**2 * m_n)) * (3 * np.pi**2 / m_n)**(5/3)
    P_n = K_n * rho**(5/3)
    return {
        'P_n': P_n,
        'equation': "P_n = K_n ρ^(5/3)",
        'description': "Neutron star support (NR limit)"
    }


def eq_89_quark_deconfinement_density() -> Dict:
    """
    Eq 89: Quark Deconfinement Density
    ρ_QCD ≈ 5-10 ρ_nuclear ≈ 1-2×10¹⁸ kg/m³
    """
    rho_nuclear = 2.3e17  # kg/m³
    return {
        'rho_nuclear': rho_nuclear,
        'rho_QCD_min': 5 * rho_nuclear,
        'rho_QCD_max': 10 * rho_nuclear,
        'equation': "ρ_QCD ≈ 5-10 ρ_nuclear",
        'description': "Density for quark-gluon phase"
    }


def eq_90_superconducting_gap(T_c: float) -> Dict:
    """
    Eq 90: BCS Superconducting Gap
    Δ(0) = 1.76 k_B T_c
    """
    k_B = C['k_B']
    Delta_0 = 1.76 * k_B * T_c
    return {
        'Delta_0': Delta_0,
        'Delta_0_eV': Delta_0 / (1.6e-19),
        'equation': "Δ(0) = 1.76 k_B T_c",
        'description': "Zero-temperature energy gap"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY 91-100: MISCELLANEOUS ASTROPHYSICS
# ═══════════════════════════════════════════════════════════════════════════════

def eq_91_einstein_ring(M: float, D_L: float, D_S: float) -> Dict:
    """
    Eq 91: Einstein Ring Radius
    θ_E = √(4GM D_LS / (c² D_L D_S))
    """
    G, c = C['G'], C['c']
    D_LS = D_S - D_L
    theta_E = np.sqrt(4 * G * M * D_LS / (c**2 * D_L * D_S))
    return {
        'theta_E': theta_E,
        'theta_E_arcsec': theta_E * 206265,
        'equation': "θ_E = √(4GM D_LS/(c² D_L D_S))",
        'description': "Gravitational lensing ring"
    }


def eq_92_magnification_point_lens(u: float) -> Dict:
    """
    Eq 92: Point Lens Magnification
    μ = (u² + 2)/(u√(u² + 4))
    """
    mu = (u**2 + 2) / (u * np.sqrt(u**2 + 4))
    return {
        'mu': mu,
        'u': u,
        'equation': "μ = (u² + 2)/(u√(u² + 4))",
        'description': "Microlensing magnification"
    }


def eq_93_shakura_sunyaev(alpha: float, H: float, c_s: float) -> Dict:
    """
    Eq 93: Shakura-Sunyaev α-Viscosity
    ν = α c_s H
    """
    nu = alpha * c_s * H
    return {
        'nu': nu,
        'alpha': alpha,
        'equation': "ν = α c_s H",
        'description': "Accretion disk viscosity prescription"
    }


def eq_94_disk_temperature(M: float, M_dot: float, r: float) -> Dict:
    """
    Eq 94: Thin Disk Temperature
    T(r) = [3GMṀ/(8πσr³)]^(1/4) × f(r)
    """
    G, sigma_sb = C['G'], C['sigma_sb']
    T = (3 * G * M * M_dot / (8 * np.pi * sigma_sb * r**3))**(1/4)
    return {
        'T': T,
        'equation': "T(r) ∝ (GMṀ/r³)^(1/4)",
        'description': "Accretion disk effective temperature"
    }


def eq_95_salpeter_imf(M: float) -> Dict:
    """
    Eq 95: Salpeter Initial Mass Function
    ξ(M) ∝ M^(-2.35)
    """
    xi = M**(-2.35)
    return {
        'xi': xi,
        'slope': -2.35,
        'equation': "ξ(M) ∝ M^(-2.35)",
        'description': "Stellar initial mass distribution"
    }


def eq_96_kennicutt_schmidt(Sigma_gas: float) -> Dict:
    """
    Eq 96: Kennicutt-Schmidt Law
    Σ_SFR = A × Σ_gas^N, N ≈ 1.4
    """
    A = 2.5e-4  # (M_☉/yr/kpc²) / (M_☉/pc²)^1.4
    N = 1.4
    Sigma_SFR = A * Sigma_gas**N
    return {
        'Sigma_SFR': Sigma_SFR,
        'N': N,
        'equation': "Σ_SFR = A × Σ_gas^1.4",
        'description': "Star formation rate vs gas surface density"
    }


def eq_97_sn_ia_peak(M_Ni: float = 0.6) -> Dict:
    """
    Eq 97: Type Ia SN Peak Luminosity (Arnett's Rule)
    L_peak ≈ L_Ni(t_rise) ≈ M_Ni × (6.45 e^(-t_r/8.8) + 1.45 e^(-t_r/111)) × 10⁴³ erg/s
    """
    # At t_rise ~ 17-20 days
    L_peak = M_Ni * 2e43  # erg/s
    return {
        'L_peak': L_peak,
        'L_peak_W': L_peak * 1e-7,
        'equation': "L_peak ≈ M_Ni × ε_Ni(t_rise)",
        'description': "Standard candle calibration"
    }


def eq_98_stroemgren_sphere(S_star: float, n: float, alpha_B: float = 2.6e-13) -> Dict:
    """
    Eq 98: Strömgren Radius
    R_S = (3 S*/(4πn² α_B))^(1/3)
    """
    R_S = (3 * S_star / (4 * np.pi * n**2 * alpha_B))**(1/3)
    return {
        'R_S': R_S,
        'R_S_pc': R_S / C['pc'],
        'equation': "R_S = (3 S*/(4πn² α_B))^(1/3)",
        'description': "HII region radius"
    }


def eq_99_cooling_function(T: float) -> Dict:
    """
    Eq 99: Astrophysical Cooling Function (approximate)
    Λ(T) ~ 10^(-22) (T/10⁶)^(-0.7) for 10⁵ < T < 10⁷ K
    """
    Lambda = 1e-22 * (T / 1e6)**(-0.7)
    return {
        'Lambda': Lambda,
        'unit': 'erg cm³/s',
        'equation': "Λ(T) ~ 10^(-22)(T/10⁶)^(-0.7)",
        'description': "Radiative cooling rate"
    }


def eq_100_ram_pressure(rho_ICM: float, v: float) -> Dict:
    """
    Eq 100: Ram Pressure Stripping
    P_ram = ρ_ICM × v²
    Stripping: P_ram > 2πGΣ_gas Σ_stars
    """
    P_ram = rho_ICM * v**2
    return {
        'P_ram': P_ram,
        'equation': "P_ram = ρ_ICM × v²",
        'description': "ICM pressure on orbiting galaxy"
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MASTER EQUATION REGISTRY
# ═══════════════════════════════════════════════════════════════════════════════

EQUATION_REGISTRY = {
    # Cosmology
    'friedmann_1': eq_01_friedmann_1,
    'friedmann_2': eq_02_friedmann_2,
    'hubble_law': eq_03_hubble_law,
    'cosmological_redshift': eq_04_cosmological_redshift,
    'luminosity_distance': eq_05_luminosity_distance,
    'age_of_universe': eq_06_age_of_universe,
    'inflation_slow_roll': eq_07_inflation_slow_roll,
    'matter_radiation_eq': eq_08_matter_radiation_equality,
    'cmb_temperature': eq_09_cmb_temperature,
    'bao_sound_horizon': eq_10_baryonic_acoustic_oscillation,
    
    # Stellar Structure
    'hydrostatic_eq': eq_11_hydrostatic_equilibrium,
    'mass_continuity': eq_12_mass_continuity,
    'radiative_transport': eq_13_radiative_transport,
    'luminosity_gen': eq_14_luminosity_generation,
    'eddington_lum': eq_15_eddington_luminosity,
    'convection_criterion': eq_16_schwarzschild_criterion,
    'ms_lifetime': eq_17_main_sequence_lifetime,
    'mass_luminosity': eq_18_mass_luminosity_relation,
    'chandrasekhar': eq_19_chandrasekhar_limit,
    'virial_theorem': eq_20_virial_theorem,
    
    # SNe & Remnants
    'sedov_taylor': eq_21_sedov_taylor,
    'sn_energy': eq_22_supernova_energy,
    'snr_transition': eq_23_snr_transition_time,
    'pulsar_spindown': eq_24_pulsar_spindown,
    'pulsar_bfield': eq_25_pulsar_magnetic_field,
    'ns_radius': eq_26_neutron_star_radius,
    'tov': eq_27_tov_equation,
    'glitch': eq_28_glitch_model,
    'magnetar_field': eq_29_magnetar_field,
    'braking_index': eq_30_braking_index,
    
    # Gravitational Waves
    'chirp_mass': eq_31_gw_chirp_mass,
    'gw_strain': eq_32_gw_strain,
    'gw_freq_evolution': eq_33_gw_frequency_evolution,
    'gw_merger_freq': eq_34_gw_merger_frequency,
    'gw_ringdown': eq_35_gw_ringdown_frequency,
    'gw_luminosity': eq_36_gw_luminosity,
    'gw_coalescence_time': eq_37_gw_coalescence_time,
    'gw_energy': eq_38_gw_energy_radiated,
    'gw_kick': eq_39_gw_kick_velocity,
    'gw_antenna': eq_40_gw_antenna_pattern,
    
    # Black Holes
    'schwarzschild': eq_41_schwarzschild_radius,
    'hawking_temp': eq_42_hawking_temperature,
    'hawking_lum': eq_43_hawking_luminosity,
    'bh_evaporation': eq_44_bh_evaporation_time,
    'isco': eq_45_isco_radius,
    'bondi_accretion': eq_46_bondi_accretion,
    'eddington_accretion': eq_47_eddington_accretion,
    'blandford_znajek': eq_48_blandford_znajek,
    'ergosphere': eq_49_kerr_metric_ergosphere,
    'penrose_efficiency': eq_50_penrose_process_efficiency,
    
    # Galaxy Dynamics
    'rotation_curve': eq_51_rotation_curve,
    'nfw_profile': eq_52_nfw_profile,
    'tully_fisher': eq_53_tully_fisher,
    'faber_jackson': eq_54_faber_jackson,
    'dynamical_mass': eq_55_dynamical_mass,
    'jeans_mass': eq_56_jeans_mass,
    'toomre_q': eq_57_toomre_Q,
    'dynamical_friction': eq_58_dynamical_friction,
    'oort_constants': eq_59_oort_constants,
    'relaxation_time': eq_60_relaxation_time,
    
    # MHD
    'alfven_velocity': eq_61_alfven_velocity,
    'plasma_beta': eq_62_plasma_beta,
    'magnetic_reynolds': eq_63_magnetic_reynolds,
    'sweet_parker': eq_64_sweet_parker_reconnection,
    'synchrotron_freq': eq_65_synchrotron_frequency,
    'synchrotron_cooling': eq_66_synchrotron_cooling,
    'faraday_rotation': eq_67_faraday_rotation,
    'parker_instability': eq_68_parker_instability,
    'mhd_induction': eq_69_ideal_mhd_induction,
    'magnetic_pressure': eq_70_magnetic_pressure_gradient,
    
    # Particle Acceleration
    'dsa_momentum': eq_71_dsa_momentum_gain,
    'dsa_spectral_index': eq_72_dsa_spectral_index,
    'cr_energy_density': eq_73_cosmic_ray_energy_density,
    'hillas': eq_74_hillas_criterion,
    'inverse_compton': eq_75_inverse_compton,
    'bremsstrahlung': eq_76_bremsstrahlung_spectrum,
    'pair_production': eq_77_photon_photon_pair_production,
    'compton_wavelength': eq_78_compton_wavelength,
    'gyrofrequency': eq_79_gyrofrequency,
    'gyroradius': eq_80_gyroradius,
    
    # Nuclear
    'pp_chain': eq_81_pp_chain_rate,
    'cno_cycle': eq_82_cno_cycle_rate,
    'triple_alpha': eq_83_triple_alpha_rate,
    'binding_energy': eq_84_nuclear_binding_energy,
    'gamow_peak': eq_85_gamow_peak,
    'neutrino_collapse': eq_86_neutrino_luminosity_collapse,
    'electron_degeneracy': eq_87_electron_degeneracy_pressure,
    'neutron_degeneracy': eq_88_neutron_degeneracy_pressure,
    'quark_deconf': eq_89_quark_deconfinement_density,
    'bcs_gap': eq_90_superconducting_gap,
    
    # Miscellaneous
    'einstein_ring': eq_91_einstein_ring,
    'point_lens_mag': eq_92_magnification_point_lens,
    'shakura_sunyaev': eq_93_shakura_sunyaev,
    'disk_temperature': eq_94_disk_temperature,
    'salpeter_imf': eq_95_salpeter_imf,
    'kennicutt_schmidt': eq_96_kennicutt_schmidt,
    'sn_ia_peak': eq_97_sn_ia_peak,
    'stroemgren': eq_98_stroemgren_sphere,
    'cooling_function': eq_99_cooling_function,
    'ram_pressure': eq_100_ram_pressure,
}


def list_all_equations() -> str:
    """List all 100 equations with descriptions"""
    lines = ["=" * 80]
    lines.append("STANDARD ASTROPHYSICS EQUATIONS (100 Equations)")
    lines.append("=" * 80)
    
    categories = [
        (1, 10, "COSMOLOGY & EXPANSION"),
        (11, 20, "STELLAR STRUCTURE & EVOLUTION"),
        (21, 30, "SUPERNOVAE & REMNANTS"),
        (31, 40, "GRAVITATIONAL WAVES"),
        (41, 50, "BLACK HOLES & COMPACT OBJECTS"),
        (51, 60, "GALAXY DYNAMICS"),
        (61, 70, "MAGNETOHYDRODYNAMICS"),
        (71, 80, "PARTICLE ACCELERATION & RADIATION"),
        (81, 90, "NUCLEAR & QUANTUM PHYSICS"),
        (91, 100, "MISCELLANEOUS ASTROPHYSICS"),
    ]
    
    for start, end, name in categories:
        lines.append(f"\n--- {name} ---")
        for i in range(start, end + 1):
            eq_name = list(EQUATION_REGISTRY.keys())[i - 1]
            lines.append(f"  {i:02d}. {eq_name}")
    
    return "\n".join(lines)


__all__ = [
    'C',
    'EQUATION_REGISTRY',
    'list_all_equations',
] + list(EQUATION_REGISTRY.keys())


if __name__ == "__main__":
    print(list_all_equations())
    
    # Demo some equations
    print("\n" + "=" * 80)
    print("DEMONSTRATION")
    print("=" * 80)
    
    # Chirp mass for GW150914-like event
    result = eq_31_gw_chirp_mass(36 * C['M_sun'], 29 * C['M_sun'])
    print(f"\nChirp Mass (GW150914-like): {result['M_c_solar']:.1f} M_☉")
    
    # Sedov-Taylor for Crab-like SNR
    result = eq_21_sedov_taylor(1e44, 1e-21, 1000 * C['yr'])
    print(f"SNR Radius after 1000 yr: {result['R']/C['pc']:.2f} pc")
    
    # Schwarzschild radius of Sgr A*
    result = eq_41_schwarzschild_radius(4e6 * C['M_sun'])
    print(f"Sgr A* Schwarzschild radius: {result['r_s']/C['AU']:.1f} AU")
