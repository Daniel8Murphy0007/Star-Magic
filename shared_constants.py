#!/usr/bin/env python3
"""
Unified Physical Constants for UQFF Star-Magic Framework
=========================================================

This module provides ALL physical constants used across:
- CondensedPhysics.py (Python UQFF engine)
- MAIN_1_CoAnQi.cpp (via shared_constants.h)
- source2.cpp (Qt6 GUI)
- index.js (JavaScript engine)

IMPORTANT: Any constant changes must be synchronized with:
  - shared_constants.h (C++)
  - index.js CONSTANTS object

Usage:
    from shared_constants import CONSTANTS
    
    # Access via dictionary
    c = CONSTANTS['c']  # 2.99792458e8
    
    # Access via Constants class
    from shared_constants import Constants
    c = Constants.c

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

from dataclasses import dataclass
from typing import Dict, Any
import math


@dataclass(frozen=True)
class Constants:
    """
    Physical constants as a frozen dataclass.
    All values synchronized with shared_constants.h.
    """
    
    # ═══════════════════════════════════════════════════════════════════════════
    # FUNDAMENTAL PHYSICAL CONSTANTS (CODATA 2018)
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Speed of light in vacuum (m/s) - exact by definition
    c: float = 2.99792458e8
    
    # Gravitational constant (m³ kg⁻¹ s⁻²)
    G: float = 6.67430e-11
    
    # Reduced Planck constant ℏ (J·s)
    hbar: float = 1.054571817e-34
    
    # Planck constant h (J·s)
    h_planck: float = 6.62607015e-34
    
    # Elementary charge (C)
    e: float = 1.602176634e-19
    
    # Boltzmann constant (J/K)
    k_B: float = 1.380649e-23
    
    # Vacuum permittivity ε₀ (F/m)
    epsilon_0: float = 8.8541878128e-12
    
    # Vacuum permeability μ₀ (H/m)
    mu_0: float = 1.25663706212e-6
    
    # Fine structure constant α ≈ 1/137
    alpha: float = 7.2973525693e-3
    
    # Avogadro constant (mol⁻¹)
    N_A: float = 6.02214076e23
    
    # Stefan-Boltzmann constant (W m⁻² K⁻⁴)
    sigma_SB: float = 5.670374419e-8
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # PARTICLE MASSES
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Electron mass (kg)
    m_e: float = 9.1093837015e-31
    
    # Proton mass (kg)
    m_p: float = 1.67262192369e-27
    
    # Neutron mass (kg)
    m_n: float = 1.67492749804e-27
    
    # Atomic mass unit (kg)
    u: float = 1.66053906660e-27
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # ASTROPHYSICAL CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Solar mass M☉ (kg)
    M_sun: float = 1.98892e30
    
    # Solar radius R☉ (m)
    R_sun: float = 6.9634e8
    
    # Solar luminosity L☉ (W)
    L_sun: float = 3.828e26
    
    # Earth mass M⊕ (kg)
    M_earth: float = 5.9722e24
    
    # Earth radius R⊕ (m)
    R_earth: float = 6.371e6
    
    # Astronomical Unit (m)
    AU: float = 1.495978707e11
    
    # Parsec (m)
    pc: float = 3.0856775814914e16
    
    # Kiloparsec (m)
    kpc: float = 3.0856775814914e19
    
    # Megaparsec (m)
    Mpc: float = 3.0856775814914e22
    
    # Light-year (m)
    ly: float = 9.4607304725808e15
    
    # Year (s)
    yr: float = 3.15576e7
    
    # Day (s)
    day: float = 86400.0
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # COSMOLOGICAL CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Hubble constant H₀ (s⁻¹) - 70 km/s/Mpc
    H0: float = 2.269e-18
    
    # Cosmological constant Λ (m⁻²)
    Lambda: float = 1.1e-52
    
    # Critical density ρ_crit (kg/m³) at z=0
    rho_crit: float = 9.47e-27
    
    # CMB temperature T_CMB (K)
    T_CMB: float = 2.7255
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UQFF-SPECIFIC CONSTANTS (Murphy Framework)
    # ═══════════════════════════════════════════════════════════════════════════
    
    # ═══════════════════════════════════════════════════════════════════════════
    # VACUUM DENSITY GRADIENT SYSTEM
    # ═══════════════════════════════════════════════════════════════════════════
    #
    # The UQFF framework uses TWO vacuum density scales that create a GRADIENT:
    #
    # 1. GRAVITATIONAL SCALE (rho_vac_UA): 7.09e-36 J/m³
    #    - Used in: Ug1-4 equations, cosmological terms, UQFF buoyancy
    #    - Represents: Cosmic vacuum energy density (dark energy scale)
    #    - Derived from: Λc²/8πG ≈ 5.96e-27 J/m³ modified by SCm factor
    #
    # 2. FIELD SCALE (rho_vac_UA_field): 1e-27 J/m³
    #    - Used in: Electric field terms, neutron production, magnetism
    #    - Represents: Local field coupling vacuum density
    #    - Derived from: Critical density ρ_c ≈ 9.47e-27 kg/m³
    #
    # GRADIENT RATIO: rho_vac_UA / rho_vac_UA_field = 7.09e-9
    #    - This ~10^9 ratio creates coupling between gravitational and field sectors
    #    - The gradient drives energy flow in DPM (Di-Pseudo-Monopole) interactions
    #    - Mathematically: ∇ρ_vac = (ρ_UA - ρ_field) / L_transition
    #
    # DO NOT "FIX" THIS BY UNIFYING - THE GRADIENT IS INTENTIONAL PHYSICS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Universal Aether vacuum density ρ_UA (J/m³) - GRAVITATIONAL SCALE
    # Used in: UQFF gravity (Ug1-4), buoyancy (Ub_i), cosmological terms
    rho_vac_UA: float = 7.09e-36
    
    # Universal Aether vacuum density ρ_UA_field (J/m³) - FIELD SCALE
    # Used in: Electric field (E = U_m/ρ_field/r), neutron production (η)
    rho_vac_UA_field: float = 1e-27
    
    # Vacuum density gradient ratio (dimensionless)
    # ∇ρ_ratio = rho_vac_UA / rho_vac_UA_field ≈ 7.09e-9
    rho_vac_gradient_ratio: float = 7.09e-9
    
    # SCm vacuum density ρ_SCm (J/m³) - Superconductive medium
    rho_vac_SCm: float = 7.09e-37
    
    # UQFF decay constant κ (day⁻¹)
    kappa: float = 0.0005
    
    # UQFF solvability [SSq]
    SSq: float = 0.57
    
    # SCm modulation factor H_SCm
    H_SCm: float = 0.99
    
    # Universal Aether correction U_UA
    U_UA: float = 0.0001
    
    # Eta coupling constant k_η
    k_eta: float = 1e-113
    
    # Beta buoyancy factor β_i
    beta_i: float = 0.603
    
    # F₀ reference force (N) - UQFF normalization
    F0: float = 1.83e71
    
    # Number of magnetic strings (compact objects)
    num_strings: float = 1e9
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UQFF COUPLING CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # LENR coupling k_LENR
    k_LENR: float = 1e-10
    
    # Activation coupling k_act
    k_act: float = 1e-14
    
    # Dark energy coupling k_DE
    k_DE: float = 1e-16
    
    # Neutron coupling k_neutron
    k_neutron: float = 1e-20
    
    # Relativistic coupling k_rel
    k_rel: float = 1e-12
    
    # Vacuum coupling k_vac
    k_vac: float = 1e-10
    
    # THz resonance coupling k_thz
    k_thz: float = 1e-15
    
    # Conduit coupling k_conduit
    k_conduit: float = 1e-18
    
    # Spooky (entanglement) coupling k_spooky
    k_spooky: float = 1e-20
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UQFF FREQUENCY CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # LENR resonance frequency ω_LENR (Hz)
    omega_LENR: float = 1.2e12
    
    # Magnetar rotation frequency Ω_g (rad/s) - SGR 1745 reference
    Omega_g: float = 7.3e-16
    
    # THz source frequency f_THz (Hz)
    f_THz: float = 1e12
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # MAGNETAR / EXTREME FIELD CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Critical (Schwinger) magnetic field B_crit (T)
    B_crit: float = 4.414e9
    
    # Critical magnetar field (T) - SGR type reference
    B_crit_magnetar: float = 4.4e13
    
    # Solar surface magnetic field (T)
    B_sun_surface: float = 1e-4
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # REFERENCE SYSTEM PARAMETERS (SGR 1745-2900)
    # ═══════════════════════════════════════════════════════════════════════════
    
    # SGR 1745 black hole mass M_bh (kg)
    M_bh_SGR1745: float = 8.155e36
    
    # SGR 1745 characteristic distance d_g (m)
    d_g_SGR1745: float = 2.55e20
    
    # SGR 1745 stellar wind density ρ_sw (kg/m³)
    rho_sw_SGR1745: float = 8e-21
    
    # SGR 1745 stellar wind efficiency ε_sw
    epsilon_sw_SGR1745: float = 0.001
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # DERIVED CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Bohr magneton μ_B = eℏ/(2m_e) (J/T)
    mu_B: float = 9.2740100783e-24
    
    # Electron Compton wavelength λ_C = h/(m_e c) (m)
    lambda_C: float = 2.42631023867e-12
    
    # Classical electron radius r_e (m)
    r_e: float = 2.8179403262e-15
    
    # Hartree energy E_H (J)
    E_H: float = 4.3597447222071e-18
    
    # Bohr radius a₀ (m)
    a_0: float = 5.29177210903e-11
    
    # Schwarzschild radius of 1 M☉ = 2GM/c² (m)
    r_s_sun: float = 2953.0
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # QUANTUM GRAVITY / BSM CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Planck mass m_P (kg)
    m_P: float = 2.176434e-8
    
    # Planck length l_P (m)
    l_P: float = 1.616255e-35
    
    # Planck time t_P (s)
    t_P: float = 5.391247e-44
    
    # Planck temperature T_P (K)
    T_P: float = 1.416784e32
    
    # Hubble time t_H = 1/H₀ (s) ≈ 13.8 Gyr
    t_H: float = 4.35e17
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # INFORMATION PARADOX CONSTANTS (Batch 21)
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Position uncertainty (26D compact dimensions)
    Delta_x_p: float = 1e-68
    
    # Wavefunction integral normalization
    integral_psi: float = 2.176e-18


# Create singleton instance
_CONSTANTS_INSTANCE = Constants()


def _to_dict(obj: Constants) -> Dict[str, float]:
    """Convert Constants dataclass to dictionary."""
    return {
        field: getattr(obj, field) 
        for field in obj.__dataclass_fields__
    }


# Dictionary interface for compatibility with existing code
CONSTANTS: Dict[str, float] = _to_dict(_CONSTANTS_INSTANCE)


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE EXPORTS
# ═══════════════════════════════════════════════════════════════════════════════

# Frequently used constants available at module level
c = CONSTANTS['c']
G = CONSTANTS['G']
hbar = CONSTANTS['hbar']
M_sun = CONSTANTS['M_sun']
epsilon_0 = CONSTANTS['epsilon_0']
mu_0 = CONSTANTS['mu_0']
k_B = CONSTANTS['k_B']
m_e = CONSTANTS['m_e']
e_charge = CONSTANTS['e']


def get_constant(name: str) -> float:
    """
    Get a constant by name with case-insensitive lookup.
    
    Args:
        name: Constant name (case-insensitive for common names)
        
    Returns:
        Constant value
        
    Raises:
        KeyError: If constant not found
    """
    # Try exact match first
    if name in CONSTANTS:
        return CONSTANTS[name]
    
    # Try lowercase
    name_lower = name.lower()
    for key in CONSTANTS:
        if key.lower() == name_lower:
            return CONSTANTS[key]
    
    raise KeyError(f"Unknown constant: {name}")


def list_constants(category: str = None) -> Dict[str, float]:
    """
    List constants, optionally filtered by category.
    
    Args:
        category: Optional filter - 'fundamental', 'astro', 'uqff', 'cosmo', etc.
        
    Returns:
        Dictionary of matching constants
    """
    category_prefixes = {
        'fundamental': ('c', 'G', 'hbar', 'h_planck', 'e', 'k_B', 'epsilon_0', 
                       'mu_0', 'alpha', 'N_A', 'sigma_SB'),
        'particle': ('m_e', 'm_p', 'm_n', 'u'),
        'astro': ('M_sun', 'R_sun', 'L_sun', 'M_earth', 'R_earth', 'AU', 
                 'pc', 'kpc', 'Mpc', 'ly', 'yr', 'day'),
        'cosmo': ('H0', 'Lambda', 'rho_crit', 'T_CMB'),
        'uqff': ('rho_vac_UA', 'rho_vac_SCm', 'kappa', 'SSq', 'H_SCm', 'U_UA',
                'k_eta', 'beta_i', 'F0', 'num_strings', 'k_LENR', 'k_act',
                'k_DE', 'k_neutron', 'k_rel', 'k_vac', 'k_thz', 'k_conduit',
                'k_spooky', 'omega_LENR', 'Omega_g', 'f_THz'),
        'magnetar': ('B_crit', 'B_crit_magnetar', 'B_sun_surface', 'M_bh_SGR1745',
                    'd_g_SGR1745', 'rho_sw_SGR1745', 'epsilon_sw_SGR1745'),
        'derived': ('mu_B', 'lambda_C', 'r_e', 'E_H', 'a_0', 'r_s_sun'),
        'qg': ('m_P', 'l_P', 't_P', 'T_P', 't_H', 'Delta_x_p', 'integral_psi'),
    }
    
    if category is None:
        return CONSTANTS.copy()
    
    if category not in category_prefixes:
        raise ValueError(f"Unknown category: {category}. "
                        f"Available: {list(category_prefixes.keys())}")
    
    return {k: CONSTANTS[k] for k in category_prefixes[category] if k in CONSTANTS}


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════

def validate_constants() -> Dict[str, Any]:
    """
    Validate physical constants for consistency.
    
    Returns:
        Dictionary with validation results
    """
    results = {'valid': True, 'checks': [], 'errors': []}
    
    # Check c² = 1/(μ₀ε₀)
    c_from_em = 1.0 / math.sqrt(CONSTANTS['mu_0'] * CONSTANTS['epsilon_0'])
    c_err = abs(c_from_em - CONSTANTS['c']) / CONSTANTS['c']
    results['checks'].append({
        'name': 'c = 1/√(μ₀ε₀)',
        'computed': c_from_em,
        'stored': CONSTANTS['c'],
        'rel_error': c_err,
        'pass': c_err < 1e-6
    })
    if c_err >= 1e-6:
        results['valid'] = False
        results['errors'].append(f"c consistency check failed: {c_err:.2e}")
    
    # Check Bohr magneton μ_B = eℏ/(2m_e)
    mu_B_computed = (CONSTANTS['e'] * CONSTANTS['hbar']) / (2 * CONSTANTS['m_e'])
    mu_B_err = abs(mu_B_computed - CONSTANTS['mu_B']) / CONSTANTS['mu_B']
    results['checks'].append({
        'name': 'μ_B = eℏ/(2m_e)',
        'computed': mu_B_computed,
        'stored': CONSTANTS['mu_B'],
        'rel_error': mu_B_err,
        'pass': mu_B_err < 1e-6
    })
    if mu_B_err >= 1e-6:
        results['valid'] = False
        results['errors'].append(f"μ_B consistency check failed: {mu_B_err:.2e}")
    
    # Check Schwarzschild radius r_s = 2GM/c²
    r_s_computed = 2 * CONSTANTS['G'] * CONSTANTS['M_sun'] / CONSTANTS['c']**2
    r_s_err = abs(r_s_computed - CONSTANTS['r_s_sun']) / CONSTANTS['r_s_sun']
    results['checks'].append({
        'name': 'r_s = 2GM☉/c²',
        'computed': r_s_computed,
        'stored': CONSTANTS['r_s_sun'],
        'rel_error': r_s_err,
        'pass': r_s_err < 1e-3
    })
    
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("UQFF Shared Constants")
    print("=" * 60)
    
    print("\nFundamental Constants:")
    for name, value in list_constants('fundamental').items():
        print(f"  {name:15s} = {value:.10e}")
    
    print("\nUQFF Constants:")
    for name, value in list_constants('uqff').items():
        print(f"  {name:20s} = {value:.6e}")
    
    print("\nValidation:")
    validation = validate_constants()
    print(f"  Overall: {'PASS' if validation['valid'] else 'FAIL'}")
    for check in validation['checks']:
        status = '✓' if check['pass'] else '✗'
        print(f"  {status} {check['name']}: error = {check['rel_error']:.2e}")
    
    print(f"\nTotal constants defined: {len(CONSTANTS)}")
