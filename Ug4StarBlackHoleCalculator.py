#!/usr/bin/env python3
"""
Ug4StarBlackHoleCalculator.py - Complete Star-Black Hole Vacuum Interactions

Implements Ug4 (Universal Gravity component 4) for modeling vacuum-mediated
forces between stars and supermassive black holes (SMBHs) at galactic scales.

Integration: Grok Thread b9a29cedc27b45dfa309ea1705721bf0 (March 5, 2026)
Source: "Star Magic: The Quest for Unity" - Daniel T. Murphy ©2025

Formula:
    Ug4 = k4 * ρ_vac * ([SCm] * M_bh) / d_g * e^(-α*t) * cos(π*t_n) * (1 + f_feedback)

Where:
    k4 = 1.0 (coupling constant for Ug4, unitless)
    ρ_vac = ~10^-9 J/m³ (vacuum energy density, cosmological constant)
    [SCm] = ~10^15 kg/m³ (Superconducting Material concentration)
    M_bh = black hole mass (kg) - e.g. 8.55×10^36 kg for Sagittarius A*
    d_g = distance from galactic center (m) - e.g. 2.55×10^20 m for Sun
    α = 0.001 day^-1 (non-linear time decay rate)
    t = time (days)
    t_n = negative time (t - t_0), enables temporal reversal
    f_feedback = feedback factor (0.1 for ΔM_BH = 1 dex)

Physics:
    - Operates at quantum levels 20-26 (galactic vacuum scales)
    - Models observable Sun-Sgr A* interactions via vacuum energy
    - Incorporates SCm modulation of vacuum fluctuations
    - Feedback loops model AGN influence on stellar orbits
    - Temporal oscillations (cos π·t_n) introduce periodicity

Applications:
    - Stellar orbit perturbations around SMBHs
    - Galactic structure formation
    - Dark energy effects at galactic scales
    - Quasar dynamics and black hole accretion feedback

Author: GitHub Copilot (implementation from thread analysis)
Date: March 5, 2026
"""

import math
import numpy as np
from typing import Dict, Any, List, Tuple, Optional

# ═══════════════════════════════════════════════════════════════
# PHYSICAL CONSTANTS (CGS + SI)
# ═══════════════════════════════════════════════════════════════

# Fundamental constants
SPEED_OF_LIGHT = 2.99792458e8  # m/s
PLANCK_CONSTANT_REDUCED = 1.054571817e-34  # J·s (ℏ)
GRAVITATIONAL_CONSTANT = 6.67430e-11  # m³/kg·s² (G)

# UQFF-specific constants
K_4 = 1.0  # Ug4 coupling constant (unitless, baseline for vacuum forces)
RHO_VAC_COSMOLOGICAL = 1e-9  # J/m³ (vacuum energy density from cosmological constant)
SCM_CONCENTRATION_STELLAR = 1e15  # kg/m³ (SCm concentration in stellar cores)
ALPHA_DECAY = 0.001  # day^-1 (non-linear time decay rate)

# Sagittarius A* (EHT 2024-2025 observations) - updated from 8.15×10^36 kg
M_BH_SGR_A_STAR = 8.55e36  # kg (4.3 × 10^6 M_sun)

# Sun to galactic center distance
D_G_SUN_TO_SGR_A = 2.55e20  # m (27,000 light-years)

# Galactic dynamics
OMEGA_GALACTIC = 7.3e-16  # rad/s (galactic spin rate at Sun's position)

# ═══════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════

def compute_feedback_factor(delta_M_BH_dex: float, base_factor: float = 0.1) -> float:
    """
    Compute feedback factor for black hole mass change.
    
    Formula: f_feedback = base_factor * delta_M_BH_dex
    
    Args:
        delta_M_BH_dex: Black hole mass change in dex (log10 scale)
                        e.g. 1.0 dex = tenfold mass increase
        base_factor: Base feedback strength (default 0.1 for ΔM_BH = 1 dex)
    
    Returns:
        Feedback factor (unitless)
    
    Physical meaning:
        Models AGN feedback regulating star formation and galactic dynamics.
        Larger accretion events (higher ΔM_BH) → stronger feedback.
    
    Example:
        >>> compute_feedback_factor(1.0)  # Tenfold mass increase
        0.1
        >>> compute_feedback_factor(0.5)  # √10 ≈ 3.16x mass increase
        0.05
    """
    return base_factor * delta_M_BH_dex


def compute_negative_time(t: float, t_0: float = 0.0) -> float:
    """
    Compute negative time t_n for temporal reversal effects.
    
    Formula: t_n = t - t_0
    
    Args:
        t: Current time (days)
        t_0: Reference time (days, default 0)
    
    Returns:
        Negative time t_n (days)
    
    Physical meaning:
        Enables temporal reversal in quasar dynamics and black hole accretion.
        When t < t_0, t_n is negative → past influences future via cos(π·t_n).
    
    Example:
        >>> compute_negative_time(100, 0)
        100.0
        >>> compute_negative_time(-50, 0)  # Past time
        -50.0
    """
    return t - t_0


def compute_time_decay_factor(t: float, alpha: float = ALPHA_DECAY) -> float:
    """
    Compute non-linear time decay factor.
    
    Formula: e^(-α*t)
    
    Args:
        t: Time (days)
        alpha: Decay rate (day^-1, default 0.001)
    
    Returns:
        Decay factor (unitless, 0 < factor ≤ 1)
    
    Physical meaning:
        Models gradual weakening of vacuum-mediated forces over time.
        α = 0.001 day^-1 → ~0.998 per day decay (slow, galactic timescales).
    
    Example:
        >>> compute_time_decay_factor(0)
        1.0
        >>> compute_time_decay_factor(1000)  # ~1000 days ≈ 2.7 years
        0.36787944117144233  # e^-1
    """
    return math.exp(-alpha * t)


def compute_temporal_cycle(t_n: float) -> float:
    """
    Compute π-cycle temporal modulation.
    
    Formula: cos(π * t_n)
    
    Args:
        t_n: Negative time (days)
    
    Returns:
        Cycle factor (-1 ≤ factor ≤ 1)
    
    Physical meaning:
        Introduces periodicity with half-period = 1 day (π cycle).
        Links to Riemann zeta function zeros (speculative).
        Oscillates between attraction (cos > 0) and repulsion (cos < 0).
    
    Example:
        >>> compute_temporal_cycle(0)
        1.0  # Maximum attraction
        >>> compute_temporal_cycle(0.5)
        6.123233995736766e-17  # ≈ 0 (neutral)
        >>> compute_temporal_cycle(1)
        -1.0  # Maximum repulsion
    """
    return math.cos(math.pi * t_n)


# ═══════════════════════════════════════════════════════════════
# MAIN UG4 CALCULATOR CLASS
# ═══════════════════════════════════════════════════════════════

class Ug4StarBlackHoleCalculator:
    """
    Calculator for Ug4 star-black hole vacuum interactions.
    
    Implements complete Ug4 formula from "Star Magic: The Quest for Unity"
    with all 8 parameters for galactic-scale vacuum-mediated forces.
    
    Capabilities:
        - Single-point Ug4 calculation
        - Time-series evolution (stellar orbit modeling)
        - Parametric studies (vary SCm, M_bh, d_g, etc.)
        - Negative time exploration (temporal reversal)
        - Feedback analysis (accretion events)
    
    Usage:
        calc = Ug4StarBlackHoleCalculator()
        result = calc.compute({
            'M_bh': 8.55e36,  # Sgr A* mass
            'd_g': 2.55e20,   # Sun-Sgr A* distance
            't': 0,           # Current time
            't_0': 0,         # Reference time
            'delta_M_BH_dex': 0.0  # No accretion
        })
    """
    
    def __init__(self):
        """Initialize Ug4 calculator with default UQFF parameters."""
        self.k4 = K_4
        self.rho_vac = RHO_VAC_COSMOLOGICAL
        self.SCm_concentration = SCM_CONCENTRATION_STELLAR
        self.alpha_decay = ALPHA_DECAY
        self.base_feedback_factor = 0.1
    
    def compute(self, dataset: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compute Ug4 star-black hole vacuum interaction.
        
        Args:
            dataset: Dictionary with keys:
                - M_bh: Black hole mass (kg) [required]
                - d_g: Distance from galactic center (m) [required]
                - t: Time (days) [default: 0.0]
                - t_0: Reference time (days) [default: 0.0]
                - delta_M_BH_dex: BH mass change in dex [default: 0.0]
                - SCm_concentration: SCm density override (kg/m³) [optional]
                - k4: Coupling constant override [optional]
                - rho_vac: Vacuum energy override (J/m³) [optional]
                - alpha_decay: Decay rate override (day^-1) [optional]
        
        Returns:
            Dictionary with keys:
                - Ug4: Energy density (J/m³)
                - Ug4_N_per_m2: Force per unit area (N/m²)
                - components: Dict of individual terms
                - parameters_used: Dict of all input parameters
                - equation: Long-form equation string
        
        Example:
            >>> calc = Ug4StarBlackHoleCalculator()
            >>> result = calc.compute({'M_bh': 8.55e36, 'd_g': 2.55e20})
            >>> print(f"Ug4 = {result['Ug4']:.3e} J/m³")
        """
        # Extract parameters with defaults
        M_bh = dataset['M_bh']  # Required
        d_g = dataset['d_g']    # Required
        t = dataset.get('t', 0.0)
        t_0 = dataset.get('t_0', 0.0)
        delta_M_BH_dex = dataset.get('delta_M_BH_dex', 0.0)
        
        # Optional overrides
        SCm_conc = dataset.get('SCm_concentration', self.SCm_concentration)
        k4 = dataset.get('k4', self.k4)
        rho_vac = dataset.get('rho_vac', self.rho_vac)
        alpha = dataset.get('alpha_decay', self.alpha_decay)
        
        # Compute temporal variables
        t_n = compute_negative_time(t, t_0)
        decay_factor = compute_time_decay_factor(t, alpha)
        cycle_factor = compute_temporal_cycle(t_n)
        f_feedback = compute_feedback_factor(delta_M_BH_dex, self.base_feedback_factor)
        
        # Ug4 formula components
        numerator = rho_vac * (SCm_conc * M_bh)
        denominator = d_g
        base_term = k4 * (numerator / denominator)
        temporal_term = decay_factor * cycle_factor
        feedback_term = 1 + f_feedback
        
        # Complete Ug4 calculation
        Ug4 = base_term * temporal_term * feedback_term
        
        # Convert to force per unit area (N/m²) for some applications
        # Using relation: Energy density [J/m³] = Pressure [N/m²]
        Ug4_N_per_m2 = Ug4
        
        # Long-form equation with solutions
        equation = (
            f"Ug4 Star-Black Hole Vacuum Interaction\n\n"
            f"Formula:\n"
            f"  Ug4 = k4 × ρ_vac × ([SCm] × M_bh)/d_g × e^(-α×t) × cos(π×t_n) × (1 + f_feedback)\n\n"
            f"Parameters:\n"
            f"  k4 = {k4} (coupling constant, unitless)\n"
            f"  ρ_vac = {rho_vac:.3e} J/m³ (vacuum energy density)\n"
            f"  [SCm] = {SCm_conc:.3e} kg/m³ (SCm concentration)\n"
            f"  M_bh = {M_bh:.3e} kg (black hole mass)\n"
            f"  d_g = {d_g:.3e} m (galactic distance)\n"
            f"  α = {alpha:.4f} day^-1 (decay rate)\n"
            f"  t = {t} days (current time)\n"
            f"  t_0 = {t_0} days (reference time)\n"
            f"  t_n = {t_n} days (negative time)\n"
            f"  ΔM_BH = {delta_M_BH_dex} dex (BH mass change)\n"
            f"  f_feedback = {f_feedback:.3f} (feedback factor)\n\n"
            f"Intermediate Calculations:\n"
            f"  e^(-α×t) = e^(-{alpha}×{t}) = {decay_factor:.6f}\n"
            f"  cos(π×t_n) = cos(π×{t_n}) = {cycle_factor:.6f}\n"
            f"  (1 + f_feedback) = 1 + {f_feedback:.3f} = {feedback_term:.3f}\n"
            f"  [SCm]×M_bh = {SCm_conc:.3e} × {M_bh:.3e} = {SCm_conc*M_bh:.3e} kg²/m³\n\n"
            f"Solution:\n"
            f"  Ug4 = {Ug4:.6e} J/m³ (energy density)\n"
            f"  Ug4 = {Ug4_N_per_m2:.6e} N/m² (force per unit area)\n"
        )
        
        return {
            'Ug4': Ug4,
            'Ug4_N_per_m2': Ug4_N_per_m2,
            'components': {
                'k4': k4,
                'rho_vac': rho_vac,
                'SCm_M_bh_product': SCm_conc * M_bh,
                'd_g': d_g,
                'decay_factor': decay_factor,
                'cycle_factor': cycle_factor,
                'feedback_term': feedback_term,
                'temporal_term': temporal_term,
                'base_term': base_term
            },
            'time_variables': {
                't': t,
                't_0': t_0,
                't_n': t_n,
                'alpha': alpha
            },
            'parameters_used': {
                'M_bh': M_bh,
                'd_g': d_g,
                'SCm_concentration': SCm_conc,
                'k4': k4,
                'rho_vac': rho_vac,
                'delta_M_BH_dex': delta_M_BH_dex,
                'f_feedback': f_feedback
            },
            'equation': equation,
            'calculator': 'Ug4StarBlackHoleCalculator',
            'thread_source': 'b9a29cedc27b45dfa309ea1705721bf0',
            'integration_date': '2026-03-05'
        }
    
    def compute_time_series(self, 
                          M_bh: float, 
                          d_g: float, 
                          t_start: float = 0.0, 
                          t_end: float = 365.0, 
                          n_points: int = 100,
                          t_0: float = 0.0,
                          delta_M_BH_dex: float = 0.0) -> Dict[str, Any]:
        """
        Compute Ug4 evolution over time (stellar orbit modeling).
        
        Args:
            M_bh: Black hole mass (kg)
            d_g: Distance from galactic center (m)
            t_start: Start time (days)
            t_end: End time (days)
            n_points: Number of time points
            t_0: Reference time for negative time (days)
            delta_M_BH_dex: Black hole mass change (dex)
        
        Returns:
            Dictionary with keys:
                - t_array: Time points (days)
                - Ug4_array: Ug4 values (J/m³)
                - mean_Ug4: Mean over time period
                - max_Ug4: Maximum value
                - min_Ug4: Minimum value
                - parameters: Input parameters
        
        Example:
            >>> calc = Ug4StarBlackHoleCalculator()
            >>> series = calc.compute_time_series(8.55e36, 2.55e20, 0, 365)
            >>> plt.plot(series['t_array'], series['Ug4_array'])
        """
        t_array = np.linspace(t_start, t_end, n_points)
        Ug4_array = []
        
        for t in t_array:
            result = self.compute({
                'M_bh': M_bh,
                'd_g': d_g,
                't': t,
                't_0': t_0,
                'delta_M_BH_dex': delta_M_BH_dex
            })
            Ug4_array.append(result['Ug4'])
        
        Ug4_array = np.array(Ug4_array)
        
        return {
            't_array': t_array,
            'Ug4_array': Ug4_array,
            'mean_Ug4': np.mean(Ug4_array),
            'max_Ug4': np.max(Ug4_array),
            'min_Ug4': np.min(Ug4_array),
            'std_Ug4': np.std(Ug4_array),
            'parameters': {
                'M_bh': M_bh,
                'd_g': d_g,
                't_start': t_start,
                't_end': t_end,
                'n_points': n_points,
                't_0': t_0,
                'delta_M_BH_dex': delta_M_BH_dex
            },
            'calculator': 'Ug4StarBlackHoleCalculator.compute_time_series'
        }


# ═══════════════════════════════════════════════════════════════
# PREDEFINED ASTROPHYSICAL SYSTEMS
# ═══════════════════════════════════════════════════════════════

# Sagittarius A* (Milky Way SMBH) - EHT 2024-2025 observations
SGR_A_STAR_SYSTEM = {
    'name': 'Sagittarius A*',
    'M_bh': M_BH_SGR_A_STAR,  # 8.55×10^36 kg (4.3×10^6 M_sun)
    'd_g': D_G_SUN_TO_SGR_A,  # 2.55×10^20 m (27,000 ly from Sun)
    'description': 'Milky Way supermassive black hole, updated EHT 2024-2025',
    'ref_time': 2025.0  # Reference epoch (years)
}

# M87* (Virgo A) - First black hole image
M87_STAR_SYSTEM = {
    'name': 'M87*',
    'M_bh': 1.3e40,  # kg (6.5×10^9 M_sun)
    'd_g': 1.65e25,  # m (~55 million ly from Earth)
    'description': 'M87 galaxy SMBH, Event Horizon Telescope first image 2019'
}

# Cygnus X-1 (stellar-mass black hole)
CYGNUS_X1_SYSTEM = {
    'name': 'Cygnus X-1',
    'M_bh': 4.3e31,  # kg (~21 M_sun)
    'd_g': 1.86e19,  # m (~6,100 ly from Earth)
    'description': 'Stellar-mass black hole in X-ray binary system'
}


# ═══════════════════════════════════════════════════════════════
# MODULE EXPORTS
# ═══════════════════════════════════════════════════════════════

__all__ = [
    'Ug4StarBlackHoleCalculator',
    'compute_feedback_factor',
    'compute_negative_time',
    'compute_time_decay_factor',
    'compute_temporal_cycle',
    'SGR_A_STAR_SYSTEM',
    'M87_STAR_SYSTEM',
    'CYGNUS_X1_SYSTEM',
    'K_4',
    'RHO_VAC_COSMOLOGICAL',
    'SCM_CONCENTRATION_STELLAR',
    'ALPHA_DECAY',
    'M_BH_SGR_A_STAR',
    'D_G_SUN_TO_SGR_A'
]


if __name__ == '__main__':
    """Run validation examples."""
    
    print("═" * 70)
    print("Ug4 Star-Black Hole Calculator Validation")
    print("═" * 70)
    
    calc = Ug4StarBlackHoleCalculator()
    
    # Example 1: Sun-Sagittarius A* at t=0
    print("\nExample 1: Sun-Sgr A* Interaction (t=0)\n")
    result1 = calc.compute(SGR_A_STAR_SYSTEM)
    print(result1['equation'])
    
    # Example 2: After 1000 days
    print("\n" + "═" * 70)
    print("\nExample 2: Sun-Sgr A* After 1000 Days\n")
    result2 = calc.compute({**SGR_A_STAR_SYSTEM, 't': 1000})
    print(f"Ug4(t=1000) = {result2['Ug4']:.6e} J/m³")
    print(f"Percent change from t=0: {100*(result2['Ug4']/result1['Ug4'] - 1):.2f}%")
    
    # Example 3: Feedback from accretion event
    print("\n" + "═" * 70)
    print("\nExample 3: Tenfold Accretion Event (ΔM_BH = 1 dex)\n")
    result3 = calc.compute({**SGR_A_STAR_SYSTEM, 'delta_M_BH_dex': 1.0})
    print(f"Ug4 with feedback = {result3['Ug4']:.6e} J/m³")
    print(f"Feedback amplification: {100*(result3['Ug4']/result1['Ug4'] - 1):.1f}%")
    
    print("\n" + "═" * 70)
    print("Validation complete. All examples executed successfully.")
    print("═" * 70)
