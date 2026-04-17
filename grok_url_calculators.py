#!/usr/bin/env python3
"""
grok_url_calculators.py - 121 Calculator Classes from Grok URL Equations
=========================================================================

Source: https://x.com/i/grok/share/683542a41e744554928bfcd8b0a19e40
Contains: 100 Standard Physics + 7 UQFF Framework + 6 MUGE System + 8 Updated UQFF

Pipeline extension of CondensedPhysics.py for the complete Grok URL equation set.
Each Calculator class follows the compute(dataset) -> dict pattern.

ARCHITECTURE (MANDATORY):
    - NO hardcoded astrophysical data
    - NO named system classes
    - ALL parameters via dataset dict
    - PURE stateless calculators

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
Created: March 1, 2026
"""

import math
from dpm_helpers import dpm_emergent_ug1, dpm_emergent_ug2

import numpy as np
from typing import Dict, Any, Optional

# ═══════════════════════════════════════════════════════════════════════════════
# FUNDAMENTAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

CONST = {
    'G': 6.674e-11,          # m³/kg/s² - Gravitational constant
    'c': 2.998e8,            # m/s - Speed of light
    'hbar': 1.055e-34,       # J·s - Reduced Planck constant
    'k_B': 1.381e-23,        # J/K - Boltzmann constant
    'e_charge': 1.602e-19,   # C - Elementary charge
    'm_p': 1.673e-27,        # kg - Proton mass
    'm_e': 9.109e-31,        # kg - Electron mass
    'sigma_T': 6.652e-29,    # m² - Thomson cross-section
    'M_sun': 1.989e30,       # kg - Solar mass
    'M_pl': 2.176e-8,        # kg - Planck mass
    't_pl': 5.391e-44,       # s - Planck time
    'l_pl': 1.616e-35,       # m - Planck length
    'mu_0': 4 * math.pi * 1e-7,  # H/m - Vacuum permeability
    'pi': math.pi,
}

# UQFF-specific constants
UQFF_CONST = {
    'F_rel': 4.30e33,           # N - Relativistic coherence force (from LEP)
    'E_LEP': 200.0,             # GeV - LEP baseline energy
    'rho_vac_SCm': 7.09e-37,    # J/m³ - Superconductive vacuum density
    'rho_vac_UA': 7.09e-36,     # J/m³ - Aether vacuum density
    'gamma_decay': 5e-5,        # day⁻¹ - Decay constant
    'omega_c': 1.585e-8,        # rad/s - Oscillation frequency
    'f_Heav': 0.01,             # Heaviside modulation
    'f_quasi': 0.01,            # Quasi-particle factor
    'P_SCm': 1.0,               # Superconductive probability
    'mu_base': 3.38e20,         # T·pm³ - Base magnetic moment
    'Q_wave': 1e12,             # THz resonance factor
}

# Cosmological parameters
COSMO = {
    'H_0': 70.0,         # km/s/Mpc
    'Omega_m': 0.3,       # Matter density parameter
    'Omega_Lambda': 0.7,  # Dark energy density parameter
    'Lambda': 1e-52,      # m⁻² - Cosmological constant
    'delta_c': 1.686,     # Critical overdensity
    'rho_c': 1e-26,       # kg/m³ - Critical density
    'eta_BBN': 6e-10,     # Baryon-to-photon ratio
    'n_gamma': 4.1e8,     # m⁻³ - CMB photon density
}


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 1-4: PROTOSTELLAR JETS AND OUTFLOWS
# ═══════════════════════════════════════════════════════════════════════════════

class AngularMomentumTransportCalculator:
    """Eq. 1: dL/dt = Ṁr²Ω - T_B (Angular Momentum Transport in accretion disks)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        M_star = dataset.get('M_star', 1.0 * CONST['M_sun'])
        M_dot = dataset.get('M_dot', 1e-7 * CONST['M_sun'] / 3.156e7)  # M_sun/yr → kg/s
        r = dataset.get('r', 10.0 * 1.496e11)  # AU → m
        B = dataset.get('B', 1e-3)  # T
        Omega = math.sqrt(G * M_star / r**3) if r > 0 else 0
        T_B = B**2 * r**3 / (4 * CONST['pi'])  # Magnetic torque
        dLdt = M_dot * r**2 * Omega - T_B
        L = M_dot * r**2 * Omega  # Angular momentum addition rate
        return {
            'dL_dt': dLdt, 'L_accretion': L, 'T_B': T_B,
            'Omega': Omega, 'r': r, 'M_dot': M_dot,
            'equation': 'dL/dt = Ṁr²Ω - T_B',
            'source': 'Grok URL Eq. 1: Angular Momentum Transport'
        }

class MHDJetVelocityCalculator:
    """Eq. 2: v_j ≈ v_K(r_A/r_0)^(1/2) (MHD Jet Velocity)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        M_star = dataset.get('M_star', 1.0 * CONST['M_sun'])
        r_0 = dataset.get('r_0', 0.1 * 1.496e11)  # AU → m
        r_A = dataset.get('r_A', 30.0 * 1.496e11)  # Alfvén radius
        v_K = math.sqrt(G * M_star / r_0) if r_0 > 0 else 0
        v_j = v_K * math.sqrt(r_A / r_0) if r_0 > 0 else 0
        return {
            'v_j': v_j, 'v_K': v_K, 'r_A': r_A, 'r_0': r_0,
            'v_j_km_s': v_j / 1e3,
            'equation': 'v_j ≈ v_K × (r_A/r_0)^(1/2)',
            'source': 'Grok URL Eq. 2: MHD Jet Velocity'
        }

class JTypeShockRankineHugoniotCalculator:
    """Eq. 3: Rankine-Hugoniot jump conditions (J-type shock)."""
    def compute(self, dataset: dict) -> dict:
        rho_1 = dataset.get('rho_1', 1e-18)  # kg/m³
        v_1 = dataset.get('v_1', 300e3)  # m/s
        P_1 = dataset.get('P_1', 1e-10)  # Pa
        gamma = dataset.get('gamma', 5.0/3.0)
        # Strong shock limit
        compression = (gamma + 1) / (gamma - 1)
        rho_2 = rho_1 * compression
        v_2 = v_1 / compression
        P_2 = 2 * rho_1 * v_1**2 / (gamma + 1)
        T_ratio = P_2 / P_1 * rho_1 / rho_2 if rho_2 > 0 and P_1 > 0 else 0
        return {
            'rho_2': rho_2, 'v_2': v_2, 'P_2': P_2,
            'compression_ratio': compression, 'T_ratio': T_ratio,
            'equation': 'ρ₁v₁=ρ₂v₂; ρ₁v₁²+P₁=ρ₂v₂²+P₂',
            'source': 'Grok URL Eq. 3: J-type Shock Rankine-Hugoniot'
        }

class CTypeShockDampingCalculator:
    """Eq. 4: v(z) ≈ v_s exp(-z/L_d) (C-type shock damping)."""
    def compute(self, dataset: dict) -> dict:
        v_s = dataset.get('v_s', 50e3)  # m/s shock velocity
        L_d = dataset.get('L_d', 1e14)  # m damping length
        z = dataset.get('z_position', 1e14)  # m position
        v_z = v_s * math.exp(-z / L_d) if L_d > 0 else 0
        return {
            'v_z': v_z, 'v_s': v_s, 'L_d': L_d, 'z': z,
            'damping_fraction': v_z / v_s if v_s > 0 else 0,
            'equation': 'v(z) ≈ v_s × exp(-z/L_d)',
            'source': 'Grok URL Eq. 4: C-type Shock Damping'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 5-7: GALAXY MERGERS AND SFR
# ═══════════════════════════════════════════════════════════════════════════════

class EPSMergerRateCalculator:
    """Eq. 5: Extended Press-Schechter merger rate."""
    def compute(self, dataset: dict) -> dict:
        sigma_M = dataset.get('sigma_M', 0.8)
        sigma_m = dataset.get('sigma_m', 1.2)
        delta_c = dataset.get('delta_c', COSMO['delta_c'])
        z = dataset.get('z', 0.0)
        delta_c_z = delta_c * (1 + z)  # Simplified growth
        dsigma = sigma_m**2 - sigma_M**2
        if dsigma > 0:
            exponent = -delta_c_z**2 / (2 * dsigma)
            rate = math.sqrt(2 / CONST['pi']) * (sigma_M / sigma_m) * abs(delta_c_z) * math.exp(exponent) / math.sqrt(dsigma)
        else:
            rate = 0.0
        return {
            'merger_rate_normalized': rate, 'delta_c_z': delta_c_z,
            'sigma_diff_sq': dsigma, 'z': z,
            'equation': 'dN/dtdM ∝ (σ_M/σ_m)|dδ_c/dz| exp(-δ_c²/2(σ_m²-σ_M²))',
            'source': 'Grok URL Eq. 5: EPS Merger Rate'
        }

class OrbitalTorqueTimeCalculator:
    """Eq. 6: t_orb = 2π√(r³/GM) (Orbital/torque timescale)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        r = dataset.get('r', 50e3 * 3.086e16)  # kpc → m
        M = dataset.get('M', 1e12 * CONST['M_sun'])
        t_orb = 2 * CONST['pi'] * math.sqrt(r**3 / (G * M)) if M > 0 else 0
        return {
            't_orb_s': t_orb, 't_orb_Myr': t_orb / 3.156e13,
            'equation': 't_orb = 2π√(r³/GM)',
            'source': 'Grok URL Eq. 6: Orbital Torque Time'
        }

class SFRDEvolutionCalculator:
    """Eq. 7: SFRD ∝ (1+z)^2.7 (Star Formation Rate Density evolution)."""
    def compute(self, dataset: dict) -> dict:
        z = dataset.get('z', 0.0)
        SFRD_0 = dataset.get('SFRD_0', 0.015)  # M_sun/yr/Mpc³
        alpha = dataset.get('alpha', 2.7)
        SFRD = SFRD_0 * (1 + z)**alpha if z <= 2 else SFRD_0 * 3**alpha * (1 + z)**(-alpha)
        return {
            'SFRD': SFRD, 'z': z, 'alpha': alpha,
            'equation': 'SFRD ∝ (1+z)^2.7 (z≤2)',
            'source': 'Grok URL Eq. 7: SFRD Evolution'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 8-9: BLACK HOLE GROWTH
# ═══════════════════════════════════════════════════════════════════════════════

class EPSBHMassFunctionCalculator:
    """Eq. 8: EPS BH mass function N(>M,z) with Harvard distribution integration."""
    
    def __init__(self):
        self.harvard_data = None
    
    def load_harvard_distribution(self, json_path: str = None) -> dict:
        """
        Load Harvard energy_distributions.json with BH mass bins.
        
        Args:
            json_path: Path to energy_distributions.json (default: data/harvard_energy_distributions.json)
            
        Returns:
            {
                'bins': array of bin centers [2.75, 3.5, ..., 8.0] in log(M_sun),
                'peaks': {
                    3.9: 'Neutron drop resonance (Kozima, 10^3-10^4 s lifetime)',
                    5.6: 'THz coherence (1.2 THz LENR, Colman-Gillespie 300 Hz)',
                    6.5: 'Vacuum feedback (Ug4 galactic, levels 20-26)',
                    7.0: 'SMBH transition (M ~ 10^7 M_sun)'
                },
                'integral': 0.064,
                'distribution': array of density values,
                'cumulative': array of cumulative distribution
            }
        """
        import json
        import os
        
        if json_path is None:
            # Default to data/ directory relative to this file
            script_dir = os.path.dirname(os.path.abspath(__file__))
            json_path = os.path.join(script_dir, 'data', 'harvard_energy_distributions.json')
        
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Extract UQFF peak interpretations
            peaks_uqff = {}
            for peak_key, peak_data in data['peaks'].items():
                log_M = float(peak_key)
                peaks_uqff[log_M] = peak_data['interpretation']
            
            # Store loaded data
            self.harvard_data = {
                'bins': data['bins']['centers'],
                'bin_edges': data['bins']['edges'],
                'peaks': peaks_uqff,
                'peak_details': data['peaks'],
                'density': data['distribution']['density'],
                'cumulative': data['distribution']['cumulative'],
                'integral': data['distribution']['integral'],
                'metadata': data['metadata'],
                'uqff_validation': data['uqff_validation'],
                'schmidt_2016': data['schmidt_2016_validation']
            }
            
            return self.harvard_data
            
        except FileNotFoundError:
            print(f"Warning: Harvard distribution file not found at {json_path}")
            return None
        except json.JSONDecodeError as e:
            print(f"Warning: Error parsing Harvard distribution JSON: {e}")
            return None
    
    def get_peak_interpretation(self, log_M_sun: float, tolerance: float = 0.3) -> str:
        """
        Get UQFF interpretation for a given log(M/M_sun) value.
        
        Args:
            log_M_sun: Logarithmic mass in solar masses
            tolerance: Tolerance for peak matching (default 0.3 dex)
            
        Returns:
            UQFF interpretation string or 'No peak match'
        """
        if self.harvard_data is None:
            self.load_harvard_distribution()
        
        if self.harvard_data is None:
            return 'Harvard data not loaded'
        
        for peak_log_M, interpretation in self.harvard_data['peaks'].items():
            if abs(log_M_sun - peak_log_M) < tolerance:
                return interpretation
        
        return 'No peak match'
    
    def compute(self, dataset: dict) -> dict:
        rho_bar = dataset.get('rho_bar', COSMO['rho_c'] * COSMO['Omega_m'])
        M_BH = dataset.get('M_BH', 1e8 * CONST['M_sun'])
        z = dataset.get('z', 0.0)
        sigma = dataset.get('sigma', 1.0)
        delta_c = COSMO['delta_c'] * (1 + z)
        from math import erfc, log10
        arg = delta_c / (math.sqrt(2) * sigma)
        N_cumulative = rho_bar / M_BH * erfc(arg) if M_BH > 0 else 0
        
        # Add Harvard distribution interpretation
        log_M_sun = log10(M_BH / CONST['M_sun']) if M_BH > 0 else 0
        peak_interpretation = self.get_peak_interpretation(log_M_sun)
        
        result = {
            'N_cumulative': N_cumulative, 'delta_c_z': delta_c,
            'erfc_arg': arg, 'M_BH': M_BH, 'log_M_sun': log_M_sun,
            'equation': 'N(>M,z) = ρ̄ ∫ dM\'/M\'² erfc(δ_c/√2σ)',
            'source': 'Grok URL Eq. 8: EPS BH Mass Function',
            'uqff_interpretation': peak_interpretation
        }
        
        # Add Harvard data if loaded
        if self.harvard_data is not None:
            result['harvard_peaks'] = self.harvard_data['peaks']
            result['harvard_integral'] = self.harvard_data['integral']
        
        return result

class EddingtonAccretionCalculator:
    """Eq. 9: Ṁ_BH = 4πGM_BH m_p / (ε_r σ_T c) (Eddington accretion rate)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']; m_p = CONST['m_p']; sigma_T = CONST['sigma_T']
        M_BH = dataset.get('M_BH', 1e8 * CONST['M_sun'])
        epsilon_r = dataset.get('epsilon_r', 0.1)
        M_dot = 4 * CONST['pi'] * G * M_BH * m_p / (epsilon_r * sigma_T * c) if epsilon_r > 0 else 0
        L_Edd = 4 * CONST['pi'] * G * M_BH * m_p * c / sigma_T
        t_Sal = M_BH * epsilon_r * sigma_T * c / (4 * CONST['pi'] * G * m_p * M_BH) if M_BH > 0 else 0
        return {
            'M_dot_kg_s': M_dot, 'M_dot_Msun_yr': M_dot * 3.156e7 / CONST['M_sun'],
            'L_Edd_W': L_Edd, 't_Salpeter_yr': t_Sal / 3.156e7,
            'equation': 'Ṁ_BH = 4πGM_BH m_p / (ε_r σ_T c)',
            'source': 'Grok URL Eq. 9: Eddington Accretion Rate'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 10-11: SUPERNOVA REMNANTS
# ═══════════════════════════════════════════════════════════════════════════════

class SedovTaylorExpansionCalculator:
    """Eq. 10: R(t) = (Et²/ρ)^(1/5) (Sedov-Taylor blast wave)."""
    def compute(self, dataset: dict) -> dict:
        E = dataset.get('E', 1e44)  # J (~10⁵¹ erg)
        t = dataset.get('t', 1e10)  # s (~300 yr)
        rho = dataset.get('rho', 1e-21)  # kg/m³
        R = (E * t**2 / rho)**(1.0/5.0) if rho > 0 else 0
        v_expansion = (2.0/5.0) * R / t if t > 0 else 0
        return {
            'R_m': R, 'R_pc': R / 3.086e16, 'v_expansion': v_expansion,
            'E': E, 't': t,
            'equation': 'R(t) = (Et²/ρ)^(1/5)',
            'source': 'Grok URL Eq. 10: Sedov-Taylor Expansion'
        }

class DSAParticleAccelerationCalculator:
    """Eq. 11: dp/dt = (4/3)(u_s²/r_d)p (Diffusive Shock Acceleration)."""
    def compute(self, dataset: dict) -> dict:
        u_s = dataset.get('u_s', 3e6)  # m/s shock speed
        r_d = dataset.get('r_d', 1e16)  # m diffusion length
        p = dataset.get('p', 1e-20)  # kg·m/s momentum
        dpdt = (4.0/3.0) * (u_s**2 / r_d) * p if r_d > 0 else 0
        t_acc = r_d / u_s if u_s > 0 else 0  # acceleration timescale
        return {
            'dp_dt': dpdt, 't_acc_s': t_acc, 'u_s': u_s,
            'equation': 'dp/dt = (4/3)(u_s²/r_d)p',
            'source': 'Grok URL Eq. 11: DSA Particle Acceleration'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 12-13: GRAVITATIONAL WAVES
# ═══════════════════════════════════════════════════════════════════════════════

class ChirpMassFormulaCalculator:
    """Eq. 12: Chirp mass M = (m₁m₂)^(3/5)/(m₁+m₂)^(1/5)."""
    def compute(self, dataset: dict) -> dict:
        m1 = dataset.get('m1', 30.0 * CONST['M_sun'])
        m2 = dataset.get('m2', 25.0 * CONST['M_sun'])
        M_chirp = (m1 * m2)**(3.0/5.0) / (m1 + m2)**(1.0/5.0) if (m1 + m2) > 0 else 0
        eta = m1 * m2 / (m1 + m2)**2 if (m1 + m2) > 0 else 0
        return {
            'M_chirp_kg': M_chirp, 'M_chirp_Msun': M_chirp / CONST['M_sun'],
            'eta': eta, 'm1': m1, 'm2': m2,
            'equation': 'M = (m₁m₂)^(3/5)/(m₁+m₂)^(1/5)',
            'source': 'Grok URL Eq. 12: Chirp Mass Formula'
        }

class QNMRingdownCalculator:
    """Eq. 13: f_QNM quasinormal mode frequency."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        M_f = dataset.get('M_f', 50.0 * CONST['M_sun'])
        a_f = dataset.get('a_f', 0.7)  # dimensionless spin
        f_QNM = (c**3 / (2 * CONST['pi'] * G * M_f)) * (0.3737 + 0.088 * a_f) if M_f > 0 else 0
        tau_QNM = 2 * G * M_f / (c**3 * 0.0889 * (1 - a_f)**0.45) if (1 - a_f) > 0 and M_f > 0 else 0
        return {
            'f_QNM_Hz': f_QNM, 'tau_QNM_s': tau_QNM,
            'M_f_Msun': M_f / CONST['M_sun'], 'a_f': a_f,
            'equation': 'f_QNM = (c³/2πGM_f)(0.3737 + 0.088a_f)',
            'source': 'Grok URL Eq. 13: QNM Ringdown Frequency'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 14-15: QUASAR JETS
# ═══════════════════════════════════════════════════════════════════════════════

class BlandfordZnajekPowerCalculator:
    """Eq. 14: P_BZ = (1/32)B²R_H⁴Ω_H²c (Blandford-Znajek jet power)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        M = dataset.get('M', 1e9 * CONST['M_sun'])
        a = dataset.get('a_spin', 0.9)  # dimensionless spin
        B = dataset.get('B', 1e4 * 1e-4)  # G → T
        R_H = dpm_emergent_ug1(M, c)  # horizon radius
        Omega_H = a * c / (2 * R_H) if R_H > 0 else 0
        P_BZ = (1.0/32.0) * B**2 * R_H**4 * Omega_H**2 * c
        return {
            'P_BZ_W': P_BZ, 'P_BZ_erg_s': P_BZ * 1e7,
            'R_H': R_H, 'Omega_H': Omega_H,
            'equation': 'P_BZ = (1/32)B²R_H⁴Ω_H²c',
            'source': 'Grok URL Eq. 14: Blandford-Znajek Jet Power'
        }

class RelativisticJetVelocityCalculator:
    """Eq. 15: Γ(z) relativistic jet Lorentz factor profile."""
    def compute(self, dataset: dict) -> dict:
        Gamma_0 = dataset.get('Gamma_0', 1.5)
        Gamma_inf = dataset.get('Gamma_inf', 15.0)
        z_acc = dataset.get('z_acc', 1e19)  # m
        z = dataset.get('z_distance', 1e19)
        Gamma = Gamma_0 + (z / z_acc) * (Gamma_inf - Gamma_0) * (1 - math.exp(-z / z_acc)) if z_acc > 0 else Gamma_0
        beta = math.sqrt(1 - 1.0/Gamma**2) if Gamma > 1 else 0
        return {
            'Gamma': Gamma, 'beta': beta, 'v_fraction_c': beta,
            'equation': 'Γ(z) = Γ₀ + (z/z_acc)(Γ_∞-Γ₀)(1-e^(-z/z_acc))',
            'source': 'Grok URL Eq. 15: Relativistic Jet Lorentz Factor'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 16-18: NEUTRON STARS
# ═══════════════════════════════════════════════════════════════════════════════

class TOVEquationCalculator:
    """Eq. 16: TOV hydrostatic equilibrium for neutron stars."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        P = dataset.get('P', 1e34)  # Pa central pressure
        rho = dataset.get('rho', 1e18)  # kg/m³ central density
        r = dataset.get('r', 1e4)  # m (10 km)
        m_r = dataset.get('m_enclosed', 1.4 * CONST['M_sun'])
        factor1 = 1 + P / (rho * c**2) if rho > 0 else 1
        factor2 = 1 + 4 * CONST['pi'] * r**3 * P / (m_r * c**2) if m_r > 0 else 1
        factor3 = 1.0 / (1 - 2 * G * m_r / (r * c**2)) if r > 0 and (1 - 2*G*m_r/(r*c**2)) > 0 else 1
        dPdr = -dpm_emergent_ug1(m_r, r) * rho * factor1 * factor2 * factor3 if r > 0 else 0
        return {
            'dP_dr': dPdr, 'P': P, 'rho': rho, 'r': r,
            'GR_correction_1': factor1, 'GR_correction_2': factor2, 'GR_correction_3': factor3,
            'equation': 'dP/dr = -Gm(r)ρ/r² × (1+P/(ρc²))(1+4πr³P/(mc²))(1-2Gm/(rc²))⁻¹',
            'source': 'Grok URL Eq. 16: TOV Equation'
        }

class PulsarSpinDownAgeCalculator:
    """Eq. 17: τ = P/(2Ṗ) (Pulsar characteristic age)."""
    def compute(self, dataset: dict) -> dict:
        P = dataset.get('P', 0.033)  # s (Crab-like)
        P_dot = dataset.get('P_dot', 4.21e-13)  # s/s
        tau = P / (2 * P_dot) if P_dot > 0 else float('inf')
        B_surface = 3.2e19 * math.sqrt(P * P_dot)  # T
        E_dot = 4 * CONST['pi']**2 * 1e45 * P_dot / P**3 if P > 0 else 0  # spin-down luminosity
        return {
            'tau_s': tau, 'tau_yr': tau / 3.156e7, 'tau_kyr': tau / 3.156e10,
            'B_surface_T': B_surface, 'E_dot_W': E_dot,
            'equation': 'τ = P/(2Ṗ)',
            'source': 'Grok URL Eq. 17: Pulsar Spin-Down Age'
        }

class GlitchRecoveryCalculator:
    """Eq. 18: Δν = (I_s/I)ν₀(1-e^(-t/τ_q)) (Pulsar glitch recovery)."""
    def compute(self, dataset: dict) -> dict:
        I_s_frac = dataset.get('I_s_over_I', 0.01)  # superfluid fraction
        nu_0 = dataset.get('nu_0', 30.0)  # Hz spin frequency
        tau_q = dataset.get('tau_q', 86400.0)  # s quench time
        t = dataset.get('t', 86400.0)  # s
        Delta_nu = I_s_frac * nu_0 * (1 - math.exp(-t / tau_q)) if tau_q > 0 else 0
        Q = 1 - math.exp(-t / tau_q) if tau_q > 0 else 1  # healing factor
        return {
            'Delta_nu': Delta_nu, 'Q_factor': Q,
            'equation': 'Δν = (I_s/I)ν₀(1-e^(-t/τ_q))',
            'source': 'Grok URL Eq. 18: Glitch Recovery'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 19-20: GAMMA-RAY BURSTS
# ═══════════════════════════════════════════════════════════════════════════════

class FireballExpansionCalculator:
    """Eq. 19: Γ(r) fireball Lorentz factor evolution."""
    def compute(self, dataset: dict) -> dict:
        R_0 = dataset.get('R_0', 1e7)  # m initial radius
        eta = dataset.get('eta', 300.0)  # baryon loading
        r = dataset.get('r', 1e12)  # m
        R_s = eta * R_0  # saturation radius
        if r < R_s:
            Gamma = r / R_0
        else:
            Gamma = eta
        beta = math.sqrt(1 - 1.0/Gamma**2) if Gamma > 1 else 0
        return {
            'Gamma': Gamma, 'beta': beta, 'R_s': R_s,
            'phase': 'acceleration' if r < R_s else 'coasting',
            'equation': 'Γ(r) = r/R₀ (r<R_s); Γ=η (r>R_s)',
            'source': 'Grok URL Eq. 19: Fireball Expansion'
        }

class AfterglowSynchrotronCalculator:
    """Eq. 20: F_ν ∝ ν^(-(p-1)/2) t^(-3(p-1)/4) (GRB afterglow)."""
    def compute(self, dataset: dict) -> dict:
        p = dataset.get('p', 2.3)  # electron power-law index
        nu = dataset.get('nu', 1e14)  # Hz
        t = dataset.get('t', 86400.0)  # s (1 day)
        nu_ref = dataset.get('nu_ref', 1e14)
        t_ref = dataset.get('t_ref', 86400.0)
        F_0 = dataset.get('F_0', 1e-26)  # W/m²/Hz reference flux
        alpha_nu = -(p - 1) / 2
        alpha_t = -3 * (p - 1) / 4
        F_nu = F_0 * (nu/nu_ref)**alpha_nu * (t/t_ref)**alpha_t if nu_ref > 0 and t_ref > 0 else 0
        return {
            'F_nu': F_nu, 'alpha_nu': alpha_nu, 'alpha_t': alpha_t,
            'equation': 'F_ν ∝ ν^(-(p-1)/2) t^(-3(p-1)/4)',
            'source': 'Grok URL Eq. 20: Afterglow Synchrotron'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 21-22: CMB
# ═══════════════════════════════════════════════════════════════════════════════

class CMBAngularPowerCalculator:
    """Eq. 21: C_ℓ angular power spectrum (simplified)."""
    def compute(self, dataset: dict) -> dict:
        ell = dataset.get('ell', 200)
        A_s = dataset.get('A_s', 2.1e-9)
        n_s = dataset.get('n_s', 0.965)
        k_pivot = dataset.get('k_pivot', 0.05)  # Mpc⁻¹
        # Simplified Sachs-Wolfe plateau + acoustic peaks
        C_ell = A_s * (ell / 200.0)**(n_s - 1) / (ell * (ell + 1)) * 4 * CONST['pi']
        D_ell = ell * (ell + 1) * C_ell / (2 * CONST['pi'])  # Standard D_ℓ
        return {
            'C_ell': C_ell, 'D_ell': D_ell, 'ell': ell,
            'equation': 'C_ℓ = (4π/(2ℓ+1))∫P(k)|Δ_ℓ^T(k)|²dk/k',
            'source': 'Grok URL Eq. 21: CMB Angular Power Spectrum'
        }

class OpticalDepthCalculator:
    """Eq. 22: τ(z) = ∫n_e σ_T c (dt/dz') dz' (CMB optical depth)."""
    def compute(self, dataset: dict) -> dict:
        z = dataset.get('z_reion', 8.0)
        n_e0 = dataset.get('n_e0', 0.25)  # m⁻³ today
        sigma_T = CONST['sigma_T']; c = CONST['c']
        H_0_si = COSMO['H_0'] * 1e3 / 3.086e22  # Hz
        # Simplified: τ ≈ n_e0 σ_T c / H_0 × 2/3 × [(1+z)^(3/2) - 1]
        tau = n_e0 * sigma_T * c / H_0_si * (2.0/3.0) * ((1 + z)**1.5 - 1) if H_0_si > 0 else 0
        return {
            'tau': tau, 'z_reion': z,
            'equation': 'τ(z) = ∫n_e σ_T c (dt/dz\') dz\'',
            'source': 'Grok URL Eq. 22: Optical Depth to Reionization'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 23-25: AGN FEEDBACK
# ═══════════════════════════════════════════════════════════════════════════════

class MomentumFeedbackCalculator:
    """Eq. 23: p_term = Ṁ_out v_out (AGN momentum feedback)."""
    def compute(self, dataset: dict) -> dict:
        L_AGN = dataset.get('L_AGN', 1e39)  # W
        c = CONST['c']
        v_out = dataset.get('v_out', 1000e3)  # m/s
        f_v = dataset.get('f_coupling', 20.0)  # momentum boost
        M_dot_out = f_v * L_AGN / (c * v_out) if v_out > 0 else 0
        p_term = M_dot_out * v_out
        return {
            'p_term': p_term, 'M_dot_out': M_dot_out,
            'M_dot_out_Msun_yr': M_dot_out * 3.156e7 / CONST['M_sun'],
            'equation': 'p_term = Ṁ_out × v_out = f(v)L_AGN/c',
            'source': 'Grok URL Eq. 23: Momentum Feedback'
        }

class BZJetPowerUpdatedCalculator:
    """Eq. 24: P_jet = (κ/16π)Φ_BH²Ω_BH²/c (updated BZ with magnetic flux)."""
    def compute(self, dataset: dict) -> dict:
        c = CONST['c']
        kappa = dataset.get('kappa', 0.5)
        Phi_BH = dataset.get('Phi_BH', 1e25)  # Wb magnetic flux
        Omega_BH = dataset.get('Omega_BH', 1e-4)  # rad/s
        P_jet = (kappa / (16 * CONST['pi'])) * Phi_BH**2 * Omega_BH**2 / c
        return {
            'P_jet_W': P_jet, 'P_jet_erg_s': P_jet * 1e7,
            'equation': 'P_jet = (κ/16π)Φ_BH²Ω_BH²/c',
            'source': 'Grok URL Eq. 24: Updated BZ Jet Power'
        }

class FeedbackDutyCycleCalculator:
    """Eq. 25: f_duty(t) = (1-exp(-t/τ_cool))(1+Ṁ_acc/Ṁ_Edd)^(-1)."""
    def compute(self, dataset: dict) -> dict:
        tau_cool = dataset.get('tau_cool', 1e8 * 3.156e7)  # yr → s
        M_dot_acc = dataset.get('M_dot_acc', 0.1)
        M_dot_Edd = dataset.get('M_dot_Edd', 1.0)
        t = dataset.get('t', 1e8 * 3.156e7)
        f_duty = (1 - math.exp(-t / tau_cool)) * 1.0 / (1 + M_dot_acc / M_dot_Edd) if tau_cool > 0 and M_dot_Edd > 0 else 0
        return {
            'f_duty': f_duty, 'tau_cool_yr': tau_cool / 3.156e7,
            'equation': 'f_duty = (1-exp(-t/τ_cool))(1+Ṁ_acc/Ṁ_Edd)⁻¹',
            'source': 'Grok URL Eq. 25: Feedback Duty Cycle'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 26-28: EXOPLANETS
# ═══════════════════════════════════════════════════════════════════════════════

class PhotoevaporationRateCalculator:
    """Eq. 26: Ṁ = ε F_X R_p³ / (GM_p/R_p) (Photoevaporation)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        epsilon = dataset.get('epsilon', 0.15)
        F_X = dataset.get('F_X', 1e3)  # W/m² X-ray flux
        R_p = dataset.get('R_p', 6.371e6 * 3.0)  # m (3 R_Earth)
        M_p = dataset.get('M_p', 5.972e24 * 10.0)  # kg (10 M_Earth)
        v_esc_sq = G * M_p / R_p if R_p > 0 else 1
        M_dot = epsilon * F_X * R_p**3 / v_esc_sq if v_esc_sq > 0 else 0
        return {
            'M_dot_kg_s': M_dot, 'M_dot_g_s': M_dot * 1e3,
            'v_esc': math.sqrt(2 * v_esc_sq),
            'equation': 'Ṁ = ε F_X R_p³/(GM_p/R_p)',
            'source': 'Grok URL Eq. 26: Photoevaporation Rate'
        }

class TypeIMigrationTorqueCalculator:
    """Eq. 27: Type-I migration torque Γ."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        M_p = dataset.get('M_p', 5.972e24 * 10.0)
        M_star = dataset.get('M_star', CONST['M_sun'])
        Sigma = dataset.get('Sigma', 1000.0)  # kg/m² surface density
        r_p = dataset.get('r_p', 1.496e11)  # m (1 AU)
        H_over_r = dataset.get('H_over_r', 0.05)
        Omega = math.sqrt(G * M_star / r_p**3) if r_p > 0 else 0
        f_coeff = dataset.get('f_coeff', -2.5)
        Gamma = f_coeff * (G * M_p)**2 * Sigma * Omega**2 * r_p**4 / (H_over_r**3) if H_over_r > 0 else 0
        t_mig = M_p * r_p**2 * Omega / abs(Gamma) if abs(Gamma) > 0 else float('inf')
        return {
            'Gamma_torque': Gamma, 't_migration_yr': t_mig / 3.156e7,
            'equation': 'Γ = f(GM_p)²Σ Ω² r_p⁴(H/r)⁻³',
            'source': 'Grok URL Eq. 27: Type-I Migration Torque'
        }

class RadialVelocitySemiAmplitudeCalculator:
    """Eq. 28: K = 2πG/P × M_p sin i / (M_*+M_p)^(2/3) / √(1-e²)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        P = dataset.get('P', 365.25 * 86400)  # s
        M_p = dataset.get('M_p', 1.898e27)  # kg (Jupiter)
        M_star = dataset.get('M_star', CONST['M_sun'])
        inc = dataset.get('i', CONST['pi']/2)  # rad
        e = dataset.get('e', 0.0)
        K = (2 * CONST['pi'] * G / P)**(1.0/3.0) * M_p * math.sin(inc) / (M_star + M_p)**(2.0/3.0) / math.sqrt(1 - e**2) if P > 0 and (1 - e**2) > 0 else 0
        return {
            'K_m_s': K, 'K_km_s': K / 1e3,
            'equation': 'K = (2πG/P)^(1/3) M_p sin i / (M_*+M_p)^(2/3) / √(1-e²)',
            'source': 'Grok URL Eq. 28: Radial Velocity Semi-Amplitude'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 29-31: DARK MATTER HALOS
# ═══════════════════════════════════════════════════════════════════════════════

class NFWDensityProfileCalculator:
    """Eq. 29: ρ(r) = ρ_s / ((r/r_s)(1+r/r_s)²) (NFW profile)."""
    def compute(self, dataset: dict) -> dict:
        rho_s = dataset.get('rho_s', 1e7 * CONST['M_sun'] / (3.086e16)**3)
        r_s = dataset.get('r_s', 20.0 * 3.086e19)  # kpc → m
        r = dataset.get('r', 8.0 * 3.086e19)  # kpc → m
        x = r / r_s if r_s > 0 else 1
        rho = rho_s / (x * (1 + x)**2) if x > 0 else 0
        M_enclosed = 4 * CONST['pi'] * rho_s * r_s**3 * (math.log(1 + x) - x/(1 + x)) if x > 0 else 0
        return {
            'rho': rho, 'rho_s': rho_s, 'x': x,
            'M_enclosed': M_enclosed, 'M_enclosed_Msun': M_enclosed / CONST['M_sun'],
            'equation': 'ρ(r) = ρ_s / ((r/r_s)(1+r/r_s)²)',
            'source': 'Grok URL Eq. 29: NFW Density Profile'
        }

class RotationCurveCalculator:
    """Eq. 30: v(r)² = GM(r)/r (Rotation curve from enclosed mass)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        M_enclosed = dataset.get('M_enclosed', 1e11 * CONST['M_sun'])
        r = dataset.get('r', 8.0 * 3.086e19)  # kpc → m
        v_sq = G * M_enclosed / r if r > 0 else 0
        v_circ = math.sqrt(v_sq) if v_sq > 0 else 0
        return {
            'v_circ': v_circ, 'v_circ_km_s': v_circ / 1e3,
            'v_sq': v_sq, 'M_enclosed': M_enclosed,
            'equation': 'v(r)² = GM(r)/r',
            'source': 'Grok URL Eq. 30: Rotation Curve'
        }

class SIDMCoreFormationCalculator:
    """Eq. 31: t_core ≈ 1/(ρ σ/m) (SIDM core formation timescale)."""
    def compute(self, dataset: dict) -> dict:
        rho = dataset.get('rho', 1e8 * CONST['M_sun'] / (3.086e19)**3)  # M_sun/kpc³
        sigma_over_m = dataset.get('sigma_over_m', 1.0)  # cm²/g → m²/kg: ×0.1
        sigma_m_si = sigma_over_m * 1e-4  # cm²/g → m²/kg
        t_core = 1.0 / (rho * sigma_m_si) if rho > 0 and sigma_m_si > 0 else float('inf')
        return {
            't_core_s': t_core, 't_core_Gyr': t_core / 3.156e16,
            'equation': 't_core ≈ 1/(ρ × σ/m)',
            'source': 'Grok URL Eq. 31: SIDM Core Formation'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 32-34: GALAXY CLUSTERS
# ═══════════════════════════════════════════════════════════════════════════════

class StrongLensingMassCalculator:
    """Eq. 32: θ_E Einstein radius for strong lensing."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        M = dataset.get('M', 1e14 * CONST['M_sun'])
        D_L = dataset.get('D_L', 1e9 * 3.086e16)  # pc → m
        D_S = dataset.get('D_S', 3e9 * 3.086e16)
        D_LS = dataset.get('D_LS', 2e9 * 3.086e16)
        theta_E = math.sqrt(4 * dpm_emergent_ug1(M, c) * D_LS / (D_L * D_S)) if D_L > 0 and D_S > 0 else 0
        return {
            'theta_E_rad': theta_E, 'theta_E_arcsec': theta_E * 206265,
            'M': M, 'M_Msun': M / CONST['M_sun'],
            'equation': 'θ_E = √(4GM/c² × D_LS/(D_L D_S))',
            'source': 'Grok URL Eq. 32: Strong Lensing Einstein Radius'
        }

class XRayMassEstimateCalculator:
    """Eq. 33: M(<r) from X-ray hydrostatic equilibrium."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; k_B = CONST['k_B']
        T = dataset.get('T', 5e7)  # K
        r = dataset.get('r', 1e6 * 3.086e16)  # pc → m
        mu = dataset.get('mu', 0.6)
        m_H = CONST['m_p']
        d_ln_rho = dataset.get('d_ln_rho_d_ln_r', -1.5)
        d_ln_T = dataset.get('d_ln_T_d_ln_r', -0.3)
        M_r = -k_B * T * r / (G * mu * m_H) * (d_ln_rho + d_ln_T)
        return {
            'M_enclosed': M_r, 'M_enclosed_Msun': M_r / CONST['M_sun'],
            'T_keV': k_B * T / CONST['e_charge'] / 1e3,
            'equation': 'M(<r) = -kTr/(Gμm_H)(d ln ρ/d ln r + d ln T/d ln r)',
            'source': 'Grok URL Eq. 33: X-ray Mass Estimate'
        }

class MergerShockMachCalculator:
    """Eq. 34: Merger shock Mach number from density jump."""
    def compute(self, dataset: dict) -> dict:
        gamma = dataset.get('gamma', 5.0/3.0)
        rho_ratio = dataset.get('rho_2_over_rho_1', 3.0)
        M_sq = ((gamma + 1) * rho_ratio + (gamma - 1)) / (2 * gamma) if gamma > 0 else 0
        M = math.sqrt(M_sq) if M_sq > 0 else 0
        return {
            'Mach': M, 'Mach_sq': M_sq, 'rho_ratio': rho_ratio,
            'equation': 'M = √(((γ+1)(ρ₂/ρ₁)+(γ-1))/(2γ))',
            'source': 'Grok URL Eq. 34: Merger Shock Mach Number'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 35-36: COSMIC VOIDS
# ═══════════════════════════════════════════════════════════════════════════════

class VoidDensityEvolutionCalculator:
    """Eq. 35: δ_v(a) void density contrast evolution."""
    def compute(self, dataset: dict) -> dict:
        delta_v0 = dataset.get('delta_v0', -0.8)
        a = dataset.get('a', 1.0)
        Omega_m = dataset.get('Omega_m', COSMO['Omega_m'])
        Omega_L = dataset.get('Omega_Lambda', COSMO['Omega_Lambda'])
        denom = (Omega_m * a + Omega_L)**(3.0/2.0) if (Omega_m * a + Omega_L) > 0 else 1
        delta_v = -(3.0/5.0) * delta_v0 / denom
        return {
            'delta_v': delta_v, 'a': a,
            'equation': 'δ_v(a) = -(3/5)(Ω_m a + Ω_Λ)^(-3/2) δ_v0',
            'source': 'Grok URL Eq. 35: Void Density Evolution'
        }

class OutflowVelocityCalculator:
    """Eq. 36: v_pec peculiar velocity from void."""
    def compute(self, dataset: dict) -> dict:
        f = dataset.get('f_growth', 0.5)
        H = dataset.get('H', COSMO['H_0'] * 1e3 / 3.086e22)  # Hz
        delta = dataset.get('delta', -0.5)
        r = dataset.get('r', 30e6 * 3.086e16)  # pc → m
        v_pec = -f * H * delta * r / 3
        return {
            'v_pec': v_pec, 'v_pec_km_s': v_pec / 1e3,
            'equation': 'v_pec = -fH/3 × δ × r',
            'source': 'Grok URL Eq. 36: Void Outflow Velocity'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 37-38: REIONIZATION
# ═══════════════════════════════════════════════════════════════════════════════

class IonizationFractionCalculator:
    """Eq. 37: dx_e/dt reionization ionization fraction evolution."""
    def compute(self, dataset: dict) -> dict:
        n_dot_gamma = dataset.get('n_dot_gamma', 1e51)  # ionizing photon rate /s/Mpc³
        epsilon_esc = dataset.get('epsilon_esc', 0.1)
        f_star = dataset.get('f_star', 0.1)
        alpha_B = dataset.get('alpha_B', 2.6e-19)  # m³/s
        n_e = dataset.get('n_e', 0.25)  # m⁻³
        C = dataset.get('C_clumping', 3.0)
        dx_e_dt = n_dot_gamma * epsilon_esc * f_star - alpha_B * n_e**2 * C
        return {
            'dx_e_dt': dx_e_dt,
            'ionization_term': n_dot_gamma * epsilon_esc * f_star,
            'recombination_term': alpha_B * n_e**2 * C,
            'equation': 'dx_e/dt = ṅ_γ ε_esc f_* - α_B n_e² C',
            'source': 'Grok URL Eq. 37: Ionization Fraction'
        }

class BubbleRadiusCalculator:
    """Eq. 38: R_b(t) = (3Ṅ_γ t / (4π n_H))^(1/3) (Strömgren bubble)."""
    def compute(self, dataset: dict) -> dict:
        N_dot_gamma = dataset.get('N_dot_gamma', 1e49)  # /s
        t = dataset.get('t', 1e14)  # s
        n_H = dataset.get('n_H', 0.25)  # m⁻³
        R_b = (3 * N_dot_gamma * t / (4 * CONST['pi'] * n_H))**(1.0/3.0) if n_H > 0 else 0
        return {
            'R_b_m': R_b, 'R_b_Mpc': R_b / 3.086e22,
            'equation': 'R_b(t) = (3Ṅ_γ t/(4πn_H))^(1/3)',
            'source': 'Grok URL Eq. 38: Reionization Bubble Radius'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 39-41: ISM
# ═══════════════════════════════════════════════════════════════════════════════

class JeansLengthCalculator:
    """Eq. 39: λ_J = √(πc_s²/(Gρ)) (Jeans length)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        c_s = dataset.get('c_s', 200.0)  # m/s
        rho = dataset.get('rho', 1e-18)  # kg/m³
        lambda_J = math.sqrt(CONST['pi'] * c_s**2 / (G * rho)) if rho > 0 else 0
        M_J = (CONST['pi'] / 6) * rho * lambda_J**3 if lambda_J > 0 else 0
        return {
            'lambda_J_m': lambda_J, 'lambda_J_pc': lambda_J / 3.086e16,
            'M_Jeans_kg': M_J, 'M_Jeans_Msun': M_J / CONST['M_sun'],
            'equation': 'λ_J = √(πc_s²/(Gρ))',
            'source': 'Grok URL Eq. 39: Jeans Length'
        }

class AlfvenVelocityCalculator:
    """Eq. 40: v_A = B/√(4πρ) (Alfvén velocity)."""
    def compute(self, dataset: dict) -> dict:
        B = dataset.get('B', 5e-10)  # T (~5 μG)
        rho = dataset.get('rho', 1e-21)  # kg/m³
        v_A = B / math.sqrt(4 * CONST['pi'] * rho) if rho > 0 else 0
        # Also compute using SI μ_0
        v_A_SI = B / math.sqrt(CONST['mu_0'] * rho) if rho > 0 else 0
        return {
            'v_A': v_A_SI, 'v_A_cgs': v_A, 'v_A_km_s': v_A_SI / 1e3,
            'equation': 'v_A = B/√(μ₀ρ)',
            'source': 'Grok URL Eq. 40: Alfvén Velocity'
        }

class TurbulentCascadeCalculator:
    """Eq. 41: ε = v_ℓ³/ℓ = const (Kolmogorov turbulent cascade)."""
    def compute(self, dataset: dict) -> dict:
        v_L = dataset.get('v_L', 10e3)  # m/s (outer scale velocity)
        L = dataset.get('L', 100 * 3.086e16)  # m (100 pc outer scale)
        ell = dataset.get('ell', 1 * 3.086e16)  # m (target scale)
        epsilon = v_L**3 / L if L > 0 else 0
        v_ell = (epsilon * ell)**(1.0/3.0) if ell > 0 else 0
        return {
            'epsilon': epsilon, 'v_ell': v_ell,
            'v_ell_km_s': v_ell / 1e3,
            'equation': 'ε = v_ℓ³/ℓ = const (Kolmogorov)',
            'source': 'Grok URL Eq. 41: Turbulent Cascade'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 42-44: STELLAR EVOLUTION
# ═══════════════════════════════════════════════════════════════════════════════

class MainSequenceLifetimeCalculator:
    """Eq. 42: τ_MS ≈ 10¹⁰(M/M_⊙)^(-2.5) yr."""
    def compute(self, dataset: dict) -> dict:
        M = dataset.get('M', CONST['M_sun'])
        M_ratio = M / CONST['M_sun']
        tau_yr = 1e10 * M_ratio**(-2.5) if M_ratio > 0 else float('inf')
        return {
            'tau_yr': tau_yr, 'tau_Gyr': tau_yr / 1e9,
            'M_Msun': M_ratio,
            'equation': 'τ_MS ≈ 10¹⁰(M/M_⊙)^(-2.5) yr',
            'source': 'Grok URL Eq. 42: Main Sequence Lifetime'
        }

class MassLuminosityRelationCalculator:
    """Eq. 43: L ∝ M^3.5 (Mass-luminosity relation)."""
    def compute(self, dataset: dict) -> dict:
        M = dataset.get('M', CONST['M_sun'])
        M_ratio = M / CONST['M_sun']
        L_sun = 3.828e26  # W
        if M_ratio < 0.43:
            alpha = 2.3
        elif M_ratio < 2:
            alpha = 4.0
        elif M_ratio < 55:
            alpha = 3.5
        else:
            alpha = 1.0
        L = L_sun * M_ratio**alpha
        return {
            'L_W': L, 'L_Lsun': L / L_sun, 'alpha': alpha,
            'M_Msun': M_ratio,
            'equation': 'L ∝ M^α (α=3.5 for 2-55 M_⊙)',
            'source': 'Grok URL Eq. 43: Mass-Luminosity Relation'
        }

class ConvectiveTurnoverCalculator:
    """Eq. 44: t_conv = H_p/v_conv (Convective turnover time)."""
    def compute(self, dataset: dict) -> dict:
        H_p = dataset.get('H_p', 1e8)  # m pressure scale height
        v_conv = dataset.get('v_conv', 100.0)  # m/s
        t_conv = H_p / v_conv if v_conv > 0 else float('inf')
        return {
            't_conv_s': t_conv, 't_conv_days': t_conv / 86400,
            'equation': 't_conv = H_p/v_conv',
            'source': 'Grok URL Eq. 44: Convective Turnover Time'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 45-46: BBN
# ═══════════════════════════════════════════════════════════════════════════════

class BaryonPhotonRatioCalculator:
    """Eq. 45: η = n_b/n_γ = 6×10⁻¹⁰ (Baryon-to-photon ratio)."""
    def compute(self, dataset: dict) -> dict:
        n_b = dataset.get('n_b', 0.25)  # m⁻³
        n_gamma = dataset.get('n_gamma', COSMO['n_gamma'])
        eta = n_b / n_gamma if n_gamma > 0 else 0
        Omega_b_h2 = eta / 2.74e-8  # Ω_b h²
        return {
            'eta': eta, 'Omega_b_h2': Omega_b_h2,
            'equation': 'η = n_b/n_γ ≈ 6×10⁻¹⁰',
            'source': 'Grok URL Eq. 45: Baryon-to-Photon Ratio'
        }

class DeuteriumBottleneckCalculator:
    """Eq. 46: t_D = (3/(32πGρ_rad))^(1/2) ≈ 180 s (D bottleneck)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        T_MeV = dataset.get('T_MeV', 0.1)
        # ρ_rad = π²/30 × g_* × (kT)⁴/(ℏc)³
        g_star = dataset.get('g_star', 10.75)
        kT = T_MeV * 1e6 * CONST['e_charge']  # MeV → J
        rho_rad = CONST['pi']**2 / 30 * g_star * kT**4 / (CONST['hbar'] * CONST['c'])**3 / CONST['c']**2
        t_D = math.sqrt(3.0 / (32 * CONST['pi'] * G * rho_rad)) if rho_rad > 0 else 0
        return {
            't_D_s': t_D, 'rho_rad': rho_rad, 'T_MeV': T_MeV,
            'equation': 't_D = √(3/(32πGρ_rad)) ≈ 180 s',
            'source': 'Grok URL Eq. 46: Deuterium Bottleneck'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 47-49: COSMOLOGY
# ═══════════════════════════════════════════════════════════════════════════════

class FirstFriedmannCalculator:
    """Eq. 47: (ȧ/a)² = 8πGρ/3 - kc²/a² + Λc²/3 (First Friedmann)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        rho = dataset.get('rho', COSMO['rho_c'])
        k = dataset.get('k', 0)  # curvature
        a = dataset.get('a', 1.0)
        Lambda = dataset.get('Lambda', COSMO['Lambda'])
        H_sq = 8 * CONST['pi'] * G * rho / 3 - k * c**2 / a**2 + Lambda * c**2 / 3 if a > 0 else 0
        H = math.sqrt(abs(H_sq))
        return {
            'H_sq': H_sq, 'H': H, 'H_km_s_Mpc': H * 3.086e22 / 1e3,
            'equation': '(ȧ/a)² = 8πGρ/3 - kc²/a² + Λc²/3',
            'source': 'Grok URL Eq. 47: First Friedmann Equation'
        }

class SecondFriedmannCalculator:
    """Eq. 48: ä/a = -4πG(ρ+3p/c²)/3 + Λc²/3 (Second Friedmann)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        rho = dataset.get('rho', COSMO['rho_c'])
        p = dataset.get('p', 0.0)  # pressure
        Lambda = dataset.get('Lambda', COSMO['Lambda'])
        a_ddot_over_a = -4 * CONST['pi'] * G * (rho + 3 * p / c**2) / 3 + Lambda * c**2 / 3
        deceleration_q = -(a_ddot_over_a) / (COSMO['H_0'] * 1e3 / 3.086e22)**2 if COSMO['H_0'] > 0 else 0
        return {
            'a_ddot_over_a': a_ddot_over_a, 'q_deceleration': deceleration_q,
            'accelerating': a_ddot_over_a > 0,
            'equation': 'ä/a = -4πG(ρ+3p/c²)/3 + Λc²/3',
            'source': 'Grok URL Eq. 48: Second Friedmann Equation'
        }

class DensityParameterCalculator:
    """Eq. 49: Ω(z) = 8πGρ(z)/(3H(z)²)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        rho = dataset.get('rho', COSMO['rho_c'])
        H = dataset.get('H', COSMO['H_0'] * 1e3 / 3.086e22)
        Omega = 8 * CONST['pi'] * G * rho / (3 * H**2) if H > 0 else 0
        rho_crit = 3 * H**2 / (8 * CONST['pi'] * G) if G > 0 else 0
        return {
            'Omega': Omega, 'rho_crit': rho_crit,
            'equation': 'Ω(z) = 8πGρ(z)/(3H(z)²)',
            'source': 'Grok URL Eq. 49: Density Parameter'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 50-52: INFLATION
# ═══════════════════════════════════════════════════════════════════════════════

class SlowRollParametersCalculator:
    """Eq. 50: ε = (1/2)(V'/V)², η = V''/V (Slow-roll parameters)."""
    def compute(self, dataset: dict) -> dict:
        V = dataset.get('V', 1e64)  # Inflaton potential (GeV⁴ in natural units)
        V_prime = dataset.get('V_prime', 1e62)  # dV/dφ
        V_double_prime = dataset.get('V_double_prime', 1e60)
        epsilon = 0.5 * (V_prime / V)**2 if V > 0 else 0
        eta = V_double_prime / V if V > 0 else 0
        n_s = 1 - 6 * epsilon + 2 * eta  # spectral index
        r = 16 * epsilon  # tensor-to-scalar ratio
        return {
            'epsilon': epsilon, 'eta': eta, 'n_s': n_s, 'r': r,
            'inflation_valid': epsilon < 1 and abs(eta) < 1,
            'equation': 'ε = (1/2)(V\'/V)²; η = V\'\'/V',
            'source': 'Grok URL Eq. 50: Slow-Roll Parameters'
        }

class CurvaturePowerSpectrumCalculator:
    """Eq. 51: P_R(k) = H²/(8π²εc_s³) (Curvature power spectrum)."""
    def compute(self, dataset: dict) -> dict:
        H = dataset.get('H_inflation', 1e13 * 1.602e-10 / CONST['hbar'])
        epsilon = dataset.get('epsilon', 0.01)
        c_s = dataset.get('c_s', 1.0)  # sound speed / c
        P_R = H**2 / (8 * CONST['pi']**2 * epsilon * c_s**3) if epsilon > 0 and c_s > 0 else 0
        A_s = P_R
        return {
            'P_R': P_R, 'A_s': A_s, 'ln_10_10_A_s': math.log(1e10 * A_s) if A_s > 0 else 0,
            'equation': 'P_R(k) = H²/(8π²εc_s³)',
            'source': 'Grok URL Eq. 51: Curvature Power Spectrum'
        }

class EFoldsCalculator:
    """Eq. 52: N = ∫H dt (Number of e-folds of inflation)."""
    def compute(self, dataset: dict) -> dict:
        H = dataset.get('H_inflation', 1e36)  # 1/s
        epsilon = dataset.get('epsilon', 0.01)
        delta_phi = dataset.get('delta_phi', 1.0)  # Planck units
        M_pl = CONST['M_pl']
        N = delta_phi / math.sqrt(2 * epsilon) if epsilon > 0 else 0
        return {
            'N_efolds': N, 'sufficient': N >= 50,
            'equation': 'N = ∫dφ/√(2ε) ≈ 50-60',
            'source': 'Grok URL Eq. 52: Number of e-Folds'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 53-54: GW BACKGROUNDS
# ═══════════════════════════════════════════════════════════════════════════════

class TensorPowerSpectrumCalculator:
    """Eq. 53: P_T(k) = 2H²/(π²M_pl²) (Primordial tensor spectrum)."""
    def compute(self, dataset: dict) -> dict:
        H = dataset.get('H_inflation', 1e13 * 1.602e-10 / CONST['hbar'])
        M_pl = CONST['M_pl']
        P_T = 2 * H**2 / (CONST['pi']**2 * M_pl**2 * CONST['c']**2) if M_pl > 0 else 0
        return {
            'P_T': P_T,
            'equation': 'P_T(k) = 2H²/(π²M_pl²c²)',
            'source': 'Grok URL Eq. 53: Tensor Power Spectrum'
        }

class StochasticGWDensityCalculator:
    """Eq. 54: Ω_GW(f) stochastic GW energy density."""
    def compute(self, dataset: dict) -> dict:
        f = dataset.get('f', 1e-8)  # Hz (nHz for PTA)
        H_0 = COSMO['H_0'] * 1e3 / 3.086e22  # Hz
        A = dataset.get('A_GW', 1e-15)  # amplitude
        alpha = dataset.get('alpha', -2.0/3.0)
        f_ref = dataset.get('f_ref', 1e-8)
        Omega_GW = (2 * CONST['pi']**2 / (3 * H_0**2)) * f**2 * A**2 * (f/f_ref)**(2*alpha) if H_0 > 0 and f_ref > 0 else 0
        return {
            'Omega_GW': Omega_GW, 'f': f, 'A_GW': A,
            'equation': 'Ω_GW(f) = (2π²/3H₀²)f²h_c²(f)',
            'source': 'Grok URL Eq. 54: Stochastic GW Energy Density'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 55-57: BINARY BH MERGERS
# ═══════════════════════════════════════════════════════════════════════════════

class InspiralFrequencyEvolutionCalculator:
    """Eq. 55: ḟ = (96/5)π^(8/3)(GM/c³)^(5/3)f^(11/3)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        M_chirp = dataset.get('M_chirp', 28.3 * CONST['M_sun'])
        f = dataset.get('f', 30.0)  # Hz
        f_dot = (96.0/5.0) * CONST['pi']**(8.0/3.0) * (G * M_chirp / c**3)**(5.0/3.0) * f**(11.0/3.0)
        return {
            'f_dot': f_dot, 'M_chirp_Msun': M_chirp / CONST['M_sun'],
            'equation': 'ḟ = (96/5)π^(8/3)(GM/c³)^(5/3)f^(11/3)',
            'source': 'Grok URL Eq. 55: Inspiral Frequency Evolution'
        }

class MergerTimeCalculator:
    """Eq. 56: t_merge = (5c⁵/(256G^(5/3)))(1/(πf_i)^(8/3))/M^(5/3)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        M_chirp = dataset.get('M_chirp', 28.3 * CONST['M_sun'])
        f_i = dataset.get('f_i', 10.0)  # Hz initial frequency
        t_merge = (5 * c**5 / (256 * G**(5.0/3.0))) / ((CONST['pi'] * f_i)**(8.0/3.0) * M_chirp**(5.0/3.0)) if f_i > 0 and M_chirp > 0 else 0
        return {
            't_merge_s': t_merge, 't_merge_yr': t_merge / 3.156e7,
            'equation': 't_merge = (5c⁵/256)(πf_i)^(-8/3)(GM)^(-5/3)',
            'source': 'Grok URL Eq. 56: Binary BH Merger Time'
        }

class RingdownDampingTimeCalculator:
    """Eq. 57: τ_ℓm = 2M_f c²/(Q_ℓm) (Ringdown damping time)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        M_f = dataset.get('M_f', 50.0 * CONST['M_sun'])
        Q = dataset.get('Q_lm', 10.0)
        a_f = dataset.get('a_f', 0.7)
        tau = 2 * G * M_f / (c**3 * (1 - a_f)**0.45 * 0.0889) if (1 - a_f) > 0 else 0
        return {
            'tau_ringdown_s': tau, 'Q_factor': Q,
            'equation': 'τ_ℓm ≈ 2GM_f/(c³ × f_decay)',
            'source': 'Grok URL Eq. 57: Ringdown Damping Time'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 58-60: SUPERNOVAE
# ═══════════════════════════════════════════════════════════════════════════════

class ArnettLawLightCurveCalculator:
    """Eq. 58: L_peak = M_Ni ε_Ni / t_d (Arnett's law)."""
    def compute(self, dataset: dict) -> dict:
        M_Ni = dataset.get('M_Ni', 0.6 * CONST['M_sun'])
        epsilon_Ni = dataset.get('epsilon_Ni', 3.9e10 * 1e-4)  # erg/g → J/kg
        t_d = dataset.get('t_d', 15 * 86400)  # days → s
        L_peak = M_Ni * epsilon_Ni / t_d if t_d > 0 else 0
        return {
            'L_peak_W': L_peak, 'L_peak_Lsun': L_peak / 3.828e26,
            'M_56_Msun': M_Ni / CONST['M_sun'],
            'equation': 'L_peak = M_Ni × ε_Ni / t_d',
            'source': 'Grok URL Eq. 58: Arnett\'s Law'
        }

class EjectaVelocityCalculator:
    """Eq. 59: v_ej = √(2E_kin/M_ej) (Supernova ejecta velocity)."""
    def compute(self, dataset: dict) -> dict:
        E_kin = dataset.get('E_kin', 1e44)  # J (~10⁵¹ erg)
        M_ej = dataset.get('M_ej', 5.0 * CONST['M_sun'])
        v_ej = math.sqrt(2 * E_kin / M_ej) if M_ej > 0 else 0
        return {
            'v_ej': v_ej, 'v_ej_km_s': v_ej / 1e3,
            'equation': 'v_ej = √(2E_kin/M_ej)',
            'source': 'Grok URL Eq. 59: Ejecta Velocity'
        }

class NucleosynthesisYieldCalculator:
    """Eq. 60: Y_i = ∫ρX_i dt (Nucleosynthesis yield)."""
    def compute(self, dataset: dict) -> dict:
        rho = dataset.get('rho', 1e9)  # kg/m³
        X_i = dataset.get('X_i', 0.01)  # mass fraction
        t = dataset.get('t', 1.0)  # s
        Y_i = rho * X_i * t
        return {
            'Y_i': Y_i, 'X_i': X_i,
            'equation': 'Y_i = ∫ρX_i dt',
            'source': 'Grok URL Eq. 60: Nucleosynthesis Yield'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 61-63: PLANETARY NEBULAE
# ═══════════════════════════════════════════════════════════════════════════════

class PNExpansionRadiusCalculator:
    """Eq. 61: R(t) = v_exp × t (PN expansion)."""
    def compute(self, dataset: dict) -> dict:
        v_exp = dataset.get('v_exp', 20e3)  # m/s (~20 km/s)
        t = dataset.get('t', 1e4 * 3.156e7)  # yr → s
        R = v_exp * t
        return {
            'R_m': R, 'R_pc': R / 3.086e16,
            'equation': 'R(t) = v_exp × t',
            'source': 'Grok URL Eq. 61: PN Expansion Radius'
        }

class IonizationFrontVelocityCalculator:
    """Eq. 62: v_IF = Ṅ_UV/(4πR²n) - α_B nR/3."""
    def compute(self, dataset: dict) -> dict:
        N_dot_UV = dataset.get('N_dot_UV', 1e47)  # /s
        R = dataset.get('R', 0.1 * 3.086e16)  # pc → m
        n = dataset.get('n', 1e4 * 1e6)  # cm⁻³ → m⁻³
        alpha_B = dataset.get('alpha_B', 2.6e-19)  # m³/s
        v_IF = N_dot_UV / (4 * CONST['pi'] * R**2 * n) - alpha_B * n * R / 3 if R > 0 and n > 0 else 0
        return {
            'v_IF': v_IF, 'v_IF_km_s': v_IF / 1e3,
            'equation': 'v_IF = Ṅ_UV/(4πR²n) - α_B nR/3',
            'source': 'Grok URL Eq. 62: Ionization Front Velocity'
        }

class AGBMassLossReimersCalculator:
    """Eq. 63: Ṁ = 4×10⁻¹³ LR/M M_⊙/yr (Reimers mass loss)."""
    def compute(self, dataset: dict) -> dict:
        L_Lsun = dataset.get('L_Lsun', 3000.0)
        R_Rsun = dataset.get('R_Rsun', 200.0)
        M_Msun = dataset.get('M_Msun', 1.5)
        eta_R = dataset.get('eta_Reimers', 1.0)
        M_dot = 4e-13 * eta_R * L_Lsun * R_Rsun / M_Msun if M_Msun > 0 else 0
        return {
            'M_dot_Msun_yr': M_dot, 'M_dot_kg_s': M_dot * CONST['M_sun'] / 3.156e7,
            'equation': 'Ṁ = 4×10⁻¹³ η LR/M (M_⊙/yr)',
            'source': 'Grok URL Eq. 63: AGB Reimers Mass Loss'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 64-69: CLUSTERS
# ═══════════════════════════════════════════════════════════════════════════════

class ShockMachNumberTempCalculator:
    """Eq. 64: Mach from temperature jump."""
    def compute(self, dataset: dict) -> dict:
        gamma = dataset.get('gamma', 5.0/3.0)
        T2_over_T1 = dataset.get('T2_over_T1', 3.0)
        M_sq = ((gamma + 1) / (gamma - 1) * T2_over_T1 - 2 / (gamma - 1)) if (gamma - 1) > 0 else 0
        M = math.sqrt(M_sq) if M_sq > 0 else 0
        return {
            'Mach': M, 'T2_over_T1': T2_over_T1,
            'equation': 'M = √((γ+1)/(γ-1) × T₂/T₁ - 2/(γ-1))',
            'source': 'Grok URL Eq. 64: Shock Mach from Temperature'
        }

class MergerTimescaleVirialCalculator:
    """Eq. 65: t_merge = r_vir/σ_v."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        r_vir = dataset.get('r_vir', 2e6 * 3.086e16)  # pc → m
        M = dataset.get('M', 1e15 * CONST['M_sun'])
        sigma_v = dataset.get('sigma_v', math.sqrt(G * M / (5 * r_vir)) * 3) if r_vir > 0 else 1e6
        t_merge = r_vir / sigma_v if sigma_v > 0 else float('inf')
        return {
            't_merge_s': t_merge, 't_merge_Gyr': t_merge / 3.156e16,
            'sigma_v_km_s': sigma_v / 1e3,
            'equation': 't_merge = r_vir/σ_v',
            'source': 'Grok URL Eq. 65: Virial Merger Timescale'
        }

class CoolCoreHeatingCalculator:
    """Eq. 66: Ė_heat = L_cool/f_duty."""
    def compute(self, dataset: dict) -> dict:
        L_cool = dataset.get('L_cool', 1e37)  # W
        f_duty = dataset.get('f_duty', 0.3)
        E_dot_heat = L_cool / f_duty if f_duty > 0 else float('inf')
        return {
            'E_dot_heat_W': E_dot_heat,
            'equation': 'Ė_heat = L_cool / f_duty',
            'source': 'Grok URL Eq. 66: Cool Core Heating Rate'
        }

class EvaporationRateCalculator:
    """Eq. 67: Ṅ = -N/t_evap (Cluster star evaporation)."""
    def compute(self, dataset: dict) -> dict:
        N = dataset.get('N', 1e5)
        t_relax = dataset.get('t_relax', 1e9 * 3.156e7)  # yr → s
        t_evap = 136 * t_relax / math.log(0.02 * N) if N > 50 else float('inf')
        N_dot = -N / t_evap if t_evap > 0 else 0
        return {
            'N_dot': N_dot, 't_evap_yr': t_evap / 3.156e7,
            'equation': 'Ṅ = -N/t_evap; t_evap = 136 t_relax/ln(0.02N)',
            'source': 'Grok URL Eq. 67: Cluster Evaporation Rate'
        }

class RelaxationTimeCalculator:
    """Eq. 68: t_relax = N/(8 ln N) × r/σ_v."""
    def compute(self, dataset: dict) -> dict:
        N = dataset.get('N', 1e5)
        r = dataset.get('r', 5 * 3.086e16)  # pc → m
        sigma_v = dataset.get('sigma_v', 5e3)  # m/s
        t_relax = N / (8 * math.log(N)) * r / sigma_v if N > 1 and sigma_v > 0 else 0
        return {
            't_relax_s': t_relax, 't_relax_Myr': t_relax / 3.156e13,
            'equation': 't_relax = N/(8 ln N) × r/σ_v',
            'source': 'Grok URL Eq. 68: Relaxation Time'
        }

class VirialMassCalculator:
    """Eq. 69: M_vir = 3σ²r_h/G (Virial mass)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        sigma = dataset.get('sigma_v', 1000e3)  # m/s
        r_h = dataset.get('r_h', 1e6 * 3.086e16)  # pc → m
        M_vir = 3 * sigma**2 * r_h / G
        return {
            'M_vir_kg': M_vir, 'M_vir_Msun': M_vir / CONST['M_sun'],
            'equation': 'M_vir = 3σ²r_h/G',
            'source': 'Grok URL Eq. 69: Virial Mass'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 70-72: QUASAR FEEDBACK
# ═══════════════════════════════════════════════════════════════════════════════

class WindTerminalVelocityCalculator:
    """Eq. 70: v_∞ = √(2GM(1-Γ)/r_launch)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        M = dataset.get('M', 1e9 * CONST['M_sun'])
        Gamma_Edd = dataset.get('Gamma_Edd', 0.9)
        r_launch = dataset.get('r_launch', 1e14)  # m
        v_inf = math.sqrt(2 * G * M * (1 - Gamma_Edd) / r_launch) if r_launch > 0 and Gamma_Edd < 1 else 0
        return {
            'v_inf': v_inf, 'v_inf_km_s': v_inf / 1e3,
            'v_inf_over_c': v_inf / CONST['c'],
            'equation': 'v_∞ = √(2GM(1-Γ)/r_launch)',
            'source': 'Grok URL Eq. 70: Wind Terminal Velocity'
        }

class IonizationParameterUCalculator:
    """Eq. 71: U = Q_H/(4πr²n_H c) (Ionization parameter)."""
    def compute(self, dataset: dict) -> dict:
        c = CONST['c']
        Q_H = dataset.get('Q_H', 1e56)  # /s ionizing photon rate
        r = dataset.get('r', 1e3 * 3.086e16)  # pc → m
        n_H = dataset.get('n_H', 1e6)  # m⁻³
        U = Q_H / (4 * CONST['pi'] * r**2 * n_H * c) if r > 0 and n_H > 0 else 0
        return {
            'U': U, 'log_U': math.log10(U) if U > 0 else -99,
            'equation': 'U = Q_H/(4πr²n_H c)',
            'source': 'Grok URL Eq. 71: Ionization Parameter'
        }

class EnergyCouplingEfficiencyCalculator:
    """Eq. 72: ε_f = Ė_kin/(Ṁ_acc c²)."""
    def compute(self, dataset: dict) -> dict:
        c = CONST['c']
        E_dot_kin = dataset.get('E_dot_kin', 1e38)  # W
        M_dot_acc = dataset.get('M_dot_acc', 1.0 * CONST['M_sun'] / 3.156e7)  # M_sun/yr → kg/s
        epsilon_f = E_dot_kin / (M_dot_acc * c**2) if M_dot_acc > 0 else 0
        return {
            'epsilon_f': epsilon_f,
            'equation': 'ε_f = Ė_kin/(Ṁ_acc c²)',
            'source': 'Grok URL Eq. 72: Energy Coupling Efficiency'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 73-75: BINARY PULSARS
# ═══════════════════════════════════════════════════════════════════════════════

class OrbitalDecayCalculator:
    """Eq. 73: Ṗ_b GW orbital decay (Peters formula)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        P_b = dataset.get('P_b', 7.75 * 3600)  # s (Hulse-Taylor)
        m1 = dataset.get('m1', 1.4 * CONST['M_sun'])
        m2 = dataset.get('m2', 1.4 * CONST['M_sun'])
        e = dataset.get('e', 0.617)
        M_total = m1 + m2
        f_e = (1 + 73.0/24.0 * e**2 + 37.0/96.0 * e**4) / (1 - e**2)**(7.0/2.0)
        P_b_dot = -(192 * CONST['pi'] / 5) * (P_b / (2 * CONST['pi']))**(-5.0/3.0) * (G * m1 * m2 / (c**3 * M_total**(1.0/3.0)))**(5.0/3.0) * f_e
        return {
            'P_b_dot': P_b_dot, 'P_b_dot_s_per_s': P_b_dot,
            'f_eccentricity': f_e,
            'equation': 'Ṗ_b = -(192π/5)(P_b/2π)^(-5/3)(Gm₁m₂/c³M^(1/3))^(5/3)f(e)',
            'source': 'Grok URL Eq. 73: GW Orbital Decay'
        }

class PeriastronAdvanceCalculator:
    """Eq. 74: ω̇ GR periastron advance."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        P_b = dataset.get('P_b', 7.75 * 3600)
        m1 = dataset.get('m1', 1.4 * CONST['M_sun'])
        m2 = dataset.get('m2', 1.4 * CONST['M_sun'])
        e = dataset.get('e', 0.617)
        M = m1 + m2
        omega_dot = 3 * (P_b / (2 * CONST['pi']))**(-5.0/3.0) * (G * M / c**3)**(2.0/3.0) / (1 - e**2)
        return {
            'omega_dot_rad_s': omega_dot,
            'omega_dot_deg_yr': omega_dot * 180 / CONST['pi'] * 3.156e7,
            'equation': 'ω̇ = 3(P_b/2π)^(-5/3)(GM/c³)^(2/3)(1-e²)⁻¹',
            'source': 'Grok URL Eq. 74: GR Periastron Advance'
        }

class KilonovaPeakLuminosityCalculator:
    """Eq. 75: L_peak ≈ 10⁴¹(M_ej/0.01)(v_ej/0.1c)(κ/1)⁻¹ erg/s."""
    def compute(self, dataset: dict) -> dict:
        c = CONST['c']
        M_ej = dataset.get('M_ej', 0.01 * CONST['M_sun'])
        v_ej = dataset.get('v_ej', 0.1 * c)
        kappa = dataset.get('kappa', 1.0)  # opacity cm²/g
        L_peak_cgs = 1e41 * (M_ej / (0.01 * CONST['M_sun'])) * (v_ej / (0.1 * c)) / kappa
        L_peak_W = L_peak_cgs * 1e-7
        return {
            'L_peak_erg_s': L_peak_cgs, 'L_peak_W': L_peak_W,
            'L_peak_Lsun': L_peak_W / 3.828e26,
            'equation': 'L_peak ≈ 10⁴¹(M_ej/0.01)(v_ej/0.1c)(κ)⁻¹ erg/s',
            'source': 'Grok URL Eq. 75: Kilonova Peak Luminosity'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 76-78: COSMIC RAYS
# ═══════════════════════════════════════════════════════════════════════════════

class FermiSecondOrderCalculator:
    """Eq. 76: dE/dt = (4/3)(v_c²/c²)(E/λ) (Fermi II acceleration)."""
    def compute(self, dataset: dict) -> dict:
        c = CONST['c']
        v_c = dataset.get('v_c', 10e3)  # m/s cloud speed
        E = dataset.get('E', 1e9 * CONST['e_charge'])  # GeV → J
        lambda_mfp = dataset.get('lambda_mfp', 3.086e16)  # pc → m
        dEdt = (4.0/3.0) * (v_c**2 / c**2) * E / lambda_mfp if lambda_mfp > 0 else 0
        t_acc = E / dEdt if dEdt > 0 else float('inf')
        return {
            'dE_dt_J_s': dEdt, 't_acc_s': t_acc, 't_acc_yr': t_acc / 3.156e7,
            'equation': 'dE/dt = (4/3)(v_c²/c²)(E/λ)',
            'source': 'Grok URL Eq. 76: Fermi Second-Order'
        }

class KneeEnergyDSACalculator:
    """Eq. 77: E_max = ZeBu_s r_g (Maximum DSA energy)."""
    def compute(self, dataset: dict) -> dict:
        e = CONST['e_charge']
        Z = dataset.get('Z', 1)
        B = dataset.get('B', 3e-10)  # T (~3 μG)
        u_s = dataset.get('u_s', 1e6)  # m/s (~1000 km/s)
        R = dataset.get('R', 10 * 3.086e16)  # pc → m
        E_max = Z * e * B * u_s * R
        return {
            'E_max_J': E_max, 'E_max_eV': E_max / CONST['e_charge'],
            'E_max_PeV': E_max / (CONST['e_charge'] * 1e15),
            'equation': 'E_max = ZeBu_s R ≈ 3×10¹⁵ Z eV',
            'source': 'Grok URL Eq. 77: Knee Energy DSA'
        }

class DiffusionCoefficientCalculator:
    """Eq. 78: D(E) = 10²⁸(E/10GeV)^δ cm²/s."""
    def compute(self, dataset: dict) -> dict:
        E_GeV = dataset.get('E_GeV', 10.0)
        delta = dataset.get('delta', 0.4)
        D_cgs = 1e28 * (E_GeV / 10.0)**delta
        D_SI = D_cgs * 1e-4  # cm²/s → m²/s
        return {
            'D_cm2_s': D_cgs, 'D_m2_s': D_SI,
            'equation': 'D(E) = 10²⁸(E/10GeV)^δ cm²/s',
            'source': 'Grok URL Eq. 78: CR Diffusion Coefficient'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 79-81: IGM
# ═══════════════════════════════════════════════════════════════════════════════

class WHIMTemperatureCalculator:
    """Eq. 79: T = μm_H GM/(2kr) (WHIM virial temperature)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; k_B = CONST['k_B']; m_H = CONST['m_p']
        mu = dataset.get('mu', 0.6)
        M = dataset.get('M', 1e12 * CONST['M_sun'])
        r = dataset.get('r', 200e3 * 3.086e16)  # kpc → m
        T = mu * m_H * G * M / (2 * k_B * r) if r > 0 else 0
        return {
            'T_K': T, 'T_keV': k_B * T / (CONST['e_charge'] * 1e3),
            'equation': 'T = μm_H GM/(2kr)',
            'source': 'Grok URL Eq. 79: WHIM Temperature'
        }

class MetalEnrichmentRateCalculator:
    """Eq. 80: Ż = Y × SFR/M_gas - Z × Ṁ_out/M_gas."""
    def compute(self, dataset: dict) -> dict:
        Y = dataset.get('Y_yield', 0.02)
        SFR = dataset.get('SFR', 1.0)  # M_sun/yr
        M_gas = dataset.get('M_gas', 1e10)  # M_sun
        Z = dataset.get('Z', 0.02)
        M_dot_out = dataset.get('M_dot_out', 0.5)  # M_sun/yr
        Z_dot = Y * SFR / M_gas - Z * M_dot_out / M_gas if M_gas > 0 else 0
        return {
            'Z_dot': Z_dot, 'Z_dot_per_Gyr': Z_dot * 1e9,
            'equation': 'Ż = Y×SFR/M_gas - Z×Ṁ_out/M_gas',
            'source': 'Grok URL Eq. 80: Metal Enrichment Rate'
        }

class CoolingTimeBremsstrahlungCalculator:
    """Eq. 81: t_cool = 3nkT/(2n_e n_i Λ(T))."""
    def compute(self, dataset: dict) -> dict:
        k_B = CONST['k_B']
        n = dataset.get('n', 1e3)  # m⁻³
        T = dataset.get('T', 1e6)  # K
        n_e = dataset.get('n_e', n * 0.88)
        n_i = dataset.get('n_i', n * 0.12)
        Lambda = dataset.get('Lambda_cool', 3e-23 * 1e-7)  # erg cm³/s → W m³
        t_cool = 3 * n * k_B * T / (2 * n_e * n_i * Lambda) if n_e > 0 and n_i > 0 and Lambda > 0 else float('inf')
        return {
            't_cool_s': t_cool, 't_cool_Gyr': t_cool / 3.156e16,
            'equation': 't_cool = 3nkT/(2n_e n_i Λ)',
            'source': 'Grok URL Eq. 81: Bremsstrahlung Cooling Time'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 82-84: FIRST GALAXIES
# ═══════════════════════════════════════════════════════════════════════════════

class PressSchechterMassFunctionCalculator:
    """Eq. 82: Press-Schechter (PS) halo mass function."""
    def compute(self, dataset: dict) -> dict:
        rho_0 = dataset.get('rho_0', COSMO['rho_c'] * COSMO['Omega_m'])
        M = dataset.get('M', 1e10 * CONST['M_sun'])
        sigma = dataset.get('sigma', 1.0)
        delta_c = dataset.get('delta_c', COSMO['delta_c'])
        d_ln_sigma = dataset.get('d_ln_sigma_d_ln_M', -0.3)
        dn_dM = math.sqrt(2 / CONST['pi']) * (rho_0 / M) * (delta_c / sigma) * abs(d_ln_sigma) * math.exp(-delta_c**2 / (2 * sigma**2)) / M if M > 0 and sigma > 0 else 0
        return {
            'dn_dM': dn_dM, 'sigma': sigma, 'delta_c': delta_c,
            'equation': 'dn/dM = √(2/π)(ρ₀/M)(δ_c/σ)|d ln σ/d ln M|exp(-δ_c²/(2σ²))',
            'source': 'Grok URL Eq. 82: Press-Schechter Mass Function'
        }

class SFEfficiencyCalculator:
    """Eq. 83: ε_* = f_b Ṁ_halo/(M_halo H(z))(1+M_halo/M_crit)⁻¹."""
    def compute(self, dataset: dict) -> dict:
        f_b = dataset.get('f_b', 0.16)
        M_dot_halo = dataset.get('M_dot_halo', 10.0)  # M_sun/yr
        M_halo = dataset.get('M_halo', 1e11)  # M_sun
        H_z = dataset.get('H_z', COSMO['H_0'])  # km/s/Mpc
        M_crit = dataset.get('M_crit', 1e12)  # M_sun
        H_z_inv_yr = 1.0 / (H_z / 978.0)  # Hubble time in Gyr → 1/yr
        epsilon_star = f_b * M_dot_halo / (M_halo * H_z_inv_yr) / (1 + M_halo / M_crit) if M_halo > 0 and M_crit > 0 else 0
        return {
            'epsilon_star': epsilon_star,
            'equation': 'ε_* = f_b Ṁ_halo/(M_halo H)(1+M_halo/M_crit)⁻¹',
            'source': 'Grok URL Eq. 83: Star Formation Efficiency'
        }

class FeedbackEnergyInjectionCalculator:
    """Eq. 84: E_fb = η Ṁ_* c² (Feedback energy injection)."""
    def compute(self, dataset: dict) -> dict:
        c = CONST['c']
        eta = dataset.get('eta_feedback', 1e-3)
        M_dot_star = dataset.get('M_dot_star', 1.0 * CONST['M_sun'] / 3.156e7)  # M_sun/yr → kg/s
        E_fb = eta * M_dot_star * c**2
        return {
            'E_fb_W': E_fb, 'E_fb_erg_s': E_fb * 1e7,
            'equation': 'E_fb = η Ṁ_* c²',
            'source': 'Grok URL Eq. 84: Feedback Energy Injection'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 85-87: QUANTUM FLUCTUATIONS
# ═══════════════════════════════════════════════════════════════════════════════

class CurvaturePerturbationAmplitudeCalculator:
    """Eq. 85: Δ_R² = H²/(8π²εM_pl²) ≈ 2.1×10⁻⁹."""
    def compute(self, dataset: dict) -> dict:
        H = dataset.get('H_inflation', 1e36)  # 1/s
        epsilon = dataset.get('epsilon', 0.01)
        M_pl = CONST['M_pl']
        Delta_R_sq = H**2 / (8 * CONST['pi']**2 * epsilon * M_pl**2 * CONST['c']**2) if epsilon > 0 else 0
        return {
            'Delta_R_sq': Delta_R_sq,
            'consistent_with_CMB': abs(Delta_R_sq - 2.1e-9) / 2.1e-9 < 0.1 if Delta_R_sq > 0 else False,
            'equation': 'Δ_R² = H²/(8π²εM_pl²) ≈ 2.1×10⁻⁹',
            'source': 'Grok URL Eq. 85: Curvature Perturbation Amplitude'
        }

class NonGaussianityFNLCalculator:
    """Eq. 86: f_NL non-Gaussianity parameter (simplified)."""
    def compute(self, dataset: dict) -> dict:
        epsilon = dataset.get('epsilon', 0.01)
        eta = dataset.get('eta_sr', 0.01)
        # Single-field slow-roll: f_NL ~ O(ε, η)
        f_NL = 5.0/12.0 * (6 * epsilon - 2 * eta)
        return {
            'f_NL': f_NL,
            'detectable': abs(f_NL) > 5,
            'equation': 'f_NL ≈ (5/12)(6ε - 2η) for single-field',
            'source': 'Grok URL Eq. 86: Non-Gaussianity f_NL'
        }

class ReheatingTemperatureCalculator:
    """Eq. 87: T_reh = (30V_end/(π²g_*))^(1/4) exp(-3N_reh/4)."""
    def compute(self, dataset: dict) -> dict:
        V_end = dataset.get('V_end', 1e64)  # energy density at end
        g_star = dataset.get('g_star', 100.0)
        N_reh = dataset.get('N_reh', 10.0)
        T_reh = (30 * V_end / (CONST['pi']**2 * g_star))**(0.25) * math.exp(-3 * N_reh / 4)
        return {
            'T_reh': T_reh, 'T_reh_GeV': T_reh * CONST['k_B'] / (1.602e-10),
            'equation': 'T_reh = (30V_end/(π²g_*))^(1/4)exp(-3N_reh/4)',
            'source': 'Grok URL Eq. 87: Reheating Temperature'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 88-90: MHD DYNAMOS
# ═══════════════════════════════════════════════════════════════════════════════

class SmallScaleDynamoGrowthCalculator:
    """Eq. 88: dE_M/dt = (3/2)(E_M/t_eddy)(Re_m^(1/2)-1)."""
    def compute(self, dataset: dict) -> dict:
        E_M = dataset.get('E_M', 1e20)  # J magnetic energy
        t_eddy = dataset.get('t_eddy', 1e14)  # s eddy turnover
        Re_m = dataset.get('Re_m', 1e10)
        growth_rate = (3.0/2.0) * (E_M / t_eddy) * (math.sqrt(Re_m) - 1)
        e_folding_time = t_eddy / (1.5 * (math.sqrt(Re_m) - 1)) if Re_m > 1 else float('inf')
        return {
            'dE_M_dt': growth_rate, 'e_folding_s': e_folding_time,
            'equation': 'dE_M/dt = (3/2)(E_M/t_eddy)(Re_m^(1/2)-1)',
            'source': 'Grok URL Eq. 88: Small-Scale Dynamo Growth'
        }

class AlfvenMachNumberCalculator:
    """Eq. 89: M_A = v_turb/v_A = √(4πρv²/B²)."""
    def compute(self, dataset: dict) -> dict:
        rho = dataset.get('rho', 1e-21)  # kg/m³
        v_turb = dataset.get('v_turb', 10e3)  # m/s
        B = dataset.get('B', 5e-10)  # T
        v_A = B / math.sqrt(CONST['mu_0'] * rho) if rho > 0 else 0
        M_A = v_turb / v_A if v_A > 0 else float('inf')
        return {
            'M_A': M_A, 'v_A': v_A,
            'sub_alfvenic': M_A < 1,
            'equation': 'M_A = v_turb/v_A',
            'source': 'Grok URL Eq. 89: Alfvén Mach Number'
        }

class FieldReversalScaleCalculator:
    """Eq. 90: ℓ_rev = (α_dynamo/η)^(1/2) ℓ_force."""
    def compute(self, dataset: dict) -> dict:
        alpha_dyn = dataset.get('alpha_dynamo', 1e3)  # m/s
        eta_diff = dataset.get('eta_diffusivity', 1e10)  # m²/s
        ell_force = dataset.get('ell_force', 100 * 3.086e16)  # pc → m
        ell_rev = math.sqrt(alpha_dyn / eta_diff) * ell_force if eta_diff > 0 else 0
        return {
            'ell_rev_m': ell_rev, 'ell_rev_pc': ell_rev / 3.086e16,
            'equation': 'ℓ_rev = √(α_dynamo/η) × ℓ_force',
            'source': 'Grok URL Eq. 90: Field Reversal Scale'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 91-93: DARK ENERGY
# ═══════════════════════════════════════════════════════════════════════════════

class EquationOfStateWCalculator:
    """Eq. 91: w = p/(ρc²) (Dark energy equation of state)."""
    def compute(self, dataset: dict) -> dict:
        c = CONST['c']
        p = dataset.get('p', -1e-10)  # Pa
        rho = dataset.get('rho', 1e-9 / c**2)  # kg/m³
        w = p / (rho * c**2) if rho > 0 else -1
        return {
            'w': w,
            'is_cosmological_constant': abs(w + 1) < 0.05,
            'is_phantom': w < -1,
            'equation': 'w = p/(ρc²)',
            'source': 'Grok URL Eq. 91: Dark Energy EoS'
        }

class CPLDensityEvolutionCalculator:
    """Eq. 92: ρ_DE(a) CPL parameterization w(a) = w_0 + w_a(1-a)."""
    def compute(self, dataset: dict) -> dict:
        rho_DE_0 = dataset.get('rho_DE_0', 1e-26 * 0.7)  # kg/m³
        w_0 = dataset.get('w_0', -1.0)
        w_a = dataset.get('w_a', 0.0)
        a = dataset.get('a', 1.0)
        w_a_val = w_0 + w_a * (1 - a)
        # ρ_DE(a) = ρ_DE,0 × a^(-3(1+w_0+w_a)) × exp(-3w_a(1-a))
        rho_DE = rho_DE_0 * a**(-3 * (1 + w_0 + w_a)) * math.exp(-3 * w_a * (1 - a)) if a > 0 else 0
        return {
            'rho_DE': rho_DE, 'w_at_a': w_a_val,
            'equation': 'ρ_DE(a) = ρ_DE,0 × a^(-3(1+w₀+w_a)) × exp(-3w_a(1-a))',
            'source': 'Grok URL Eq. 92: CPL Dark Energy Evolution'
        }

class GrowthSuppressionCalculator:
    """Eq. 93: D(a) growth factor with dark energy suppression."""
    def compute(self, dataset: dict) -> dict:
        Omega_m = dataset.get('Omega_m', COSMO['Omega_m'])
        a = dataset.get('a', 1.0)
        # Approximate growth factor D(a) ~ a for matter domination
        # With Λ: D(a) ≈ (5Ω_m/2) × H(a)/H_0 × ∫...
        # Simplified Carroll et al. 1992 approximation
        Omega_m_a = Omega_m / (Omega_m + (1 - Omega_m) * a**3) if a > 0 else Omega_m
        D = a * Omega_m_a**0.6 / Omega_m**0.6 if Omega_m > 0 else a
        f_growth = Omega_m_a**0.55  # growth rate
        return {
            'D_growth': D, 'f_growth': f_growth, 'Omega_m_a': Omega_m_a,
            'equation': 'D(a) ≈ a × (Ω_m(a)/Ω_m,0)^0.6',
            'source': 'Grok URL Eq. 93: Growth Suppression'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 94-96: BLACK HOLE THERMODYNAMICS
# ═══════════════════════════════════════════════════════════════════════════════

class HawkingTemperatureCalculator:
    """Eq. 94: T_H = ℏc³/(8πGMk_B) (Hawking temperature)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']; hbar = CONST['hbar']; k_B = CONST['k_B']
        M = dataset.get('M', CONST['M_sun'])
        T_H = hbar * c**3 / (8 * CONST['pi'] * G * M * k_B) if M > 0 else 0
        return {
            'T_H_K': T_H, 'M_Msun': M / CONST['M_sun'],
            'equation': 'T_H = ℏc³/(8πGMk_B)',
            'source': 'Grok URL Eq. 94: Hawking Temperature'
        }

class BekensteinHawkingEntropyCalculator:
    """Eq. 95: S = k_B c³ A/(4Gℏ) = 4πk_B GM²/(ℏc)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']; hbar = CONST['hbar']; k_B = CONST['k_B']
        M = dataset.get('M', CONST['M_sun'])
        r_s = 2 * dpm_emergent_ug1(M, c)
        A = 4 * CONST['pi'] * r_s**2
        S = k_B * c**3 * A / (4 * G * hbar)
        return {
            'S_J_K': S, 'S_over_kB': S / k_B, 'A_m2': A, 'r_s_m': r_s,
            'equation': 'S = k_B c³ A/(4Gℏ)',
            'source': 'Grok URL Eq. 95: Bekenstein-Hawking Entropy'
        }

class EvaporationLifetimeCalculator:
    """Eq. 96: τ_evap = 5120πG²M³/(ℏc⁴)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']; hbar = CONST['hbar']
        M = dataset.get('M', CONST['M_sun'])
        tau = 5120 * CONST['pi'] * G**2 * M**3 / (hbar * c**4)
        return {
            'tau_evap_s': tau, 'tau_evap_yr': tau / 3.156e7,
            'equation': 'τ_evap = 5120πG²M³/(ℏc⁴)',
            'source': 'Grok URL Eq. 96: BH Evaporation Lifetime'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 97-99: LQC BOUNCE
# ═══════════════════════════════════════════════════════════════════════════════

class LQCEffectiveFriedmannCalculator:
    """Eq. 97: H² = (8πGρ/3)(1-ρ/ρ_crit) (LQC modified Friedmann)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        rho = dataset.get('rho', 1e90)  # kg/m³
        rho_crit_LQC = dataset.get('rho_crit_LQC', 0.41 * 5.155e96)  # ρ_Planck
        H_sq = (8 * CONST['pi'] * G * rho / 3) * (1 - rho / rho_crit_LQC)
        bounce = rho >= rho_crit_LQC
        return {
            'H_sq': H_sq, 'H': math.sqrt(abs(H_sq)),
            'bounce': bounce, 'rho_over_rho_crit': rho / rho_crit_LQC,
            'equation': 'H² = (8πGρ/3)(1-ρ/ρ_crit)',
            'source': 'Grok URL Eq. 97: LQC Effective Friedmann'
        }

class LQCPerturbationSpectrumCalculator:
    """Eq. 98: P(k) ∝ k^(n_s-1)(1+k/k_*)^(-α) (LQC perturbation)."""
    def compute(self, dataset: dict) -> dict:
        k = dataset.get('k', 0.05)  # Mpc⁻¹
        n_s = dataset.get('n_s', 0.965)
        k_star = dataset.get('k_star', 0.001)  # bounce scale
        alpha = dataset.get('alpha', 2.0)
        A_s = dataset.get('A_s', 2.1e-9)
        P = A_s * k**(n_s - 1) * (1 + k / k_star)**(-alpha)
        return {
            'P_k': P, 'k': k, 'LQC_suppression': (1 + k / k_star)**(-alpha),
            'equation': 'P(k) ∝ k^(n_s-1)(1+k/k_*)^(-α)',
            'source': 'Grok URL Eq. 98: LQC Perturbation Spectrum'
        }

class BounceTimescaleCalculator:
    """Eq. 99: t_b ~ √(3/(8πGρ_crit)) ~ t_Pl."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        rho_crit_LQC = dataset.get('rho_crit_LQC', 0.41 * 5.155e96)
        t_b = math.sqrt(3.0 / (8 * CONST['pi'] * G * rho_crit_LQC)) if rho_crit_LQC > 0 else 0
        return {
            't_bounce_s': t_b, 't_bounce_over_t_Pl': t_b / CONST['t_pl'],
            'equation': 't_b ~ √(3/(8πGρ_crit)) ~ t_Pl',
            'source': 'Grok URL Eq. 99: Bounce Timescale'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR 100: ROCHE LOBE
# ═══════════════════════════════════════════════════════════════════════════════

class RocheLobeRadiusCalculator:
    """Eq. 100: R_L = 0.49q^(2/3)a / (0.6q^(2/3)+ln(1+q^(1/3))) (Eggleton)."""
    def compute(self, dataset: dict) -> dict:
        q = dataset.get('q', 0.001)  # M_p/M_*
        a = dataset.get('a', 1.496e11)  # m (1 AU)
        q23 = q**(2.0/3.0)
        q13 = q**(1.0/3.0)
        R_L = 0.49 * q23 * a / (0.6 * q23 + math.log(1 + q13)) if (0.6 * q23 + math.log(1 + q13)) > 0 else 0
        return {
            'R_L_m': R_L, 'R_L_AU': R_L / 1.496e11,
            'R_L_over_a': R_L / a if a > 0 else 0,
            'equation': 'R_L = 0.49q^(2/3)a/(0.6q^(2/3)+ln(1+q^(1/3)))',
            'source': 'Grok URL Eq. 100: Roche Lobe Radius (Eggleton)'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 101-107: UQFF FRAMEWORK
# ═══════════════════════════════════════════════════════════════════════════════

class UQFFBuoyancyCalculator:
    """Eq. 101: F_UBii = F_rel × (E_cm/E_LEP) × Q_wave × g(r,t)."""
    def compute(self, dataset: dict) -> dict:
        F_rel = UQFF_CONST['F_rel']
        E_cm = dataset.get('E_cm', 1e3)  # GeV
        E_LEP = UQFF_CONST['E_LEP']
        Q_wave = dataset.get('Q_wave', UQFF_CONST['Q_wave'])
        g = dataset.get('g', 9.81)  # m/s² local gravity
        r = dataset.get('r', 1.0)  # m
        t = dataset.get('t', 1.0)  # s
        rho_vac = dataset.get('rho_vac', UQFF_CONST['rho_vac_SCm'])
        F_UBii = F_rel * (E_cm / E_LEP) * Q_wave * g
        return {
            'F_UBii': F_UBii, 'F_rel': F_rel,
            'E_ratio': E_cm / E_LEP, 'Q_wave': Q_wave,
            'equation': 'F_UBii = F_rel × (E_cm/E_LEP) × Q_wave × g(r,t)',
            'source': 'Grok URL Eq. 101: UQFF Universal Buoyancy'
        }

class UQFFMagnetismCalculator:
    """Eq. 102: Um(t,r,n) = Σ_j [μ_j(t)/r_j × (1-e^(-γt cos(πtn))) × φ_j] × P_SCm × E_react × (1+10¹³f_Heav) × (1+f_quasi)."""
    def compute(self, dataset: dict) -> dict:
        mu_base = UQFF_CONST['mu_base']
        omega_c = UQFF_CONST['omega_c']
        gamma = UQFF_CONST['gamma_decay']
        P_SCm = UQFF_CONST['P_SCm']
        f_Heav = UQFF_CONST['f_Heav']
        f_quasi = UQFF_CONST['f_quasi']
        t = dataset.get('t', 1.0)  # days
        r = dataset.get('r', 1.0)  # m
        n = dataset.get('n', 1)  # quantum level
        phi = dataset.get('phi', 1.0)
        mu_j = (1e3 + 0.4 * math.sin(omega_c * t)) * mu_base
        temporal = 1 - math.exp(-gamma * t * math.cos(CONST['pi'] * t * n))
        E_react = 1e46 * math.exp(-0.0005 * t)
        Um = (mu_j / r) * temporal * phi * P_SCm * E_react * (1 + 1e13 * f_Heav) * (1 + f_quasi) if r > 0 else 0
        return {
            'Um': Um, 'mu_j': mu_j, 'temporal_factor': temporal,
            'E_react': E_react,
            'equation': 'Um(t,r,n) = Σ[μ_j/r_j × (1-e^(-γt cos(πtn)))] × P_SCm × E_react × Heav × quasi',
            'source': 'Grok URL Eq. 102: UQFF Universal Magnetism'
        }

class UQFFElectricFieldCalculator:
    """Eq. 103: E = Um × ρ_vac,[UA] / r."""
    def compute(self, dataset: dict) -> dict:
        Um = dataset.get('Um', 1e20)
        rho_vac_UA = dataset.get('rho_vac_UA', UQFF_CONST['rho_vac_UA'])
        r = dataset.get('r', 1.0)
        E_field = Um * rho_vac_UA / r if r > 0 else 0
        return {
            'E_field': E_field, 'Um': Um, 'rho_vac_UA': rho_vac_UA,
            'equation': 'E = Um × ρ_vac,[UA] / r',
            'source': 'Grok URL Eq. 103: UQFF Electric Field'
        }

class UQFFNeutronProductionCalculator:
    """Eq. 104: η = k_η × e^(-[SSq]^n/26) × e^(-π-t) × Um / ρ_vac."""
    def compute(self, dataset: dict) -> dict:
        k_eta = dataset.get('k_eta', 1e-113)
        SSq = dataset.get('SSq', 0.57)
        n = dataset.get('n', 1)
        t = dataset.get('t', 1.0)
        Um = dataset.get('Um', 1e20)
        rho_vac = dataset.get('rho_vac', UQFF_CONST['rho_vac_UA'])
        exponent_SSq = -SSq**n / 26
        exponent_pi = -(CONST['pi'] + t)
        eta = k_eta * math.exp(exponent_SSq) * math.exp(exponent_pi) * Um / rho_vac if rho_vac > 0 else 0
        return {
            'eta': eta, 'k_eta': k_eta, 'SSq': SSq, 'n': n,
            'equation': 'η = k_η × e^(-[SSq]^n/26) × e^(-π-t) × Um/ρ_vac',
            'source': 'Grok URL Eq. 104: UQFF Neutron Production'
        }

class UQFFPseudoMonopoleCalculator:
    """Eq. 105: δ_n = (2π)^n / 6 (Pseudo-monopole states)."""
    def compute(self, dataset: dict) -> dict:
        n = dataset.get('n', 1)
        delta_n = (2 * CONST['pi'])**n / 6
        return {
            'delta_n': delta_n, 'n': n,
            'equation': 'δ_n = (2π)^n / 6',
            'source': 'Grok URL Eq. 105: UQFF Pseudo-Monopole States'
        }

class UQFFVacuumDensityCalculator:
    """Eq. 106: ρ_vac,[UA']:SCm(n,t) = 10⁻²³ × (0.1)^n × e^(-[SSq]^n/26) × e^(-π-t)."""
    def compute(self, dataset: dict) -> dict:
        n = dataset.get('n', 1)
        t = dataset.get('t', 1.0)
        SSq = dataset.get('SSq', 0.57)
        rho_vac = 1e-23 * (0.1)**n * math.exp(-SSq**n / 26) * math.exp(-(CONST['pi'] + t))
        # Also compute per-level density
        rho_level = 10**(-n * 2)  # J/m³ approximation for level n
        return {
            'rho_vac': rho_vac, 'rho_level': rho_level, 'n': n,
            'equation': 'ρ_vac,[UA\']:SCm(n,t) = 10⁻²³ × (0.1)^n × e^(-[SSq]^n/26) × e^(-π-t)',
            'source': 'Grok URL Eq. 106: UQFF Vacuum Density'
        }

class UQFFGinzburgLandauCalculator:
    """Eq. 107: L_Ug = |∇ψ|² - (m²/2)|ψ|² + (λ/4)|ψ|⁴ (Ginzburg-Landau)."""
    def compute(self, dataset: dict) -> dict:
        psi = dataset.get('psi', 1.0)  # order parameter magnitude
        grad_psi = dataset.get('grad_psi', 0.1)
        m_sq = dataset.get('m_squared', -1.0)  # negative for SSB
        lambda_4 = dataset.get('lambda_coupling', 0.1)
        L_Ug = grad_psi**2 - (m_sq / 2) * psi**2 + (lambda_4 / 4) * psi**4
        # Find minimum: d/d|ψ| = 0 → |ψ|² = m²/λ (when m² < 0)
        psi_min = math.sqrt(-m_sq / lambda_4) if m_sq < 0 and lambda_4 > 0 else 0
        return {
            'L_Ug': L_Ug, 'psi_min': psi_min,
            'SSB': m_sq < 0,
            'equation': 'L_Ug = |∇ψ|² - (m²/2)|ψ|² + (λ/4)|ψ|⁴',
            'source': 'Grok URL Eq. 107: UQFF Ginzburg-Landau'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 108-113: MUGE SYSTEMS
# ═══════════════════════════════════════════════════════════════════════════════

class MUGEHydrogenAtomCalculator:
    """Eq. 108: MUGE for Hydrogen atom with superconductive species."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        m_eff = dataset.get('m_eff', CONST['m_p'])
        m_p = CONST['m_p']
        r = dataset.get('r', 5.29e-11)  # Bohr radius
        M_Z = dataset.get('M_Z', 1.673e-27)  # nuclear mass
        r_Z = dataset.get('r_Z', 1e-15)  # nuclear radius
        f_sc = dataset.get('f_sc', 0.2)
        H_0 = dataset.get('H_0', COSMO['H_0'] * 1e3 / 3.086e22)
        t = dataset.get('t', 4.35e17)  # s (age of universe)
        g_Newton = dpm_emergent_ug1(m_eff, r) * m_p if r > 0 else 0
        g_Z = dpm_emergent_ug1(M_Z, r_Z) * (1 + f_sc) if r_Z > 0 else 0
        evolution = math.exp(H_0 * t / CONST['c']) if CONST['c'] > 0 else 1
        g_MUGE = g_Newton + g_Z * evolution
        return {
            'g_MUGE': g_MUGE, 'g_Newton': g_Newton,
            'g_nuclear': g_Z, 'evolution_factor': evolution,
            'UQFF_contribution_pct': (g_MUGE - g_Newton) / g_MUGE * 100 if g_MUGE > 0 else 0,
            'equation': 'g_MUGE = Gm_eff m_p/r² + Σ(GM_Z/r_Z²)(1+f_sc)e^(H₀t/c)',
            'source': 'Grok URL Eq. 108: MUGE Hydrogen Atom'
        }

class MUGERingsOfRelativityCalculator:
    """Eq. 109: MUGE for Einstein Rings (Rings of Relativity)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']; c = CONST['c']
        M = dataset.get('M', 1e14 * CONST['M_sun'])
        r = dataset.get('r', 1e23)  # m
        z = dataset.get('z', 0.5)
        B = dataset.get('B', 1e-9)  # T
        B_crit = dataset.get('B_crit', 4.4e13)  # T
        H_z = dataset.get('H_z', COSMO['H_0'] * 1e3 / 3.086e22)
        t = dataset.get('t', 4.35e17)
        g_GR = dpm_emergent_ug1(M, r) if r > 0 else 0
        H_term = 1 + H_z * t
        B_term = 1 - B / B_crit
        L_t = dataset.get('L_t', 0.01)  # luminosity evolution
        g_Rings = g_GR * H_term * B_term * (1 + L_t)
        # Einstein radius
        D_L = dataset.get('D_L', 1e25)
        D_S = dataset.get('D_S', 3e25)
        D_LS = dataset.get('D_LS', 2e25)
        R_E = math.sqrt(4 * dpm_emergent_ug1(M, c) * D_LS * D_L / D_S) if D_S > 0 else 0
        return {
            'g_Rings': g_Rings, 'g_GR': g_GR, 'R_Einstein': R_E,
            'UQFF_contribution_pct': 30.0,
            'equation': 'g_Rings = (GM/r²)(1+Ht)(1-B/B_crit)(1+L(t))',
            'source': 'Grok URL Eq. 109: MUGE Rings of Relativity'
        }

class MUGEMagnetarCalculator:
    """Eq. 110: MUGE for Magnetar B-field decay."""
    def compute(self, dataset: dict) -> dict:
        B_0 = dataset.get('B_0', 1e14)  # T (converted from 10¹⁰ G)
        t = dataset.get('t', 1000.0)  # years
        tau_decay = dataset.get('tau_decay', 4000.0)  # years
        B_t = B_0 * math.exp(-t / tau_decay) if tau_decay > 0 else B_0
        # Spin angular velocity
        P = dataset.get('P', 5.0)  # s (period)
        Omega = 2 * CONST['pi'] / P if P > 0 else 0
        # Surface gravity (billions × Earth)
        M = dataset.get('M', 1.4 * CONST['M_sun'])
        R = dataset.get('R', 1e4)  # m (10 km)
        g_surface = CONST['G'] * M / R**2 if R > 0 else 0
        return {
            'B_t': B_t, 'B_t_over_B0': B_t / B_0 if B_0 > 0 else 0,
            'Omega': Omega, 'g_surface': g_surface,
            'g_over_Earth': g_surface / 9.81,
            'UQFF_contribution_pct': 50.0,
            'equation': 'B(t) = B₀ exp(-t/τ); UQFF via B_crit coupling',
            'source': 'Grok URL Eq. 110: MUGE Magnetar'
        }

class MUGEGlobularClusterCalculator:
    """Eq. 111: MUGE for Globular Clusters (core collapse + BH likelihood)."""
    def compute(self, dataset: dict) -> dict:
        t = dataset.get('t', 10.0)  # Gyr
        t_cc = dataset.get('t_cc', 12.0)  # Gyr core collapse time
        alpha = dataset.get('alpha', 0.75)
        M = dataset.get('M', 1e6 * CONST['M_sun'])
        Fe_H = dataset.get('Fe_H', -1.5)  # metallicity
        Y = dataset.get('Y', 0.30)  # He abundance
        f_core = (1 - t / t_cc)**alpha if t < t_cc else 0
        # BH likelihood for massive clusters
        f_BH = 0.70 + 0.20 * min(M / (1e6 * CONST['M_sun']), 1.0)  # 70-90%
        return {
            'f_core': f_core, 'f_BH': f_BH,
            'core_collapsed': t >= t_cc,
            'UQFF_contribution_pct': 20.0,
            'equation': 'f_core = (1-t/t_cc)^α; f_BH = 0.70-0.90',
            'source': 'Grok URL Eq. 111: MUGE Globular Cluster'
        }

class MUGESagittariusAStarCalculator:
    """Eq. 112: MUGE for Sgr A* SMBH mass growth."""
    def compute(self, dataset: dict) -> dict:
        M_0 = dataset.get('M_0', 4.3e6 * CONST['M_sun'])
        M_dot = dataset.get('M_dot', 1e-4)  # M_sun/yr (dimensionless for growth)
        t = dataset.get('t', 1e10)  # yr
        tau = dataset.get('tau', 9e9)  # yr growth timescale
        spin_misalign = dataset.get('spin_misalignment_deg', 30.0)
        M_t = M_0 * (1 + M_dot * math.exp(-t / tau))
        # Include Λ term
        Lambda = COSMO['Lambda']
        g_at_horizon = CONST['G'] * M_t / (2 * CONST['G'] * M_t / CONST['c']**2)**2 if M_t > 0 else 0
        return {
            'M_t': M_t, 'M_t_Msun': M_t / CONST['M_sun'],
            'g_horizon': g_at_horizon,
            'spin_misalignment': spin_misalign,
            'UQFF_contribution_pct': 35.0,
            'equation': 'M(t) = M₀(1+Ṁ exp(-t/τ)); Kerr+Λ+UQFF',
            'source': 'Grok URL Eq. 112: MUGE Sagittarius A*'
        }

class MUGESunPlanetarySystemCalculator:
    """Eq. 113: MUGE for Sun Planetary System with resonance."""
    def compute(self, dataset: dict) -> dict:
        A_1 = dataset.get('A_1', 1.0)  # resonance amplitude
        A_2 = dataset.get('A_2', 0.5)
        f = dataset.get('f', 1.0)  # Hz
        t = dataset.get('t', 0.0)  # s
        phi = dataset.get('phi', 0.0)  # rad
        k = dataset.get('k', 1.0)
        f_dp = dataset.get('f_dp', 2.0)  # di-pseudo-monopole freq
        phi_dp = dataset.get('phi_dp', 0.0)
        SC_m = dataset.get('SC_m', 1.0)
        # Resonance potential
        U_r = A_1 * math.sin(2 * CONST['pi'] * f * t) + A_2 * math.sin(2 * CONST['pi'] * f * t + phi)
        # Di-pseudo-monopole reciprocation
        U_dp = k * (A_1 * A_2 / f_dp**2) * math.cos(phi_dp) if f_dp > 0 else 0
        # 26 quantum level densities
        levels = {i: 10**(-i * 2) for i in range(1, 27)}
        return {
            'U_resonance': U_r, 'U_dipseudo': U_dp, 'SC_m': SC_m,
            'quantum_levels_sample': {1: levels[1], 10: levels[10], 13: levels[13], 26: levels[26]},
            'UQFF_contribution_pct': 25.0,
            'equation': 'U_r = AΣsin(2πft + φ); U_dp = k(A₁A₂/f²)cos(φ)',
            'source': 'Grok URL Eq. 113: MUGE Sun Planetary System'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATORS 114-121: UPDATED UQFF APPLICATIONS
# ═══════════════════════════════════════════════════════════════════════════════

class UQFFProtostellarJetCalculator:
    """Eq. 114: F_UBii,jet and Um_jet for protostellar jets."""
    def compute(self, dataset: dict) -> dict:
        F_rel = UQFF_CONST['F_rel']
        E_shock = dataset.get('E_shock', 1e3)  # GeV
        E_LEP = UQFF_CONST['E_LEP']
        Q_wave = dataset.get('Q_wave', UQFF_CONST['Q_wave'])
        g_disk = dataset.get('g_disk', 1.0)  # m/s²
        tau_damp = dataset.get('tau_damp', 1e6)  # s
        t = dataset.get('t', 1e5)  # s
        T_B = dataset.get('T_B', 1e20)  # N·m
        B = dataset.get('B', 1e-3)  # T
        r_A = dataset.get('r_A', 1e12)  # m
        gamma = UQFF_CONST['gamma_decay']
        n = dataset.get('n', 1)
        F_UBii_jet = -F_rel * (E_shock / E_LEP) * Q_wave * g_disk * math.exp(-t / tau_damp) if tau_damp > 0 else 0
        mu_j = UQFF_CONST['mu_base'] * 1e3
        r_j = dataset.get('r', 1e12)
        temporal = 1 - math.exp(-gamma * t * math.cos(CONST['pi'] * t * n))
        Um_jet = mu_j / r_j * temporal * T_B / (B * r_A**2) if r_j > 0 and B > 0 and r_A > 0 else 0
        return {
            'F_UBii_jet': F_UBii_jet, 'Um_jet': Um_jet,
            'equation': 'F_UBii,jet = -F_rel(E_shock/E_LEP)Q_wave × g × e^(-t/τ)',
            'source': 'Grok URL Eq. 114: UQFF Protostellar Jet'
        }

class UQFFGalaxyMergerCalculator:
    """Eq. 115: F_UBii,merger and Um_merger for galaxy mergers."""
    def compute(self, dataset: dict) -> dict:
        F_rel = UQFF_CONST['F_rel']
        E_burst = dataset.get('E_burst', 1e5)  # GeV
        E_LEP = UQFF_CONST['E_LEP']
        Q_wave_z = dataset.get('Q_wave_z', UQFF_CONST['Q_wave'])
        g_halo = dataset.get('g_halo', 1e-10)  # m/s²
        z = dataset.get('z', 1.0)
        m = dataset.get('m_exponent', 2.0)
        sigma_M = dataset.get('sigma_M', 0.8)
        sigma_m = dataset.get('sigma_m', 1.2)
        delta_c = dataset.get('delta_c', COSMO['delta_c'])
        F_UBii_merger = F_rel * (E_burst / E_LEP) * Q_wave_z * g_halo * (1 + z)**m
        dsigma = sigma_m**2 - sigma_M**2
        exponent = -delta_c**2 / (2 * dsigma) if dsigma > 0 else 0
        mu = UQFF_CONST['mu_base']
        rho_vac = UQFF_CONST['rho_vac_SCm']
        gamma = UQFF_CONST['gamma_decay']
        t = dataset.get('t', 1.0)
        Um_merger = mu * (1e3 * rho_vac) * math.exp(exponent) * (1 + z)**(m/2) * (1 - math.exp(-gamma * t))
        return {
            'F_UBii_merger': F_UBii_merger, 'Um_merger': Um_merger,
            'equation': 'F_UBii,merger = +F_rel(E_burst/E_LEP)Q_z × g_halo × (1+z)^m',
            'source': 'Grok URL Eq. 115: UQFF Galaxy Merger'
        }

class UQFFBlackHoleGrowthCalculator:
    """Eq. 116: F_UBii,BH and Um_BH for black hole accretion growth."""
    def compute(self, dataset: dict) -> dict:
        F_rel = UQFF_CONST['F_rel']; c = CONST['c']; G = CONST['G']
        M_dot_BH = dataset.get('M_dot_BH', 1.0 * CONST['M_sun'] / 3.156e7)  # kg/s
        E_LEP = UQFF_CONST['E_LEP']
        Q_wave = dataset.get('Q_wave', UQFF_CONST['Q_wave'])
        M_BH = dataset.get('M_BH', 1e8 * CONST['M_sun'])
        r = dataset.get('r', 1e14)  # m
        sigma = dataset.get('sigma', 1.0)
        delta_c = dataset.get('delta_c', COSMO['delta_c'])
        a = dataset.get('a_spin', 0.9)
        gamma = UQFF_CONST['gamma_decay']
        t = dataset.get('t', 1.0)
        rho_vac = UQFF_CONST['rho_vac_SCm']
        E_acc = M_dot_BH * c**2  # GeV equivalent handled by ratio
        F_UBii_BH = -F_rel * (E_acc / (E_LEP * 1.602e-10)) * Q_wave * (4 * CONST['pi'] * G * M_BH / (c**2 * r))
        erfc_val = math.erfc(delta_c / (math.sqrt(2) * sigma))
        R_H = dpm_emergent_ug1(M_BH, c)
        B_sq_term = dataset.get('B', 1e-1)**2 * R_H**4 / (4 * CONST['pi'] * c)
        mu = UQFF_CONST['mu_base'] * rho_vac
        Um_BH = mu * (1 - math.exp(-gamma * t)) * (a * c**3 / (2 * G * M_BH)) * B_sq_term if M_BH > 0 else 0
        return {
            'F_UBii_BH': F_UBii_BH, 'Um_BH': Um_BH,
            'erfc_collapse': erfc_val,
            'equation': 'F_UBii,BH = -F_rel(Ṁc²/E_LEP)Q × (4πGM/(c²r)) × erfc',
            'source': 'Grok URL Eq. 116: UQFF Black Hole Growth'
        }

class UQFFIssingAnyonCalculator:
    """Eq. 117: F_UBii,anyons with Ising braiding (2025 insight)."""
    def compute(self, dataset: dict) -> dict:
        F_rel = UQFF_CONST['F_rel']
        E_anyons = dataset.get('E_anyons', 1e-3)  # GeV (meV scale)
        E_LEP = UQFF_CONST['E_LEP']
        Q_wave = dataset.get('Q_wave', UQFF_CONST['Q_wave'])
        g = dataset.get('g', 9.81)
        delta_c = dataset.get('delta_c', COSMO['delta_c'])
        sigma = dataset.get('sigma', 1.0)
        F_UBii_anyons = -F_rel * (E_anyons / E_LEP) * Q_wave * g * math.exp(-delta_c**2 / (2 * sigma**2))
        # Topological protection factor
        nu = dataset.get('nu_topological', 5.0/2.0)  # FQHE filling fraction
        return {
            'F_UBii_anyons': F_UBii_anyons, 'nu': nu,
            'equation': 'F_UBii,anyons = -F_rel(E_anyons/E_LEP)Q × g × exp(-δ_c²/2σ²)',
            'source': 'Grok URL Eq. 117: UQFF Ising Anyons (2025)'
        }

class UQFFPolaritonQFTCalculator:
    """Eq. 118: Um_polariton with curved spacetime QFT."""
    def compute(self, dataset: dict) -> dict:
        mu_base = UQFF_CONST['mu_base']
        gamma = UQFF_CONST['gamma_decay']
        t = dataset.get('t', 1.0)
        r = dataset.get('r', 1e-6)  # m (microcavity)
        n = dataset.get('n', 1)
        v_sound = dataset.get('v_sound', 1e6)  # m/s (polariton sound)
        c = CONST['c']
        Delta_T = dataset.get('Delta_T', 0.01)  # K (Unruh-like temp analog)
        T = dataset.get('T', 300.0)  # K
        mu_j = (1e3 + 0.4 * math.sin(UQFF_CONST['omega_c'] * t)) * mu_base
        temporal = 1 - math.exp(-gamma * t)
        Um_pol = mu_j / r * temporal * (v_sound**2 / c**2) * (1 + Delta_T / T) if r > 0 and c > 0 and T > 0 else 0
        return {
            'Um_polariton': Um_pol,
            'v_sound_over_c': v_sound / c,
            'equation': 'Um_pol = Σ μ_j/r × (1-e^(-γt))(v²/c²)(1+ΔT/T)',
            'source': 'Grok URL Eq. 118: UQFF Polariton QFT (2025)'
        }

class UQFFUTe2TopologicalCalculator:
    """Eq. 119: δ_n,UTe2 with topological factor for UTe2 superconductor."""
    def compute(self, dataset: dict) -> dict:
        n = dataset.get('n', 1)
        SSq = dataset.get('SSq', 0.57)
        f_topo = dataset.get('f_topo', 0.1)  # topological factor
        t = dataset.get('t', 1.0)
        delta_n = (2 * CONST['pi'])**n / 6 * math.exp(-SSq**n / 26) * (1 + f_topo) * math.exp(-(CONST['pi'] + t))
        return {
            'delta_n_UTe2': delta_n, 'n': n, 'f_topo': f_topo,
            'equation': 'δ_n,UTe2 = (2π)^n/6 × e^(-[SSq]^n/26)(1+f_topo)e^(-π-t)',
            'source': 'Grok URL Eq. 119: UQFF UTe2 Topological (2025)'
        }

class UQFFElectricUniverseRatioCalculator:
    """Eq. 120: R = F_EM / F_g (Electric Universe proof via UQFF)."""
    def compute(self, dataset: dict) -> dict:
        G = CONST['G']
        q = dataset.get('q', CONST['e_charge'])
        Um = dataset.get('Um', 1e20)
        rho_vac = dataset.get('rho_vac', UQFF_CONST['rho_vac_UA'])
        v = dataset.get('v', 1e6)  # m/s
        r = dataset.get('r', 1e-15)  # m (nuclear)
        M = dataset.get('M', CONST['m_p'])
        m = dataset.get('m', CONST['m_p'])
        F_EM = q * Um * rho_vac * v / r if r > 0 else 0
        F_g = dpm_emergent_ug1(M, r) * m if r > 0 else 0
        R = F_EM / F_g if F_g > 0 else float('inf')
        return {
            'R_ratio': R, 'log10_R': math.log10(R) if R > 0 else 0,
            'F_EM': F_EM, 'F_g': F_g,
            'exceeds_EU_claim': R > 1e39,
            'equation': 'R = (q×Um×ρ_vac×v/r) / (G×M×m/r²)',
            'source': 'Grok URL Eq. 120: UQFF Electric Universe Ratio'
        }

class UQFFGyroTorqueNullificationCalculator:
    """Eq. 121: τ + Ui = 0 (Gyroscopic torque nullification by UQFF)."""
    def compute(self, dataset: dict) -> dict:
        I = dataset.get('I', 1.0)  # kg·m² moment of inertia
        omega = dataset.get('omega', 100.0)  # rad/s spin rate
        alpha = dataset.get('alpha', 1.0)  # rad/s² precession
        F_UBii = dataset.get('F_UBii', 1e10)  # N
        r = dataset.get('r', 0.1)  # m lever arm
        theta = dataset.get('theta', CONST['pi']/4)  # rad
        tau_gyro = I * omega * alpha
        Ui = -F_UBii * r * math.sin(theta)
        residual = tau_gyro + Ui
        nullified = abs(residual) < abs(tau_gyro) * 0.01 if tau_gyro != 0 else False
        return {
            'tau_gyro': tau_gyro, 'Ui_UQFF': Ui,
            'residual': residual, 'nullified': nullified,
            'nullification_pct': (1 - abs(residual) / abs(tau_gyro)) * 100 if tau_gyro != 0 else 100,
            'equation': 'τ + Ui = 0; τ=Iωα, Ui=-F_UBii×r×sin(θ)',
            'source': 'Grok URL Eq. 121: UQFF Gyro Torque Nullification'
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MASTER REGISTRY - ALL 121 CALCULATOR CLASSES
# ═══════════════════════════════════════════════════════════════════════════════

GROK_URL_CALCULATORS = {
    # Protostellar Jets (1-4)
    'AngularMomentumTransportCalculator': AngularMomentumTransportCalculator(),
    'MHDJetVelocityCalculator': MHDJetVelocityCalculator(),
    'JTypeShockRankineHugoniotCalculator': JTypeShockRankineHugoniotCalculator(),
    'CTypeShockDampingCalculator': CTypeShockDampingCalculator(),
    # Galaxy Mergers (5-7)
    'EPSMergerRateCalculator': EPSMergerRateCalculator(),
    'OrbitalTorqueTimeCalculator': OrbitalTorqueTimeCalculator(),
    'SFRDEvolutionCalculator': SFRDEvolutionCalculator(),
    # BH Growth (8-9)
    'EPSBHMassFunctionCalculator': EPSBHMassFunctionCalculator(),
    'EddingtonAccretionCalculator': EddingtonAccretionCalculator(),
    # SN Remnants (10-11)
    'SedovTaylorExpansionCalculator': SedovTaylorExpansionCalculator(),
    'DSAParticleAccelerationCalculator': DSAParticleAccelerationCalculator(),
    # GW Mergers (12-13)
    'ChirpMassFormulaCalculator': ChirpMassFormulaCalculator(),
    'QNMRingdownCalculator': QNMRingdownCalculator(),
    # Quasar Jets (14-15)
    'BlandfordZnajekPowerCalculator': BlandfordZnajekPowerCalculator(),
    'RelativisticJetVelocityCalculator': RelativisticJetVelocityCalculator(),
    # Neutron Stars (16-18)
    'TOVEquationCalculator': TOVEquationCalculator(),
    'PulsarSpinDownAgeCalculator': PulsarSpinDownAgeCalculator(),
    'GlitchRecoveryCalculator': GlitchRecoveryCalculator(),
    # GRBs (19-20)
    'FireballExpansionCalculator': FireballExpansionCalculator(),
    'AfterglowSynchrotronCalculator': AfterglowSynchrotronCalculator(),
    # CMB (21-22)
    'CMBAngularPowerCalculator': CMBAngularPowerCalculator(),
    'OpticalDepthCalculator': OpticalDepthCalculator(),
    # AGN Feedback (23-25)
    'MomentumFeedbackCalculator': MomentumFeedbackCalculator(),
    'BZJetPowerUpdatedCalculator': BZJetPowerUpdatedCalculator(),
    'FeedbackDutyCycleCalculator': FeedbackDutyCycleCalculator(),
    # Exoplanets (26-28)
    'PhotoevaporationRateCalculator': PhotoevaporationRateCalculator(),
    'TypeIMigrationTorqueCalculator': TypeIMigrationTorqueCalculator(),
    'RadialVelocitySemiAmplitudeCalculator': RadialVelocitySemiAmplitudeCalculator(),
    # Dark Matter (29-31)
    'NFWDensityProfileCalculator': NFWDensityProfileCalculator(),
    'RotationCurveCalculator': RotationCurveCalculator(),
    'SIDMCoreFormationCalculator': SIDMCoreFormationCalculator(),
    # Galaxy Clusters (32-34)
    'StrongLensingMassCalculator': StrongLensingMassCalculator(),
    'XRayMassEstimateCalculator': XRayMassEstimateCalculator(),
    'MergerShockMachCalculator': MergerShockMachCalculator(),
    # Cosmic Voids (35-36)
    'VoidDensityEvolutionCalculator': VoidDensityEvolutionCalculator(),
    'OutflowVelocityCalculator': OutflowVelocityCalculator(),
    # Reionization (37-38)
    'IonizationFractionCalculator': IonizationFractionCalculator(),
    'BubbleRadiusCalculator': BubbleRadiusCalculator(),
    # ISM (39-41)
    'JeansLengthCalculator': JeansLengthCalculator(),
    'AlfvenVelocityCalculator': AlfvenVelocityCalculator(),
    'TurbulentCascadeCalculator': TurbulentCascadeCalculator(),
    # Stellar Evolution (42-44)
    'MainSequenceLifetimeCalculator': MainSequenceLifetimeCalculator(),
    'MassLuminosityRelationCalculator': MassLuminosityRelationCalculator(),
    'ConvectiveTurnoverCalculator': ConvectiveTurnoverCalculator(),
    # BBN (45-46)
    'BaryonPhotonRatioCalculator': BaryonPhotonRatioCalculator(),
    'DeuteriumBottleneckCalculator': DeuteriumBottleneckCalculator(),
    # Cosmology (47-49)
    'FirstFriedmannCalculator': FirstFriedmannCalculator(),
    'SecondFriedmannCalculator': SecondFriedmannCalculator(),
    'DensityParameterCalculator': DensityParameterCalculator(),
    # Inflation (50-52)
    'SlowRollParametersCalculator': SlowRollParametersCalculator(),
    'CurvaturePowerSpectrumCalculator': CurvaturePowerSpectrumCalculator(),
    'EFoldsCalculator': EFoldsCalculator(),
    # GW Backgrounds (53-54)
    'TensorPowerSpectrumCalculator': TensorPowerSpectrumCalculator(),
    'StochasticGWDensityCalculator': StochasticGWDensityCalculator(),
    # Binary BH (55-57)
    'InspiralFrequencyEvolutionCalculator': InspiralFrequencyEvolutionCalculator(),
    'MergerTimeCalculator': MergerTimeCalculator(),
    'RingdownDampingTimeCalculator': RingdownDampingTimeCalculator(),
    # Supernovae (58-60)
    'ArnettLawLightCurveCalculator': ArnettLawLightCurveCalculator(),
    'EjectaVelocityCalculator': EjectaVelocityCalculator(),
    'NucleosynthesisYieldCalculator': NucleosynthesisYieldCalculator(),
    # Planetary Nebulae (61-63)
    'PNExpansionRadiusCalculator': PNExpansionRadiusCalculator(),
    'IonizationFrontVelocityCalculator': IonizationFrontVelocityCalculator(),
    'AGBMassLossReimersCalculator': AGBMassLossReimersCalculator(),
    # Clusters (64-69)
    'ShockMachNumberTempCalculator': ShockMachNumberTempCalculator(),
    'MergerTimescaleVirialCalculator': MergerTimescaleVirialCalculator(),
    'CoolCoreHeatingCalculator': CoolCoreHeatingCalculator(),
    'EvaporationRateCalculator': EvaporationRateCalculator(),
    'RelaxationTimeCalculator': RelaxationTimeCalculator(),
    'VirialMassCalculator': VirialMassCalculator(),
    # Quasar Feedback (70-72)
    'WindTerminalVelocityCalculator': WindTerminalVelocityCalculator(),
    'IonizationParameterUCalculator': IonizationParameterUCalculator(),
    'EnergyCouplingEfficiencyCalculator': EnergyCouplingEfficiencyCalculator(),
    # Binary Pulsars (73-75)
    'OrbitalDecayCalculator': OrbitalDecayCalculator(),
    'PeriastronAdvanceCalculator': PeriastronAdvanceCalculator(),
    'KilonovaPeakLuminosityCalculator': KilonovaPeakLuminosityCalculator(),
    # Cosmic Rays (76-78)
    'FermiSecondOrderCalculator': FermiSecondOrderCalculator(),
    'KneeEnergyDSACalculator': KneeEnergyDSACalculator(),
    'DiffusionCoefficientCalculator': DiffusionCoefficientCalculator(),
    # IGM (79-81)
    'WHIMTemperatureCalculator': WHIMTemperatureCalculator(),
    'MetalEnrichmentRateCalculator': MetalEnrichmentRateCalculator(),
    'CoolingTimeBremsstrahlungCalculator': CoolingTimeBremsstrahlungCalculator(),
    # First Galaxies (82-84)
    'PressSchechterMassFunctionCalculator': PressSchechterMassFunctionCalculator(),
    'SFEfficiencyCalculator': SFEfficiencyCalculator(),
    'FeedbackEnergyInjectionCalculator': FeedbackEnergyInjectionCalculator(),
    # Quantum Fluctuations (85-87)
    'CurvaturePerturbationAmplitudeCalculator': CurvaturePerturbationAmplitudeCalculator(),
    'NonGaussianityFNLCalculator': NonGaussianityFNLCalculator(),
    'ReheatingTemperatureCalculator': ReheatingTemperatureCalculator(),
    # MHD Dynamos (88-90)
    'SmallScaleDynamoGrowthCalculator': SmallScaleDynamoGrowthCalculator(),
    'AlfvenMachNumberCalculator': AlfvenMachNumberCalculator(),
    'FieldReversalScaleCalculator': FieldReversalScaleCalculator(),
    # Dark Energy (91-93)
    'EquationOfStateWCalculator': EquationOfStateWCalculator(),
    'CPLDensityEvolutionCalculator': CPLDensityEvolutionCalculator(),
    'GrowthSuppressionCalculator': GrowthSuppressionCalculator(),
    # BH Thermodynamics (94-96)
    'HawkingTemperatureCalculator': HawkingTemperatureCalculator(),
    'BekensteinHawkingEntropyCalculator': BekensteinHawkingEntropyCalculator(),
    'EvaporationLifetimeCalculator': EvaporationLifetimeCalculator(),
    # LQC Bounce (97-99)
    'LQCEffectiveFriedmannCalculator': LQCEffectiveFriedmannCalculator(),
    'LQCPerturbationSpectrumCalculator': LQCPerturbationSpectrumCalculator(),
    'BounceTimescaleCalculator': BounceTimescaleCalculator(),
    # Exoplanets (100)
    'RocheLobeRadiusCalculator': RocheLobeRadiusCalculator(),
    # UQFF Framework (101-107)
    'UQFFBuoyancyCalculator': UQFFBuoyancyCalculator(),
    'UQFFMagnetismCalculator': UQFFMagnetismCalculator(),
    'UQFFElectricFieldCalculator': UQFFElectricFieldCalculator(),
    'UQFFNeutronProductionCalculator': UQFFNeutronProductionCalculator(),
    'UQFFPseudoMonopoleCalculator': UQFFPseudoMonopoleCalculator(),
    'UQFFVacuumDensityCalculator': UQFFVacuumDensityCalculator(),
    'UQFFGinzburgLandauCalculator': UQFFGinzburgLandauCalculator(),
    # MUGE Systems (108-113)
    'MUGEHydrogenAtomCalculator': MUGEHydrogenAtomCalculator(),
    'MUGERingsOfRelativityCalculator': MUGERingsOfRelativityCalculator(),
    'MUGEMagnetarCalculator': MUGEMagnetarCalculator(),
    'MUGEGlobularClusterCalculator': MUGEGlobularClusterCalculator(),
    'MUGESagittariusAStarCalculator': MUGESagittariusAStarCalculator(),
    'MUGESunPlanetarySystemCalculator': MUGESunPlanetarySystemCalculator(),
    # Updated UQFF (114-121)
    'UQFFProtostellarJetCalculator': UQFFProtostellarJetCalculator(),
    'UQFFGalaxyMergerCalculator': UQFFGalaxyMergerCalculator(),
    'UQFFBlackHoleGrowthCalculator': UQFFBlackHoleGrowthCalculator(),
    'UQFFIssingAnyonCalculator': UQFFIssingAnyonCalculator(),
    'UQFFPolaritonQFTCalculator': UQFFPolaritonQFTCalculator(),
    'UQFFUTe2TopologicalCalculator': UQFFUTe2TopologicalCalculator(),
    'UQFFElectricUniverseRatioCalculator': UQFFElectricUniverseRatioCalculator(),
    'UQFFGyroTorqueNullificationCalculator': UQFFGyroTorqueNullificationCalculator(),
}


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE __all__ EXPORTS
# ═══════════════════════════════════════════════════════════════════════════════

__all__ = [
    # Constants
    'CONST', 'UQFF_CONST', 'COSMO',
    # Registry
    'GROK_URL_CALCULATORS',
    # Standard Physics (1-100)
    'AngularMomentumTransportCalculator', 'MHDJetVelocityCalculator',
    'JTypeShockRankineHugoniotCalculator', 'CTypeShockDampingCalculator',
    'EPSMergerRateCalculator', 'OrbitalTorqueTimeCalculator', 'SFRDEvolutionCalculator',
    'EPSBHMassFunctionCalculator', 'EddingtonAccretionCalculator',
    'SedovTaylorExpansionCalculator', 'DSAParticleAccelerationCalculator',
    'ChirpMassFormulaCalculator', 'QNMRingdownCalculator',
    'BlandfordZnajekPowerCalculator', 'RelativisticJetVelocityCalculator',
    'TOVEquationCalculator', 'PulsarSpinDownAgeCalculator', 'GlitchRecoveryCalculator',
    'FireballExpansionCalculator', 'AfterglowSynchrotronCalculator',
    'CMBAngularPowerCalculator', 'OpticalDepthCalculator',
    'MomentumFeedbackCalculator', 'BZJetPowerUpdatedCalculator', 'FeedbackDutyCycleCalculator',
    'PhotoevaporationRateCalculator', 'TypeIMigrationTorqueCalculator', 'RadialVelocitySemiAmplitudeCalculator',
    'NFWDensityProfileCalculator', 'RotationCurveCalculator', 'SIDMCoreFormationCalculator',
    'StrongLensingMassCalculator', 'XRayMassEstimateCalculator', 'MergerShockMachCalculator',
    'VoidDensityEvolutionCalculator', 'OutflowVelocityCalculator',
    'IonizationFractionCalculator', 'BubbleRadiusCalculator',
    'JeansLengthCalculator', 'AlfvenVelocityCalculator', 'TurbulentCascadeCalculator',
    'MainSequenceLifetimeCalculator', 'MassLuminosityRelationCalculator', 'ConvectiveTurnoverCalculator',
    'BaryonPhotonRatioCalculator', 'DeuteriumBottleneckCalculator',
    'FirstFriedmannCalculator', 'SecondFriedmannCalculator', 'DensityParameterCalculator',
    'SlowRollParametersCalculator', 'CurvaturePowerSpectrumCalculator', 'EFoldsCalculator',
    'TensorPowerSpectrumCalculator', 'StochasticGWDensityCalculator',
    'InspiralFrequencyEvolutionCalculator', 'MergerTimeCalculator', 'RingdownDampingTimeCalculator',
    'ArnettLawLightCurveCalculator', 'EjectaVelocityCalculator', 'NucleosynthesisYieldCalculator',
    'PNExpansionRadiusCalculator', 'IonizationFrontVelocityCalculator', 'AGBMassLossReimersCalculator',
    'ShockMachNumberTempCalculator', 'MergerTimescaleVirialCalculator', 'CoolCoreHeatingCalculator',
    'EvaporationRateCalculator', 'RelaxationTimeCalculator', 'VirialMassCalculator',
    'WindTerminalVelocityCalculator', 'IonizationParameterUCalculator', 'EnergyCouplingEfficiencyCalculator',
    'OrbitalDecayCalculator', 'PeriastronAdvanceCalculator', 'KilonovaPeakLuminosityCalculator',
    'FermiSecondOrderCalculator', 'KneeEnergyDSACalculator', 'DiffusionCoefficientCalculator',
    'WHIMTemperatureCalculator', 'MetalEnrichmentRateCalculator', 'CoolingTimeBremsstrahlungCalculator',
    'PressSchechterMassFunctionCalculator', 'SFEfficiencyCalculator', 'FeedbackEnergyInjectionCalculator',
    'CurvaturePerturbationAmplitudeCalculator', 'NonGaussianityFNLCalculator', 'ReheatingTemperatureCalculator',
    'SmallScaleDynamoGrowthCalculator', 'AlfvenMachNumberCalculator', 'FieldReversalScaleCalculator',
    'EquationOfStateWCalculator', 'CPLDensityEvolutionCalculator', 'GrowthSuppressionCalculator',
    'HawkingTemperatureCalculator', 'BekensteinHawkingEntropyCalculator', 'EvaporationLifetimeCalculator',
    'LQCEffectiveFriedmannCalculator', 'LQCPerturbationSpectrumCalculator', 'BounceTimescaleCalculator',
    'RocheLobeRadiusCalculator',
    # UQFF Framework (101-107)
    'UQFFBuoyancyCalculator', 'UQFFMagnetismCalculator', 'UQFFElectricFieldCalculator',
    'UQFFNeutronProductionCalculator', 'UQFFPseudoMonopoleCalculator',
    'UQFFVacuumDensityCalculator', 'UQFFGinzburgLandauCalculator',
    # MUGE Systems (108-113)
    'MUGEHydrogenAtomCalculator', 'MUGERingsOfRelativityCalculator',
    'MUGEMagnetarCalculator', 'MUGEGlobularClusterCalculator',
    'MUGESagittariusAStarCalculator', 'MUGESunPlanetarySystemCalculator',
    # Updated UQFF (114-121)
    'UQFFProtostellarJetCalculator', 'UQFFGalaxyMergerCalculator',
    'UQFFBlackHoleGrowthCalculator', 'UQFFIssingAnyonCalculator',
    'UQFFPolaritonQFTCalculator', 'UQFFUTe2TopologicalCalculator',
    'UQFFElectricUniverseRatioCalculator', 'UQFFGyroTorqueNullificationCalculator',
]
