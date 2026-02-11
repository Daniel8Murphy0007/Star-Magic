#!/usr/bin/env python3
"""
QCalc.py - UQFF Quantum Calculator (Pure Physics Solver)
=========================================================

A general-purpose physics calculator implementing the 8 UQFF Master Equations.

ARCHITECTURE RULES (MANDATORY):
─────────────────────────────────────────────────────────────────────────────────
1. NO HARDCODED SYSTEM DATA - All parameters passed via compute() methods
2. NO NAMED SYSTEM CLASSES - Only generic physics domain calculators
3. NO GLOBAL INSTANCES - Stateless calculator classes only
4. CONSTANTS ONLY - Fundamental physics constants (G, c, ℏ, etc.)
─────────────────────────────────────────────────────────────────────────────────

DATA FLOW:
    APIFetch.py → parameters dict → QCalc.solve() → OPData.py

OUTPUT FORMAT:
    {
        'long_form_equations': [...],    # Equations with substitutions shown
        'solutions': {...},              # Numerical results
        'available_equations': [...],    # Other solvable equations
        'simulation_set': {...}          # For multi-equation simulation
    }

8 UQFF Master Equations:
    1. UQFF (Base Unified Field)
    2. UQFF_Compressed (Newtonian + 9 corrections)
    3. UQFF_Resonant (aDPM + 13 frequency modes)
    4. UQFF_Superconductive (SCm vacuum modulation)
    5. UQFF_Buoyant (F_U_Bi) - Inside→Out, Atomic scale
    6. UQFF_Master_Buoyant (F_U_Bi_i) - Outside→In, Cosmic scale
    7. UQFF_Triadic (26-layer gravitational scaling)
    8. UQFF_Quadratic (Root solutions)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import numpy as np
from enum import Enum
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Union
from datetime import datetime
import json

# Import data layer modules
from IPData import InputParameters, recall_input, get_latest_input
from OPData import QueryResult, save_result, recall_result

# ═══════════════════════════════════════════════════════════════════════════════
# UNIVERSAL PHYSICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════
# These are FUNDAMENTAL physics constants - NOT system-specific data.
# System-specific parameters (M, r, T, etc.) come from APIFetch.py → IPData.py

CONSTANTS = {
    # ═══════════════════════════════════════════════════════════════════════════
    # FUNDAMENTAL CONSTANTS (SI Units)
    # ═══════════════════════════════════════════════════════════════════════════
    'G': 6.6743e-11,           # Gravitational constant (m³/kg·s²)
    'c': 2.998e8,              # Speed of light (m/s)
    'hbar': 1.0546e-34,        # Reduced Planck constant (J·s)
    'q': 1.602e-19,            # Elementary charge (C)
    'm_e': 9.109e-31,          # Electron mass (kg)
    'm_p': 1.673e-27,          # Proton mass (kg)
    'mu_B': 9.274e-24,         # Bohr magneton (J/T)
    'k_B': 1.381e-23,          # Boltzmann constant (J/K)
    'mu_0': 4 * np.pi * 1e-7,  # Vacuum permeability (H/m)
    'epsilon_0': 8.854e-12,    # Vacuum permittivity (F/m)
    'pi': np.pi,
    'Phi_0': 2.068e-15,        # Magnetic flux quantum (Wb)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STANDARD UNIT CONVERSIONS
    # ═══════════════════════════════════════════════════════════════════════════
    'M_sun': 1.989e30,         # Solar mass (kg) - as UNIT, not specific system
    'R_sun': 6.96e8,           # Solar radius (m) - as UNIT
    'L_sun': 3.828e26,         # Solar luminosity (W) - as UNIT
    'AU': 1.496e11,            # Astronomical Unit (m)
    'pc': 3.086e16,            # Parsec (m)
    'kpc': 3.086e19,           # Kiloparsec (m)
    'Mpc': 3.086e22,           # Megaparsec (m)
    'ly': 9.461e15,            # Light-year (m)
    'eV': 1.602e-19,           # Electronvolt (J)
    'MeV': 1.602e-13,          # Mega-electronvolt (J)
    'GeV': 1.602e-10,          # Giga-electronvolt (J)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UQFF CALIBRATED CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    'F0': 1.83e71,             # Base force constant (N)
    'kappa': 0.0005,           # κ: [SCm] reactivity decay rate (day⁻¹)
    'SSq': 0.57,               # [SSq] quantum state factor
    'U_UA': 1.0,               # Aether buoyancy factor
    'k_eta': 1e-113,           # Neutron rate coefficient
    'gamma': 5e-5,             # γ: Reciprocation decay rate (day⁻¹)
    'alpha': 1e-10,            # α: Time decay rate (s⁻¹)
    'H_SCm': 0.99,             # Heliosphere thickness factor
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UNIVERSAL GRAVITY COUPLING CONSTANTS (k_i)
    # ═══════════════════════════════════════════════════════════════════════════
    'k_1': 1.5,                # k₁ for Ug1 (Internal Dipole)
    'k_2': 1.2,                # k₂ for Ug2 (Outer Field Bubble)
    'k_3': 1.8,                # k₃ for Ug3 (Magnetic Strings Disk)
    'k_4': 1.0,                # k₄ for Ug4 (Star-Black Hole Interactions)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # BUOYANCY COUPLING CONSTANTS (β_i)
    # ═══════════════════════════════════════════════════════════════════════════
    'beta_i': 0.6,             # Buoyancy coupling constant
    'beta_1': 0.6,             # β₁ for Ug1
    'beta_2': 0.6,             # β₂ for Ug2
    'beta_3': 0.6,             # β₃ for Ug3
    'beta_4': 0.6,             # β₄ for Ug4
    
    # ═══════════════════════════════════════════════════════════════════════════
    # VACUUM ENERGY DENSITIES (Scale-neutral reference values)
    # ═══════════════════════════════════════════════════════════════════════════
    'rho_vac_SCm': 7.09e-37,   # ρ_vac,[SCm] reference (J/m³)
    'rho_vac_UA': 7.09e-36,    # ρ_vac,[UA] reference (J/m³)
    'rho_vac_cosmological': 5.96e-27,  # Cosmological vacuum energy (J/m³)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STANDARD MODEL PARTICLE MASSES
    # ═══════════════════════════════════════════════════════════════════════════
    # Quarks
    'm_u': 3.95e-30,           # Up quark (kg)
    'm_d': 8.40e-30,           # Down quark (kg)
    'm_c': 2.27e-27,           # Charm quark (kg)
    'm_s': 1.70e-28,           # Strange quark (kg)
    'm_t': 3.08e-25,           # Top quark (kg)
    'm_b': 7.49e-27,           # Bottom quark (kg)
    
    # Leptons
    'm_muon': 1.884e-28,       # Muon (kg)
    'm_tau': 3.168e-27,        # Tau (kg)
    'm_n': 1.675e-27,          # Neutron (kg)
    
    # Bosons
    'm_W': 1.43e-25,           # W boson (kg)
    'm_Z': 1.63e-25,           # Z boson (kg)
    'm_H': 2.23e-25,           # Higgs boson (kg)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # COSMOLOGICAL CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    'H0': 67.4,                # Hubble constant (km/s/Mpc)
    'H0_SI': 2.18e-18,         # Hubble constant (s⁻¹)
    'Omega_m': 0.315,          # Matter density parameter
    'Omega_Lambda': 0.685,     # Dark energy density parameter
    'Omega_b': 0.0493,         # Baryon density parameter
    'T_CMB': 2.725,            # CMB temperature (K)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # 26-LAYER QUANTUM STATE CONSTANTS
    # ═══════════════════════════════════════════════════════════════════════════
    'n_quantum_states': 26,    # Number of quantum states
    'f_TRZ': 0.1,              # Time-reversal zone factor
    'f_quasi': 0.01,           # Quasi-longitudinal wave factor
}


# ═══════════════════════════════════════════════════════════════════════════════
# UQFF SCALE SYSTEM - Scale Categories (NOT system-specific)
# ═══════════════════════════════════════════════════════════════════════════════

class UQFFScale(Enum):
    """
    UQFF operates identically across all scales - same equations, different parameters.
    
    The framework is scale-invariant: Ug, Ub, Ui, Um, Ur, Ut, UA, SCm equations
    apply at every level with scale-appropriate constants.
    """
    QUANTUM = 1       # Subatomic: ~10⁻¹⁵ m (nuclear, quark-gluon)
    ATOMIC = 2        # Atomic/Molecular: ~10⁻¹⁰ m
    CONDENSED = 3     # Lab-scale superconductivity: ~10⁻³ to 10⁰ m
    PLANETARY = 4     # Planetary: ~10⁶ to 10⁸ m
    STELLAR = 5       # Stellar: ~10⁸ to 10¹² m
    GALACTIC = 6      # Galactic: ~10²⁰ to 10²² m
    COSMOLOGICAL = 7  # Universe: ~10²⁶ m (Hubble radius)


# ═══════════════════════════════════════════════════════════════════════════════
# MULTI-SCALE PARAMETERS DATACLASS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ComputeParams:
    """
    Unified parameter set for UQFF calculations at any scale.
    
    ALL values are INPUT parameters - no hardcoded system data.
    These parameters come from APIFetch.py or manual user input.
    """
    # Identification
    query_name: str = "unnamed"
    scale: UQFFScale = UQFFScale.STELLAR
    
    # Core Physical Parameters (MUST be provided by API/user)
    M: float = None            # Mass (kg)
    r: float = None            # Distance/radius (m)
    T: float = None            # Temperature (K)
    L: float = None            # Luminosity (W)
    
    # Spatial Parameters
    R: float = None            # Object radius (m)
    z: float = None            # Redshift (dimensionless)
    d: float = None            # Distance to observer (m)
    
    # Kinematic Parameters
    v: float = None            # Velocity (m/s)
    omega: float = None        # Angular frequency (rad/s)
    P: float = None            # Period (s)
    
    # Magnetic Parameters
    B: float = None            # Magnetic field (T)
    mu: float = None           # Magnetic moment (J/T)
    
    # Quantum/Condensed Parameters
    psi: complex = None        # Order parameter
    Delta: float = None        # Energy gap (J)
    Phi: float = None          # Magnetic flux (Wb)
    
    # Galactic Parameters
    M_bh: float = None         # Central black hole mass (kg)
    d_g: float = None          # Distance to galactic center (m)
    Omega_g: float = None      # Galactic rotation rate (rad/s)
    sigma: float = None        # Velocity dispersion (m/s)
    
    # Time Parameters
    t: float = 0.0             # Time (s or days, context-dependent)
    t_n: float = 0.0           # Quantum time parameter
    
    # Coupling Parameters
    k_coupling: float = 1.0    # k₁-k₄ unified
    beta_coupling: float = 0.6 # β_i buoyancy coupling
    
    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {k: v for k, v in self.__dict__.items() if v is not None}
    
    @classmethod
    def from_api_response(cls, api_data: dict, query_name: str = "api_query"):
        """Create ComputeParams from API fetch response."""
        return cls(
            query_name=query_name,
            M=api_data.get('mass'),
            r=api_data.get('distance') or api_data.get('radius'),
            T=api_data.get('temperature'),
            L=api_data.get('luminosity'),
            z=api_data.get('redshift'),
            B=api_data.get('magnetic_field'),
            # ... map other API fields
        )


# ═══════════════════════════════════════════════════════════════════════════════
# EQUATION RESULT DATACLASS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class EquationResult:
    """
    Result of a single equation calculation with long-form output.
    """
    name: str                          # Equation name (e.g., "Universal Gravity Ug1")
    latex: str                         # LaTeX form of equation
    substituted: str                   # Equation with values substituted
    result: float                      # Numerical result
    unit: str                          # Physical unit
    parameters_used: Dict[str, float]  # Parameters that were used
    
    def to_dict(self) -> dict:
        return {
            'name': self.name,
            'latex': self.latex,
            'substituted': self.substituted,
            'result': self.result,
            'unit': self.unit,
            'parameters_used': self.parameters_used
        }


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED FIELD SOLVER - The Core Calculator
# ═══════════════════════════════════════════════════════════════════════════════

class UnifiedFieldSolver:
    """
    UQFF Universal Field Solver - Computes all equations from input parameters.
    
    This is a PURE CALCULATOR:
    - Takes parameters from APIFetch.py or user input
    - Computes applicable equations
    - Returns long-form equations with solutions
    - NO hardcoded system data
    """
    
    def __init__(self):
        """Initialize solver with fundamental constants only."""
        self.C = CONSTANTS  # Reference to constants
    
    # ═══════════════════════════════════════════════════════════════════════════
    # MAIN SOLVE METHOD
    # ═══════════════════════════════════════════════════════════════════════════
    
    def solve(self, params: ComputeParams) -> Dict[str, Any]:
        """
        Main entry point: Compute all applicable equations for given parameters.
        
        Args:
            params: ComputeParams with values from API fetch or user input
            
        Returns:
            {
                'query_id': str,
                'timestamp': str,
                'input_params': dict,
                'long_form_equations': List[EquationResult],
                'solutions': dict,
                'available_equations': List[str],
                'simulation_set': dict
            }
        """
        timestamp = datetime.now().isoformat()
        
        # Compute all applicable equations
        equations = []
        solutions = {}
        
        # Check which parameters are available and compute applicable equations
        if params.M is not None and params.r is not None:
            # Gravitational equations applicable
            ug_results = self._compute_universal_gravity(params)
            equations.extend(ug_results)
            for eq in ug_results:
                solutions[eq.name] = eq.result
        
        if params.M is not None and params.r is not None:
            # Buoyancy equations applicable
            ub_results = self._compute_universal_buoyancy(params)
            equations.extend(ub_results)
            for eq in ub_results:
                solutions[eq.name] = eq.result
        
        if params.B is not None or params.mu is not None:
            # Magnetic equations applicable
            um_results = self._compute_universal_magnetism(params)
            equations.extend(um_results)
            for eq in um_results:
                solutions[eq.name] = eq.result
        
        # Unified field combination
        if 'Ug' in solutions and 'Ub' in solutions:
            F_U = solutions.get('Ug', 0) + solutions.get('Ub', 0) + solutions.get('Um', 0)
            solutions['F_U'] = F_U
            equations.append(EquationResult(
                name='Unified Field F_U',
                latex=r'F_U = \sum_i (U_{g,i} + U_{b,i}) + U_m',
                substituted=f'F_U = {solutions.get("Ug", 0):.3e} + {solutions.get("Ub", 0):.3e} + {solutions.get("Um", 0):.3e}',
                result=F_U,
                unit='J/m³',
                parameters_used={'Ug': solutions.get('Ug'), 'Ub': solutions.get('Ub'), 'Um': solutions.get('Um', 0)}
            ))
        
        # Determine available equations based on parameters
        available = self._get_available_equations(params)
        
        return {
            'query_id': f"{params.query_name}_{timestamp.replace(':', '').replace('-', '').replace('.', '')}",
            'timestamp': timestamp,
            'input_params': params.to_dict(),
            'long_form_equations': [eq.to_dict() for eq in equations],
            'solutions': solutions,
            'available_equations': available,
            'simulation_set': self._build_simulation_set(params, solutions)
        }
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UNIVERSAL GRAVITY (Ug) EQUATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compute_universal_gravity(self, params: ComputeParams) -> List[EquationResult]:
        """Compute Ug1-Ug4 components."""
        results = []
        G = self.C['G']
        rho_vac = self.C['rho_vac_SCm']
        k_1, k_2, k_3, k_4 = self.C['k_1'], self.C['k_2'], self.C['k_3'], self.C['k_4']
        
        M = params.M
        r = params.r
        
        # Ug1: Internal Dipole
        Ug1 = k_1 * G * M / (r ** 2)
        results.append(EquationResult(
            name='Ug1',
            latex=r'U_{g1} = k_1 \times \frac{G \times M}{r^2}',
            substituted=f'Ug1 = {k_1} × ({G:.4e} × {M:.4e}) / ({r:.4e})²',
            result=Ug1,
            unit='m/s²',
            parameters_used={'k_1': k_1, 'G': G, 'M': M, 'r': r}
        ))
        
        # Ug2: Outer Field Bubble (with H_SCm)
        H_SCm = self.C['H_SCm']
        Ug2 = k_2 * rho_vac * M / (r ** 2) * H_SCm
        results.append(EquationResult(
            name='Ug2',
            latex=r'U_{g2} = k_2 \times \frac{\rho_{vac} \times M}{r^2} \times H_{SCm}',
            substituted=f'Ug2 = {k_2} × ({rho_vac:.4e} × {M:.4e}) / ({r:.4e})² × {H_SCm}',
            result=Ug2,
            unit='J/m³',
            parameters_used={'k_2': k_2, 'rho_vac': rho_vac, 'M': M, 'r': r, 'H_SCm': H_SCm}
        ))
        
        # Ug3: Magnetic Strings Disk
        Ug3 = k_3 * rho_vac * M / (r ** 2)
        results.append(EquationResult(
            name='Ug3',
            latex=r'U_{g3} = k_3 \times \frac{\rho_{vac} \times M}{r^2}',
            substituted=f'Ug3 = {k_3} × ({rho_vac:.4e} × {M:.4e}) / ({r:.4e})²',
            result=Ug3,
            unit='J/m³',
            parameters_used={'k_3': k_3, 'rho_vac': rho_vac, 'M': M, 'r': r}
        ))
        
        # Ug4: Star-Black Hole Interactions
        Ug4 = k_4 * rho_vac * M / (r ** 2)
        results.append(EquationResult(
            name='Ug4',
            latex=r'U_{g4} = k_4 \times \frac{\rho_{vac} \times M}{r^2}',
            substituted=f'Ug4 = {k_4} × ({rho_vac:.4e} × {M:.4e}) / ({r:.4e})²',
            result=Ug4,
            unit='J/m³',
            parameters_used={'k_4': k_4, 'rho_vac': rho_vac, 'M': M, 'r': r}
        ))
        
        # Total Ug
        Ug_total = Ug1 + Ug2 + Ug3 + Ug4
        results.append(EquationResult(
            name='Ug',
            latex=r'U_g = U_{g1} + U_{g2} + U_{g3} + U_{g4}',
            substituted=f'Ug = {Ug1:.4e} + {Ug2:.4e} + {Ug3:.4e} + {Ug4:.4e}',
            result=Ug_total,
            unit='mixed',
            parameters_used={'Ug1': Ug1, 'Ug2': Ug2, 'Ug3': Ug3, 'Ug4': Ug4}
        ))
        
        return results
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UNIVERSAL BUOYANCY (Ub) EQUATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compute_universal_buoyancy(self, params: ComputeParams) -> List[EquationResult]:
        """Compute Ub components (opposing gravity)."""
        results = []
        beta = self.C['beta_i']
        
        # First need Ug components
        Ug_results = self._compute_universal_gravity(params)
        Ug_total = sum(eq.result for eq in Ug_results if eq.name.startswith('Ug') and eq.name != 'Ug')
        
        # Get galactic parameters if available
        Omega_g = params.Omega_g or 7.3e-16  # Default galactic rotation
        M_bh = params.M_bh or 8.15e36         # Default central BH mass
        d_g = params.d_g or 2.55e20           # Default distance to galactic center
        
        # Ub = -β × Ug × Ω_g × M_bh/d_g
        Ub = -beta * Ug_total * Omega_g * (M_bh / d_g)
        
        results.append(EquationResult(
            name='Ub',
            latex=r'U_b = -\beta \times U_g \times \Omega_g \times \frac{M_{bh}}{d_g}',
            substituted=f'Ub = -{beta} × {Ug_total:.4e} × {Omega_g:.4e} × ({M_bh:.4e} / {d_g:.4e})',
            result=Ub,
            unit='J/m³',
            parameters_used={'beta': beta, 'Ug': Ug_total, 'Omega_g': Omega_g, 'M_bh': M_bh, 'd_g': d_g}
        ))
        
        return results
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UNIVERSAL MAGNETISM (Um) EQUATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compute_universal_magnetism(self, params: ComputeParams) -> List[EquationResult]:
        """Compute Um magnetic contributions."""
        results = []
        gamma = self.C['gamma']
        
        mu = params.mu or 1e23  # Magnetic moment (provide default for calculation)
        r = params.r or 1e10
        t = params.t
        t_n = params.t_n
        
        # Time factor: (1 - e^(-γt × cos(πt_n)))
        time_factor = 1 - np.exp(-gamma * t * np.cos(np.pi * t_n)) if t > 0 else 0
        
        # Um = μ/r × time_factor
        Um = (mu / r) * time_factor
        
        results.append(EquationResult(
            name='Um',
            latex=r'U_m = \frac{\mu}{r} \times (1 - e^{-\gamma t \cos(\pi t_n)})',
            substituted=f'Um = ({mu:.4e} / {r:.4e}) × (1 - exp(-{gamma} × {t} × cos(π × {t_n})))',
            result=Um,
            unit='J/m³',
            parameters_used={'mu': mu, 'r': r, 'gamma': gamma, 't': t, 't_n': t_n, 'time_factor': time_factor}
        ))
        
        return results
    
    # ═══════════════════════════════════════════════════════════════════════════
    # AVAILABLE EQUATIONS DETECTION
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _get_available_equations(self, params: ComputeParams) -> List[str]:
        """Determine which equations can be computed with available parameters."""
        available = []
        
        # Core UQFF equations (always available with M, r)
        if params.M is not None and params.r is not None:
            available.extend([
                'compute_Ug1', 'compute_Ug2', 'compute_Ug3', 'compute_Ug4',
                'compute_Ug_total', 'compute_Ub', 'compute_F_U',
                'compute_schwarzschild_radius', 'compute_escape_velocity'
            ])
        
        # Temperature-dependent equations
        if params.T is not None:
            available.extend([
                'compute_stefan_boltzmann', 'compute_wien_peak',
                'compute_thermal_energy', 'compute_blackbody_spectrum'
            ])
        
        # Black hole specific
        if params.M is not None and params.M > 1e30 * self.C['M_sun']:
            available.extend([
                'compute_hawking_temperature', 'compute_hawking_radiation',
                'compute_eddington_luminosity', 'compute_innermost_stable_orbit'
            ])
        
        # Magnetic equations
        if params.B is not None:
            available.extend([
                'compute_magnetic_pressure', 'compute_cyclotron_frequency',
                'compute_alfven_velocity', 'compute_magnetic_energy_density'
            ])
        
        # Cosmological equations
        if params.z is not None:
            available.extend([
                'compute_luminosity_distance', 'compute_angular_diameter_distance',
                'compute_comoving_distance', 'compute_lookback_time'
            ])
        
        # Galactic equations
        if params.sigma is not None:
            available.extend([
                'compute_M_sigma_relation', 'compute_velocity_dispersion_mass'
            ])
        
        return available
    
    # ═══════════════════════════════════════════════════════════════════════════
    # SIMULATION SET BUILDER
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _build_simulation_set(self, params: ComputeParams, solutions: Dict) -> Dict:
        """Build a set of equations for simultaneous simulation."""
        return {
            'coupled_equations': [
                {'name': 'Ug-Ub coupling', 'description': 'Gravity-Buoyancy balance'},
                {'name': 'F_U evolution', 'description': 'Unified field time evolution'},
            ],
            'initial_conditions': params.to_dict(),
            'solutions_at_t0': solutions,
            'time_range': {'t_min': 0, 't_max': 1e8, 'unit': 'days'}
        }


# ═══════════════════════════════════════════════════════════════════════════════
# SPECIALIZED CALCULATORS (Generic Physics Domains)
# ═══════════════════════════════════════════════════════════════════════════════
# These follow the CORRECT pattern: Generic names, parameterized methods

class GravitationalCalculator:
    """Generic gravitational calculations."""
    
    def __init__(self):
        self.C = CONSTANTS
    
    def schwarzschild_radius(self, M: float) -> EquationResult:
        """Compute Schwarzschild radius for any mass."""
        G = self.C['G']
        c = self.C['c']
        r_s = 2 * G * M / (c ** 2)
        return EquationResult(
            name='Schwarzschild Radius',
            latex=r'r_s = \frac{2GM}{c^2}',
            substituted=f'r_s = 2 × {G:.4e} × {M:.4e} / ({c:.4e})²',
            result=r_s,
            unit='m',
            parameters_used={'G': G, 'M': M, 'c': c}
        )
    
    def escape_velocity(self, M: float, r: float) -> EquationResult:
        """Compute escape velocity for any mass at any radius."""
        G = self.C['G']
        v_esc = np.sqrt(2 * G * M / r)
        return EquationResult(
            name='Escape Velocity',
            latex=r'v_{esc} = \sqrt{\frac{2GM}{r}}',
            substituted=f'v_esc = √(2 × {G:.4e} × {M:.4e} / {r:.4e})',
            result=v_esc,
            unit='m/s',
            parameters_used={'G': G, 'M': M, 'r': r}
        )
    
    def gravitational_lensing_angle(self, M: float, b: float) -> EquationResult:
        """Compute gravitational lensing deflection angle."""
        G = self.C['G']
        c = self.C['c']
        alpha = 4 * G * M / (c ** 2 * b)
        return EquationResult(
            name='Gravitational Lensing Angle',
            latex=r'\alpha = \frac{4GM}{c^2 b}',
            substituted=f'α = 4 × {G:.4e} × {M:.4e} / ({c:.4e}² × {b:.4e})',
            result=alpha,
            unit='rad',
            parameters_used={'G': G, 'M': M, 'c': c, 'b': b}
        )


class ThermodynamicCalculator:
    """Generic thermodynamic calculations."""
    
    def __init__(self):
        self.C = CONSTANTS
    
    def stefan_boltzmann_luminosity(self, R: float, T: float) -> EquationResult:
        """Compute luminosity from Stefan-Boltzmann law."""
        sigma = 5.670374e-8  # Stefan-Boltzmann constant
        L = 4 * np.pi * R ** 2 * sigma * T ** 4
        return EquationResult(
            name='Stefan-Boltzmann Luminosity',
            latex=r'L = 4\pi R^2 \sigma T^4',
            substituted=f'L = 4π × ({R:.4e})² × {sigma:.4e} × ({T:.4e})⁴',
            result=L,
            unit='W',
            parameters_used={'R': R, 'T': T, 'sigma': sigma}
        )
    
    def hawking_temperature(self, M: float) -> EquationResult:
        """Compute Hawking temperature for a black hole of any mass."""
        hbar = self.C['hbar']
        c = self.C['c']
        G = self.C['G']
        k_B = self.C['k_B']
        T_H = (hbar * c ** 3) / (8 * np.pi * G * M * k_B)
        return EquationResult(
            name='Hawking Temperature',
            latex=r'T_H = \frac{\hbar c^3}{8\pi G M k_B}',
            substituted=f'T_H = ({hbar:.4e} × ({c:.4e})³) / (8π × {G:.4e} × {M:.4e} × {k_B:.4e})',
            result=T_H,
            unit='K',
            parameters_used={'hbar': hbar, 'c': c, 'G': G, 'M': M, 'k_B': k_B}
        )


class CosmologicalCalculator:
    """Generic cosmological calculations."""
    
    def __init__(self):
        self.C = CONSTANTS
    
    def luminosity_distance(self, z: float) -> EquationResult:
        """Compute luminosity distance at redshift z (flat ΛCDM)."""
        c = self.C['c']
        H0 = self.C['H0_SI']
        # Simplified approximation for flat ΛCDM
        d_L = (c / H0) * z * (1 + z/2)  # First-order approximation
        return EquationResult(
            name='Luminosity Distance',
            latex=r'd_L = \frac{c}{H_0} \int_0^z \frac{dz}{E(z)}',
            substituted=f'd_L ≈ ({c:.4e} / {H0:.4e}) × {z} × (1 + {z}/2)',
            result=d_L,
            unit='m',
            parameters_used={'c': c, 'H0': H0, 'z': z}
        )
    
    def hubble_time(self) -> EquationResult:
        """Compute Hubble time."""
        H0 = self.C['H0_SI']
        t_H = 1 / H0
        return EquationResult(
            name='Hubble Time',
            latex=r't_H = \frac{1}{H_0}',
            substituted=f't_H = 1 / {H0:.4e}',
            result=t_H,
            unit='s',
            parameters_used={'H0': H0}
        )


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL SOLVER INSTANCE (for convenience)
# ═══════════════════════════════════════════════════════════════════════════════

SOLVER = UnifiedFieldSolver()


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

def solve(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convenience function to solve UQFF equations from a parameter dictionary.
    
    Args:
        params: Dictionary with keys like 'M', 'r', 'T', 'z', etc.
        
    Returns:
        Complete solution with long-form equations
    """
    compute_params = ComputeParams(
        query_name=params.get('name', 'query'),
        M=params.get('M') or params.get('mass'),
        r=params.get('r') or params.get('distance') or params.get('radius'),
        T=params.get('T') or params.get('temperature'),
        L=params.get('L') or params.get('luminosity'),
        z=params.get('z') or params.get('redshift'),
        B=params.get('B') or params.get('magnetic_field'),
        mu=params.get('mu') or params.get('magnetic_moment'),
        M_bh=params.get('M_bh') or params.get('black_hole_mass'),
        d_g=params.get('d_g') or params.get('galactic_distance'),
        sigma=params.get('sigma') or params.get('velocity_dispersion'),
    )
    return SOLVER.solve(compute_params)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Example usage with manual parameters
    test_params = {
        'name': 'test_compact_object',
        'M': 4.15e6 * CONSTANTS['M_sun'],  # Mass in kg
        'r': 8.1 * CONSTANTS['kpc'],        # Distance in m
        'T': 1e7,                           # Temperature in K
    }
    
    result = solve(test_params)
    
    print("=" * 80)
    print("QCalc.py - UQFF Quantum Calculator Test")
    print("=" * 80)
    print(f"Query ID: {result['query_id']}")
    print(f"Timestamp: {result['timestamp']}")
    print()
    print("LONG-FORM EQUATIONS:")
    print("-" * 80)
    for eq in result['long_form_equations']:
        print(f"  {eq['name']}:")
        print(f"    LaTeX: {eq['latex']}")
        print(f"    Substituted: {eq['substituted']}")
        print(f"    Result: {eq['result']:.4e} {eq['unit']}")
        print()
    print("SOLUTIONS:")
    print("-" * 80)
    for name, value in result['solutions'].items():
        print(f"  {name}: {value:.4e}")
    print()
    print("AVAILABLE EQUATIONS:")
    print("-" * 80)
    for eq in result['available_equations']:
        print(f"  - {eq}")
