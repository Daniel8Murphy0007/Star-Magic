"""
updated_uqff_2025_module.py
UQFF 2025 Extensions Module

Updated UQFF framework incorporating 2025 physics insights:
- Ising Anyons (non-Abelian statistics)
- Polariton QFT (curved spacetime analogs)
- UTe2 Nodal Superconductivity
- Enhanced F_UBii with error-resistant stabilization
- Updated Um for polariton systems

Based on Grok AI conversation (August 2025) with references from:
- Universal Buoyancy_08April2025.docx
- Universal Magnetism_17Mar2025.docx
- Universal Quantum Framework_01May2025.docx

© 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple, List, Callable
from dataclasses import dataclass, field
from enum import Enum
import math
from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell

from abc import ABC, abstractmethod


# ═══════════════════════════════════════════════════════════════════════════════
# PHYSICAL CONSTANTS (UQFF 2025 Enhanced)
# ═══════════════════════════════════════════════════════════════════════════════

CONSTANTS_2025 = {
    # Fundamental
    'G': 6.6743e-11,           # Gravitational constant (m³/kg·s²)
    'c': 2.998e8,              # Speed of light (m/s)
    'hbar': 1.0546e-34,        # Reduced Planck constant (J·s)
    'k_B': 1.381e-23,          # Boltzmann constant (J/K)
    'e': 1.602e-19,            # Elementary charge (C)
    'mu_0': 4 * np.pi * 1e-7,  # Vacuum permeability (H/m)
    'epsilon_0': 8.854e-12,    # Vacuum permittivity (F/m)
    
    # UQFF Calibration (2025 refined)
    'kappa_UQFF': 0.0005,      # κ decay rate (1/day)
    'SSq': 0.57,               # Superconductive state factor
    'H_SCm': 0.99,             # SCm Higgs coupling
    'U_UA': 0.0001,            # UA coupling
    'k_eta': 1e-113,           # Neutron production calibration
    'beta_i': 0.603,           # Buoyancy index
    'rho_vac_UA': 7.09e-36,    # Vacuum density UA (J/m³)
    'rho_vac_SCm': 7.09e-37,   # Vacuum density SCm (J/m³)
    'F_rel': 4.30e33,          # Relativistic coherence force (N)
    
    # Anyon Parameters (2025)
    'theta_fibonacci': 3 * np.pi / 5,  # Fibonacci anyon topological angle
    'theta_ising': np.pi / 8,          # Ising anyon angle
    'd_fibonacci': (1 + np.sqrt(5)) / 2,  # Golden ratio (quantum dimension)
    'd_ising': np.sqrt(2),             # Ising anyon quantum dimension
    
    # Polariton Parameters (2025)
    'Omega_Rabi_typical': 1e12,        # Typical Rabi splitting (rad/s)
    'gamma_polariton': 1e10,           # Polariton decay rate (1/s)
    'n_cavity': 1e4,                   # Photons per cavity mode
    'g_coupling': 1e11,                # Light-matter coupling (rad/s)
    
    # UTe2 Parameters (2025)
    'T_c_UTe2': 1.6,                   # Critical temperature (K)
    'Delta_0_UTe2': 0.18e-3 * 1.602e-19,  # Gap magnitude (J), ~0.18 meV
    'mu_B': 9.274e-24,                 # Bohr magneton (J/T)
    'H_c2_UTe2': 60,                   # Upper critical field (T)
}


# ═══════════════════════════════════════════════════════════════════════════════
# ANYON FRAMEWORK (2025)
# ═══════════════════════════════════════════════════════════════════════════════

class AnyonType(Enum):
    """Enumeration of anyon types for topological computation."""
    ABELIAN = "abelian"
    FIBONACCI = "fibonacci"
    ISING = "ising"
    MAJORANA = "majorana"


@dataclass
class AnyonState:
    """Representation of an anyon quantum state."""
    type: AnyonType
    fusion_channel: str
    topological_charge: float
    braid_phase: complex
    error_rate: float


class IsingAnyonCalculator:
    """
    Ising Anyon Calculator for UQFF Integration.
    
    Ising anyons are non-Abelian anyons emerging from:
    - Majorana zero modes in topological superconductors
    - ν = 5/2 fractional quantum Hall states
    - Kitaev chain edge states
    
    Key properties:
    - Quantum dimension d = √2
    - Topological angle θ = π/8
    - Fusion rules: σ × σ = 1 + ψ (non-deterministic)
    
    UQFF extension:
    F_UBii,anyons = F_UBii × (1 + d_ising × exp(-Γ_error × t))
    """
    
    def __init__(self):
        self.constants = CONSTANTS_2025.copy()
        self.d_ising = self.constants['d_ising']
        self.theta_ising = self.constants['theta_ising']
        
    def compute_fusion_matrix(self) -> np.ndarray:
        """
        Compute Ising anyon F-matrix (Fusion matrix).
        
        The 6j-symbol structure for σ × σ × σ → σ fusion.
        """
        # Ising F-matrix (2×2 for σσσ→σ channels)
        F = np.array([
            [1/np.sqrt(2), 1/np.sqrt(2)],
            [1/np.sqrt(2), -1/np.sqrt(2)]
        ])
        return F
    
    def compute_braid_matrix(self) -> np.ndarray:
        """
        Compute Ising anyon R-matrix (Braiding matrix).
        
        R^{σσ}_ψ encodes the phase acquired under σ-σ exchange.
        """
        # Ising R-matrix
        R = np.array([
            [np.exp(-1j * np.pi / 8), 0],
            [0, np.exp(3j * np.pi / 8)]
        ])
        return R
    
    def compute_topological_protection(self, gap: float, T: float) -> float:
        """
        Compute topological protection factor.
        
        P_prot = exp(-Δ / k_B T) gives thermal excitation probability.
        True topological protection requires Δ >> k_B T.
        
        Args:
            gap: Energy gap Δ (J)
            T: Temperature (K)
            
        Returns:
            Protection factor (0-1, higher is better)
        """
        k_B = self.constants['k_B']
        if T == 0:
            return 1.0
        ratio = gap / (k_B * T)
        # Protection factor: 1 when gap >> k_B T
        return 1 - np.exp(-ratio)
    
    def compute_error_rate(self, t: float, T: float, gap: float,
                           quasiparticle_poisoning: float = 1e-6) -> float:
        """
        Compute quantum error rate for Ising anyon qubit.
        
        Γ_error = Γ_thermal + Γ_poisoning + Γ_dephasing
        
        Args:
            t: Time evolution (s)
            T: Temperature (K)
            gap: Energy gap (J)
            quasiparticle_poisoning: QP poisoning rate (1/s)
            
        Returns:
            Total error rate
        """
        k_B = self.constants['k_B']
        hbar = self.constants['hbar']
        
        # Thermal error rate
        Gamma_thermal = gap / hbar * np.exp(-gap / (k_B * T)) if T > 0 else 0
        
        # Quasiparticle poisoning
        Gamma_qp = quasiparticle_poisoning
        
        # Dephasing from environment (simplified)
        Gamma_dephasing = 1e-8  # Low for topological systems
        
        return Gamma_thermal + Gamma_qp + Gamma_dephasing
    
    def F_UBii_anyons(self, F_UBii_base: float, t: float, T: float,
                      gap: float) -> float:
        """
        Extended UQFF buoyancy with Ising anyon stabilization.
        
        F_UBii,anyons = F_UBii × (1 + d_ising × exp(-Γ_error × t))
        
        This represents error-resistant stabilization through
        topological protection.
        
        Args:
            F_UBii_base: Base buoyancy force (N)
            t: Time (s)
            T: Temperature (K)
            gap: Topological gap (J)
            
        Returns:
            Enhanced buoyancy force
        """
        Gamma_error = self.compute_error_rate(t, T, gap)
        stabilization = 1 + self.d_ising * np.exp(-Gamma_error * t)
        
        return F_UBii_base * stabilization
    
    def compute_braiding_unitary(self, n_exchanges: int) -> np.ndarray:
        """
        Compute the unitary from n successive σ-σ exchanges.
        
        U_braid = R^n
        
        Args:
            n_exchanges: Number of exchanges
            
        Returns:
            Unitary transformation matrix
        """
        R = self.compute_braid_matrix()
        return np.linalg.matrix_power(R, n_exchanges)
    
    def compute_quantum_dimension(self) -> float:
        """Return quantum dimension d_σ = √2."""
        return self.d_ising
    
    def to_dict(self) -> Dict[str, Any]:
        """Export calculator state to dictionary."""
        return {
            'type': 'IsingAnyonCalculator',
            'd_ising': self.d_ising,
            'theta_ising': self.theta_ising,
            'F_matrix': self.compute_fusion_matrix().tolist(),
            'R_matrix': [[str(c) for c in row] for row in self.compute_braid_matrix()]
        }


# ═══════════════════════════════════════════════════════════════════════════════
# POLARITON QFT FRAMEWORK (2025)
# ═══════════════════════════════════════════════════════════════════════════════

class PolaritonQFTCalculator:
    """
    Polariton QFT Calculator for curved spacetime analogs.
    
    Exciton-polaritons in microcavities provide:
    - Analog curved spacetime (effective metric)
    - Hawking-like radiation
    - Superfluid behavior
    - Nonlinear quantum optics
    
    UQFF extension:
    Um,polariton = Um × (1 + g|ψ|² / Ω_Rabi) × curved spacetime factor
    
    The effective metric emerges from the polariton dispersion:
    g_μν ∝ ∂²ω/∂k_μ∂k_ν
    """
    
    def __init__(self):
        self.constants = CONSTANTS_2025.copy()
        self.Omega_Rabi = self.constants['Omega_Rabi_typical']
        self.gamma = self.constants['gamma_polariton']
        self.g_coupling = self.constants['g_coupling']
        
    def compute_polariton_dispersion(self, k: np.ndarray,
                                      cavity_detuning: float = 0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute upper and lower polariton dispersion.
        
        E_± = (E_C + E_X)/2 ± √[(E_C - E_X)²/4 + (ħΩ_Rabi/2)²]
        
        Args:
            k: Wavevector array (1/m)
            cavity_detuning: δ = E_C - E_X at k=0 (J)
            
        Returns:
            (E_LP, E_UP): Lower and upper polariton energies
        """
        hbar = self.constants['hbar']
        c = self.constants['c']
        
        # Cavity photon dispersion (parabolic)
        m_photon = 1e-35  # Effective cavity photon mass (kg)
        E_C = hbar * c * np.abs(k) + hbar**2 * k**2 / (2 * m_photon) + cavity_detuning / 2
        
        # Exciton energy (nominally flat)
        E_X = np.zeros_like(k) - cavity_detuning / 2
        
        # Polariton energies
        term1 = (E_C + E_X) / 2
        term2 = np.sqrt((E_C - E_X)**2 / 4 + (hbar * self.Omega_Rabi / 2)**2)
        
        E_LP = term1 - term2  # Lower polariton
        E_UP = term1 + term2  # Upper polariton
        
        return E_LP, E_UP
    
    def compute_hopfield_coefficients(self, E_C: float, E_X: float) -> Tuple[float, float]:
        """
        Compute Hopfield coefficients (exciton and photon weights).
        
        |X|² + |C|² = 1
        
        Args:
            E_C: Cavity energy (J)
            E_X: Exciton energy (J)
            
        Returns:
            (|X|², |C|²): Exciton and cavity photon fractions
        """
        hbar = self.constants['hbar']
        delta = E_C - E_X
        Omega = hbar * self.Omega_Rabi
        
        denominator = np.sqrt(delta**2 + Omega**2)
        
        X_sq = 0.5 * (1 + delta / denominator)  # Exciton fraction
        C_sq = 0.5 * (1 - delta / denominator)  # Photon fraction
        
        return X_sq, C_sq
    
    def compute_effective_metric(self, psi: complex, velocity_field: np.ndarray) -> np.ndarray:
        """
        Compute effective spacetime metric from polariton superfluid.
        
        For a flowing superfluid with velocity v, the acoustic metric is:
        g_μν = (ρ/c_s) × [-(c_s² - v²), -v_i; -v_j, δ_ij]
        
        This allows simulation of:
        - Black hole horizons (where v = c_s)
        - Hawking radiation analogs
        - de Sitter spacetime
        
        Args:
            psi: Polariton wavefunction
            velocity_field: 3D velocity field (m/s)
            
        Returns:
            4×4 effective metric tensor
        """
        # Polariton density and speed of sound
        rho = np.abs(psi)**2
        g_int = 1e-23  # Interaction strength (J·m²)
        m_eff = 1e-35  # Effective mass (kg)
        c_s = np.sqrt(g_int * rho / m_eff)  # Speed of sound
        
        v = np.linalg.norm(velocity_field)
        
        # Construct acoustic metric (3+1 dimensions)
        g = np.zeros((4, 4))
        g[0, 0] = -(c_s**2 - v**2)  # g_tt
        
        # Space-time mixing
        for i in range(3):
            g[0, i+1] = -velocity_field[i]
            g[i+1, 0] = -velocity_field[i]
            
        # Spatial components
        for i in range(3):
            g[i+1, i+1] = 1.0
            
        # Normalization factor
        if rho > 0 and c_s > 0:
            g *= rho / c_s
            
        return g
    
    def compute_hawking_temperature(self, surface_gravity: float) -> float:
        """
        Compute analog Hawking temperature.
        
        T_H = ħκ / (2π k_B c_s)
        
        where κ is the surface gravity at the acoustic horizon.
        
        Args:
            surface_gravity: κ (m/s²)
            
        Returns:
            Hawking temperature (K)
        """
        hbar = self.constants['hbar']
        k_B = self.constants['k_B']
        c_s = 1e4  # Typical sound speed in polariton superfluid (m/s)
        
        return hbar * surface_gravity / (2 * np.pi * k_B * c_s)
    
    def Um_polariton(self, Um_base: float, psi: complex,
                     curvature: float = 0) -> float:
        """
        Extended UQFF magnetism with polariton corrections.
        
        Um,polariton = Um × (1 + g|ψ|² / Ω_Rabi) × (1 + R_curv)
        
        where:
        - g|ψ|² gives polariton-polariton interaction energy
        - R_curv is the effective spacetime curvature
        
        Args:
            Um_base: Base magnetism term
            psi: Polariton wavefunction
            curvature: Effective Ricci scalar (1/m²)
            
        Returns:
            Enhanced magnetism term
        """
        # Interaction enhancement
        g_int = 1e-23  # Interaction constant
        interaction_factor = 1 + g_int * np.abs(psi)**2 / (self.constants['hbar'] * self.Omega_Rabi)
        
        # Curvature correction
        curvature_factor = 1 + curvature * 1e-20  # Normalized
        
        return Um_base * interaction_factor * curvature_factor
    
    def compute_parametric_amplification(self, pump_power: float,
                                          signal_detuning: float) -> float:
        """
        Compute parametric down-conversion gain.
        
        For polariton parametric oscillator (OPO).
        
        Args:
            pump_power: Input pump power (W)
            signal_detuning: Signal-pump detuning (rad/s)
            
        Returns:
            Amplification gain
        """
        # Simplified parametric gain
        threshold_power = 1e-3  # mW typical
        gain = np.sqrt(pump_power / threshold_power) if pump_power < threshold_power * 100 else 10
        
        # Detuning-dependent phase matching
        bandwidth = self.Omega_Rabi / 10
        detuning_factor = np.exp(-signal_detuning**2 / bandwidth**2)
        
        return gain * detuning_factor
    
    def to_dict(self) -> Dict[str, Any]:
        """Export calculator state to dictionary."""
        return {
            'type': 'PolaritonQFTCalculator',
            'Omega_Rabi': self.Omega_Rabi,
            'gamma': self.gamma,
            'g_coupling': self.g_coupling
        }


# ═══════════════════════════════════════════════════════════════════════════════
# UTe2 SUPERCONDUCTIVITY FRAMEWORK (2025)
# ═══════════════════════════════════════════════════════════════════════════════

class UTe2SuperconductivityCalculator:
    """
    UTe2 Nodal Superconductivity Calculator.
    
    UTe2 is a spin-triplet topological superconductor with:
    - T_c ≈ 1.6 K (recently discovered 2019)
    - Spin-triplet pairing (likely p-wave)
    - Nodal gap structure
    - Re-entrant superconductivity at high fields
    - Multiple phases under pressure
    
    UQFF relevance:
    - Provides laboratory for SCm (superconductive mass modulation)
    - Demonstrates H_SCm coupling at high fields
    - Validates UQFF predictions for exotic pairing
    
    Extension:
    [SSq]_UTe2 = 0.57 × (1 - T/T_c) × f_nodal(θ) × (1 + H/H_c2)^α
    """
    
    def __init__(self):
        self.constants = CONSTANTS_2025.copy()
        self.T_c = self.constants['T_c_UTe2']
        self.Delta_0 = self.constants['Delta_0_UTe2']
        self.H_c2 = self.constants['H_c2_UTe2']
        
    def compute_gap_function(self, theta: float, phi: float, T: float) -> float:
        """
        Compute UTe2 nodal gap function.
        
        For triplet p-wave pairing:
        Δ(k) = Δ_0 × d(k) × (1 - (T/T_c)^2)
        
        where d(k) is the d-vector characterizing spin-triplet pairing.
        
        Args:
            theta: Polar angle on Fermi surface
            phi: Azimuthal angle on Fermi surface
            T: Temperature (K)
            
        Returns:
            Gap magnitude (J)
        """
        if T >= self.T_c:
            return 0.0
            
        # BCS-like temperature dependence
        temp_factor = np.sqrt(1 - (T / self.T_c)**2)
        
        # Nodal structure (simplified A_u representation)
        # d(k) ≈ z × sin(k_z) for typical UTe2 models
        nodal_factor = np.abs(np.sin(theta) * np.cos(phi))
        
        return self.Delta_0 * temp_factor * nodal_factor
    
    def compute_d_vector(self, k: np.ndarray) -> np.ndarray:
        """
        Compute spin-triplet d-vector.
        
        The d-vector d = (d_x, d_y, d_z) characterizes the pairing:
        Δ = (d · σ)(iσ_y) where σ are Pauli matrices.
        
        For UTe2 in A_u representation:
        d(k) ≈ z × sin(k_z c) where c is the c-axis lattice parameter.
        
        Args:
            k: Momentum vector (k_x, k_y, k_z) in 1/Å
            
        Returns:
            d-vector as numpy array
        """
        c_lattice = 4.16e-10  # c-axis lattice parameter (m)
        k_z = k[2] if len(k) > 2 else 0
        
        d_z = np.sin(k_z * c_lattice)
        
        return np.array([0, 0, d_z])
    
    def compute_specific_heat(self, T: float, gamma_n: float = 0.1) -> float:
        """
        Compute electronic specific heat C_e(T).
        
        For nodal superconductors:
        C/T ~ T^2 at low T (nodal quasiparticles)
        
        Args:
            T: Temperature (K)
            gamma_n: Normal state Sommerfeld coefficient (J/mol·K²)
            
        Returns:
            C/T (J/mol·K²)
        """
        if T >= self.T_c:
            return gamma_n
            
        # BCS jump at T_c
        jump_ratio = 1.43  # Weak coupling value
        
        # Nodal suppression (power law)
        x = T / self.T_c
        if x < 0.1:
            return gamma_n * jump_ratio * (x**2)  # C/T ~ T² for nodes
        else:
            return gamma_n * jump_ratio * np.exp(-self.Delta_0 / (self.constants['k_B'] * T))
    
    def compute_penetration_depth(self, T: float) -> float:
        """
        Compute London penetration depth λ(T).
        
        For nodal superconductors:
        λ(T) - λ(0) ~ T at low T
        
        Args:
            T: Temperature (K)
            
        Returns:
            Penetration depth (m)
        """
        lambda_0 = 5e-7  # λ(0) ~ 500 nm typical
        
        if T >= self.T_c:
            return float('inf')
            
        x = T / self.T_c
        
        # Nodal: Δλ/λ_0 ~ T/T_c
        delta_lambda = lambda_0 * x
        
        return lambda_0 + delta_lambda
    
    def compute_upper_critical_field(self, T: float, angle: float = 0) -> float:
        """
        Compute upper critical field H_c2(T).
        
        UTe2 shows highly anisotropic H_c2 with re-entrance.
        
        Args:
            T: Temperature (K)
            angle: Field angle from c-axis (rad)
            
        Returns:
            H_c2 (T)
        """
        if T >= self.T_c:
            return 0.0
            
        # Basic Ginzburg-Landau behavior
        H_c2_0 = self.H_c2
        H_c2_T = H_c2_0 * (1 - (T / self.T_c)**2)
        
        # Angular anisotropy (simplified)
        gamma_aniso = 3.0  # H_c2 anisotropy ratio
        angular_factor = np.sqrt(np.cos(angle)**2 + np.sin(angle)**2 / gamma_aniso**2)
        
        return H_c2_T / angular_factor
    
    def compute_SSq_UTe2(self, T: float, H: float, theta: float = 0) -> float:
        """
        Compute UQFF superconductive state factor for UTe2.
        
        [SSq]_UTe2 = 0.57 × (1 - T/T_c) × f_nodal(θ) × (1 + H/H_c2)^α
        
        The re-entrant behavior at high fields is encoded in the
        non-monotonic H dependence.
        
        Args:
            T: Temperature (K)
            H: Applied magnetic field (T)
            theta: Gap node angle
            
        Returns:
            [SSq] factor for UTe2
        """
        SSq_base = self.constants['SSq']  # 0.57
        
        if T >= self.T_c:
            return 0.0
            
        # Temperature factor
        temp_factor = 1 - T / self.T_c
        
        # Nodal structure
        f_nodal = np.abs(np.sin(2 * theta))  # Nodes at θ = 0, π/2, ...
        
        # Field factor with re-entrance
        H_c2 = self.compute_upper_critical_field(T)
        alpha = 0.3  # Re-entrance exponent
        
        if H < H_c2:
            field_factor = (1 - H / H_c2)**(1 - alpha)
        else:
            # Re-entrant regime
            field_factor = 0.1 * np.exp(-(H - H_c2) / 10)
            
        return SSq_base * temp_factor * (1 + 0.1 * f_nodal) * field_factor
    
    def compute_Knight_shift(self, T: float, direction: str = 'c') -> float:
        """
        Compute Knight shift reduction below T_c.
        
        For spin-triplet superconductors:
        - No reduction for d ⊥ H
        - Full reduction for d || H
        
        Args:
            T: Temperature (K)
            direction: Field direction ('a', 'b', or 'c')
            
        Returns:
            Knight shift reduction factor (0 to 1)
        """
        if T >= self.T_c:
            return 1.0
            
        # d-vector along c-axis for A_u representation
        # K_c shows reduction, K_ab doesn't
        if direction == 'c':
            reduction = 1 - 0.8 * (1 - T / self.T_c)
        else:
            reduction = 1.0  # No reduction for H ⊥ d
            
        return reduction
    
    def to_dict(self) -> Dict[str, Any]:
        """Export calculator state to dictionary."""
        return {
            'type': 'UTe2SuperconductivityCalculator',
            'T_c': self.T_c,
            'Delta_0': self.Delta_0,
            'H_c2': self.H_c2
        }


# ═══════════════════════════════════════════════════════════════════════════════
# ENHANCED UQFF 2025 UNIFIED CALCULATOR
# ═══════════════════════════════════════════════════════════════════════════════

class UQFF2025Calculator:
    """
    Unified UQFF 2025 Calculator.
    
    Integrates all 2025 extensions:
    - Ising anyon error-resistant buoyancy
    - Polariton curved spacetime magnetism
    - UTe2 superconductive state modulation
    
    Master equation (2025):
    F_U(r,t) = Σ_i [Ug1_i + Ug2_i + Ug3_i + Ug4_i] × [SSq]
              + F_UBii,anyons + Um,polariton + [SSq]_UTe2 corrections
    """
    
    def __init__(self):
        self.constants = CONSTANTS_2025.copy()
        self.ising = IsingAnyonCalculator()
        self.polariton = PolaritonQFTCalculator()
        self.ute2 = UTe2SuperconductivityCalculator()
        
    def compute_F_U_2025(self, r: float, t: float, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compute full UQFF 2025 unified field.
        
        Args:
            r: Distance (m)
            t: Time (s)
            params: Dictionary with:
                - M: Mass (kg)
                - mu: Magnetic moment (A·m²)
                - Q: Charge (C)
                - T: Temperature (K)
                - H: Magnetic field (T)
                - psi_polariton: Polariton wavefunction (optional)
                - gap_anyon: Anyon topological gap (J, optional)
                
        Returns:
            Dictionary with all field components
        """
        # Extract parameters
        M = params.get('M', 1e30)
        mu = params.get('mu', 1e20)
        Q = params.get('Q', 1e10)
        T = params.get('T', 1.0)
        H = params.get('H', 0.0)
        psi_polariton = params.get('psi_polariton', 0.0)
        gap_anyon = params.get('gap_anyon', 1e-21)
        
        G = self.constants['G']
        c = self.constants['c']
        
        # Base UQFF components
        Ug1 = G * mu / r**3  # Magnetic dipole
        Ug2 = G * Q**2 / r**4  # Charge-reactivity
        Ug3 = 0.01 * dpm_ug1_seed(M, r) * np.sin(2 * np.pi * t / 86400)  # String rotation
        Ug4 = c**2 * self.constants['rho_vac_UA'] / r  # Vacuum concentration
        
        # Base buoyancy
        delta_rho = 1e-10  # Density differential
        V = (4/3) * np.pi * r**3
        g_local = dpm_ug1_seed(M, r)
        F_UBii_base = -delta_rho * V * g_local
        
        # Base magnetism
        gamma_Um = 5e-5 / 86400
        Um_base = (mu / r) * (1 - np.exp(-gamma_Um * t))
        
        # === 2025 EXTENSIONS ===
        
        # Ising anyon enhanced buoyancy
        F_UBii_anyons = self.ising.F_UBii_anyons(F_UBii_base, t, T, gap_anyon)
        
        # Polariton enhanced magnetism
        Um_polariton = self.polariton.Um_polariton(Um_base, psi_polariton)
        
        # UTe2 superconductive correction
        SSq_UTe2 = self.ute2.compute_SSq_UTe2(T, H)
        
        # Total unified field
        F_U_total = (Ug1 + Ug2 + Ug3 + Ug4) * (1 + SSq_UTe2) + F_UBii_anyons
        
        return {
            'F_U_total': F_U_total,
            'components': {
                'Ug1': Ug1,
                'Ug2': Ug2,
                'Ug3': Ug3,
                'Ug4': Ug4,
                'F_UBii_base': F_UBii_base,
                'F_UBii_anyons': F_UBii_anyons,
                'Um_base': Um_base,
                'Um_polariton': Um_polariton,
            },
            'uqff_2025': {
                'SSq_UTe2': SSq_UTe2,
                'anyon_stabilization': F_UBii_anyons / F_UBii_base if F_UBii_base != 0 else 1,
                'polariton_enhancement': Um_polariton / Um_base if Um_base != 0 else 1,
            },
            'parameters': params,
            'equations': [
                r"F_U(r,t) = \sum_i [Ug_{1-4,i}] \times [SSq] + F_{UBii,anyons} + U_{m,polariton}",
                r"F_{UBii,anyons} = F_{UBii} \times (1 + d_{ising} e^{-\Gamma_{error} t})",
                r"U_{m,polariton} = U_m \times (1 + g|\psi|^2 / \Omega_{Rabi})",
                r"[SSq]_{UTe2} = 0.57 \times (1-T/T_c) \times f_{nodal} \times (1+H/H_{c2})^\alpha"
            ]
        }
    
    def compute_topological_quantum_number(self, winding: int, filling: float) -> float:
        """
        Compute topological invariant for UQFF system.
        
        ν = n × σ_xy / (e²/h)
        
        where n is the winding number and σ_xy is Hall conductance.
        
        Args:
            winding: Winding number (integer)
            filling: Fractional filling factor (e.g., 0.5 for ν = 5/2)
            
        Returns:
            Topological quantum number
        """
        h = 2 * np.pi * self.constants['hbar']
        e = self.constants['e']
        
        sigma_xy = filling * e**2 / h  # Hall conductance
        
        return winding * sigma_xy * h / e**2
    
    def to_dict(self) -> Dict[str, Any]:
        """Export full calculator state."""
        return {
            'type': 'UQFF2025Calculator',
            'constants': self.constants,
            'ising': self.ising.to_dict(),
            'polariton': self.polariton.to_dict(),
            'ute2': self.ute2.to_dict()
        }


# ═══════════════════════════════════════════════════════════════════════════════
# REGISTRY EXPORT
# ═══════════════════════════════════════════════════════════════════════════════

# Global calculator instances
ISING_ANYON = IsingAnyonCalculator()
POLARITON_QFT = PolaritonQFTCalculator()
UTE2_SC = UTe2SuperconductivityCalculator()
UQFF_2025 = UQFF2025Calculator()

# Calculator registry
UQFF_2025_CALCULATORS = {
    'IsingAnyonCalculator': ISING_ANYON,
    'PolaritonQFTCalculator': POLARITON_QFT,
    'UTe2SuperconductivityCalculator': UTE2_SC,
    'UQFF2025Calculator': UQFF_2025,
}


# ═══════════════════════════════════════════════════════════════════════════════
# TEST HARNESS
# ═══════════════════════════════════════════════════════════════════════════════

def test_uqff_2025():
    """Test all UQFF 2025 calculators."""
    print("=" * 80)
    print("UQFF 2025 Extensions Test Suite")
    print("=" * 80)
    
    # Test 1: Ising Anyons
    print("\n1. ISING ANYON CALCULATOR")
    print("-" * 40)
    F_matrix = ISING_ANYON.compute_fusion_matrix()
    R_matrix = ISING_ANYON.compute_braid_matrix()
    d_ising = ISING_ANYON.compute_quantum_dimension()
    
    print(f"   Quantum dimension d_σ = {d_ising:.4f} (√2 = {np.sqrt(2):.4f})")
    print(f"   F-matrix:\n{F_matrix}")
    
    # Test anyon-enhanced buoyancy
    F_base = 1e-10  # Base buoyancy force
    T = 0.1  # 100 mK
    gap = 1e-21  # Topological gap
    F_anyons = ISING_ANYON.F_UBii_anyons(F_base, t=1e-3, T=T, gap=gap)
    print(f"   F_UBii enhancement factor: {F_anyons/F_base:.4f}")
    
    # Test 2: Polariton QFT
    print("\n2. POLARITON QFT CALCULATOR")
    print("-" * 40)
    k = np.linspace(0, 1e7, 100)  # 0 to 10^7 m^-1
    E_LP, E_UP = POLARITON_QFT.compute_polariton_dispersion(k)
    
    print(f"   Rabi splitting: {POLARITON_QFT.Omega_Rabi:.2e} rad/s ({POLARITON_QFT.Omega_Rabi/2e12:.1f} meV)")
    print(f"   E_LP(k=0): {E_LP[0]:.2e} J")
    print(f"   E_UP(k=0): {E_UP[0]:.2e} J")
    
    # Test effective metric
    psi = 1e5 + 0j  # Polariton density
    v = np.array([1e4, 0, 0])  # Flow velocity
    metric = POLARITON_QFT.compute_effective_metric(psi, v)
    print(f"   Effective metric g_00: {metric[0,0]:.2e}")
    
    # Test polariton-enhanced magnetism
    Um_base = 1e-20
    Um_pol = POLARITON_QFT.Um_polariton(Um_base, psi)
    print(f"   Um enhancement factor: {Um_pol/Um_base:.4f}")
    
    # Test 3: UTe2 Superconductivity
    print("\n3. UTe2 SUPERCONDUCTIVITY CALCULATOR")
    print("-" * 40)
    
    print(f"   T_c = {UTE2_SC.T_c:.2f} K")
    print(f"   Δ_0 = {UTE2_SC.Delta_0/1.602e-19*1000:.3f} meV")
    print(f"   H_c2(T=0) = {UTE2_SC.H_c2:.1f} T")
    
    # Test gap function
    gap_node = UTE2_SC.compute_gap_function(theta=0, phi=0, T=0.8)
    gap_max = UTE2_SC.compute_gap_function(theta=np.pi/4, phi=np.pi/4, T=0.8)
    print(f"   Δ(node, T=0.8K) = {gap_node/1.602e-19*1000:.4f} meV")
    print(f"   Δ(antinode, T=0.8K) = {gap_max/1.602e-19*1000:.4f} meV")
    
    # Test SSq_UTe2
    SSq = UTE2_SC.compute_SSq_UTe2(T=0.8, H=5.0, theta=np.pi/4)
    print(f"   [SSq]_UTe2(T=0.8K, H=5T) = {SSq:.4f}")
    
    # Test 4: Unified UQFF 2025
    print("\n4. UNIFIED UQFF 2025 CALCULATOR")
    print("-" * 40)
    
    result = UQFF_2025.compute_F_U_2025(
        r=1e6,  # 1000 km
        t=1000,  # 1000 s
        params={
            'M': 1e24,  # Earth-like
            'mu': 8e22,
            'Q': 1e8,
            'T': 0.5,
            'H': 10,
            'psi_polariton': 1e4,
            'gap_anyon': 1e-21
        }
    )
    
    print(f"   F_U_total = {result['F_U_total']:.6e}")
    print(f"   Ug1 (magnetic dipole) = {result['components']['Ug1']:.6e}")
    print(f"   Ug4 (vacuum) = {result['components']['Ug4']:.6e}")
    print(f"   Anyon stabilization = {result['uqff_2025']['anyon_stabilization']:.4f}")
    print(f"   Polariton enhancement = {result['uqff_2025']['polariton_enhancement']:.4f}")
    print(f"   [SSq]_UTe2 = {result['uqff_2025']['SSq_UTe2']:.4f}")
    
    print("\n" + "=" * 80)
    print("All UQFF 2025 tests completed successfully!")
    print("=" * 80)


if __name__ == "__main__":
    test_uqff_2025()
