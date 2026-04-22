#!/usr/bin/env python3
"""
QCalc_Advanced.py - Advanced Physics Extensions
===============================================

Advanced theoretical extensions beyond Phase 4:
- Wormhole metrics (Morris-Thorne, Einstein-Rosen)
- Higher-order General Relativity corrections
- Spatial variation of stress-energy tensor
- Full Christoffel symbols and geodesics
- Black hole thermodynamics with aether corrections
- Cosmological evolution with UQFF

DESIGN PRINCIPLES:
- Extends QCalc.py without modifying core
- Experimental/speculative features
- Research-grade implementations
- Full tensor calculus support

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable
import sys
import os

# Import core components
sys.path.insert(0, os.path.dirname(__file__))
from QCalc import CONSTANTS, ComputeParams, EquationResult


# ═══════════════════════════════════════════════════════════════════════════════
# WORMHOLE METRICS
# ═══════════════════════════════════════════════════════════════════════════════

class MorrisThorneWormhole:
    """
    Morris-Thorne traversable wormhole metric.
    
    Metric:
        ds² = -e^(2Φ(r)) dt² + (1 - b(r)/r)^(-1) dr² + r² dΩ²
    
    Where:
        Φ(r) = shape function (redshift)
        b(r) = throat function (radius)
    
    Requires exotic matter: ρ + P < 0
    
    Reference: Morris & Thorne, Am. J. Phys. 56, 395 (1988)
    """
    
    def __init__(self, b0: float = 1e3, rho_scale: float = 1e-23):
        """
        Args:
            b0: Throat radius (m)
            rho_scale: Exotic matter density scale (kg/m³)
        """
        self.b0 = b0
        self.rho_scale = rho_scale
        self.c = CONSTANTS['c']
        self.G = CONSTANTS['G']
    
    def shape_function(self, r: float) -> float:
        """
        Φ(r): Redshift function.
        
        Must satisfy: dΦ/dr → 0 as r → ∞ (asymptotically flat)
        
        Args:
            r: Radial coordinate (m)
            
        Returns:
            Φ(r): Redshift factor
        """
        # Simple form: Φ = 0 (zero tidal forces at throat)
        # Can be modified for specific scenarios
        return 0.0
    
    def throat_function(self, r: float) -> float:
        """
        b(r): Throat radius function.
        
        Must satisfy:
            b(b0) = b0 (throat condition)
            b(r) < r for r > b0 (no horizon)
            db/dr < 1 at r = b0 (flare-out condition)
        
        Args:
            r: Radial coordinate (m)
            
        Returns:
            b(r): Throat radius
        """
        # Ellis drainhole form
        return self.b0 / np.sqrt(1 + (r/self.b0)**2)
    
    def metric_components(self, r: float, theta: float = np.pi/2) -> Dict[str, float]:
        """
        Compute metric tensor components.
        
        Args:
            r: Radial coordinate (m)
            theta: Polar angle (equatorial plane: π/2)
            
        Returns:
            Dictionary with g_tt, g_rr, g_theta_theta, g_phi_phi
        """
        Phi = self.shape_function(r)
        b = self.throat_function(r)
        
        # Check validity (no horizon)
        if r <= b:
            raise ValueError(f"r = {r} m is inside throat radius b = {b} m!")
        
        g_tt = -np.exp(2 * Phi)
        g_rr = 1.0 / (1.0 - b / r)
        g_theta_theta = r**2
        g_phi_phi = r**2 * np.sin(theta)**2
        
        return {
            'g_tt': g_tt,
            'g_rr': g_rr,
            'g_theta_theta': g_theta_theta,
            'g_phi_phi': g_phi_phi,
            'Phi': Phi,
            'b': b
        }
    
    def exotic_matter_density(self, r: float) -> float:
        """
        Compute exotic matter density required to support wormhole.
        
        From Einstein equations:
            rho + P = -(b - r*b') / (8πG r² (1 - b/r)^1.5)
        
        Args:
            r: Radial coordinate (m)
            
        Returns:
            rho + P (kg/m³) - MUST BE NEGATIVE for exotic matter
        """
        b = self.throat_function(r)
        
        # Numerical derivative db/dr
        dr = 1e-6  # Small step
        b_plus = self.throat_function(r + dr)
        b_prime = (b_plus - b) / dr
        
        # Weak energy condition violation
        numerator = -(b - r * b_prime)
        denominator = 8 * np.pi * self.G * r**2 * (1 - b/r)**1.5
        
        rho_plus_P = numerator / denominator
        
        return rho_plus_P
    
    def is_traversable(self, r: float) -> bool:
        """
        Check traversability conditions at radius r.
        
        Conditions:
            1. r > b(r) (outside throat)
            2. db/dr < 1 at throat (flare-out)
            3. ρ + P < 0 (exotic matter present)
        
        Args:
            r: Radial coordinate (m)
            
        Returns:
            True if traversable at this radius
        """
        try:
            b = self.throat_function(r)
            
            # Condition 1: Outside throat
            if r <= b:
                return False
            
            # Condition 2: Flare-out (check numerically)
            dr = 1e-6
            b_plus = self.throat_function(r + dr)
            db_dr = (b_plus - b) / dr
            
            if db_dr >= 1.0:
                return False
            
            # Condition 3: Exotic matter
            rho_plus_P = self.exotic_matter_density(r)
            if rho_plus_P >= 0:
                return False
            
            return True
            
        except Exception:
            return False


# ═══════════════════════════════════════════════════════════════════════════════
# HIGHER-ORDER GR CORRECTIONS
# ═══════════════════════════════════════════════════════════════════════════════

class HigherOrderGR:
    """
    Second-order general relativity corrections.
    
    Perturbation expansion:
        g_μν = g_μν^(0) + ε g_μν^(1) + ε² g_μν^(2) + ...
    
    Where:
        g^(0) = Minkowski (flat)
        g^(1) = First-order (aether metric from Phase 4)
        g^(2) = Second-order (this class)
    """
    
    def __init__(self, eta: float = 1e-22):
        """
        Args:
            eta: Aether coupling constant (from Phase 4)
        """
        self.eta = eta
        self.eta2 = eta ** 2  # Second-order coupling
    
    def second_order_perturbation(
        self,
        T_s_mu_nu: np.ndarray,
        delta_g_mu_nu_1: np.ndarray
    ) -> np.ndarray:
        """
        Compute second-order metric perturbation.
        
        Formula:
            δ²g_μν ~ η² × (∂_α T_s^αβ)² + η² × (δg^(1))²
        
        Args:
            T_s_mu_nu: Stress-energy tensor (4×4)
            delta_g_mu_nu_1: First-order perturbation (4×4)
            
        Returns:
            δ²g_μν: Second-order correction (4×4)
        """
        # Term 1: Quadratic in stress-energy (simplified)
        T_squared = np.einsum('ij,ij->', T_s_mu_nu, T_s_mu_nu)
        term1 = self.eta2 * T_squared * np.eye(4)
        
        # Term 2: Quadratic in first-order perturbation
        delta_g_squared = np.einsum('ij,ij->', delta_g_mu_nu_1, delta_g_mu_nu_1)
        term2 = self.eta2 * delta_g_squared * np.eye(4)
        
        # Total second-order correction
        delta_g_2 = term1 + term2
        
        return delta_g_2
    
    def full_metric_to_second_order(
        self,
        g_0: np.ndarray,
        delta_g_1: np.ndarray,
        delta_g_2: np.ndarray
    ) -> np.ndarray:
        """
        Full metric including second-order corrections.
        
        g = g^(0) + δg^(1) + δ²g^(2)
        
        Args:
            g_0: Minkowski metric (4×4)
            delta_g_1: First-order perturbation (4×4)
            delta_g_2: Second-order perturbation (4×4)
            
        Returns:
            g_full: Complete metric to O(η²)
        """
        return g_0 + delta_g_1 + delta_g_2


# ═══════════════════════════════════════════════════════════════════════════════
# SPATIAL VARIATION OF STRESS-ENERGY
# ═══════════════════════════════════════════════════════════════════════════════

class SpatiallyVaryingStressEnergy:
    """
    Stress-energy tensor with spatial gradients: T_s^μν(x, t).
    
    Conserved: ∇_μ T_s^μν = 0
    
    Applications:
    - Cosmological evolution (scale factor a(t))
    - Matter distribution (density profiles)
    - Inhomogeneous vacuum energy
    """
    
    def __init__(self, lambda_UA_0: float = 7.09e-36, scale_radius: float = 1e26):
        """
        Args:
            lambda_UA_0: Central vacuum density (J/m³)
            scale_radius: Characteristic length scale (m)
        """
        self.lambda_UA_0 = lambda_UA_0
        self.R_scale = scale_radius
        self.c = CONSTANTS['c']
    
    def vacuum_density_profile(self, r: float, profile_type: str = 'exponential') -> float:
        """
        Vacuum energy density as function of radius.
        
        Profiles:
            'uniform': λ(r) = λ_0 (constant)
            'exponential': λ(r) = λ_0 × exp(-r/R_scale)
            'power_law': λ(r) = λ_0 × (1 + r/R_scale)^(-3)
            'cored': λ(r) = λ_0 × (1 + (r/R_scale)²)^(-1)
        
        Args:
            r: Radial coordinate (m)
            profile_type: Profile name
            
        Returns:
            λ(r): Vacuum density at radius r
        """
        if profile_type == 'uniform':
            return self.lambda_UA_0
        
        elif profile_type == 'exponential':
            return self.lambda_UA_0 * np.exp(-r / self.R_scale)
        
        elif profile_type == 'power_law':
            return self.lambda_UA_0 / (1 + r / self.R_scale)**3
        
        elif profile_type == 'cored':
            return self.lambda_UA_0 / (1 + (r / self.R_scale)**2)
        
        else:
            raise ValueError(f"Unknown profile type: {profile_type}")
    
    def gradient_vacuum_density(
        self,
        r: float,
        profile_type: str = 'exponential'
    ) -> float:
        """
        Radial gradient: dλ/dr
        
        Args:
            r: Radial coordinate (m)
            profile_type: Profile name
            
        Returns:
            dλ/dr: Radial derivative
        """
        if profile_type == 'uniform':
            return 0.0
        
        elif profile_type == 'exponential':
            return -self.lambda_UA_0 / self.R_scale * np.exp(-r / self.R_scale)
        
        elif profile_type == 'power_law':
            return -3 * self.lambda_UA_0 / (self.R_scale * (1 + r / self.R_scale)**4)
        
        elif profile_type == 'cored':
            return -2 * self.lambda_UA_0 * r / (self.R_scale**2 * (1 + (r / self.R_scale)**2)**2)
        
        else:
            raise ValueError(f"Unknown profile type: {profile_type}")
    
    def stress_energy_tensor_spatial(
        self,
        r: float,
        profile_type: str = 'exponential'
    ) -> np.ndarray:
        """
        Stress-energy tensor with spatial variation.
        
        Perfect fluid form:
            T^00 = ρ(r) c²
            T^ii = -P(r) = -ρ(r) c² / 3
        
        Args:
            r: Radial coordinate (m)
            profile_type: Density profile type
            
        Returns:
            T_s^μν(r): 4×4 tensor
        """
        lambda_r = self.vacuum_density_profile(r, profile_type)
        rho = lambda_r / self.c**2  # Energy density → mass density
        P = -rho / 3.0  # Relativistic pressure
        
        T_s = np.zeros((4, 4))
        T_s[0, 0] = rho * self.c**2  # Energy density
        T_s[1, 1] = -P               # Radial pressure
        T_s[2, 2] = -P               # Tangential pressure
        T_s[3, 3] = -P               # Tangential pressure
        
        return T_s
    
    def check_conservation(
        self,
        r: float,
        dr: float = 1e6,
        profile_type: str = 'exponential'
    ) -> float:
        """
        Check energy-momentum conservation: ∇_μ T_s^μν ≈ 0
        
        In spherical symmetry:
            ∇_μ T^μr = ∂_r (r² T^rr) / r² + T^θθ / r + ...
        
        Args:
            r: Radial coordinate (m)
            dr: Step size for numerical derivative (m)
            profile_type: Profile type
            
        Returns:
            |∇_μ T^μν|: Conservation violation magnitude
        """
        # Compute T at r and r+dr
        T_r = self.stress_energy_tensor_spatial(r, profile_type)
        T_r_plus = self.stress_energy_tensor_spatial(r + dr, profile_type)
        
        # Radial component of divergence (approximation)
        div_T = (T_r_plus[1, 1] - T_r[1, 1]) / dr + 2 * T_r[1, 1] / r
        
        return abs(div_T)


# ═══════════════════════════════════════════════════════════════════════════════
# FULL CHRISTOFFEL SYMBOLS
# ═══════════════════════════════════════════════════════════════════════════════

class ChristoffelCalculator:
    """
    Compute full Christoffel symbols from metric.
    
    Formula:
        Γ^λ_μν = ½ g^λσ (∂_μ g_σν + ∂_ν g_μσ - ∂_σ g_μν)
    
    Applications:
    - Geodesic equations (particle trajectories)
    - Parallel transport
    - Curvature tensors (Riemann, Ricci, Weyl)
    """
    
    def __init__(self):
        pass
    
    def numerical_derivative_metric(
        self,
        metric_func: Callable[[np.ndarray], np.ndarray],
        x: np.ndarray,
        mu: int,
        dx: float = 1e-6
    ) -> np.ndarray:
        """
        Numerical derivative of metric: ∂_μ g_αβ
        
        Args:
            metric_func: Function(x) -> g_αβ (4×4)
            x: Position 4-vector [t, x, y, z]
            mu: Index for derivative (0=t, 1=x, 2=y, 3=z)
            dx: Step size
            
        Returns:
            ∂_μ g_αβ: 4×4 array
        """
        x_plus = x.copy()
        x_plus[mu] += dx
        
        g_plus = metric_func(x_plus)
        g = metric_func(x)
        
        return (g_plus - g) / dx
    
    def christoffel_symbols(
        self,
        metric_func: Callable[[np.ndarray], np.ndarray],
        x: np.ndarray,
        dx: float = 1e-6
    ) -> np.ndarray:
        """
        Compute all Christoffel symbols Γ^λ_μν
        
        Args:
            metric_func: Function(x) -> g_αβ (4×4)
            x: Position 4-vector [t, x, y, z]
            dx: Step size for numerical derivatives
            
        Returns:
            Γ: Array of shape (4, 4, 4) containing Γ^λ_μν
        """
        # Get metric and inverse at position x
        g = metric_func(x)
        g_inv = np.linalg.inv(g)
        
        # Compute all metric derivatives
        dg = np.zeros((4, 4, 4))  # dg[μ, α, β] = ∂_μ g_αβ
        for mu in range(4):
            dg[mu] = self.numerical_derivative_metric(metric_func, x, mu, dx)
        
        # Compute Christoffel symbols
        Gamma = np.zeros((4, 4, 4))  # Gamma[λ, μ, ν] = Γ^λ_μν
        
        for lam in range(4):
            for mu in range(4):
                for nu in range(4):
                    # Γ^λ_μν = ½ g^λσ (∂_μ g_σν + ∂_ν g_μσ - ∂_σ g_μν)
                    sum_term = 0.0
                    for sigma in range(4):
                        sum_term += g_inv[lam, sigma] * (
                            dg[mu, sigma, nu] + dg[nu, mu, sigma] - dg[sigma, mu, nu]
                        )
                    Gamma[lam, mu, nu] = 0.5 * sum_term
        
        return Gamma
    
    def geodesic_equation_rhs(
        self,
        x: np.ndarray,
        v: np.ndarray,
        Gamma: np.ndarray
    ) -> np.ndarray:
        """
        Right-hand side of geodesic equation:
            d²x^λ/dτ² = -Γ^λ_μν dx^μ/dτ dx^ν/dτ
        
        Args:
            x: Position 4-vector [t, x, y, z]
            v: Velocity 4-vector [dt/dτ, dx/dτ, dy/dτ, dz/dτ]
            Gamma: Christoffel symbols (4, 4, 4)
            
        Returns:
            a: Acceleration 4-vector d²x^λ/dτ²
        """
        a = np.zeros(4)
        
        for lam in range(4):
            for mu in range(4):
                for nu in range(4):
                    a[lam] -= Gamma[lam, mu, nu] * v[mu] * v[nu]
        
        return a


# ═══════════════════════════════════════════════════════════════════════════════
# BLACK HOLE THERMODYNAMICS WITH AETHER
# ═══════════════════════════════════════════════════════════════════════════════

class AetherBlackHoleThermodynamics:
    """
    Black hole thermodynamics with aether corrections.
    
    Hawking temperature:
        T_H = ℏ c³ / (8π k_B G M)
    
    With aether modification:
        T_H' = T_H × (1 + α × η × T_s)
    
    Where α is dimensionless coupling (~O(1))
    """
    
    def __init__(self, alpha: float = 1.0, eta: float = 1e-22):
        """
        Args:
            alpha: Aether-temperature coupling (dimensionless)
            eta: Aether coupling constant
        """
        self.alpha = alpha
        self.eta = eta
        self.c = CONSTANTS['c']
        self.G = CONSTANTS['G']
        self.h_bar = CONSTANTS['hbar']
        self.k_B = 1.380649e-23  # Boltzmann constant (J/K)
    
    def hawking_temperature(self, M: float) -> float:
        """
        Standard Hawking temperature.
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            T_H: Temperature (K)
        """
        return self.h_bar * self.c**3 / (8 * np.pi * self.k_B * self.G * M)
    
    def aether_corrected_temperature(
        self,
        M: float,
        T_s_00: float
    ) -> float:
        """
        Hawking temperature with aether correction.
        
        Args:
            M: Black hole mass (kg)
            T_s_00: Stress-energy component (kg/m³ c²)
            
        Returns:
            T_H': Modified temperature (K)
        """
        T_H = self.hawking_temperature(M)
        correction = 1.0 + self.alpha * self.eta * T_s_00
        return T_H * correction
    
    def evaporation_time(self, M: float) -> float:
        """
        Black hole evaporation time (standard).
        
        τ = (5120 π G² M³) / (ℏ c⁴)
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            τ: Evaporation time (s)
        """
        return (5120 * np.pi * self.G**2 * M**3) / (self.h_bar * self.c**4)
    
    def schwarzschild_radius(self, M: float) -> float:
        """
        Schwarzschild radius: r_s = 2GM/c²
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            r_s: Event horizon radius (m)
        """
        return 2 * self.G * M / self.c**2


# ═══════════════════════════════════════════════════════════════════════════════
# COSMOLOGICAL EVOLUTION
# ═══════════════════════════════════════════════════════════════════════════════

class CosmologicalEvolution:
    """
    Cosmological expansion with UQFF vacuum energy.
    
    Friedmann equations with aether:
        H² = (8πG/3) × (ρ_m + ρ_aether) - k/a² + Λ/3
        
    Where:
        H = ȧ/a (Hubble parameter)
        ρ_aether = T_s^00 from Phase 4
        Λ = Cosmological constant
    """
    
    def __init__(self, H0: float = 67.4, Omega_m: float = 0.315):
        """
        Args:
            H0: Hubble constant (km/s/Mpc)
            Omega_m: Matter density parameter (current)
        """
        self.H0_km = H0
        self.H0_SI = H0 * 1e3 / CONSTANTS['Mpc']  # Convert to 1/s
        self.Omega_m = Omega_m
        self.Omega_Lambda = 1.0 - Omega_m  # Flat universe
        self.c = CONSTANTS['c']
        self.G = CONSTANTS['G']
        
        # Critical density
        self.rho_crit_0 = 3 * self.H0_SI**2 / (8 * np.pi * self.G)
    
    def hubble_parameter(self, z: float, Omega_aether: float = 0.0) -> float:
        """
        Hubble parameter as function of redshift.
        
        H(z) = H0 × sqrt[Ω_m (1+z)³ + Ω_Λ + Ω_aether (1+z)^n]
        
        Where n depends on aether equation of state: P = w ρ
        For w = -1/3 (radiation-like): n = 4
        For w = -1 (dark energy-like): n = 0
        
        Args:
            z: Redshift
            Omega_aether: Aether density parameter
            
        Returns:
            H(z): Hubble parameter (SI units: 1/s)
        """
        # Standard ΛCDM
        E_z_squared = (
            self.Omega_m * (1 + z)**3 +
            self.Omega_Lambda +
            Omega_aether * (1 + z)**4  # Assuming w = -1/3
        )
        
        return self.H0_SI * np.sqrt(E_z_squared)
    
    def comoving_distance(self, z: float, Omega_aether: float = 0.0, n_points: int = 1000) -> float:
        """
        Comoving distance to redshift z.
        
        d_c = c ∫_0^z dz' / H(z')
        
        Args:
            z: Redshift
            Omega_aether: Aether density parameter
            n_points: Integration points
            
        Returns:
            d_c: Comoving distance (m)
        """
        z_arr = np.linspace(0, z, n_points)
        H_arr = np.array([self.hubble_parameter(z_i, Omega_aether) for z_i in z_arr])
        
        # Trapezoidal integration
        integrand = 1.0 / H_arr
        d_c = self.c * np.trapz(integrand, z_arr)
        
        return d_c


# ═══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLES
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("="*80)
    print("QCalc Advanced Extensions - Test Suite")
    print("="*80)
    
    # Test 1: Morris-Thorne Wormhole
    print("\n[TEST 1] Morris-Thorne Wormhole")
    wormhole = MorrisThorneWormhole(b0=1e3, rho_scale=1e-23)
    
    r_test = 2e3  # 2× throat radius
    metric = wormhole.metric_components(r_test)
    rho_plus_P = wormhole.exotic_matter_density(r_test)
    traversable = wormhole.is_traversable(r_test)
    
    print(f"  Throat radius: {wormhole.b0:.2e} m")
    print(f"  Test radius: {r_test:.2e} m")
    print(f"  g_tt = {metric['g_tt']:.6f}")
    print(f"  g_rr = {metric['g_rr']:.6f}")
    print(f"  ρ + P = {rho_plus_P:.4e} kg/m³ ({'EXOTIC' if rho_plus_P < 0 else 'NORMAL'})")
    print(f"  Traversable: {'YES' if traversable else 'NO'}")
    
    # Test 2: Higher-Order GR
    print("\n[TEST 2] Higher-Order GR Corrections")
    ho_gr = HigherOrderGR(eta=1e-22)
    
    # Mock stress-energy and first-order perturbation
    T_s = np.eye(4) * 1e8  # Simplified
    delta_g_1 = np.eye(4) * 1e-14
    
    delta_g_2 = ho_gr.second_order_perturbation(T_s, delta_g_1)
    
    print(f"  η = {ho_gr.eta:.2e}")
    print(f"  η² = {ho_gr.eta2:.2e}")
    print(f"  |δg^(1)| ~ {np.max(np.abs(delta_g_1)):.2e}")
    print(f"  |δ²g^(2)| ~ {np.max(np.abs(delta_g_2)):.2e}")
    print(f"  Ratio: |δ²g|/|δg| = {np.max(np.abs(delta_g_2))/np.max(np.abs(delta_g_1)):.2e}")
    
    # Test 3: Spatially Varying Stress-Energy
    print("\n[TEST 3] Spatially Varying Stress-Energy")
    spatial_T = SpatiallyVaryingStressEnergy(lambda_UA_0=7.09e-36, scale_radius=1e26)
    
    radii = [1e25, 1e26, 1e27]  # m
    for r in radii:
        lambda_r = spatial_T.vacuum_density_profile(r, 'exponential')
        dlambda_dr = spatial_T.gradient_vacuum_density(r, 'exponential')
        T_s_r = spatial_T.stress_energy_tensor_spatial(r, 'exponential')
        
        print(f"\n  r = {r:.2e} m:")
        print(f"    λ(r) = {lambda_r:.4e} J/m³")
        print(f"    dλ/dr = {dlambda_dr:.4e} J/m⁴")
        print(f"    T^00 = {T_s_r[0,0]:.4e} kg/m³ c²")
    
    # Test 4: Black Hole Thermodynamics
    print("\n[TEST 4] Black Hole Thermodynamics with Aether")
    bh_thermo = AetherBlackHoleThermodynamics(alpha=1.0, eta=1e-22)
    
    M_bh = 10 * CONSTANTS['M_sun']  # 10 solar masses
    T_s_00 = 1.0975e8  # From Phase 4
    
    T_H_standard = bh_thermo.hawking_temperature(M_bh)
    T_H_aether = bh_thermo.aether_corrected_temperature(M_bh, T_s_00)
    tau_evap = bh_thermo.evaporation_time(M_bh)
    r_s = bh_thermo.schwarzschild_radius(M_bh)
    
    print(f"  M_BH = {M_bh/CONSTANTS['M_sun']:.1f} M_sun")
    print(f"  r_s = {r_s:.2e} m ({r_s/1e3:.2f} km)")
    print(f"  T_H (standard) = {T_H_standard:.4e} K")
    print(f"  T_H (aether) = {T_H_aether:.4e} K")
    print(f"  Correction: {(T_H_aether/T_H_standard - 1)*100:.2e}%")
    print(f"  τ_evap = {tau_evap:.4e} s ({tau_evap/(365.25*86400):.4e} years)")
    
    # Test 5: Cosmological Evolution
    print("\n[TEST 5] Cosmological Evolution with Aether")
    cosmo = CosmologicalEvolution(H0=67.4, Omega_m=0.315)
    
    print(f"  H0 = {cosmo.H0_km:.2f} km/s/Mpc")
    print(f"  Ω_m = {cosmo.Omega_m:.3f}")
    print(f"  Ω_Λ = {cosmo.Omega_Lambda:.3f}")
    print(f"  ρ_crit,0 = {cosmo.rho_crit_0:.4e} kg/m³")
    
    # Compute for different redshifts
    redshifts = [0.0, 0.5, 1.0, 2.0, 5.0]
    Omega_aether = 1e-10  # Small aether contribution
    
    print(f"\n  With Ω_aether = {Omega_aether:.2e}:")
    for z in redshifts:
        H_z = cosmo.hubble_parameter(z, Omega_aether)
        d_c = cosmo.comoving_distance(z, Omega_aether)
        
        print(f"    z = {z:.1f}: H(z) = {H_z*1e3/CONSTANTS['Mpc']:.2f} km/s/Mpc, "
              f"d_c = {d_c/CONSTANTS['Mpc']:.1f} Mpc")
    
    print("\n" + "="*80)
    print("All advanced extensions tests complete!")
    print("="*80)
