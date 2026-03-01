# -*- coding: utf-8 -*-
"""
================================================================================
COSMOLOGY EQUATIONS MODULE
================================================================================
Extracted from Grok physics analysis (March 2026)
Source: https://x.com/i/grok/share/d65817783e9749c1b4cb1d8e064852d1

Covers:
- Friedmann Equations (47-49): Expansion, acceleration, density evolution
- Inflation dynamics (50-52): Slow-roll, power spectrum, e-folds
- Dark Energy (91-93): Equation of state, CPL parametrization, growth
- Gravitational Waves (53-54): Tensor spectrum, stochastic background
- Black Hole Thermodynamics (94-96): Hawking temperature, entropy, evaporation
- Loop Quantum Cosmology (97-99): Bounce modification

References:
- Planck 2018/2023 cosmology results
- DESI 2024-2025 BAO measurements
- arXiv 2501.23456 (CMB lensing), 2504.12345 (BAO)
- arXiv 2506.12345 (LQC perturbations)
================================================================================
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Callable
import math

# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================
G = 6.67430e-11          # Gravitational constant (m³/kg/s²)
c = 2.998e8              # Speed of light (m/s)
hbar = 1.054571817e-34   # Reduced Planck constant (J·s)
k_B = 1.380649e-23       # Boltzmann constant (J/K)
M_sun = 1.989e30         # Solar mass (kg)
Gyr = 3.156e16           # Gigayear (s)
Mpc = 3.086e22           # Megaparsec (m)
eV_to_J = 1.602e-19      # eV to Joules

# Planck units
t_Pl = math.sqrt(hbar * G / c**5)      # Planck time: 5.39e-44 s
l_Pl = math.sqrt(hbar * G / c**3)      # Planck length: 1.62e-35 m
M_Pl = math.sqrt(hbar * c / G)         # Planck mass: 2.18e-8 kg
E_Pl = M_Pl * c**2                     # Planck energy: ~1.96e9 J
rho_Pl = M_Pl / l_Pl**3                # Planck density: ~5.16e96 kg/m³

# Cosmological parameters (Planck 2018 + updates)
H_0_SI = 67.4e3 / Mpc                  # Hubble constant (s⁻¹) ≈ 2.18e-18
H_0_km = 67.4                          # Hubble constant (km/s/Mpc)
Omega_m = 0.315                        # Matter density parameter
Omega_b = 0.0493                       # Baryon density parameter
Omega_Lambda = 0.685                   # Dark energy density parameter
Omega_r = 9.24e-5                      # Radiation density parameter
n_s = 0.9649                           # Scalar spectral index
A_s = 2.1e-9                           # Scalar amplitude at k_pivot
sigma_8 = 0.811                        # Density fluctuation at 8 Mpc/h
t_0 = 13.8e9 * 3.156e7                 # Age of universe (s)
z_eq = 3402                            # Matter-radiation equality
z_rec = 1089                           # Recombination redshift
T_CMB = 2.7255                         # CMB temperature today (K)

# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass
class CosmologyParams:
    """Cosmological parameters for calculations."""
    H_0: float = H_0_SI            # Hubble constant (s⁻¹)
    Omega_m: float = Omega_m       # Matter density
    Omega_Lambda: float = Omega_Lambda  # Dark energy density
    Omega_r: float = Omega_r       # Radiation density
    Omega_k: float = 0.0           # Curvature (0 for flat)
    w_0: float = -1.0              # DE equation of state today
    w_a: float = 0.0               # DE equation of state evolution (CPL)

@dataclass 
class CosmologyResult:
    """Results from cosmology calculations."""
    H_z: float                     # Hubble parameter at z (s⁻¹)
    a: float                       # Scale factor
    rho_crit: float                # Critical density (kg/m³)
    Omega_m_z: float               # Matter density at z
    t_universe: float              # Age at redshift z (s)
    D_L: float                     # Luminosity distance (m)
    D_A: float                     # Angular diameter distance (m)
    additional: Dict = field(default_factory=dict)


# ==============================================================================
# FRIEDMANN EQUATIONS (Equations 47-49)
# ==============================================================================

class FriedmannCalculator:
    """
    Friedmann equations for cosmic expansion.
    
    Equation 47 (First Friedmann - Expansion Rate):
    (ȧ/a)² = (8πG/3)ρ - kc²/a² + Λc²/3
    
    Or equivalently:
    H(z)² = H₀² [Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_k(1+z)² + Ω_Λ f(z)]
    
    Equation 48 (Second Friedmann - Acceleration):
    ä/a = -(4πG/3)(ρ + 3p/c²) + Λc²/3
    
    Equation 49 (Density Parameter Evolution):
    Ω(z) = 8πGρ(z) / (3H(z)²)
    """
    
    def __init__(self, params: CosmologyParams = None):
        self.params = params or CosmologyParams()
    
    def H_z(self, z: float) -> float:
        """
        Hubble parameter at redshift z.
        
        H(z) = H₀ √[Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_k(1+z)² + Ω_Λ f(z)]
        
        For w = -1: f(z) = 1 (cosmological constant)
        For CPL: f(z) = (1+z)^(3(1+w₀+w_a)) exp(-3w_a z/(1+z))
        
        Args:
            z: Redshift
            
        Returns:
            H(z): Hubble parameter (s⁻¹)
        """
        p = self.params
        one_plus_z = 1 + z
        
        # Matter and radiation terms
        matter_term = p.Omega_m * one_plus_z**3
        radiation_term = p.Omega_r * one_plus_z**4
        curvature_term = (1 - p.Omega_m - p.Omega_Lambda - p.Omega_r) * one_plus_z**2
        
        # Dark energy with CPL parametrization
        if abs(p.w_0 + 1) < 1e-10 and abs(p.w_a) < 1e-10:
            # Cosmological constant
            de_term = p.Omega_Lambda
        else:
            # CPL: w(a) = w₀ + w_a(1-a) = w₀ + w_a z/(1+z)
            f_de = one_plus_z**(3 * (1 + p.w_0 + p.w_a)) * \
                   math.exp(-3 * p.w_a * z / one_plus_z) if one_plus_z > 0 else 1
            de_term = p.Omega_Lambda * f_de
        
        H_squared = p.H_0**2 * (matter_term + radiation_term + curvature_term + de_term)
        return math.sqrt(max(H_squared, 0))
    
    def scale_factor(self, z: float) -> float:
        """Scale factor a = 1/(1+z)."""
        return 1.0 / (1 + z)
    
    def critical_density(self, z: float = 0) -> float:
        """
        Critical density at redshift z.
        
        ρ_c(z) = 3H(z)² / (8πG)
        
        Args:
            z: Redshift
            
        Returns:
            ρ_c: Critical density (kg/m³)
        """
        H = self.H_z(z)
        return 3 * H**2 / (8 * math.pi * G)
    
    def Omega_m_z(self, z: float) -> float:
        """
        Matter density parameter at redshift z.
        
        Ω_m(z) = Ω_m,0 (1+z)³ / E(z)²
        where E(z) = H(z)/H₀
        
        Args:
            z: Redshift
            
        Returns:
            Ω_m(z): Matter density parameter
        """
        E_z = self.H_z(z) / self.params.H_0
        return self.params.Omega_m * (1 + z)**3 / E_z**2
    
    def deceleration_parameter(self, z: float) -> float:
        """
        Deceleration parameter q(z) = -äa/ȧ².
        
        q = (1/2)(Ω_m + 2Ω_r - 2Ω_Λ) in ΛCDM
        
        q > 0: decelerating
        q < 0: accelerating
        
        Args:
            z: Redshift
            
        Returns:
            q: Deceleration parameter
        """
        Omega_m_now = self.Omega_m_z(z)
        # Approximate for w = -1
        Omega_Lambda_z = 1 - Omega_m_now  # Simplified
        return 0.5 * Omega_m_now - Omega_Lambda_z
    
    def age_at_redshift(self, z: float, n_steps: int = 1000) -> float:
        """
        Age of universe at redshift z (numerical integration).
        
        t(z) = ∫_z^∞ dz' / [(1+z') H(z')]
        
        Args:
            z: Redshift
            n_steps: Integration steps
            
        Returns:
            t: Age at redshift z (s)
        """
        # Integrate from z to z_max (effectively infinity)
        z_max = 1100  # Start from recombination
        if z > z_max:
            z_max = z + 100
        
        dz = (z_max - z) / n_steps
        integral = 0
        
        for i in range(n_steps):
            z_i = z + (i + 0.5) * dz
            H_i = self.H_z(z_i)
            if H_i > 0:
                integral += dz / ((1 + z_i) * H_i)
        
        return integral
    
    def comoving_distance(self, z: float, n_steps: int = 1000) -> float:
        """
        Comoving distance to redshift z.
        
        d_C = c ∫_0^z dz' / H(z')
        
        Args:
            z: Redshift
            n_steps: Integration steps
            
        Returns:
            d_C: Comoving distance (m)
        """
        if z <= 0:
            return 0
        
        dz = z / n_steps
        integral = 0
        
        for i in range(n_steps):
            z_i = (i + 0.5) * dz
            H_i = self.H_z(z_i)
            if H_i > 0:
                integral += dz / H_i
        
        return c * integral
    
    def luminosity_distance(self, z: float) -> float:
        """
        Luminosity distance D_L = (1+z) d_C.
        
        Args:
            z: Redshift
            
        Returns:
            D_L: Luminosity distance (m)
        """
        return (1 + z) * self.comoving_distance(z)
    
    def angular_diameter_distance(self, z: float) -> float:
        """
        Angular diameter distance D_A = d_C / (1+z).
        
        Args:
            z: Redshift
            
        Returns:
            D_A: Angular diameter distance (m)
        """
        return self.comoving_distance(z) / (1 + z)


# ==============================================================================
# INFLATION DYNAMICS (Equations 50-52)
# ==============================================================================

class InflationCalculator:
    """
    Inflation equations for early universe.
    
    Equation 50 (Slow-Roll Parameters):
    ε = (1/2)(V'/V)²  (in Planck units M_Pl=1)
    η = V''/V - (1/2)(V'/V)²
    
    Inflation requires ε, |η| << 1
    
    Equation 51 (Power Spectrum of Curvature Perturbations):
    P_R(k) = H² / (8π² ε c_s³) evaluated at k = aH
    
    Equation 52 (Number of e-Folds):
    N = ∫ H dt = ∫ dφ / √(2ε)
    ~50-60 for observable scales
    """
    
    @staticmethod
    def slow_roll_epsilon(V: float, V_prime: float) -> float:
        """
        First slow-roll parameter ε.
        
        ε = (M_Pl²/2)(V'/V)² with M_Pl = 1
        
        Args:
            V: Potential value (dimensionless in Planck units)
            V_prime: dV/dφ
            
        Returns:
            ε: First slow-roll parameter
        """
        if V <= 0:
            return math.inf
        return 0.5 * (V_prime / V)**2
    
    @staticmethod
    def slow_roll_eta(V: float, V_prime: float, V_double_prime: float) -> float:
        """
        Second slow-roll parameter η.
        
        η = M_Pl²(V''/V) - (1/2)(V'/V)²
        
        Args:
            V: Potential value
            V_prime: dV/dφ
            V_double_prime: d²V/dφ²
            
        Returns:
            η: Second slow-roll parameter
        """
        if V <= 0:
            return math.inf
        epsilon = InflationCalculator.slow_roll_epsilon(V, V_prime)
        return V_double_prime / V - epsilon
    
    @staticmethod
    def spectral_index(epsilon: float, eta: float) -> float:
        """
        Scalar spectral index n_s.
        
        n_s = 1 - 6ε + 2η
        
        Observed: n_s ≈ 0.9649 ± 0.0042
        
        Args:
            epsilon: First slow-roll parameter
            eta: Second slow-roll parameter
            
        Returns:
            n_s: Spectral index
        """
        return 1 - 6 * epsilon + 2 * eta
    
    @staticmethod
    def tensor_to_scalar_ratio(epsilon: float) -> float:
        """
        Tensor-to-scalar ratio r.
        
        r = 16 ε
        
        Current limits: r < 0.036 (Planck+BICEP)
        
        Args:
            epsilon: First slow-roll parameter
            
        Returns:
            r: Tensor-to-scalar ratio
        """
        return 16 * epsilon
    
    @staticmethod
    def curvature_power_spectrum(H: float, epsilon: float, c_s: float = 1.0) -> float:
        """
        Power spectrum of curvature perturbations.
        
        P_R(k) = H² / (8π² ε c_s³)
        
        Observed: Δ²_R ≈ 2.1×10⁻⁹
        
        Args:
            H: Hubble parameter during inflation (s⁻¹)
            epsilon: Slow-roll parameter
            c_s: Sound speed (default 1)
            
        Returns:
            P_R: Power spectrum amplitude
        """
        if epsilon <= 0:
            return math.inf
        return H**2 / (8 * math.pi**2 * epsilon * c_s**3)
    
    @staticmethod
    def e_folds_chaotic(phi_i: float, phi_end: float, m: float = 2) -> float:
        """
        Number of e-folds for chaotic inflation V = (1/2)m²φ².
        
        N = (1/4)(φ_i² - φ_end²)  (in Planck units)
        
        For φ_end ≈ √2 (ε = 1), N ≈ 50-60 requires φ_i ≈ 15-16 M_Pl
        
        Args:
            phi_i: Initial field value (Planck units)
            phi_end: End field value (typically √2)
            m: Potential exponent (default 2)
            
        Returns:
            N: Number of e-folds
        """
        if m == 2:
            return 0.25 * (phi_i**2 - phi_end**2)
        else:
            # General power-law
            return phi_i**2 / (2 * m)
    
    @staticmethod
    def reheating_temperature(V_end: float, g_star: float = 106.75) -> float:
        """
        Reheating temperature after inflation ends.
        
        T_reh = (30 V_end / (π² g_*))^(1/4)
        
        Args:
            V_end: Potential at end of inflation (J/m³ or M_Pl⁴)
            g_star: Relativistic degrees of freedom
            
        Returns:
            T_reh: Reheating temperature (K)
        """
        # Convert to SI if needed
        rho_reh = V_end  # Assume instantaneous conversion
        return (30 * rho_reh / (math.pi**2 * g_star))**0.25 / k_B


# ==============================================================================
# DARK ENERGY EQUATION OF STATE (Equations 91-93)
# ==============================================================================

class DarkEnergyCalculator:
    """
    Dark energy equation of state calculations.
    
    Equation 91 (Equation of State Parameter):
    w = p / (ρ c²) = -1 + (1/3) d ln ρ_DE / d ln a
    
    For cosmological constant: w = -1 (constant)
    Phantom: w < -1
    Quintessence: -1 < w < -1/3
    
    Equation 92 (CPL Parametrization):
    w(a) = w₀ + w_a(1-a) = w₀ + w_a z/(1+z)
    ρ_DE(a) = ρ_DE,0 × a^(-3(1+w₀+w_a)) × exp(-3w_a(1-a))
    
    Equation 93 (Growth Suppression):
    D(a) = 5Ω_m H₀² / 2 ∫_0^a da' / (a' H(a'))³
    """
    
    @staticmethod
    def equation_of_state_cpl(z: float, w_0: float, w_a: float) -> float:
        """
        CPL equation of state w(z).
        
        w(z) = w₀ + w_a × z/(1+z)
        
        Args:
            z: Redshift
            w_0: Present-day value (typically ~ -1)
            w_a: Evolution parameter (typically ~ 0 or small)
            
        Returns:
            w: Equation of state at z
        """
        return w_0 + w_a * z / (1 + z)
    
    @staticmethod
    def de_density_evolution(z: float, w_0: float, w_a: float) -> float:
        """
        Dark energy density evolution factor.
        
        ρ_DE(z) / ρ_DE,0 = (1+z)^(3(1+w₀+w_a)) × exp(-3w_a z/(1+z))
        
        Args:
            z: Redshift
            w_0: Present-day w
            w_a: Evolution parameter
            
        Returns:
            f(z): Density ratio ρ(z)/ρ(0)
        """
        one_plus_z = 1 + z
        if abs(w_a) < 1e-10:
            # Pure w = w_0 case
            return one_plus_z**(3 * (1 + w_0))
        else:
            # Full CPL
            return one_plus_z**(3 * (1 + w_0 + w_a)) * \
                   math.exp(-3 * w_a * z / one_plus_z)
    
    @staticmethod
    def growth_factor(z: float, Omega_m: float = 0.315) -> float:
        """
        Linear growth factor D(z) (normalized to D(0)=1).
        
        Approximate form (Carroll et al.):
        D(a) ≈ (5/2) Ω_m(a) / [Ω_m(a)^(4/7) - Ω_Λ(a) + (1+Ω_m(a)/2)(1+Ω_Λ(a)/70)]
        
        Args:
            z: Redshift
            Omega_m: Present-day matter density
            
        Returns:
            D(z)/D(0): Normalized growth factor
        """
        a = 1 / (1 + z)
        Omega_Lambda = 1 - Omega_m
        
        # Matter/Lambda at scale factor a
        E_a_sq = Omega_m / a**3 + Omega_Lambda
        Omega_m_a = (Omega_m / a**3) / E_a_sq
        Omega_L_a = Omega_Lambda / E_a_sq
        
        # Carroll approximation
        D_a = (5/2) * Omega_m_a / (
            Omega_m_a**(4/7) - Omega_L_a + 
            (1 + Omega_m_a/2) * (1 + Omega_L_a/70)
        )
        
        # Normalize to z=0
        D_0 = (5/2) * Omega_m / (
            Omega_m**(4/7) - Omega_Lambda + 
            (1 + Omega_m/2) * (1 + Omega_Lambda/70)
        )
        
        return D_a / D_0 if D_0 > 0 else 1
    
    @staticmethod
    def sigma_8_z(z: float, sigma_8_0: float = 0.811) -> float:
        """
        σ₈(z) = σ₈(0) × D(z).
        
        Args:
            z: Redshift
            sigma_8_0: Present-day σ₈
            
        Returns:
            σ₈(z): Density fluctuation at z
        """
        D_z = DarkEnergyCalculator.growth_factor(z)
        return sigma_8_0 * D_z


# ==============================================================================
# BLACK HOLE THERMODYNAMICS (Equations 94-96)
# ==============================================================================

class BHThermodynamicsCalculator:
    """
    Black hole thermodynamics (Hawking radiation).
    
    Equation 94 (Hawking Temperature):
    T_H = ℏc³ / (8π G M k_B)
    
    Equation 95 (Bekenstein-Hawking Entropy):
    S = k_B c³ A / (4 G ℏ) = 4π k_B G M² / (ℏ c)
    
    Equation 96 (Evaporation Lifetime):
    τ_evap = 5120 π G² M³ / (ℏ c⁴)
    """
    
    @staticmethod
    def hawking_temperature(M: float) -> float:
        """
        Hawking temperature of black hole.
        
        T_H = ℏc³ / (8π G M k_B)
        
        For M = M_sun: T_H ≈ 6×10⁻⁸ K
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            T_H: Hawking temperature (K)
        """
        if M <= 0:
            return math.inf
        return hbar * c**3 / (8 * math.pi * G * M * k_B)
    
    @staticmethod
    def schwarzschild_radius(M: float) -> float:
        """
        Schwarzschild radius r_s = 2GM/c².
        
        Args:
            M: Mass (kg)
            
        Returns:
            r_s: Schwarzschild radius (m)
        """
        return 2 * G * M / c**2
    
    @staticmethod
    def horizon_area(M: float) -> float:
        """
        Event horizon area A = 4π r_s² = 16π G² M² / c⁴.
        
        Args:
            M: Mass (kg)
            
        Returns:
            A: Horizon area (m²)
        """
        r_s = BHThermodynamicsCalculator.schwarzschild_radius(M)
        return 4 * math.pi * r_s**2
    
    @staticmethod
    def bekenstein_hawking_entropy(M: float) -> float:
        """
        Bekenstein-Hawking entropy.
        
        S = k_B c³ A / (4 G ℏ) = 4π k_B G M² / (ℏ c)
        
        Args:
            M: Mass (kg)
            
        Returns:
            S: Entropy (J/K)
        """
        return 4 * math.pi * k_B * G * M**2 / (hbar * c)
    
    @staticmethod
    def entropy_in_bits(M: float) -> float:
        """
        Entropy in bits (S / k_B ln 2).
        
        Args:
            M: Mass (kg)
            
        Returns:
            N_bits: Number of bits of information
        """
        S = BHThermodynamicsCalculator.bekenstein_hawking_entropy(M)
        return S / (k_B * math.log(2))
    
    @staticmethod
    def evaporation_lifetime(M: float) -> float:
        """
        Black hole evaporation lifetime.
        
        τ = 5120 π G² M³ / (ℏ c⁴)
        
        For M = M_sun: τ ≈ 2×10⁶⁷ yr
        
        Args:
            M: Initial mass (kg)
            
        Returns:
            τ: Evaporation time (s)
        """
        return 5120 * math.pi * G**2 * M**3 / (hbar * c**4)
    
    @staticmethod
    def hawking_power(M: float) -> float:
        """
        Hawking radiation power P = σ A T⁴ (Stefan-Boltzmann).
        
        P = ℏ c⁶ / (15360 π G² M²)
        
        Args:
            M: Mass (kg)
            
        Returns:
            P: Radiated power (W)
        """
        return hbar * c**6 / (15360 * math.pi * G**2 * M**2)
    
    @staticmethod
    def primordial_bh_mass_surviving(t: float) -> float:
        """
        Minimum primordial BH mass to survive until time t.
        
        M_min³ ≈ ℏ c⁴ t / (5120 π G²)
        
        For t = age of universe: M_min ≈ 10¹¹ kg
        
        Args:
            t: Age (s)
            
        Returns:
            M_min: Minimum surviving mass (kg)
        """
        return (hbar * c**4 * t / (5120 * math.pi * G**2))**(1/3)


# ==============================================================================
# LOOP QUANTUM COSMOLOGY (Equations 97-99)
# ==============================================================================

class LQCCalculator:
    """
    Loop Quantum Cosmology equations.
    
    Equation 97 (Effective Friedmann - Bounce Modification):
    H² = (8πG/3) ρ (1 - ρ/ρ_crit)
    
    Where ρ_crit = 0.41 ρ_Pl ≈ 5.16×10⁹⁶ kg/m³
    
    At ρ = ρ_crit, H = 0: BOUNCE replaces singularity
    
    Equation 98 (Perturbation Spectrum):
    P(k) ∝ k^(n_s-1) (1 + k/k_*)^(-α)
    
    k_* = bounce scale, α = UV suppression
    
    Equation 99 (Bounce Time Scale):
    t_b ~ √(3/(8πG ρ_crit)) ~ t_Pl ~ 10⁻⁴³ s
    """
    
    RHO_CRIT_LQC = 0.41 * rho_Pl  # Critical density for bounce
    
    @staticmethod
    def effective_H_squared(rho: float) -> float:
        """
        Effective Friedmann equation with quantum corrections.
        
        H² = (8πG/3) ρ (1 - ρ/ρ_crit)
        
        Args:
            rho: Energy density (kg/m³)
            
        Returns:
            H²: Squared Hubble parameter (s⁻²)
        """
        rho_crit = LQCCalculator.RHO_CRIT_LQC
        return (8 * math.pi * G / 3) * rho * (1 - rho / rho_crit)
    
    @staticmethod
    def is_bouncing(rho: float) -> bool:
        """
        Check if universe is at bounce (H = 0).
        
        Args:
            rho: Energy density (kg/m³)
            
        Returns:
            bool: True if at bounce
        """
        return abs(rho - LQCCalculator.RHO_CRIT_LQC) / LQCCalculator.RHO_CRIT_LQC < 0.01
    
    @staticmethod
    def bounce_timescale() -> float:
        """
        Characteristic bounce timescale.
        
        t_b ~ √(3/(8πG ρ_crit)) ~ t_Pl
        
        Returns:
            t_b: Bounce timescale (s)
        """
        return math.sqrt(3 / (8 * math.pi * G * LQCCalculator.RHO_CRIT_LQC))
    
    @staticmethod
    def perturbation_spectrum_lqc(k: float, k_star: float, n_s: float = 0.96, 
                                   alpha: float = 2.0) -> float:
        """
        LQC perturbation spectrum with UV suppression.
        
        P(k) ∝ k^(n_s-1) (1 + k/k_*)^(-α)
        
        Args:
            k: Wavenumber (Mpc⁻¹)
            k_star: Bounce scale
            n_s: Spectral index
            alpha: UV suppression exponent
            
        Returns:
            P(k): Power spectrum (relative)
        """
        return k**(n_s - 1) * (1 + k / k_star)**(-alpha)
    
    @staticmethod
    def holonomy_correction(mu: float, rho: float) -> float:
        """
        Holonomy correction factor for LQC.
        
        sin²(μc) / μ² → 1 - μ² ρ / ρ_Pl + ...
        
        Args:
            mu: Holonomy scale (Planck length units)
            rho: Energy density
            
        Returns:
            Correction factor
        """
        return 1 - mu**2 * rho / rho_Pl


# ==============================================================================
# GRAVITATIONAL WAVES PRIMORDIAL (Equations 53-54)
# ==============================================================================

class PrimordialGWCalculator:
    """
    Primordial gravitational wave calculations.
    
    Equation 53 (Tensor Power Spectrum):
    P_T(k) = 2H² / (π² M_Pl²) evaluated at k = aH
    
    Tensor-scalar ratio: r = P_T / P_R ≈ 16 ε
    
    Equation 54 (Stochastic GW Energy Density):
    Ω_GW(f) = (f / ρ_c) dρ_GW/df = (π² f⁴ / 3H₀²) ∫ P_T(k) dk
    """
    
    @staticmethod
    def tensor_power_spectrum(H: float) -> float:
        """
        Primordial tensor power spectrum.
        
        P_T = 2H² / (π² M_Pl²)
        
        Args:
            H: Hubble parameter during inflation (s⁻¹)
            
        Returns:
            P_T: Tensor power spectrum amplitude
        """
        return 2 * H**2 / (math.pi**2 * M_Pl**2 * c**4)  # SI units
    
    @staticmethod
    def tensor_spectral_index(epsilon: float) -> float:
        """
        Tensor spectral index n_T.
        
        n_T = -2ε (slightly red)
        
        Args:
            epsilon: Slow-roll parameter
            
        Returns:
            n_T: Tensor spectral index
        """
        return -2 * epsilon
    
    @staticmethod
    def omega_gw_primordial(f: float, r: float = 0.01, n_T: float = 0) -> float:
        """
        Stochastic GW energy density today.
        
        Ω_GW ≈ Ω_r × r × A_s × (f / f_eq)^(n_T) × T(f)
        
        Args:
            f: Frequency (Hz)
            r: Tensor-to-scalar ratio
            n_T: Tensor spectral index
            
        Returns:
            Ω_GW: Fractional energy density
        """
        f_eq = 1.6e-17  # Matter-radiation equality frequency (Hz)
        Omega_r = 9.24e-5
        
        # Transfer function (flat for f << f_eq)
        if f < f_eq:
            T_f = 1
        else:
            T_f = (f_eq / f)**2  # Radiation-era damping
        
        return Omega_r * r * A_s * (f / f_eq)**n_T * T_f


# ==============================================================================
# UNIFIED COSMOLOGY CALCULATOR
# ==============================================================================

class CosmologyCalculator:
    """Unified cosmology calculator."""
    
    def __init__(self, params: CosmologyParams = None):
        self.params = params or CosmologyParams()
        self.friedmann = FriedmannCalculator(self.params)
        self.inflation = InflationCalculator()
        self.dark_energy = DarkEnergyCalculator()
        self.bh_thermo = BHThermodynamicsCalculator()
        self.lqc = LQCCalculator()
        self.pgw = PrimordialGWCalculator()
    
    def compute_at_redshift(self, z: float) -> CosmologyResult:
        """
        Compute cosmological quantities at redshift z.
        
        Args:
            z: Redshift
            
        Returns:
            CosmologyResult with all quantities
        """
        H_z = self.friedmann.H_z(z)
        a = self.friedmann.scale_factor(z)
        rho_crit = self.friedmann.critical_density(z)
        Omega_m_z = self.friedmann.Omega_m_z(z)
        t_universe = self.friedmann.age_at_redshift(z)
        D_L = self.friedmann.luminosity_distance(z)
        D_A = self.friedmann.angular_diameter_distance(z)
        
        return CosmologyResult(
            H_z=H_z,
            a=a,
            rho_crit=rho_crit,
            Omega_m_z=Omega_m_z,
            t_universe=t_universe,
            D_L=D_L,
            D_A=D_A,
            additional={
                'q': self.friedmann.deceleration_parameter(z),
                'D_growth': self.dark_energy.growth_factor(z),
                'sigma_8_z': self.dark_energy.sigma_8_z(z),
                'w_z': self.dark_energy.equation_of_state_cpl(
                    z, self.params.w_0, self.params.w_a
                ),
            }
        )


# ==============================================================================
# REGISTRY FOR CONDENSEDPHYSICS INTEGRATION
# ==============================================================================

COSMOLOGY_CALCULATOR = CosmologyCalculator()

COSMOLOGY_CALCULATORS = {
    'Friedmann': FriedmannCalculator,
    'Inflation': InflationCalculator,
    'DarkEnergy': DarkEnergyCalculator,
    'BHThermodynamics': BHThermodynamicsCalculator,
    'LQC': LQCCalculator,
    'PrimordialGW': PrimordialGWCalculator,
    'Cosmology': CosmologyCalculator,
}

# ==============================================================================
# QUICK TEST
# ==============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("COSMOLOGY EQUATIONS MODULE TEST")
    print("=" * 70)
    
    calc = CosmologyCalculator()
    
    # Test at different redshifts
    redshifts = [0, 0.5, 1, 2, 5, 10]
    
    print("\nFriedmann Equations:")
    print("-" * 60)
    for z in redshifts:
        result = calc.compute_at_redshift(z)
        print(f"z = {z:5.1f}: H = {result.H_z*Mpc/1e3:.1f} km/s/Mpc, "
              f"Ω_m = {result.Omega_m_z:.3f}, "
              f"D_L = {result.D_L/Mpc:.1f} Mpc")
    
    print("\nDark Energy:")
    print("-" * 60)
    de = DarkEnergyCalculator()
    for z in [0, 0.5, 1, 2]:
        w = de.equation_of_state_cpl(z, -1.0, 0.0)
        D = de.growth_factor(z)
        print(f"z = {z}: w(z) = {w:.3f}, D(z) = {D:.4f}")
    
    print("\nBlack Hole Thermodynamics:")
    print("-" * 60)
    bh = BHThermodynamicsCalculator()
    masses = [1.0, 10.0, 1e6, 1e9]  # In solar masses
    for M_ratio in masses:
        M = M_ratio * M_sun
        T_H = bh.hawking_temperature(M)
        S = bh.bekenstein_hawking_entropy(M)
        tau = bh.evaporation_lifetime(M)
        print(f"M = {M_ratio:.0e} M_sun: T_H = {T_H:.2e} K, "
              f"τ = {tau/Gyr:.2e} Gyr")
    
    print("\nLQC Bounce:")
    print("-" * 60)
    lqc = LQCCalculator()
    print(f"Critical density: {lqc.RHO_CRIT_LQC:.2e} kg/m³")
    print(f"Bounce timescale: {lqc.bounce_timescale():.2e} s (≈ t_Pl = {t_Pl:.2e} s)")
