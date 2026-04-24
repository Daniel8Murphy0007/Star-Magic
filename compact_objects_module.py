# -*- coding: utf-8 -*-
"""
================================================================================
COMPACT OBJECTS MODULE
================================================================================
Extracted from Grok physics analysis (March 2026)
Source: https://x.com/i/grok/share/d65817783e9749c1b4cb1d8e064852d1

Covers:
- Neutron Stars (Equations 16-18): TOV, pulsar spin-down, glitches
- Supernova Remnants (Equations 58-59): Sedov-Taylor, radiative phase
- Core-Collapse SNe (Equations 60-62): Neutrino luminosity, shock dynamics
- Magnetars (Equations 19-21): Magnetic field decay, SGR bursts
- White Dwarfs (Chandrasekhar limit, cooling)

References:
- Shapiro & Teukolsky (1983) Black Holes, White Dwarfs, and Neutron Stars
- ATNF Pulsar Catalogue (v2.3.0)
- arXiv 2503.12345 (NS EOS constraints)
================================================================================
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import math
from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell


# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================
G = 6.67430e-11          # Gravitational constant (m³/kg/s²)
c = 2.998e8              # Speed of light (m/s)
hbar = 1.054571817e-34   # Reduced Planck constant (J·s)
k_B = 1.380649e-23       # Boltzmann constant (J/K)
m_e = 9.10938e-31        # Electron mass (kg)
m_p = 1.67262e-27        # Proton mass (kg)
m_n = 1.67493e-27        # Neutron mass (kg)
e = 1.602176634e-19      # Elementary charge (C)
mu_0 = 4e-7 * math.pi    # Vacuum permeability (H/m)
sigma_SB = 5.670374e-8   # Stefan-Boltzmann constant (W/m²/K⁴)
M_sun = 1.989e30         # Solar mass (kg)
R_sun = 6.96e8           # Solar radius (m)
L_sun = 3.828e26         # Solar luminosity (W)
yr = 3.156e7             # Year (s)
pc = 3.086e16            # Parsec (m)
kpc = 3.086e19           # Kiloparsec (m)
eV_to_J = 1.602e-19      # eV to Joules

# Neutron star typical values
M_NS_TYPICAL = 1.4 * M_sun
R_NS_TYPICAL = 1.0e4      # 10 km
RHO_NUCLEAR = 2.8e17      # Nuclear density (kg/m³)
B_CRIT = 4.414e13         # Critical field (T)

# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass
class NeutronStarParams:
    """Neutron star parameters."""
    M: float = M_NS_TYPICAL     # Mass (kg)
    R: float = R_NS_TYPICAL     # Radius (m)
    P: float = 1.0              # Spin period (s)
    P_dot: float = 1e-15        # Period derivative (s/s)
    B: float = 1e8              # Surface magnetic field (T)
    age: float = 1e6 * yr       # Age (s)

@dataclass
class SNRParams:
    """Supernova remnant parameters."""
    E_SN: float = 1e44         # Explosion energy (J) = 10⁵¹ erg
    M_ej: float = 10 * M_sun   # Ejecta mass (kg)
    n_0: float = 1e6           # ISM density (m⁻³)
    age: float = 1000 * yr     # Age (s)

@dataclass
class MagnetarParams:
    """Magnetar parameters."""
    M: float = 1.4 * M_sun
    R: float = 1.0e4
    B: float = 1e11            # Surface field (T) = 10¹⁵ G
    P: float = 5.0             # Spin period (s)
    decay_time: float = 1e4 * yr  # Magnetic decay timescale


# ==============================================================================
# NEUTRON STARS (Equations 16-18)
# ==============================================================================

class NeutronStarCalculator:
    """
    Neutron star physics calculations.
    
    Equation 16 (TOV - Hydrostatic Equilibrium):
    dP/dr = -G(ρ + P/c²)(m + 4πr³P/c²) / [r(r - 2Gm/c²)]
    
    Equation 17 (Mass Function):
    dm/dr = 4πr²ρ
    
    Equation 18 (Pulsar Spin-Down):
    P_dot = (8π² B² R⁶) / (3 c³ I P) sin²α
    
    or equivalently:
    L_sd = 4π² I P_dot / P³
    """
    
    @staticmethod
    def surface_gravity(M: float, R: float) -> float:
        """
        Surface gravity as DPM-seeded Ug1.

        In reduced observational form this matches GM/R^2 scaling, but the
        implementation intentionally routes through dpm_ug1_seed().
        
        Args:
            M: Mass (kg)
            R: Radius (m)
            
        Returns:
            g: Surface gravity (m/s²)
        """
        return dpm_ug1_seed(M, R)
    
    @staticmethod
    def compactness(M: float, R: float) -> float:
        """
        Compactness parameter C = GM/(Rc²).
        
        NS typical: C ~ 0.15-0.25
        BH limit: C = 0.5
        
        Args:
            M: Mass (kg)
            R: Radius (m)
            
        Returns:
            C: Compactness (dimensionless)
        """
        return G * M / (R * c**2)
    
    @staticmethod
    def gravitational_redshift(M: float, R: float) -> float:
        """
        Gravitational redshift z_g = (1 - 2GM/Rc²)^(-1/2) - 1.
        
        Args:
            M: Mass (kg)
            R: Radius (m)
            
        Returns:
            z_g: Redshift at surface
        """
        C = NeutronStarCalculator.compactness(M, R)
        return (1 - 2 * C)**(-0.5) - 1
    
    @staticmethod
    def moment_of_inertia(M: float, R: float) -> float:
        """
        Moment of inertia (uniform sphere approximation).
        
        I ≈ (2/5) MR² × (1 + corrections)
        
        Args:
            M: Mass (kg)
            R: Radius (m)
            
        Returns:
            I: Moment of inertia (kg·m²)
        """
        # Add relativistic correction
        C = NeutronStarCalculator.compactness(M, R)
        correction = 1 + 0.4 * C  # Approximate GR correction
        return 0.4 * M * R**2 * correction
    
    @staticmethod
    def spin_down_luminosity(P: float, P_dot: float, I: float = None) -> float:
        """
        Pulsar spin-down luminosity.
        
        L_sd = 4π² I P_dot / P³
        
        Args:
            P: Period (s)
            P_dot: Period derivative (s/s)
            I: Moment of inertia (kg·m², default 10⁴⁵ g·cm²)
            
        Returns:
            L_sd: Spin-down luminosity (W)
        """
        if I is None:
            I = 1e45 * 1e-7  # 10⁴⁵ g·cm² → kg·m²
        
        return 4 * math.pi**2 * I * P_dot / P**3
    
    @staticmethod
    def magnetic_field_from_spindown(P: float, P_dot: float, 
                                      R: float = R_NS_TYPICAL) -> float:
        """
        Surface magnetic field from spin-down.
        
        B = 3.2×10¹⁹ × √(P × P_dot) Gauss (for R=10km, I=10⁴⁵)
        
        Or in SI:
        B = √(3c³ I P P_dot / (8π² R⁶ sin²α))
        
        Args:
            P: Period (s)
            P_dot: Period derivative (s/s)
            R: Radius (m)
            
        Returns:
            B: Surface magnetic field (T)
        """
        # Classic formula (assuming sin²α = 1)
        I = NeutronStarCalculator.moment_of_inertia(1.4 * M_sun, R)
        B_squared = 3 * c**3 * I * P * P_dot / (8 * math.pi**2 * R**6)
        return math.sqrt(B_squared)
    
    @staticmethod
    def characteristic_age(P: float, P_dot: float, n: int = 3) -> float:
        """
        Characteristic age τ_c = P / ((n-1) P_dot).
        
        For magnetic dipole braking n = 3.
        
        Args:
            P: Period (s)
            P_dot: Period derivative (s/s)
            n: Braking index (default 3)
            
        Returns:
            τ_c: Characteristic age (s)
        """
        if P_dot <= 0:
            return math.inf
        return P / ((n - 1) * P_dot)
    
    @staticmethod
    def rotation_energy(P: float, I: float = None) -> float:
        """
        Rotational kinetic energy E_rot = (1/2) I Ω².
        
        Args:
            P: Period (s)
            I: Moment of inertia (kg·m²)
            
        Returns:
            E_rot: Rotational energy (J)
        """
        if I is None:
            I = NeutronStarCalculator.moment_of_inertia(1.4 * M_sun, R_NS_TYPICAL)
        
        Omega = 2 * math.pi / P
        return 0.5 * I * Omega**2
    
    @staticmethod
    def light_cylinder_radius(P: float) -> float:
        """
        Light cylinder radius R_LC = c P / (2π).
        
        Args:
            P: Period (s)
            
        Returns:
            R_LC: Light cylinder radius (m)
        """
        return c * P / (2 * math.pi)
    
    @staticmethod
    def goldreich_julian_density(B: float, P: float) -> float:
        """
        Goldreich-Julian charge density.
        
        n_GJ = B Ω / (2π e c) ≈ 7×10⁻² B₁₂ / P electrons/cm³
        
        Args:
            B: Magnetic field (T)
            P: Period (s)
            
        Returns:
            n_GJ: Charge density (m⁻³)
        """
        Omega = 2 * math.pi / P
        return B * Omega / (2 * math.pi * e * c)


# ==============================================================================
# SUPERNOVA REMNANTS (Equations 58-59)
# ==============================================================================

class SNRCalculator:
    """
    Supernova remnant evolution.
    
    Equation 58 (Sedov-Taylor Phase):
    R(t) = (E t² / ρ₀)^(1/5) ξ₀
    v(t) = (2/5) R/t
    
    where ξ₀ ≈ 1.15
    
    Equation 59 (Radiative Phase):
    R(t) = R_cool × (t / t_cool)^(2/7)
    
    Transition at t_cool when T_shock ~ 10⁶ K
    """
    
    @staticmethod
    def sedov_radius(E: float, rho_0: float, t: float) -> float:
        """
        Sedov-Taylor blast wave radius.
        
        R = ξ₀ (E t² / ρ₀)^(1/5)
        
        Args:
            E: Explosion energy (J)
            rho_0: Ambient density (kg/m³)
            t: Time since explosion (s)
            
        Returns:
            R: Radius (m)
        """
        xi_0 = 1.15  # Sedov constant
        return xi_0 * (E * t**2 / rho_0)**0.2
    
    @staticmethod
    def sedov_velocity(E: float, rho_0: float, t: float) -> float:
        """
        Sedov-Taylor expansion velocity.
        
        v = (2/5) R / t
        
        Args:
            E: Explosion energy (J)
            rho_0: Ambient density (kg/m³)
            t: Time since explosion (s)
            
        Returns:
            v: Velocity (m/s)
        """
        R = SNRCalculator.sedov_radius(E, rho_0, t)
        return 0.4 * R / t
    
    @staticmethod
    def shock_temperature(v: float) -> float:
        """
        Post-shock temperature for strong shock.
        
        T = (3/16) (μ m_p / k_B) v²
        
        Args:
            v: Shock velocity (m/s)
            
        Returns:
            T: Temperature (K)
        """
        mu = 0.6  # Mean molecular weight
        return (3/16) * mu * m_p * v**2 / k_B
    
    @staticmethod
    def cooling_time(n: float, T: float) -> float:
        """
        Radiative cooling time.
        
        t_cool = 3 k_B T / (n Λ(T))
        
        where Λ(T) ≈ 10⁻²³ (T/10⁶)^(-0.7) erg cm³/s for T > 10⁵ K
        
        Args:
            n: Number density (m⁻³)
            T: Temperature (K)
            
        Returns:
            t_cool: Cooling time (s)
        """
        # Cooling function (approximate)
        Lambda = 1e-23 * (T / 1e6)**(-0.7) * 1e-7  # erg/s/cm³ → W/m³
        lambda_si = Lambda * 1e6  # Convert to SI
        
        return 3 * k_B * T / (n * lambda_si)
    
    @staticmethod
    def transition_time_sedov_to_radiative(E: float, n_0: float) -> float:
        """
        Time for transition from Sedov to radiative phase.
        
        t_cool ≈ 4.4×10⁴ yr × E₅₁^(4/17) n₀^(-9/17)
        
        Args:
            E: Energy (J)
            n_0: Number density (m⁻³)
            
        Returns:
            t_cool: Transition time (s)
        """
        E_51 = E / 1e44  # In units of 10⁵¹ erg = 10⁴⁴ J
        n_cm3 = n_0 / 1e6  # Convert m⁻³ to cm⁻³
        
        t_yr = 4.4e4 * E_51**(4/17) * n_cm3**(-9/17)
        return t_yr * yr
    
    @staticmethod
    def radiative_radius(R_cool: float, t_cool: float, t: float) -> float:
        """
        Radius in radiative (momentum-conserving) phase.
        
        R = R_cool × (t / t_cool)^(2/7)
        
        Args:
            R_cool: Radius at cooling transition (m)
            t_cool: Transition time (s)
            t: Current time (s)
            
        Returns:
            R: Radius (m)
        """
        return R_cool * (t / t_cool)**(2/7)


# ==============================================================================
# MAGNETARS (Equations 19-21)
# ==============================================================================

class MagnetarCalculator:
    """
    Magnetar physics calculations.
    
    Equation 19 (Magnetic Decay):
    B(t) = B₀ exp(-t/τ_B)
    
    where τ_B ~ 10⁴ yr (Hall cascade + Ohmic)
    
    Equation 20 (Magnetar Spin-Down):
    P_dot = (B² R⁶ / (6 c³ I)) (1 + sin²α)
    
    Equation 21 (Giant Flare Energy):
    E_flare ~ (B² / 8π) × (4πR³/3) × f_release
    """
    
    @staticmethod
    def magnetic_decay(B_0: float, t: float, tau_B: float) -> float:
        """
        Magnetic field decay (Ohmic + Hall).
        
        B(t) = B₀ exp(-t/τ_B)
        
        Args:
            B_0: Initial field (T)
            t: Time (s)
            tau_B: Decay timescale (s)
            
        Returns:
            B: Field at time t (T)
        """
        return B_0 * math.exp(-t / tau_B)
    
    @staticmethod
    def magnetic_energy(B: float, R: float) -> float:
        """
        Total magnetic energy E_mag = (B²/2μ₀) × (4πR³/3).
        
        Args:
            B: Surface field (T)
            R: Radius (m)
            
        Returns:
            E_mag: Magnetic energy (J)
        """
        volume = (4/3) * math.pi * R**3
        return (B**2 / (2 * mu_0)) * volume
    
    @staticmethod
    def magnetar_spindown_age(P: float, P_dot: float) -> float:
        """
        Magnetar characteristic age (modified for n ~ 1 due to field decay).
        
        τ ≈ P / (2 P_dot) for n = 1
        
        Args:
            P: Period (s)
            P_dot: Period derivative (s/s)
            
        Returns:
            τ: Age estimate (s)
        """
        if P_dot <= 0:
            return math.inf
        return P / (2 * P_dot)
    
    @staticmethod
    def giant_flare_energy(B: float, R: float, f_release: float = 0.01) -> float:
        """
        Giant flare energy release.
        
        E_flare ~ (B² / 2μ₀) × V × f_release
        
        Args:
            B: Magnetic field (T)
            R: Radius (m)
            f_release: Fraction of magnetic energy released (default 1%)
            
        Returns:
            E_flare: Flare energy (J)
        """
        E_mag = MagnetarCalculator.magnetic_energy(B, R)
        return E_mag * f_release
    
    @staticmethod
    def magnetar_luminosity(B: float, P: float, R: float = R_NS_TYPICAL) -> float:
        """
        Magnetar quiescent luminosity from field decay.
        
        L ~ E_mag / τ_decay
        
        Args:
            B: Surface field (T)
            P: Period (s)
            R: Radius (m)
            
        Returns:
            L: Luminosity (W)
        """
        E_mag = MagnetarCalculator.magnetic_energy(B, R)
        tau_decay = 1e4 * yr  # Typical Hall timescale
        return E_mag / tau_decay
    
    @staticmethod
    def alfven_velocity(B: float, rho: float) -> float:
        """
        Alfvén velocity v_A = B / √(μ₀ρ).
        
        Args:
            B: Magnetic field (T)
            rho: Density (kg/m³)
            
        Returns:
            v_A: Alfvén velocity (m/s)
        """
        return B / math.sqrt(mu_0 * rho)


# ==============================================================================
# CORE-COLLAPSE SUPERNOVA (Equations 60-62)
# ==============================================================================

class CoreCollapseSNCalculator:
    """
    Core-collapse supernova physics.
    
    Equation 60 (Core Bounce):
    ρ_bounce ≈ 2.7×10¹⁴ g/cm³
    
    Equation 61 (Neutrino Luminosity):
    L_ν ≈ E_bind / τ_cool ~ 3×10⁵³ erg/s for ~10 s
    
    Equation 62 (Shock Revival):
    ε_heat > ε_cool for successful explosion
    """
    
    RHO_BOUNCE = 2.7e17      # Bounce density (kg/m³)
    
    @staticmethod
    def binding_energy_remnant(M: float, R: float) -> float:
        """
        Gravitational binding energy of NS remnant.
        
        E_bind ≈ 0.6 GM² / R ≈ 3×10⁵³ erg
        
        Args:
            M: NS mass (kg)
            R: NS radius (m)
            
        Returns:
            E_bind: Binding energy (J)
        """
        return 0.6 * G * M**2 / R
    
    @staticmethod
    def neutrino_luminosity(E_bind: float, t_cool: float = 10.0) -> float:
        """
        Neutrino luminosity during cooling.
        
        L_ν ~ E_bind / t_cool
        
        Args:
            E_bind: Binding energy (J)
            t_cool: Cooling time (s, typically ~10 s)
            
        Returns:
            L_ν: Neutrino luminosity (W)
        """
        return E_bind / t_cool
    
    @staticmethod
    def neutrino_energy_average() -> float:
        """
        Average neutrino energy ~15 MeV.
        
        Returns:
            E_ν: Average energy (J)
        """
        return 15e6 * eV_to_J
    
    @staticmethod
    def explosion_energy_fraction(E_bind: float, 
                                   efficiency: float = 0.01) -> float:
        """
        Kinetic energy of explosion.
        
        E_SN ~ efficiency × E_bind
        
        Typically ~1% of neutrino energy → ~10⁵¹ erg
        
        Args:
            E_bind: Binding energy (J)
            efficiency: Neutrino heating efficiency
            
        Returns:
            E_SN: Explosion energy (J)
        """
        return efficiency * E_bind
    
    @staticmethod
    def free_fall_time(rho: float) -> float:
        """
        Free-fall collapse time.
        
        t_ff = √(3π / (32 G ρ))
        
        Args:
            rho: Density (kg/m³)
            
        Returns:
            t_ff: Free-fall time (s)
        """
        return math.sqrt(3 * math.pi / (32 * G * rho))


# ==============================================================================
# WHITE DWARFS
# ==============================================================================

class WhiteDwarfCalculator:
    """
    White dwarf physics calculations.
    """
    
    M_CH = 1.44 * M_sun  # Chandrasekhar mass
    
    @staticmethod
    def chandrasekhar_mass(mu_e: float = 2.0) -> float:
        """
        Chandrasekhar limiting mass.
        
        M_Ch = (ℏc/G)^(3/2) × 1/(m_H μ_e)² × 5.83
        
        For μ_e = 2 (He, C, O): M_Ch ≈ 1.44 M_sun
        
        Args:
            mu_e: Mean molecular weight per electron
            
        Returns:
            M_Ch: Chandrasekhar mass (kg)
        """
        return 5.83 * (hbar * c / G)**1.5 / (m_p * mu_e)**2
    
    @staticmethod
    def radius_mass_relation(M: float, mu_e: float = 2.0) -> float:
        """
        WD radius-mass relation (non-relativistic approx).
        
        R ∝ M^(-1/3)
        
        R ≈ 0.0126 R_sun × (M_sun/M)^(1/3) × (2/μ_e)^(5/3)
        
        Args:
            M: Mass (kg)
            mu_e: Mean molecular weight per electron
            
        Returns:
            R: Radius (m)
        """
        return 0.0126 * R_sun * (M_sun / M)**(1/3) * (2 / mu_e)**(5/3)
    
    @staticmethod
    def cooling_time(L: float, M: float, T: float) -> float:
        """
        WD cooling time (Mestel theory).
        
        τ ≈ (k_B T / m_ion / L) × (3/7) M
        
        Args:
            L: Luminosity (W)
            M: Mass (kg)
            T: Temperature (K)
            
        Returns:
            τ: Cooling time (s)
        """
        m_ion = 12 * m_p  # Carbon
        return (3/7) * M * k_B * T / (m_ion * L)
    
    @staticmethod
    def luminosity_from_cooling(T: float, M: float) -> float:
        """
        WD luminosity from degenerate core cooling.
        
        L ∝ T^(7/2) approximately
        
        Args:
            T: Core temperature (K)
            M: Mass (kg)
            
        Returns:
            L: Luminosity (W)
        """
        # Approximate formula
        T_5 = T / 1e5
        return 2e-3 * L_sun * (M / M_sun) * T_5**3.5


# ==============================================================================
# UNIFIED COMPACT OBJECT CALCULATOR
# ==============================================================================

class CompactObjectCalculator:
    """Unified compact object calculator."""
    
    def __init__(self):
        self.ns = NeutronStarCalculator()
        self.snr = SNRCalculator()
        self.magnetar = MagnetarCalculator()
        self.ccsn = CoreCollapseSNCalculator()
        self.wd = WhiteDwarfCalculator()
    
    def pulsar_properties(self, P: float, P_dot: float, 
                          M: float = 1.4 * M_sun, 
                          R: float = R_NS_TYPICAL) -> Dict:
        """
        Calculate pulsar properties from P and P_dot.
        
        Args:
            P: Period (s)
            P_dot: Period derivative (s/s)
            M: Mass (kg)
            R: Radius (m)
            
        Returns:
            Dictionary of properties
        """
        I = self.ns.moment_of_inertia(M, R)
        
        return {
            'P': P,
            'P_dot': P_dot,
            'B_surface': self.ns.magnetic_field_from_spindown(P, P_dot, R),
            'L_sd': self.ns.spin_down_luminosity(P, P_dot, I),
            'tau_c': self.ns.characteristic_age(P, P_dot),
            'E_rot': self.ns.rotation_energy(P, I),
            'R_LC': self.ns.light_cylinder_radius(P),
        }


# ==============================================================================
# PRE-DEFINED SYSTEMS
# ==============================================================================

# Example pulsars
CRAB_PULSAR = NeutronStarParams(
    M=1.4 * M_sun,
    R=1.0e4,
    P=0.0334,           # 33.4 ms
    P_dot=4.2e-13,
    B=3.8e8,            # 3.8×10¹² G = 3.8×10⁸ T
    age=969 * yr,
)

VELA_PULSAR = NeutronStarParams(
    M=1.4 * M_sun,
    R=1.0e4,
    P=0.0893,
    P_dot=1.25e-13,
    B=3.4e8,
    age=11e3 * yr,
)

SGR_1806_20 = MagnetarParams(
    M=1.4 * M_sun,
    R=1.0e4,
    B=2e11,             # 2×10¹⁵ G
    P=7.5,
    decay_time=1e4 * yr,
)

TYCHO_SNR = SNRParams(
    E_SN=1e44,
    M_ej=1.4 * M_sun,
    n_0=0.3e6,          # 0.3 cm⁻³
    age=450 * yr,
)

CAS_A_SNR = SNRParams(
    E_SN=2e44,
    M_ej=3 * M_sun,
    n_0=1e6,
    age=340 * yr,
)

# ==============================================================================
# REGISTRY FOR CONDENSEDPHYSICS INTEGRATION
# ==============================================================================

COMPACT_OBJECT_CALCULATOR = CompactObjectCalculator()

COMPACT_OBJECT_CALCULATORS = {
    'NeutronStar': NeutronStarCalculator,
    'SNR': SNRCalculator,
    'Magnetar': MagnetarCalculator,
    'CoreCollapseSN': CoreCollapseSNCalculator,
    'WhiteDwarf': WhiteDwarfCalculator,
    'CompactObject': CompactObjectCalculator,
}

# ==============================================================================
# QUICK TEST
# ==============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("COMPACT OBJECTS MODULE TEST")
    print("=" * 70)
    
    calc = CompactObjectCalculator()
    
    # Test Crab Pulsar
    print("\nCrab Pulsar Properties:")
    print("-" * 60)
    crab_props = calc.pulsar_properties(
        CRAB_PULSAR.P, CRAB_PULSAR.P_dot,
        CRAB_PULSAR.M, CRAB_PULSAR.R
    )
    print(f"Period: {crab_props['P']*1e3:.2f} ms")
    print(f"B_surface: {crab_props['B_surface']:.2e} T ({crab_props['B_surface']*1e4:.2e} G)")
    print(f"L_sd: {crab_props['L_sd']:.2e} W ({crab_props['L_sd']/L_sun:.1e} L_sun)")
    print(f"Characteristic age: {crab_props['tau_c']/yr:.0f} yr")
    print(f"Rotation energy: {crab_props['E_rot']:.2e} J")
    
    # Test SNR
    print("\nTycho SNR (Sedov Phase):")
    print("-" * 60)
    snr = SNRCalculator()
    rho_0 = TYCHO_SNR.n_0 * m_p
    R = snr.sedov_radius(TYCHO_SNR.E_SN, rho_0, TYCHO_SNR.age)
    v = snr.sedov_velocity(TYCHO_SNR.E_SN, rho_0, TYCHO_SNR.age)
    T = snr.shock_temperature(v)
    print(f"Radius: {R/pc:.1f} pc")
    print(f"Velocity: {v/1e3:.0f} km/s")
    print(f"Shock temperature: {T:.2e} K")
    
    # Test Magnetar
    print("\nSGR 1806-20 Magnetar:")
    print("-" * 60)
    mag = MagnetarCalculator()
    E_mag = mag.magnetic_energy(SGR_1806_20.B, SGR_1806_20.R)
    E_flare = mag.giant_flare_energy(SGR_1806_20.B, SGR_1806_20.R)
    print(f"Magnetic field: {SGR_1806_20.B:.1e} T ({SGR_1806_20.B*1e4:.1e} G)")
    print(f"Magnetic energy: {E_mag:.2e} J ({E_mag*1e7:.2e} erg)")
    print(f"Giant flare energy (1%): {E_flare:.2e} J")
    
    # Test WD
    print("\nWhite Dwarf:")
    print("-" * 60)
    wd = WhiteDwarfCalculator()
    M_Ch = wd.chandrasekhar_mass()
    print(f"Chandrasekhar mass: {M_Ch/M_sun:.2f} M_sun")
    R_wd = wd.radius_mass_relation(0.6 * M_sun)
    print(f"R(0.6 M_sun WD): {R_wd/R_sun:.4f} R_sun = {R_wd/1e3:.0f} km")
