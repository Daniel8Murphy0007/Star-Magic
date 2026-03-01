# -*- coding: utf-8 -*-
"""
================================================================================
GALAXY MERGERS AND STAR FORMATION RATE MODULE
================================================================================
Extracted from Grok physics analysis (March 2026)
Source: https://x.com/i/grok/share/d65817783e9749c1b4cb1d8e064852d1

Covers:
- Extended Press-Schechter (EPS) merger rates (Equations 5, 8)
- Star Formation Rate functions (quiescent + starburst, Equations 6-7)
- Black hole mass functions and accretion (Equations 8-9)
- Gravitational wave chirp mass (Equations 12-13)
- UQFF-tailored merger buoyancy and magnetism

References:
- Bond et al. 1991 (EPS formalism)
- Lacey & Cole 1993 (merger trees)
- arXiv 2403.11428 (merger rates z=4-9)
- arXiv 2209.03983 (galaxy-scale simulations)
- GWTC catalogs (LIGO/Virgo/KAGRA)
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
k_B = 1.380649e-23       # Boltzmann constant (J/K)
m_p = 1.6726219e-27      # Proton mass (kg)
sigma_T = 6.6524e-29     # Thomson cross-section (m²)
M_sun = 1.989e30         # Solar mass (kg)
Gyr = 3.156e16           # Gigayear (s)
Mpc = 3.086e22           # Megaparsec (m)
eV_to_J = 1.602e-19      # eV to Joules

# Cosmological parameters (Planck 2018)
H_0 = 67.4e3 / Mpc       # Hubble constant (s⁻¹)
Omega_m = 0.315          # Matter density parameter
Omega_Lambda = 0.685     # Dark energy density parameter
delta_c = 1.686          # Critical overdensity for collapse
sigma_8 = 0.811          # Power spectrum normalization at 8 Mpc/h

# UQFF Constants
F_REL = 4.30e33          # Relativistic coherence baseline (N)
E_LEP = 200e9 * eV_to_J  # LEP 1998 baseline energy (200 GeV in J)
Q_WAVE_BASE = 1e12       # Base wave resonance factor

# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass
class GalaxySystem:
    """Galaxy/merger system parameters."""
    name: str
    M_halo: float              # Halo mass (kg)
    M_stellar: float           # Stellar mass (kg)
    M_gas: float               # Gas mass (kg)
    z: float                   # Redshift
    SFR: float                 # Star formation rate (M_sun/yr → kg/s)
    M_BH: float                # Central black hole mass (kg)
    r_virial: float            # Virial radius (m)
    sigma_v: float             # Velocity dispersion (m/s)
    description: str = ""

@dataclass
class MergerResult:
    """Results from merger/SFR calculations."""
    merger_rate: float         # Mergers per Gyr (per volume)
    SFR_quiescent: float       # Quiescent SFR (kg/s)
    SFR_burst: float           # Merger-induced starburst SFR (kg/s)
    t_dyn: float               # Dynamical time (s)
    N_BH_cumulative: float     # Cumulative BH number above threshold
    M_dot_BH: float            # BH accretion rate (kg/s)
    t_Sal: float               # Salpeter time (s)
    F_UBii: float              # UQFF buoyancy force
    additional: Dict = field(default_factory=dict)


# ==============================================================================
# EXTENDED PRESS-SCHECHTER (EPS) CALCULATORS (Equations 5, 8)
# ==============================================================================

class EPSMergerCalculator:
    """
    Extended Press-Schechter formalism for halo merger rates.
    
    Equation 5: dN/dt dM = (2/π)^(1/2) (ρ̄/M) (σ_m/σ_M) |dδ_c/dz| × 
                          exp(-δ_c²/(2(σ_m² - σ_M²))) |dσ_M/dM|
    
    Where:
    - σ_M² = ∫ P(k) W²(kR) d³k/(2π)³ (variance on mass scale M)
    - δ_c ≈ 1.686 (critical overdensity)
    - P(k) ∝ k^(n_s-4) (CDM power spectrum, n_s ≈ 0.96)
    
    Merger rate ∝ (1+z)^m where m ~ 0.7-2.5 depending on mass.
    """
    
    @staticmethod
    def sigma_M(M: float, z: float = 0) -> float:
        """
        RMS density fluctuation at mass scale M.
        
        Approximate scaling: σ(M) ≈ σ_8 × (M / M_8)^(-α) × D(z)
        where M_8 ≈ 6×10^14 M_sun (8 Mpc/h sphere), α ≈ 0.2 for CDM
        
        Args:
            M: Mass (kg)
            z: Redshift
            
        Returns:
            σ_M: Dimensionless variance
        """
        M_8 = 6e14 * M_sun  # 8 Mpc/h sphere
        alpha = 0.2
        D_z = 1 / (1 + z)  # Approximate growth factor
        return sigma_8 * (M / M_8)**(-alpha) * D_z
    
    @staticmethod
    def critical_overdensity(z: float) -> float:
        """
        Critical overdensity for collapse (matter era approximation).
        
        δ_c(z) = 1.686 × (1 + z) in matter era
        
        Args:
            z: Redshift
            
        Returns:
            δ_c: Critical overdensity
        """
        # More accurate approximation including Lambda
        Omega_m_z = Omega_m * (1 + z)**3 / (Omega_m * (1 + z)**3 + Omega_Lambda)
        return delta_c * Omega_m_z**0.0055  # Weak z-dependence
    
    @staticmethod
    def erfc(x: float) -> float:
        """Complementary error function approximation."""
        # Abramowitz & Stegun approximation
        a1, a2, a3, a4 = 0.254829592, -0.284496736, 1.421413741, -1.453152027
        a5, p = 1.061405429, 0.3275911
        t = 1.0 / (1.0 + p * abs(x))
        erf_approx = 1 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1) * t * math.exp(-x*x)
        if x < 0:
            erf_approx = -erf_approx
        return 1 - erf_approx
    
    @staticmethod
    def merger_rate_density(M: float, z: float, m_exponent: float = 2.0) -> float:
        """
        Halo merger rate scaling with redshift.
        
        Rate ∝ (1+z)^m where m ~ 0.7-2.5
        
        Args:
            M: Halo mass (kg)
            z: Redshift
            m_exponent: Redshift exponent (mass-dependent, default 2.0)
            
        Returns:
            Merger rate scaling (dimensionless, normalized to z=0)
        """
        return (1 + z)**m_exponent
    
    @staticmethod
    def cumulative_BH_number(M_threshold: float, z: float, rho_bar: float) -> float:
        """
        Equation 8: Cumulative BH number above mass threshold.
        
        N(>M,z) = ρ̄ ∫_M^∞ (dM'/M'²) erfc(δ_c(z) / (√2 σ(M',z)))
        
        Args:
            M_threshold: Minimum BH mass (kg)
            z: Redshift
            rho_bar: Mean cosmic density (kg/m³)
            
        Returns:
            N: Cumulative number per unit volume (m⁻³)
        """
        delta_c_z = EPSMergerCalculator.critical_overdensity(z)
        sigma = EPSMergerCalculator.sigma_M(M_threshold, z)
        
        if sigma <= 0:
            return 0
        
        arg = delta_c_z / (math.sqrt(2) * sigma)
        erfc_val = EPSMergerCalculator.erfc(arg)
        
        # Simplified integral (leading term)
        N = rho_bar * erfc_val / M_threshold
        return N


# ==============================================================================
# STAR FORMATION RATE CALCULATORS (Equations 6-7)
# ==============================================================================

class StarFormationCalculator:
    """
    Star formation rate calculations.
    
    Equation 6 (Quiescent SFR):
    Ṁ_* = ε M_gas / t_dyn
    
    Where:
    - ε ~ 0.01-0.1 (efficiency)
    - t_dyn = √(3π / (32 G ρ)) (dynamical time)
    - For Toomre Q < 1, disk instability triggers bursts
    
    Equation 7 (Merger-Induced Starburst):
    Ṁ_burst = Ṁ_gas,inflow × ε_burst
    Ṁ_gas,inflow ≈ M_gas / t_orb
    t_orb = 2π √(r³ / GM)
    
    Enhancement factor: 10-100× at low z, less at high z (gas-rich)
    """
    
    @staticmethod
    def dynamical_time(rho: float) -> float:
        """
        Dynamical (free-fall) time.
        
        t_dyn = √(3π / (32 G ρ))
        
        Args:
            rho: Gas density (kg/m³)
            
        Returns:
            t_dyn: Dynamical time (s)
        """
        if rho <= 0:
            return math.inf
        return math.sqrt(3 * math.pi / (32 * G * rho))
    
    @staticmethod
    def orbital_time(r: float, M: float) -> float:
        """
        Orbital time at radius r.
        
        t_orb = 2π √(r³ / GM)
        
        Args:
            r: Orbital radius (m)
            M: Enclosed mass (kg)
            
        Returns:
            t_orb: Orbital time (s)
        """
        if M <= 0:
            return math.inf
        return 2 * math.pi * math.sqrt(r**3 / (G * M))
    
    @staticmethod
    def quiescent_SFR(M_gas: float, t_dyn: float, epsilon: float = 0.02) -> float:
        """
        Quiescent mode star formation rate.
        
        Ṁ_* = ε M_gas / t_dyn
        
        Args:
            M_gas: Gas mass (kg)
            t_dyn: Dynamical time (s)
            epsilon: Star formation efficiency (0.01-0.1)
            
        Returns:
            SFR: Star formation rate (kg/s)
        """
        if t_dyn <= 0:
            return 0
        return epsilon * M_gas / t_dyn
    
    @staticmethod
    def starburst_SFR(
        M_gas: float, 
        t_orb: float, 
        epsilon_burst: float = 0.5,
        enhancement: float = 10.0
    ) -> float:
        """
        Merger-induced starburst SFR.
        
        Ṁ_burst = (M_gas / t_orb) × ε_burst × enhancement
        
        Args:
            M_gas: Gas mass (kg)
            t_orb: Orbital time (s)
            epsilon_burst: Burst efficiency (higher than quiescent)
            enhancement: Enhancement factor (10-100)
            
        Returns:
            SFR_burst: Starburst rate (kg/s)
        """
        if t_orb <= 0:
            return 0
        M_dot_inflow = M_gas / t_orb
        return M_dot_inflow * epsilon_burst * enhancement
    
    @staticmethod
    def cosmic_SFRD(z: float) -> float:
        """
        Cosmic star formation rate density evolution.
        
        SFRD ∝ (1+z)^2.7 to z~2, then flat/decline
        
        Args:
            z: Redshift
            
        Returns:
            SFRD: Relative SFRD (normalized to z=0)
        """
        if z <= 2:
            return (1 + z)**2.7
        else:
            # Plateau/decline at z > 2
            return (1 + 2)**2.7 * math.exp(-(z - 2) / 2)
    
    @staticmethod
    def toomre_Q(sigma_v: float, Omega: float, Sigma: float) -> float:
        """
        Toomre stability parameter.
        
        Q = σ_v × Ω / (π G Σ)
        Q < 1 → disk instability → bursts
        
        Args:
            sigma_v: Velocity dispersion (m/s)
            Omega: Angular frequency (rad/s)
            Sigma: Surface density (kg/m²)
            
        Returns:
            Q: Toomre parameter
        """
        if Sigma <= 0:
            return math.inf
        return sigma_v * Omega / (math.pi * G * Sigma)


# ==============================================================================
# BLACK HOLE ACCRETION CALCULATOR (Equation 9)
# ==============================================================================

class BlackHoleAccretionCalculator:
    """
    Black hole accretion rate calculations.
    
    Equation 9: Eddington-limited accretion
    Ṁ_BH = 4π G M_BH m_p / (ε_r σ_T c)
    
    Where:
    - ε_r ~ 0.1 (radiative efficiency)
    - σ_T = Thomson cross-section
    - Growth timescale t_Sal ≈ 45 Myr
    """
    
    EPSILON_R = 0.1  # Radiative efficiency
    
    @staticmethod
    def eddington_luminosity(M_BH: float) -> float:
        """
        Eddington luminosity.
        
        L_Edd = 4π G M_BH m_p c / σ_T
        
        Args:
            M_BH: Black hole mass (kg)
            
        Returns:
            L_Edd: Eddington luminosity (W)
        """
        return 4 * math.pi * G * M_BH * m_p * c / sigma_T
    
    @staticmethod
    def eddington_accretion_rate(M_BH: float, epsilon_r: float = 0.1) -> float:
        """
        Eddington-limited accretion rate.
        
        Ṁ_Edd = L_Edd / (ε_r c²) = 4π G M_BH m_p / (ε_r σ_T c)
        
        Args:
            M_BH: Black hole mass (kg)
            epsilon_r: Radiative efficiency
            
        Returns:
            M_dot: Accretion rate (kg/s)
        """
        return 4 * math.pi * G * M_BH * m_p / (epsilon_r * sigma_T * c)
    
    @staticmethod
    def salpeter_time(epsilon_r: float = 0.1) -> float:
        """
        Salpeter (e-folding) time for BH growth.
        
        t_Sal = ε_r σ_T c / (4π G m_p) ≈ 45 Myr
        
        Args:
            epsilon_r: Radiative efficiency
            
        Returns:
            t_Sal: Salpeter time (s)
        """
        return epsilon_r * sigma_T * c / (4 * math.pi * G * m_p)
    
    @staticmethod
    def mass_growth(M_0: float, t: float, t_Sal: float) -> float:
        """
        Exponential mass growth at Eddington rate.
        
        M(t) = M_0 × exp(t / t_Sal)
        
        Args:
            M_0: Initial mass (kg)
            t: Time (s)
            t_Sal: Salpeter time (s)
            
        Returns:
            M: Final mass (kg)
        """
        if t_Sal <= 0:
            return M_0
        return M_0 * math.exp(t / t_Sal)


# ==============================================================================
# GRAVITATIONAL WAVE CHIRP MASS (Equations 12-13)
# ==============================================================================

class GWChirpMassCalculator:
    """
    Gravitational wave chirp mass and inspiral calculations.
    
    Equation 12: Chirp mass from inspiral
    M = (m1 m2)^(3/5) / (m1 + m2)^(1/5) 
      = (c³/G) × (5/96 π^(-8/3) f^(-11/3) ḟ)^(3/5)
    
    Equation 13: Ringdown QNM frequency
    f_QNM = (c³/2πGM_f) × (0.3737 + 0.088 a_f + ...)
    """
    
    @staticmethod
    def chirp_mass(m1: float, m2: float) -> float:
        """
        Chirp mass from component masses.
        
        M = (m1 m2)^(3/5) / (m1 + m2)^(1/5)
        
        Args:
            m1, m2: Component masses (kg)
            
        Returns:
            M_chirp: Chirp mass (kg)
        """
        return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    
    @staticmethod
    def chirp_mass_from_GW(f: float, f_dot: float) -> float:
        """
        Chirp mass from GW frequency and derivative.
        
        M = (c³/G) × (5/96 π^(-8/3) f^(-11/3) ḟ)^(3/5)
        
        Args:
            f: GW frequency (Hz)
            f_dot: Frequency derivative (Hz/s)
            
        Returns:
            M_chirp: Chirp mass (kg)
        """
        factor = (5/96) * math.pi**(-8/3) * f**(-11/3) * f_dot
        return (c**3 / G) * factor**(3/5)
    
    @staticmethod
    def ringdown_frequency(M_f: float, a_f: float) -> float:
        """
        Ringdown quasi-normal mode frequency.
        
        f_QNM = (c³/2πGM_f) × (0.3737 + 0.088 a_f + ...)
        
        Args:
            M_f: Final (remnant) BH mass (kg)
            a_f: Dimensionless spin (0 to ~0.998)
            
        Returns:
            f_QNM: QNM frequency (Hz)
        """
        f_factor = 0.3737 + 0.088 * a_f  # Leading terms
        return (c**3 / (2 * math.pi * G * M_f)) * f_factor
    
    @staticmethod
    def inspiral_frequency_evolution(M_chirp: float, f: float) -> float:
        """
        GW frequency evolution rate.
        
        ḟ = (96/5) π^(8/3) (G M/c³)^(5/3) f^(11/3)
        
        Args:
            M_chirp: Chirp mass (kg)
            f: GW frequency (Hz)
            
        Returns:
            f_dot: Frequency derivative (Hz/s)
        """
        return (96/5) * math.pi**(8/3) * (G * M_chirp / c**3)**(5/3) * f**(11/3)
    
    @staticmethod
    def merger_time(M_chirp: float, f_initial: float) -> float:
        """
        Time to merger from initial frequency.
        
        t_merge = (5/256) (c⁵/G⁵/³) (1/(π f_i)^(8/3)) M^(-5/3)
        
        Args:
            M_chirp: Chirp mass (kg)
            f_initial: Initial GW frequency (Hz)
            
        Returns:
            t_merge: Time to merger (s)
        """
        return (5/256) * (c**5 / G**(5/3)) * \
               (1 / (math.pi * f_initial)**(8/3)) * M_chirp**(-5/3)


# ==============================================================================
# UQFF-TAILORED CALCULATORS
# ==============================================================================

class UQFFMergerBuoyancyCalculator:
    """
    UQFF Buoyancy in Merger-Induced Bursts (F_UBii,merger).
    
    Positive for SFR enhancement:
    
    F_UBii,merger = +F_rel × (E_burst/E_LEP) × Q_wave,z × g_halo(z) × (1+z)^m
    
    Where:
    - F_rel = 4.30×10³³ N
    - E_burst ~ 10-100 × quiescent SFR energy
    - Q_wave,z = redshift-modulated resonance
    - g_halo = G M / r² (virial)
    - m ~ 2.5 (EPS exponent)
    
    Positive sign → enhancement/instability (drives bursts)
    """
    
    @staticmethod
    def burst_energy(SFR_burst: float, t_burst: float = 1e8 * 3.156e7) -> float:
        """
        Starburst energy output.
        
        E_burst = SFR × ε_nuc × c² × t_burst (ε_nuc ~ 0.007 for H fusion)
        
        Args:
            SFR_burst: Star formation rate (kg/s)
            t_burst: Burst duration (s), default 100 Myr
            
        Returns:
            E_burst: Total energy (J)
        """
        epsilon_nuc = 0.007  # Nuclear efficiency
        return SFR_burst * epsilon_nuc * c**2 * t_burst
    
    @staticmethod
    def redshift_wave_factor(z: float) -> float:
        """
        Redshift-modulated resonance factor.
        
        Q_wave,z = 10¹² × (1 + z)^0.5
        
        Args:
            z: Redshift
            
        Returns:
            Q_wave: Resonance factor
        """
        return Q_WAVE_BASE * (1 + z)**0.5
    
    @staticmethod
    def halo_gravity(M_halo: float, r_virial: float) -> float:
        """
        Virial halo gravity.
        
        g_halo = G M / r²
        
        Args:
            M_halo: Halo mass (kg)
            r_virial: Virial radius (m)
            
        Returns:
            g: Gravitational acceleration (m/s²)
        """
        return G * M_halo / r_virial**2
    
    @staticmethod
    def compute_F_UBii_merger(
        M_halo: float,
        r_virial: float,
        SFR_burst: float,
        z: float,
        m_exponent: float = 2.5
    ) -> float:
        """
        Compute UQFF buoyancy for merger enhancement.
        
        F_UBii,merger = +F_rel × (E_burst/E_LEP) × Q_wave,z × g_halo × (1+z)^m
        
        Positive → drives instability/bursts
        
        Args:
            M_halo: Halo mass (kg)
            r_virial: Virial radius (m)
            SFR_burst: Starburst SFR (kg/s)
            z: Redshift
            m_exponent: Redshift evolution exponent
            
        Returns:
            F_UBii: Buoyancy force (N), positive for enhancement
        """
        E_burst = UQFFMergerBuoyancyCalculator.burst_energy(SFR_burst)
        Q_wave = UQFFMergerBuoyancyCalculator.redshift_wave_factor(z)
        g_halo = UQFFMergerBuoyancyCalculator.halo_gravity(M_halo, r_virial)
        
        energy_ratio = E_burst / E_LEP if E_LEP > 0 else 0
        z_factor = (1 + z)**m_exponent
        
        F_UBii = F_REL * energy_ratio * Q_wave * g_halo * z_factor
        return F_UBii  # Positive for enhancement


class UQFFBHAccretionBuoyancyCalculator:
    """
    UQFF Buoyancy in BH Accretion (F_UBii,BH).
    
    Negative for feedback regulation:
    
    F_UBii,BH = -F_rel × (Ṁ_BH c²/E_LEP) × Q_wave × (4πGM/c²r) × erfc(δ_c/(√2 σ))
    
    Negative sign → feedback regulation (limits growth)
    """
    
    @staticmethod
    def compute_F_UBii_BH(
        M_BH: float,
        r: float,
        M_dot: float,
        z: float
    ) -> float:
        """
        Compute UQFF buoyancy for BH accretion feedback.
        
        F_UBii,BH = -F_rel × (Ṁc²/E_LEP) × Q_wave × (4πGM/c²r) × erfc(...)
        
        Negative → limits accretion via jet feedback
        
        Args:
            M_BH: Black hole mass (kg)
            r: Characteristic radius (m)
            M_dot: Accretion rate (kg/s)
            z: Redshift
            
        Returns:
            F_UBii: Buoyancy force (N), negative for regulation
        """
        # Energy ratio
        E_acc = M_dot * c**2  # Accretion luminosity
        energy_ratio = E_acc / E_LEP if E_LEP > 0 else 0
        
        # Wave factor
        Q_wave = Q_WAVE_BASE * (1 + z)**0.3
        
        # Gravitational factor (Schwarzschild-like)
        r_S = 2 * G * M_BH / c**2  # Schwarzschild radius
        g_factor = 4 * math.pi * G * M_BH / (c**2 * r) if r > 0 else 0
        
        # EPS probability factor
        delta_c_z = EPSMergerCalculator.critical_overdensity(z)
        sigma = EPSMergerCalculator.sigma_M(M_BH, z)
        erfc_factor = EPSMergerCalculator.erfc(delta_c_z / (math.sqrt(2) * sigma)) if sigma > 0 else 1
        
        # Negative for feedback regulation
        F_UBii = -F_REL * energy_ratio * Q_wave * g_factor * erfc_factor
        return F_UBii


# ==============================================================================
# UNIFIED GALAXY/MERGER CALCULATOR
# ==============================================================================

class GalaxyMergerCalculator:
    """Unified calculator for galaxy mergers and SFR."""
    
    def __init__(self):
        self.eps = EPSMergerCalculator()
        self.sfr = StarFormationCalculator()
        self.bh = BlackHoleAccretionCalculator()
        self.gw = GWChirpMassCalculator()
        self.uqff_merger = UQFFMergerBuoyancyCalculator()
        self.uqff_bh = UQFFBHAccretionBuoyancyCalculator()
    
    def compute(self, system: GalaxySystem) -> MergerResult:
        """
        Complete galaxy/merger physics calculation.
        
        Args:
            system: GalaxySystem parameters
            
        Returns:
            MergerResult with all computed values
        """
        # Dynamical times
        rho_mean = 3 * system.M_halo / (4 * math.pi * system.r_virial**3)
        t_dyn = self.sfr.dynamical_time(rho_mean)
        t_orb = self.sfr.orbital_time(system.r_virial, system.M_halo)
        
        # Star formation rates
        SFR_quiescent = self.sfr.quiescent_SFR(system.M_gas, t_dyn)
        SFR_burst = self.sfr.starburst_SFR(system.M_gas, t_orb)
        
        # Merger rate scaling
        merger_rate = self.eps.merger_rate_density(system.M_halo, system.z)
        
        # BH properties
        M_dot_BH = self.bh.eddington_accretion_rate(system.M_BH)
        t_Sal = self.bh.salpeter_time()
        
        # EPS BH number
        rho_bar = Omega_m * 3 * H_0**2 / (8 * math.pi * G)  # Mean cosmic density
        N_BH = self.eps.cumulative_BH_number(system.M_BH, system.z, rho_bar)
        
        # UQFF terms
        F_UBii = self.uqff_merger.compute_F_UBii_merger(
            system.M_halo, system.r_virial, SFR_burst, system.z
        )
        
        return MergerResult(
            merger_rate=merger_rate,
            SFR_quiescent=SFR_quiescent,
            SFR_burst=SFR_burst,
            t_dyn=t_dyn,
            N_BH_cumulative=N_BH,
            M_dot_BH=M_dot_BH,
            t_Sal=t_Sal,
            F_UBii=F_UBii,
            additional={
                't_orbital': t_orb,
                'L_Edd': self.bh.eddington_luminosity(system.M_BH),
                'cosmic_SFRD': self.sfr.cosmic_SFRD(system.z),
                'sigma_M': self.eps.sigma_M(system.M_halo, system.z)
            }
        )


# ==============================================================================
# PRE-DEFINED GALAXY SYSTEMS
# ==============================================================================

# Milky Way-like galaxy
MILKY_WAY_SYSTEM = GalaxySystem(
    name="Milky Way (z=0)",
    M_halo=1e12 * M_sun,
    M_stellar=6e10 * M_sun,
    M_gas=1e10 * M_sun,
    z=0.0,
    SFR=2.0 * M_sun / (3.156e7),  # 2 M_sun/yr
    M_BH=4e6 * M_sun,             # Sgr A*
    r_virial=200e3 * 3.086e16,    # 200 kpc
    sigma_v=200e3,                # 200 km/s
    description="Local Universe Milky Way-type disk galaxy"
)

# High-z starburst (z~2 cosmic noon)
HIGH_Z_STARBURST = GalaxySystem(
    name="z=2 Starburst",
    M_halo=5e11 * M_sun,
    M_stellar=1e10 * M_sun,
    M_gas=5e10 * M_sun,           # Gas-rich
    z=2.0,
    SFR=100 * M_sun / (3.156e7),  # 100 M_sun/yr
    M_BH=1e8 * M_sun,
    r_virial=50e3 * 3.086e16,     # 50 kpc
    sigma_v=250e3,
    description="Cosmic noon (z=2) gas-rich starburst"
)

# Merger remnant (recent major merger)
MERGER_REMNANT = GalaxySystem(
    name="Recent Merger Remnant",
    M_halo=2e12 * M_sun,
    M_stellar=2e11 * M_sun,
    M_gas=3e10 * M_sun,
    z=0.5,
    SFR=50 * M_sun / (3.156e7),   # Post-starburst
    M_BH=1e9 * M_sun,
    r_virial=300e3 * 3.086e16,
    sigma_v=300e3,
    description="Post-merger system with enhanced BH growth"
)

# High-z quasar host (z~5)
HIGH_Z_QUASAR = GalaxySystem(
    name="z=5 Quasar Host",
    M_halo=1e13 * M_sun,
    M_stellar=1e11 * M_sun,
    M_gas=1e11 * M_sun,
    z=5.0,
    SFR=500 * M_sun / (3.156e7),  # Extreme SFR
    M_BH=1e9 * M_sun,             # Already massive BH
    r_virial=100e3 * 3.086e16,
    sigma_v=400e3,
    description="Early universe massive quasar host"
)

# Binary BH merger progenitor
BINARY_BH_SYSTEM = GalaxySystem(
    name="Binary BH Host (GW Source)",
    M_halo=5e11 * M_sun,
    M_stellar=5e10 * M_sun,
    M_gas=5e9 * M_sun,
    z=0.2,
    SFR=5 * M_sun / (3.156e7),
    M_BH=30 * M_sun,              # Stellar-mass BH
    r_virial=100e3 * 3.086e16,
    sigma_v=180e3,
    description="Host of stellar-mass BBH merger"
)

# ==============================================================================
# REGISTRY FOR CONDENSEDPHYSICS INTEGRATION
# ==============================================================================

GALAXY_MERGER_CALCULATOR = GalaxyMergerCalculator()

GALAXY_SYSTEMS = {
    'MILKY_WAY': MILKY_WAY_SYSTEM,
    'HIGH_Z_STARBURST': HIGH_Z_STARBURST,
    'MERGER_REMNANT': MERGER_REMNANT,
    'HIGH_Z_QUASAR': HIGH_Z_QUASAR,
    'BINARY_BH': BINARY_BH_SYSTEM,
}

MERGER_SFR_CALCULATORS = {
    'EPSMerger': EPSMergerCalculator,
    'StarFormation': StarFormationCalculator,
    'BlackHoleAccretion': BlackHoleAccretionCalculator,
    'GWChirpMass': GWChirpMassCalculator,
    'UQFFMergerBuoyancy': UQFFMergerBuoyancyCalculator,
    'UQFFBHAccretion': UQFFBHAccretionBuoyancyCalculator,
    'GalaxyMerger': GalaxyMergerCalculator,
}

# ==============================================================================
# QUICK TEST
# ==============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("GALAXY MERGERS AND SFR MODULE TEST")
    print("=" * 70)
    
    calc = GalaxyMergerCalculator()
    
    for name, system in GALAXY_SYSTEMS.items():
        result = calc.compute(system)
        SFR_solar = result.SFR_quiescent * 3.156e7 / M_sun  # Convert to M_sun/yr
        SFR_burst_solar = result.SFR_burst * 3.156e7 / M_sun
        
        print(f"\n{system.name} (z={system.z}):")
        print(f"  SFR_quiescent = {SFR_solar:.2f} M_sun/yr")
        print(f"  SFR_burst = {SFR_burst_solar:.1f} M_sun/yr")
        print(f"  t_dyn = {result.t_dyn/Gyr:.3f} Gyr")
        print(f"  t_Sal = {result.t_Sal/Gyr*1e3:.1f} Myr")
        print(f"  Merger rate factor = {result.merger_rate:.2f}")
        print(f"  F_UBii = {result.F_UBii:.2e} N (positive=enhancing)")
