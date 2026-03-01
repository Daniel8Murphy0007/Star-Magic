"""
cosmic_domains_calculators.py
Cosmic Domains Calculator Module

Comprehensive calculators for astrophysical domains extracted from Grok AI conversation:
- Interstellar Medium (ISM)
- Stellar Evolution
- Big Bang Nucleosynthesis (BBN)
- Cosmic Voids
- Reionization
- AGN Outflows
- Binary Pulsars
- Cosmic Rays
- Intergalactic Medium (IGM)
- Dark Energy
- Black Hole Thermodynamics
- Loop Quantum Cosmology (LQC)

All calculators follow CondensedPhysics.py architecture pattern.

© 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple, List
from dataclasses import dataclass, field
from enum import Enum
import math


# ═══════════════════════════════════════════════════════════════════════════════
# PHYSICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

CONSTANTS = {
    # Fundamental
    'G': 6.6743e-11,           # Gravitational constant (m³/kg·s²)
    'c': 2.998e8,              # Speed of light (m/s)
    'hbar': 1.0546e-34,        # Reduced Planck constant (J·s)
    'k_B': 1.381e-23,          # Boltzmann constant (J/K)
    'e': 1.602e-19,            # Elementary charge (C)
    'm_p': 1.673e-27,          # Proton mass (kg)
    'm_e': 9.109e-31,          # Electron mass (kg)
    'sigma_T': 6.652e-29,      # Thomson cross-section (m²)
    'sigma_SB': 5.67e-8,       # Stefan-Boltzmann constant (W/m²/K⁴)
    'a_rad': 7.5657e-16,       # Radiation constant (J/m³/K⁴)
    
    # Astrophysical
    'M_sun': 1.989e30,         # Solar mass (kg)
    'R_sun': 6.96e8,           # Solar radius (m)
    'L_sun': 3.828e26,         # Solar luminosity (W)
    'AU': 1.496e11,            # Astronomical unit (m)
    'pc': 3.086e16,            # Parsec (m)
    'Mpc': 3.086e22,           # Megaparsec (m)
    'ly': 9.461e15,            # Light year (m)
    
    # Cosmological
    'H_0': 2.2e-18,            # Hubble constant (1/s) ≈ 70 km/s/Mpc
    'Omega_m': 0.3,            # Matter density parameter
    'Omega_Lambda': 0.7,       # Dark energy density parameter
    'Omega_b': 0.05,           # Baryon density parameter
    'T_CMB': 2.725,            # CMB temperature (K)
    'rho_crit': 9.47e-27,      # Critical density (kg/m³)
    
    # UQFF Calibration
    'kappa_UQFF': 0.0005,      # κ decay rate (1/day)
    'SSq': 0.57,               # Superconductive state factor
    'rho_vac_UA': 7.09e-36,    # Vacuum density UA (J/m³)
    
    # Nuclear
    'Q_pn': 1.293e6 * 1.602e-19,  # n-p mass diff (J)
    'tau_n': 879.4,            # Neutron lifetime (s)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. INTERSTELLAR MEDIUM (ISM) CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class ISMPhaseCalculator:
    """
    Interstellar Medium phase structure calculator.
    
    ISM phases:
    - Cold Neutral Medium (CNM): T ~ 100 K, n ~ 30 cm⁻³
    - Warm Neutral Medium (WNM): T ~ 8000 K, n ~ 0.5 cm⁻³
    - Warm Ionized Medium (WIM): T ~ 8000 K, n ~ 0.3 cm⁻³
    - Hot Ionized Medium (HIM): T ~ 10⁶ K, n ~ 0.003 cm⁻³
    
    Pressure equilibrium: P = nkT ≈ 3000 k_B K/cm³
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        self.phases = {
            'CNM': {'T': 100, 'n': 30},
            'WNM': {'T': 8000, 'n': 0.5},
            'WIM': {'T': 8000, 'n': 0.3},
            'HIM': {'T': 1e6, 'n': 0.003}
        }
        
    def compute_pressure(self, T: float, n: float) -> float:
        """
        Compute thermal pressure.
        
        P = nkT
        
        Args:
            T: Temperature (K)
            n: Number density (cm⁻³)
            
        Returns:
            Pressure (Pa)
        """
        k_B = self.constants['k_B']
        n_SI = n * 1e6  # Convert cm⁻³ to m⁻³
        return n_SI * k_B * T
    
    def compute_cooling_function(self, T: float, Z: float = 1.0) -> float:
        """
        Compute radiative cooling function Λ(T).
        
        Λ(T) ∝ T^α where α depends on dominant process:
        - T < 10⁴ K: Metal line cooling, α ≈ -0.7
        - 10⁴ < T < 10⁷ K: Bremsstrahlung + lines, α ≈ -0.5
        - T > 10⁷ K: Bremsstrahlung, α ≈ 0.5
        
        Args:
            T: Temperature (K)
            Z: Metallicity relative to solar
            
        Returns:
            Cooling rate (erg cm³/s)
        """
        if T < 1e4:
            Lambda = 1e-26 * Z * (T / 1e4)**(-0.7)
        elif T < 1e7:
            Lambda = 1e-23 * Z * (T / 1e6)**(-0.5)
        else:
            Lambda = 1e-23 * (T / 1e6)**0.5
            
        return Lambda
    
    def compute_jeans_mass(self, T: float, n: float, mu: float = 1.0) -> float:
        """
        Compute Jeans mass for gravitational collapse.
        
        M_J = (5kT / Gμm_p)^(3/2) × (3 / 4πρ)^(1/2)
        
        Args:
            T: Temperature (K)
            n: Number density (cm⁻³)
            mu: Mean molecular weight
            
        Returns:
            Jeans mass (kg)
        """
        G = self.constants['G']
        k_B = self.constants['k_B']
        m_p = self.constants['m_p']
        
        rho = n * 1e6 * mu * m_p  # Mass density (kg/m³)
        
        M_J = (5 * k_B * T / (G * mu * m_p))**1.5 * np.sqrt(3 / (4 * np.pi * rho))
        
        return M_J
    
    def compute_sound_speed(self, T: float, gamma: float = 5/3,
                            mu: float = 1.0) -> float:
        """
        Compute isothermal/adiabatic sound speed.
        
        c_s = √(γkT / μm_p)
        
        Args:
            T: Temperature (K)
            gamma: Adiabatic index
            mu: Mean molecular weight
            
        Returns:
            Sound speed (m/s)
        """
        k_B = self.constants['k_B']
        m_p = self.constants['m_p']
        
        return np.sqrt(gamma * k_B * T / (mu * m_p))


# ═══════════════════════════════════════════════════════════════════════════════
# 2. STELLAR EVOLUTION CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class StellarEvolutionCalculator:
    """
    Stellar evolution equations calculator.
    
    Key relations:
    - Main sequence lifetime: t_MS ~ 10¹⁰ (M/M_⊙)^(-2.5) yr
    - Mass-luminosity: L ~ M^3.5 (main sequence)
    - Chandrasekhar limit: M_Ch ~ 1.44 M_⊙
    - Eddington luminosity: L_Edd ~ 3.2×10⁴ (M/M_⊙) L_⊙
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_main_sequence_lifetime(self, M: float) -> float:
        """
        Compute main sequence lifetime.
        
        t_MS ≈ 10¹⁰ × (M_⊙/M)^2.5 yr
        
        Args:
            M: Stellar mass (kg)
            
        Returns:
            MS lifetime (s)
        """
        M_sun = self.constants['M_sun']
        t_sun = 10e9 * 365.25 * 86400  # 10 Gyr in seconds
        
        return t_sun * (M_sun / M)**2.5
    
    def compute_luminosity(self, M: float, phase: str = 'MS') -> float:
        """
        Compute stellar luminosity.
        
        Main sequence: L ∝ M^3.5
        RGB: L ∝ M^6
        
        Args:
            M: Stellar mass (kg)
            phase: Evolution phase ('MS', 'RGB', 'AGB')
            
        Returns:
            Luminosity (W)
        """
        L_sun = self.constants['L_sun']
        M_sun = self.constants['M_sun']
        
        if phase == 'MS':
            return L_sun * (M / M_sun)**3.5
        elif phase == 'RGB':
            return L_sun * (M / M_sun)**6.0
        else:
            return L_sun * (M / M_sun)**4.0
    
    def compute_radius(self, M: float, phase: str = 'MS') -> float:
        """
        Compute stellar radius.
        
        Main sequence: R ∝ M^0.8
        
        Args:
            M: Stellar mass (kg)
            phase: Evolution phase
            
        Returns:
            Radius (m)
        """
        R_sun = self.constants['R_sun']
        M_sun = self.constants['M_sun']
        
        if phase == 'MS':
            return R_sun * (M / M_sun)**0.8
        elif phase == 'RGB':
            return R_sun * 100 * (M / M_sun)**0.5
        else:
            return R_sun * (M / M_sun)**0.6
    
    def compute_effective_temperature(self, L: float, R: float) -> float:
        """
        Compute effective temperature from Stefan-Boltzmann.
        
        L = 4πR²σT_eff⁴
        
        Args:
            L: Luminosity (W)
            R: Radius (m)
            
        Returns:
            Effective temperature (K)
        """
        sigma_SB = self.constants['sigma_SB']
        
        return (L / (4 * np.pi * R**2 * sigma_SB))**0.25
    
    def compute_eddington_luminosity(self, M: float, kappa: float = 0.34) -> float:
        """
        Compute Eddington luminosity.
        
        L_Edd = 4πGMc/κ
        
        Args:
            M: Stellar mass (kg)
            kappa: Opacity (m²/kg), default electron scattering
            
        Returns:
            Eddington luminosity (W)
        """
        G = self.constants['G']
        c = self.constants['c']
        
        return 4 * np.pi * G * M * c / kappa
    
    def compute_chandrasekhar_mass(self, mu_e: float = 2.0) -> float:
        """
        Compute Chandrasekhar mass limit.
        
        M_Ch ≈ 5.83/μ_e² M_⊙
        
        Args:
            mu_e: Mean molecular weight per electron
            
        Returns:
            Chandrasekhar mass (kg)
        """
        M_sun = self.constants['M_sun']
        return 5.83 / mu_e**2 * M_sun


# ═══════════════════════════════════════════════════════════════════════════════
# 3. BIG BANG NUCLEOSYNTHESIS (BBN) CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class BBNCalculator:
    """
    Big Bang Nucleosynthesis calculator.
    
    Key reactions:
    - n ↔ p + e⁻ + ν̄_e (freeze-out at T ~ 0.8 MeV)
    - p + n → D + γ
    - D + D → ³He + n, ³H + p
    - ³He + D → ⁴He + p
    
    Primordial abundances:
    - Y_p (⁴He mass fraction) ~ 0.247
    - D/H ~ 2.5 × 10⁻⁵
    - ⁷Li/H ~ 1.6 × 10⁻¹⁰ (lithium problem!)
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_freeze_out_temperature(self) -> float:
        """
        Compute weak interaction freeze-out temperature.
        
        Freeze-out when Γ_weak ~ H
        T_f ~ 0.8 MeV
        
        Returns:
            Freeze-out temperature (K)
        """
        # T_f ≈ 0.8 MeV ≈ 9.3 × 10⁹ K
        return 0.8e6 * self.constants['e'] / self.constants['k_B']
    
    def compute_neutron_to_proton_ratio(self, T: float) -> float:
        """
        Compute equilibrium n/p ratio.
        
        n/p = exp(-Q/kT) where Q = (m_n - m_p)c²
        
        Args:
            T: Temperature (K)
            
        Returns:
            n/p ratio
        """
        Q = self.constants['Q_pn']
        k_B = self.constants['k_B']
        
        return np.exp(-Q / (k_B * T))
    
    def compute_helium_mass_fraction(self, n_p_ratio: float,
                                      delta_n: float = 0.0) -> float:
        """
        Compute primordial helium-4 mass fraction.
        
        Y_p = 2(n/p) / (1 + n/p) × (1 - decay losses)
        
        Assuming all neutrons end up in ⁴He.
        
        Args:
            n_p_ratio: Neutron to proton ratio at nucleosynthesis
            delta_n: Neutron decay correction
            
        Returns:
            Y_p mass fraction
        """
        # Account for neutron decay before nucleosynthesis
        n_p_eff = n_p_ratio * (1 - delta_n)
        
        return 2 * n_p_eff / (1 + n_p_eff)
    
    def compute_deuterium_abundance(self, eta: float) -> float:
        """
        Compute primordial deuterium abundance.
        
        D/H depends sensitively on baryon-to-photon ratio η.
        D/H ≈ 2.6 × 10⁻⁵ × (η / 6×10⁻¹⁰)^(-1.6)
        
        Args:
            eta: Baryon-to-photon ratio
            
        Returns:
            D/H ratio
        """
        eta_0 = 6e-10  # Standard BBN value
        return 2.6e-5 * (eta / eta_0)**(-1.6)
    
    def compute_lithium_abundance(self, eta: float) -> float:
        """
        Compute primordial lithium-7 abundance.
        
        ⁷Li/H ≈ 5 × 10⁻¹⁰ (observed ~ 1.6 × 10⁻¹⁰, lithium problem)
        
        Args:
            eta: Baryon-to-photon ratio
            
        Returns:
            ⁷Li/H ratio (theoretical)
        """
        eta_0 = 6e-10
        return 5e-10 * (eta / eta_0)**2.0


# ═══════════════════════════════════════════════════════════════════════════════
# 4. COSMIC VOID CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class CosmicVoidCalculator:
    """
    Cosmic void structure and evolution calculator.
    
    Typical voids:
    - Radius: 10-50 Mpc
    - Density contrast: δ ≈ -0.8 to -0.95
    - Expansion rate: ~10% faster than Hubble
    
    Dynamics: ṙ = H(z)r × (1 - δ/3)
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_void_profile(self, r: float, R_void: float,
                             delta_c: float = -0.8) -> float:
        """
        Compute void density profile.
        
        Commonly modeled as:
        δ(r) = δ_c × (1 - (r/R_v)²) for r < R_v
        
        Args:
            r: Radius (Mpc)
            R_void: Void radius (Mpc)
            delta_c: Central density contrast
            
        Returns:
            Density contrast δ
        """
        if r >= R_void:
            return 0.0
        return delta_c * (1 - (r / R_void)**2)
    
    def compute_void_expansion_rate(self, H: float, delta: float) -> float:
        """
        Compute local expansion rate in void.
        
        H_local = H × (1 + f(Ω_m)δ/3)
        where f ≈ Ω_m^0.55 is the growth rate
        
        Args:
            H: Background Hubble rate (1/s)
            delta: Local density contrast
            
        Returns:
            Local Hubble rate (1/s)
        """
        Omega_m = self.constants['Omega_m']
        f = Omega_m**0.55
        
        return H * (1 + f * delta / 3)
    
    def compute_integrated_sachs_wolfe(self, delta: float,
                                        R_void: float) -> float:
        """
        Compute ISW effect from void.
        
        ΔT/T ~ 2 ∫ d(Φ)/dt × dl/c
        
        For a void: ΔT/T ~ (Ω_Λ/Ω_m) × δ × (R/c/H_0)²
        
        Args:
            delta: Void density contrast
            R_void: Void radius (Mpc)
            
        Returns:
            Temperature fluctuation ΔT/T
        """
        Omega_Lambda = self.constants['Omega_Lambda']
        Omega_m = self.constants['Omega_m']
        H_0 = self.constants['H_0']
        c = self.constants['c']
        Mpc = self.constants['Mpc']
        
        R = R_void * Mpc
        
        return (Omega_Lambda / Omega_m) * abs(delta) * (R * H_0 / c)**2


# ═══════════════════════════════════════════════════════════════════════════════
# 5. REIONIZATION CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class ReionizationCalculator:
    """
    Cosmic reionization calculator.
    
    Epoch: z ~ 6-12
    Sources: First stars (Pop III) and galaxies
    
    Key physics:
    - Ionization rate: Γ_HI ~ n_γ σ_HI c
    - Recombination time: t_rec ~ 1 / (α_B n_e)
    - Stromgren sphere: R_S = (3N_γ / 4π n²α_B)^(1/3)
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        self.alpha_B = 2.6e-19  # Case B recombination coefficient (m³/s)
        
    def compute_ionization_fraction(self, n_gamma: float, n_H: float,
                                     t: float) -> float:
        """
        Compute ionization fraction evolution.
        
        dx_e/dt = (1-x_e)Γ - x_e²n_H α_B
        
        Equilibrium: x_e ~ √(Γ / n_H α_B)
        
        Args:
            n_gamma: Ionizing photon density (m⁻³)
            n_H: Hydrogen density (m⁻³)
            t: Time (s)
            
        Returns:
            Ionization fraction x_e
        """
        sigma_HI = 6.3e-22  # HI photoionization cross-section at 13.6 eV (m²)
        c = self.constants['c']
        
        Gamma = n_gamma * sigma_HI * c
        
        x_e_eq = np.sqrt(Gamma / (n_H * self.alpha_B))
        
        return min(x_e_eq, 1.0)
    
    def compute_stromgren_radius(self, N_dot: float, n_H: float) -> float:
        """
        Compute Strömgren sphere radius.
        
        R_S = (3Ṅ / 4π n²α_B)^(1/3)
        
        Args:
            N_dot: Ionizing photon rate (photons/s)
            n_H: Hydrogen density (m⁻³)
            
        Returns:
            Strömgren radius (m)
        """
        return (3 * N_dot / (4 * np.pi * n_H**2 * self.alpha_B))**(1/3)
    
    def compute_optical_depth_reionization(self, z_re: float,
                                            Delta_z: float = 1.0) -> float:
        """
        Compute Thomson optical depth to reionization.
        
        τ_e ~ ∫ n_e σ_T c dt
        
        Simplified: τ ≈ 0.05 × ((1+z_re)/10)^1.5
        
        Args:
            z_re: Reionization redshift
            Delta_z: Reionization duration
            
        Returns:
            Optical depth τ
        """
        return 0.05 * ((1 + z_re) / 10)**1.5
    
    def compute_gunn_peterson_trough(self, x_HI: float, z: float) -> float:
        """
        Compute Gunn-Peterson optical depth.
        
        τ_GP = (π e² f λ_α / m_e c H) × n_HI
        
        Complete absorption for x_HI > 10⁻⁴
        
        Args:
            x_HI: Neutral fraction
            z: Redshift
            
        Returns:
            GP optical depth
        """
        # Simplified estimate
        tau_GP_full = 6e5 * ((1 + z) / 7)**1.5  # Full neutral
        return tau_GP_full * x_HI


# ═══════════════════════════════════════════════════════════════════════════════
# 6. AGN OUTFLOWS CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class AGNOutflowCalculator:
    """
    AGN-driven outflows and feedback calculator.
    
    Key physics:
    - Momentum-driven: Ṗ ~ L/c
    - Energy-driven: Ė ~ 0.05 L_AGN
    - Mass loading: η = Ṁ_out / SFR
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_momentum_rate(self, L_AGN: float, boost: float = 1.0) -> float:
        """
        Compute radiation momentum rate.
        
        Ṗ_rad = τ_IR L/c (momentum boost from IR trapping)
        
        Args:
            L_AGN: AGN luminosity (W)
            boost: Momentum boost factor (τ_IR)
            
        Returns:
            Momentum rate (kg m/s²)
        """
        c = self.constants['c']
        return boost * L_AGN / c
    
    def compute_outflow_velocity(self, L_AGN: float, M_out: float,
                                  R: float) -> float:
        """
        Compute momentum-driven outflow velocity.
        
        v = (2L R / Ṁc)^(1/2) for momentum-driven
        
        Args:
            L_AGN: AGN luminosity (W)
            M_out: Outflowing mass (kg)
            R: Radius (m)
            
        Returns:
            Outflow velocity (m/s)
        """
        c = self.constants['c']
        M_dot = M_out / (R / 1e3)  # Rough flow time estimate
        
        return np.sqrt(2 * L_AGN * R / (M_dot * c))
    
    def compute_energy_coupling(self, L_AGN: float,
                                 epsilon: float = 0.05) -> float:
        """
        Compute energy coupling to ISM.
        
        Ė_out ~ ε L_AGN
        
        Args:
            L_AGN: AGN luminosity (W)
            epsilon: Coupling efficiency (typically 0.01-0.1)
            
        Returns:
            Energy injection rate (W)
        """
        return epsilon * L_AGN
    
    def compute_sphere_of_influence(self, M_BH: float,
                                     sigma: float) -> float:
        """
        Compute BH sphere of influence.
        
        r_SOI = GM_BH / σ²
        
        Args:
            M_BH: Black hole mass (kg)
            sigma: Velocity dispersion (m/s)
            
        Returns:
            Sphere of influence radius (m)
        """
        G = self.constants['G']
        return G * M_BH / sigma**2


# ═══════════════════════════════════════════════════════════════════════════════
# 7. BINARY PULSAR CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class BinaryPulsarCalculator:
    """
    Binary pulsar orbital dynamics and GW emission calculator.
    
    Key effects:
    - Periastron advance: ω̇ (relativistic)
    - Gravitational redshift + time dilation: γ
    - Orbital decay: Ṗ_b (GW emission)
    - Shapiro delay: r, s
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_periastron_advance(self, M_total: float, a: float,
                                    e: float) -> float:
        """
        Compute relativistic periastron advance.
        
        ω̇ = 3(2πG M_total / c² P_b)^(2/3) / (1-e²)
        
        Args:
            M_total: Total system mass (kg)
            a: Semi-major axis (m)
            e: Eccentricity
            
        Returns:
            Periastron advance rate (rad/s)
        """
        G = self.constants['G']
        c = self.constants['c']
        
        P_b = 2 * np.pi * np.sqrt(a**3 / (G * M_total))
        
        n = 2 * np.pi / P_b  # Mean motion
        
        return 3 * (G * M_total * n / c**3)**(2/3) * n / (1 - e**2)
    
    def compute_gravitational_wave_luminosity(self, M1: float, M2: float,
                                               a: float, e: float) -> float:
        """
        Compute GW luminosity (Peters & Mathews 1963).
        
        L_GW = (32/5)(G⁴/c⁵) M₁²M₂²(M₁+M₂)/a⁵ f(e)
        
        Args:
            M1, M2: Component masses (kg)
            a: Semi-major axis (m)
            e: Eccentricity
            
        Returns:
            GW luminosity (W)
        """
        G = self.constants['G']
        c = self.constants['c']
        
        # Eccentricity enhancement factor
        f_e = (1 + 73/24 * e**2 + 37/96 * e**4) / (1 - e**2)**3.5
        
        M_total = M1 + M2
        mu = M1 * M2 / M_total  # Reduced mass
        
        return (32/5) * (G**4 / c**5) * mu**2 * M_total**3 / a**5 * f_e
    
    def compute_orbital_decay_rate(self, M1: float, M2: float,
                                    a: float, e: float) -> float:
        """
        Compute orbital decay rate from GW emission.
        
        ȧ = -(64/5)(G³/c⁵) M₁M₂(M₁+M₂)/a³ f(e)
        
        Args:
            M1, M2: Component masses (kg)
            a: Semi-major axis (m)
            e: Eccentricity
            
        Returns:
            Orbital decay rate (m/s)
        """
        G = self.constants['G']
        c = self.constants['c']
        
        M_total = M1 + M2
        
        f_e = (1 + 73/24 * e**2 + 37/96 * e**4) / (1 - e**2)**3.5
        
        return -(64/5) * (G**3 / c**5) * M1 * M2 * M_total / a**3 * f_e
    
    def compute_merger_time(self, M1: float, M2: float,
                            a: float, e: float) -> float:
        """
        Compute time to merger.
        
        t_merge ~ (5/256)(c⁵/G³)(a⁴/(M₁M₂M_total)) × g(e)
        
        Args:
            M1, M2: Component masses (kg)
            a: Initial semi-major axis (m)
            e: Initial eccentricity
            
        Returns:
            Merger time (s)
        """
        G = self.constants['G']
        c = self.constants['c']
        
        M_total = M1 + M2
        
        # Eccentricity factor (circular approximation for simplicity)
        g_e = (1 - e**2)**3.5 / (1 + 73/24 * e**2 + 37/96 * e**4)
        
        return (5/256) * (c**5 / G**3) * a**4 / (M1 * M2 * M_total) * g_e


# ═══════════════════════════════════════════════════════════════════════════════
# 8. COSMIC RAY CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class CosmicRayCalculator:
    """
    Cosmic ray spectrum and propagation calculator.
    
    Spectrum: dN/dE ∝ E^(-γ)
    - Below knee (3 PeV): γ ≈ 2.7
    - Above knee: γ ≈ 3.1
    - Above ankle (3 EeV): γ ≈ 2.6 (extragalactic)
    
    GZK cutoff at ~50 EeV
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_spectrum(self, E: float, gamma: float = 2.7,
                          norm: float = 1.0) -> float:
        """
        Compute cosmic ray flux.
        
        dN/dE = A × E^(-γ)
        
        Args:
            E: Energy (eV)
            gamma: Spectral index
            norm: Normalization
            
        Returns:
            Flux (particles/m²/s/sr/eV)
        """
        return norm * E**(-gamma)
    
    def compute_larmor_radius(self, E: float, Z: int, B: float) -> float:
        """
        Compute Larmor radius.
        
        r_L = E / (ZeB)
        
        Args:
            E: Energy (J)
            Z: Charge number
            B: Magnetic field (T)
            
        Returns:
            Larmor radius (m)
        """
        e = self.constants['e']
        return E / (Z * e * B)
    
    def compute_hillas_limit(self, B: float, R: float, Z: int,
                             beta: float = 1.0) -> float:
        """
        Compute Hillas maximum energy.
        
        E_max = ZeBRβc
        
        Args:
            B: Magnetic field (T)
            R: Source size (m)
            Z: Charge number
            beta: Velocity/c
            
        Returns:
            Maximum energy (J)
        """
        e = self.constants['e']
        c = self.constants['c']
        return Z * e * B * R * beta * c
    
    def compute_dsa_knee_energy(self, Z: int = 1, B_uG: float = 3.0,
                                 u_s_km_s: float = 3000.0, R_pc: float = 10.0,
                                 verbose: bool = False) -> dict:
        """
        Compute DSA (Diffusive Shock Acceleration) maximum energy - the Knee.
        
        Full formula:
            E_max = Z × e × B × u_s × r_g
            E_max ≈ 3×10^15 × Z × (B/3μG) × (u_s/10^3 km/s) × (R/10 pc) eV
        
        This is the Hillas criterion applied to SNR shock acceleration,
        which explains the "knee" in the cosmic ray spectrum at ~3 PeV.
        
        Args:
            Z: Atomic number (charge), default 1 (proton)
            B_uG: Magnetic field strength in microgauss (default 3 μG for ISM)
            u_s_km_s: Shock velocity in km/s (default 3000 km/s for young SNR)
            R_pc: Shock radius in parsecs (default 10 pc)
            verbose: Print detailed output
            
        Returns:
            dict with:
                E_max_eV: Maximum energy in eV
                E_max_PeV: Maximum energy in PeV
                E_max_J: Maximum energy in Joules
                equation: LaTeX-style equation string
                equation_unicode: Unicode equation string
                parameters: Input parameters used
        """
        # Convert to SI units
        B_T = B_uG * 1e-10  # μG to Tesla (1 G = 1e-4 T, 1 μG = 1e-10 T)
        u_s_m_s = u_s_km_s * 1e3  # km/s to m/s
        R_m = R_pc * 3.086e16  # pc to meters
        
        e = self.constants['e']  # 1.602e-19 C
        c = self.constants['c']  # 3e8 m/s
        
        # Full DSA formula: E_max = Z × e × B × R × u_s
        E_max_J = Z * e * B_T * R_m * u_s_m_s
        
        # Convert to eV
        E_max_eV = E_max_J / e  # J / (J/eV) = eV
        E_max_PeV = E_max_eV / 1e15
        
        # Empirical scaling formula (calibrated to observations):
        # E ≈ 3e15 × Z × (B/3μG) × (u_s/10^3 km/s) × (R/10 pc) eV
        E_scaling_eV = 3e15 * Z * (B_uG / 3.0) * (u_s_km_s / 1000.0) * (R_pc / 10.0)
        
        result = {
            'E_max_eV': E_max_eV,
            'E_max_PeV': E_max_PeV,
            'E_max_J': E_max_J,
            
            # Empirical scaling result
            'E_scaling_eV': E_scaling_eV,
            'E_scaling_PeV': E_scaling_eV / 1e15,
            
            # Equations in different formats
            'equation_latex': r"E_{\max} = Ze B u_s r_g \approx 3\times10^{15} Z \left(\frac{B}{3\mu G}\right) \left(\frac{u_s}{10^3\,\mathrm{km/s}}\right) \left(\frac{R}{10\,\mathrm{pc}}\right) \mathrm{eV}",
            'equation_unicode': "E_max = Ze B u_s r_g ≈ 3×10¹⁵ Z(B/3μG)(u_s/10³km/s)(R/10pc) eV",
            'equation_ascii': "E_max = Ze B u_s r_g ~ 3e15 * Z * (B/3uG) * (u_s/1000km/s) * (R/10pc) eV",
            
            'parameters': {
                'Z': Z,
                'B_uG': B_uG,
                'u_s_km_s': u_s_km_s,
                'R_pc': R_pc,
            },
            
            'description': "DSA maximum energy (cosmic ray knee) from SNR shock acceleration",
            'physics': "Diffusive Shock Acceleration at supernova remnant shocks"
        }
        
        if verbose:
            print("=" * 80)
            print("DSA KNEE ENERGY CALCULATION (Galactic SNR Maximum)")
            print("=" * 80)
            print(f"\nEquation: {result['equation_unicode']}")
            print(f"\nParameters:")
            print(f"  Z (atomic number)     = {Z}")
            print(f"  B (magnetic field)    = {B_uG} μG = {B_T:.2e} T")
            print(f"  u_s (shock velocity)  = {u_s_km_s} km/s = {u_s_m_s:.2e} m/s")
            print(f"  R (shock radius)      = {R_pc} pc = {R_m:.2e} m")
            print(f"\nResults:")
            print(f"  Exact:   E_max = {E_max_eV:.3e} eV = {E_max_PeV:.3f} PeV")
            print(f"  Scaling: E_max ≈ {E_scaling_eV:.3e} eV = {E_scaling_eV/1e15:.3f} PeV")
            print(f"\nNote: Scaling formula is empirically calibrated to observations.")
            print(f"      Exact formula uses E = ZeBRu_s directly.")
            print("=" * 80)
        
        return result
    
    def compute_gzk_horizon(self, E: float) -> float:
        """
        Compute GZK horizon distance.
        
        Above ~50 EeV, CR interact with CMB:
        p + γ_CMB → Δ⁺ → p + π⁰ or n + π⁺
        
        Args:
            E: Cosmic ray energy (eV)
            
        Returns:
            Horizon distance (Mpc)
        """
        E_GZK = 5e19  # GZK threshold (eV)
        
        if E < E_GZK:
            return float('inf')
            
        # Approximate horizon
        return 100 * (E_GZK / E)**0.5  # Mpc


# ═══════════════════════════════════════════════════════════════════════════════
# 9. INTERGALACTIC MEDIUM (IGM) CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class IGMCalculator:
    """
    Intergalactic medium properties calculator.
    
    IGM contains most baryons in the universe:
    - Lyα forest (T ~ 10⁴ K)
    - WHIM (T ~ 10⁵-10⁷ K, Warm-Hot Intergalactic Medium)
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_baryon_density(self, z: float) -> float:
        """
        Compute mean baryon density at redshift z.
        
        ρ_b(z) = ρ_b,0 × (1+z)³
        
        Args:
            z: Redshift
            
        Returns:
            Baryon density (kg/m³)
        """
        Omega_b = self.constants['Omega_b']
        rho_crit = self.constants['rho_crit']
        
        return Omega_b * rho_crit * (1 + z)**3
    
    def compute_jeans_filtering_mass(self, z: float, T: float = 1e4) -> float:
        """
        Compute Jeans filtering mass in IGM.
        
        M_F captures integrated thermal history effect on structure.
        
        Args:
            z: Redshift
            T: IGM temperature (K)
            
        Returns:
            Filtering mass (kg)
        """
        G = self.constants['G']
        k_B = self.constants['k_B']
        m_p = self.constants['m_p']
        H_0 = self.constants['H_0']
        
        H = H_0 * np.sqrt(self.constants['Omega_m'] * (1 + z)**3 + 
                          self.constants['Omega_Lambda'])
        
        c_s = np.sqrt(k_B * T / m_p)
        
        lambda_J = c_s / H
        rho = self.compute_baryon_density(z)
        
        return (4 * np.pi / 3) * rho * (lambda_J / 2)**3
    
    def compute_lyman_alpha_optical_depth(self, N_HI: float) -> float:
        """
        Compute Lyα optical depth for a given column density.
        
        τ = N_HI × σ_Lyα
        
        Args:
            N_HI: HI column density (m⁻²)
            
        Returns:
            Optical depth
        """
        sigma_Lya = 4.5e-22  # Lyα cross-section at line center (m²)
        return N_HI * sigma_Lya


# ═══════════════════════════════════════════════════════════════════════════════
# 10. DARK ENERGY CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class DarkEnergyCalculator:
    """
    Dark energy equation of state calculator.
    
    Parameterizations:
    - Cosmological constant: w = -1
    - CPL: w(a) = w_0 + w_a(1-a)
    - Quintessence: time-varying w
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
    def compute_equation_of_state_CPL(self, a: float, w_0: float = -1.0,
                                       w_a: float = 0.0) -> float:
        """
        Compute CPL equation of state.
        
        w(a) = w_0 + w_a(1-a)
        
        Args:
            a: Scale factor
            w_0: Present-day equation of state
            w_a: Evolution parameter
            
        Returns:
            w(a)
        """
        return w_0 + w_a * (1 - a)
    
    def compute_dark_energy_density(self, a: float, rho_de_0: float,
                                     w_0: float = -1.0,
                                     w_a: float = 0.0) -> float:
        """
        Compute dark energy density evolution.
        
        ρ_DE(a) = ρ_DE,0 × exp(3∫(1+w)/a da)
        
        Args:
            a: Scale factor
            rho_de_0: Present dark energy density (kg/m³)
            w_0, w_a: CPL parameters
            
        Returns:
            Dark energy density (kg/m³)
        """
        # For CPL: ln(ρ/ρ_0) = -3(1+w_0+w_a)ln(a) + 3w_a(a-1)
        exponent = -3 * (1 + w_0 + w_a) * np.log(a) + 3 * w_a * (a - 1)
        return rho_de_0 * np.exp(exponent)
    
    def compute_hubble_rate(self, z: float, H_0: float = None,
                            Omega_m: float = None,
                            Omega_de: float = None,
                            w: float = -1.0) -> float:
        """
        Compute Hubble rate including dark energy.
        
        H(z) = H_0 √[Ω_m(1+z)³ + Ω_DE(1+z)^(3(1+w))]
        
        Args:
            z: Redshift
            H_0: Hubble constant (1/s)
            Omega_m: Matter density parameter
            Omega_de: Dark energy density parameter
            w: Dark energy equation of state
            
        Returns:
            Hubble rate (1/s)
        """
        H_0 = H_0 or self.constants['H_0']
        Omega_m = Omega_m or self.constants['Omega_m']
        Omega_de = Omega_de or self.constants['Omega_Lambda']
        
        matter_term = Omega_m * (1 + z)**3
        de_term = Omega_de * (1 + z)**(3 * (1 + w))
        
        return H_0 * np.sqrt(matter_term + de_term)


# ═══════════════════════════════════════════════════════════════════════════════
# 11. BLACK HOLE THERMODYNAMICS CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class BlackHoleThermodynamicsCalculator:
    """
    Black hole thermodynamics calculator.
    
    Laws:
    - 0th: κ (surface gravity) constant on horizon
    - 1st: dM = κdA/(8πG) + ΩdJ + ΦdQ
    - 2nd: δA ≥ 0
    - 3rd: κ → 0 impossible in finite processes
    
    Bekenstein-Hawking entropy: S = kA/(4ℓ_P²)
    Hawking temperature: T_H = ℏκ/(2πk_B c)
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        self.l_P = np.sqrt(self.constants['hbar'] * self.constants['G'] / 
                          self.constants['c']**3)  # Planck length
        
    def compute_hawking_temperature(self, M: float) -> float:
        """
        Compute Hawking temperature.
        
        T_H = ℏc³ / (8πGMk_B)
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            Hawking temperature (K)
        """
        hbar = self.constants['hbar']
        c = self.constants['c']
        G = self.constants['G']
        k_B = self.constants['k_B']
        
        return hbar * c**3 / (8 * np.pi * G * M * k_B)
    
    def compute_bekenstein_hawking_entropy(self, M: float) -> float:
        """
        Compute Bekenstein-Hawking entropy.
        
        S = πkc³A / (2Gℏ) = 4π(GM/c²)²k / ℓ_P²
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            Entropy (J/K)
        """
        G = self.constants['G']
        c = self.constants['c']
        k_B = self.constants['k_B']
        
        r_s = 2 * G * M / c**2
        A = 4 * np.pi * r_s**2
        
        return k_B * A / (4 * self.l_P**2)
    
    def compute_hawking_luminosity(self, M: float) -> float:
        """
        Compute Hawking radiation power.
        
        L = ℏc⁶ / (15360πG²M²) (for massless emission)
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            Luminosity (W)
        """
        hbar = self.constants['hbar']
        c = self.constants['c']
        G = self.constants['G']
        
        return hbar * c**6 / (15360 * np.pi * G**2 * M**2)
    
    def compute_evaporation_time(self, M: float) -> float:
        """
        Compute black hole evaporation time.
        
        t_evap = 5120πG²M³ / (ℏc⁴)
        
        Args:
            M: Initial mass (kg)
            
        Returns:
            Evaporation time (s)
        """
        hbar = self.constants['hbar']
        c = self.constants['c']
        G = self.constants['G']
        
        return 5120 * np.pi * G**2 * M**3 / (hbar * c**4)
    
    def compute_specific_heat(self, M: float) -> float:
        """
        Compute black hole specific heat.
        
        C = dE/dT = -8π(GM)²k_B / (ℏc)
        
        Negative specific heat → thermodynamically unstable!
        
        Args:
            M: Black hole mass (kg)
            
        Returns:
            Specific heat (J/K)
        """
        G = self.constants['G']
        k_B = self.constants['k_B']
        hbar = self.constants['hbar']
        c = self.constants['c']
        
        return -8 * np.pi * (G * M)**2 * k_B / (hbar * c)


# ═══════════════════════════════════════════════════════════════════════════════
# 12. LOOP QUANTUM COSMOLOGY (LQC) CALCULATORS
# ═══════════════════════════════════════════════════════════════════════════════

class LoopQuantumCosmologyCalculator:
    """
    Loop Quantum Cosmology bounce calculator.
    
    LQC replaces Big Bang singularity with quantum bounce.
    Critical density: ρ_c ~ 0.41 ρ_P (Planck density)
    
    Modified Friedmann: H² = (8πG/3)ρ(1 - ρ/ρ_c)
    """
    
    def __init__(self):
        self.constants = CONSTANTS.copy()
        
        # Planck density
        G = self.constants['G']
        c = self.constants['c']
        hbar = self.constants['hbar']
        
        self.rho_P = c**5 / (hbar * G**2)  # Planck density
        self.rho_c = 0.41 * self.rho_P     # LQC critical density
        
    def compute_modified_friedmann_H(self, rho: float, k: float = 0) -> float:
        """
        Compute modified Friedmann equation in LQC.
        
        H² = (8πG/3)ρ(1 - ρ/ρ_c) - k/a²
        
        Bounce occurs when ρ = ρ_c (H → 0)
        
        Args:
            rho: Energy density (kg/m³)
            k: Curvature (0 for flat)
            
        Returns:
            Hubble rate squared (1/s²)
        """
        G = self.constants['G']
        
        H_sq = (8 * np.pi * G / 3) * rho * (1 - rho / self.rho_c)
        
        return max(H_sq, 0)  # H² ≥ 0
    
    def compute_bounce_scale_factor(self, rho_0: float, a_0: float = 1) -> float:
        """
        Compute scale factor at bounce.
        
        At bounce: ρ = ρ_c, so using ρ ∝ a^{-3}:
        a_bounce/a_0 = (ρ_0/ρ_c)^{1/3}
        
        Args:
            rho_0: Current matter density (kg/m³)
            a_0: Current scale factor
            
        Returns:
            Scale factor at bounce
        """
        return a_0 * (rho_0 / self.rho_c)**(1/3)
    
    def compute_holonomy_correction(self, a: float, gamma: float = 0.2375) -> float:
        """
        Compute holonomy correction factor.
        
        sin(λβ)/(λβ) where λ ~ √Δ/a
        
        Args:
            a: Scale factor
            gamma: Barbero-Immirzi parameter
            
        Returns:
            Correction factor
        """
        G = self.constants['G']
        hbar = self.constants['hbar']
        c = self.constants['c']
        
        # Area gap
        Delta = 4 * np.sqrt(3) * np.pi * gamma * G * hbar / c**3
        
        lambda_param = np.sqrt(Delta) / a
        
        if lambda_param < 1e-10:
            return 1.0
            
        x = lambda_param * 1.0  # β ~ 1 normalization
        return np.sin(x) / x if x != 0 else 1.0
    
    def compute_pre_bounce_evolution(self, t_before_bounce: float,
                                      rho_bounce: float = None) -> dict:
        """
        Compute universe state before bounce.
        
        In contracting phase, ρ increases until ρ_c.
        
        Args:
            t_before_bounce: Time before bounce (positive, in Planck times)
            rho_bounce: Density at bounce (default ρ_c)
            
        Returns:
            Dictionary with pre-bounce state
        """
        rho_bounce = rho_bounce or self.rho_c
        
        # Simplified pre-bounce: a ∝ |t|^{2/3}
        t_P = np.sqrt(self.constants['hbar'] * self.constants['G'] / 
                      self.constants['c']**5)  # Planck time
        
        t = t_before_bounce * t_P
        
        # Contracting phase
        a_rel = (t_before_bounce / 1)**(-2/3)  # Relative to bounce
        rho = rho_bounce * a_rel**3
        
        H_sq = self.compute_modified_friedmann_H(rho)
        
        return {
            't_before_bounce': t,
            'a_relative': a_rel,
            'rho': rho,
            'H': -np.sqrt(H_sq),  # Negative for contraction
            'contracting': True
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR REGISTRY
# ═══════════════════════════════════════════════════════════════════════════════

# Global calculator instances
ISM_PHASE = ISMPhaseCalculator()
STELLAR_EVOLUTION = StellarEvolutionCalculator()
BBN = BBNCalculator()
COSMIC_VOID = CosmicVoidCalculator()
REIONIZATION = ReionizationCalculator()
AGN_OUTFLOW = AGNOutflowCalculator()
BINARY_PULSAR = BinaryPulsarCalculator()
COSMIC_RAY = CosmicRayCalculator()
IGM = IGMCalculator()
DARK_ENERGY = DarkEnergyCalculator()
BH_THERMODYNAMICS = BlackHoleThermodynamicsCalculator()
LQC = LoopQuantumCosmologyCalculator()

# Registry
COSMIC_DOMAIN_CALCULATORS = {
    'ISMPhaseCalculator': ISM_PHASE,
    'StellarEvolutionCalculator': STELLAR_EVOLUTION,
    'BBNCalculator': BBN,
    'CosmicVoidCalculator': COSMIC_VOID,
    'ReionizationCalculator': REIONIZATION,
    'AGNOutflowCalculator': AGN_OUTFLOW,
    'BinaryPulsarCalculator': BINARY_PULSAR,
    'CosmicRayCalculator': COSMIC_RAY,
    'IGMCalculator': IGM,
    'DarkEnergyCalculator': DARK_ENERGY,
    'BlackHoleThermodynamicsCalculator': BH_THERMODYNAMICS,
    'LoopQuantumCosmologyCalculator': LQC,
}


# ═══════════════════════════════════════════════════════════════════════════════
# TEST HARNESS
# ═══════════════════════════════════════════════════════════════════════════════

def test_cosmic_domain_calculators():
    """Test all cosmic domain calculators."""
    print("=" * 80)
    print("Cosmic Domain Calculators Test Suite")
    print("=" * 80)
    
    M_sun = CONSTANTS['M_sun']
    
    # Test 1: ISM
    print("\n1. ISM PHASE CALCULATOR")
    for phase, props in ISM_PHASE.phases.items():
        P = ISM_PHASE.compute_pressure(props['T'], props['n'])
        print(f"   {phase}: T={props['T']:.0e}K, n={props['n']}cm⁻³, P={P:.3e}Pa")
    
    # Test 2: Stellar Evolution
    print("\n2. STELLAR EVOLUTION CALCULATOR")
    for M_ratio in [0.5, 1.0, 10.0]:
        M = M_ratio * M_sun
        t_MS = STELLAR_EVOLUTION.compute_main_sequence_lifetime(M)
        L = STELLAR_EVOLUTION.compute_luminosity(M)
        print(f"   {M_ratio:.1f}M☉: t_MS={t_MS/(1e9*3.15e7):.2f}Gyr, L={L/CONSTANTS['L_sun']:.1e}L☉")
    
    # Test 3: BBN
    print("\n3. BBN CALCULATOR")
    T_freeze = BBN.compute_freeze_out_temperature()
    n_p = BBN.compute_neutron_to_proton_ratio(T_freeze)
    Y_p = BBN.compute_helium_mass_fraction(n_p * np.exp(-180/879.4))
    print(f"   T_freeze = {T_freeze:.2e}K")
    print(f"   n/p at freeze-out = {n_p:.4f}")
    print(f"   Y_p (predicted) = {Y_p:.3f}")
    
    # Test 4: Binary Pulsar
    print("\n4. BINARY PULSAR CALCULATOR")
    M1 = 1.4 * M_sun
    M2 = 1.4 * M_sun
    a = 1e9  # ~1 million km
    e = 0.6
    L_GW = BINARY_PULSAR.compute_gravitational_wave_luminosity(M1, M2, a, e)
    t_merge = BINARY_PULSAR.compute_merger_time(M1, M2, a, e)
    print(f"   L_GW = {L_GW:.3e}W ({L_GW/CONSTANTS['L_sun']:.3e}L☉)")
    print(f"   t_merge = {t_merge/(1e6*3.15e7):.2f}Myr")
    
    # Test 5: BH Thermodynamics
    print("\n5. BLACK HOLE THERMODYNAMICS CALCULATOR")
    for M_ratio in [1, 1e6, 1e9]:
        M = M_ratio * M_sun
        T_H = BH_THERMODYNAMICS.compute_hawking_temperature(M)
        S = BH_THERMODYNAMICS.compute_bekenstein_hawking_entropy(M)
        t_evap = BH_THERMODYNAMICS.compute_evaporation_time(M)
        label = f"{M_ratio:.0e}" if M_ratio > 1 else "1"
        print(f"   {label}M☉: T_H={T_H:.2e}K, S/k={S/CONSTANTS['k_B']:.2e}, t_evap={t_evap:.2e}s")
    
    # Test 6: LQC
    print("\n6. LOOP QUANTUM COSMOLOGY CALCULATOR")
    print(f"   ρ_Planck = {LQC.rho_P:.3e}kg/m³")
    print(f"   ρ_critical (LQC) = {LQC.rho_c:.3e}kg/m³")
    rho_test = 0.5 * LQC.rho_c
    H_sq = LQC.compute_modified_friedmann_H(rho_test)
    print(f"   H(ρ=0.5ρ_c) = {np.sqrt(H_sq):.3e}/s")
    H_sq_bounce = LQC.compute_modified_friedmann_H(LQC.rho_c)
    print(f"   H(ρ=ρ_c) = {np.sqrt(H_sq_bounce):.3e}/s (BOUNCE!)")
    
    print("\n" + "=" * 80)
    print("All cosmic domain calculator tests completed!")
    print("=" * 80)


if __name__ == "__main__":
    test_cosmic_domain_calculators()
