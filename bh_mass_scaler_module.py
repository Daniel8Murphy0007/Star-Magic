"""
Black Hole Mass Scaler and UQFF Energy Distribution Module
==========================================================

Extracted from Grok UQFF conversation (August 2025).
Implements Harvard BH mass distribution integration with UQFF energy scaling.

Key Physics:
- BH mass function φ(log M) from Harvard/EPS models
- Energy scaling: E_BH = M c² (relativistic)
- UQFF scaler: S = E_BH / E_LEP × Q_res
- Buoyancy modulation for merger hierarchies
- Information paradox and Hawking temperature integration

Source: Harvard BH mass distributions, GWTC-3/4 merger data
Theoretical: Extended Press-Schechter (EPS), Bekenstein-Hawking

Author: Daniel T. Murphy
Framework: Universal Quantum Field Superconductive Framework (UQFF)
"""

import math
from typing import Dict, Any, List, Tuple, Optional
from dataclasses import dataclass

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

# Fundamental
G = 6.674e-11        # Gravitational constant (m³/kg/s²)
c = 2.998e8          # Speed of light (m/s)
hbar = 1.055e-34     # Reduced Planck constant (J·s)
k_B = 1.381e-23      # Boltzmann constant (J/K)
e = 1.602e-19        # Elementary charge (C)

# Astronomical
M_sun = 1.989e30     # Solar mass (kg)
pc = 3.086e16        # Parsec (m)
Mpc = pc * 1e6       # Megaparsec (m)
Gyr = 3.156e16       # Gigayear (s)

# UQFF Constants
F_rel = 4.30e33              # Relativistic coherence force (N)
E_LEP = 200e9 * e            # LEP baseline energy (J) = 200 GeV
rho_vac_UA = 7.09e-36        # Aether vacuum density (J/m³)
kappa_UQFF = 0.0005          # κ calibration (day^-1)
SSq_factor = 0.57            # [SSq] factor


@dataclass
class BHMassDistribution:
    """Black hole mass distribution data structure."""
    bins: List[float]         # log10(M/M_sun) bin edges
    values: List[float]       # φ(log M) values (probability density)
    source: str = "Harvard"   # Data source
    redshift: float = 0.0     # Mean redshift of sample


@dataclass 
class BHSystem:
    """Individual black hole system parameters."""
    name: str
    mass_solar: float         # M/M_sun
    spin: float               # Dimensionless spin a/a_max (0 to 1)
    redshift: float           # Cosmological redshift
    accretion_rate: float = 0.0  # Eddington ratio


# =============================================================================
# HARVARD BH MASS DISTRIBUTION DATA
# =============================================================================

# Extracted from Grok conversation: Harvard energy distribution JSON
# Peaks at log M ≈ 5.2-5.35 (~1.6-2.2×10^5 M_⊙)
HARVARD_DISTRIBUTION = BHMassDistribution(
    bins=[
        3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
        4.0, 4.1, 4.2, 4.25, 4.5, 5.0, 5.1, 5.2, 5.25, 5.3, 5.35,
        5.5, 5.6, 6.0, 6.5, 7.0, 7.5, 8.0
    ],
    values=[
        1e-6, 2e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 0.01,
        0.036, 0.035, 0.034, 0.033, 0.025, 0.030, 0.035, 0.040, 0.042, 0.040, 0.038,
        0.030, 0.025, 0.015, 0.010, 0.005, 0.002, 0.001
    ],
    source="Harvard/EPS",
    redshift=0.5  # Mean z of distribution
)

# GWTC merger peaks from LIGO O3/O4
GWTC_MERGER_PEAKS = {
    "primary_peak": 8.0,      # log10(M) for ~10^8 M_sun
    "chirp_mass_peak": 1.5,   # log10(M_chirp/M_sun) ~ 30 M_sun
    "gap_lower": 2.3,         # Lower mass gap: 2-5 M_sun
    "gap_upper": 4.7,         # Upper mass gap: 50-130 M_sun (PISN)
}


# =============================================================================
# BH MASS → ENERGY SCALER
# =============================================================================

class BHMassEnergyScaler:
    """
    Scale BH mass to UQFF energy domain.
    
    Key Relations:
    - E_BH = M c² (total rest mass energy)
    - In natural units: 1 M_⊙ c² ≈ 1.12×10^54 J ≈ 7×10^65 GeV
    - UQFF scaler S = E_BH / E_LEP × Q_res
    - Buoyancy from scaler: F_U_Bi_i = F_rel × S × g_local
    """
    
    def __init__(self, distribution: BHMassDistribution = HARVARD_DISTRIBUTION):
        self.distribution = distribution
    
    @staticmethod
    def mass_to_energy_joules(M_solar: float) -> float:
        """
        Convert solar masses to Joules.
        
        E = M c²
        
        Parameters:
            M_solar: Mass in solar masses
            
        Returns:
            Energy in Joules
        """
        return M_solar * M_sun * c**2
    
    @staticmethod
    def mass_to_energy_GeV(M_solar: float) -> float:
        """
        Convert solar masses to GeV.
        
        1 J = 6.242×10^9 GeV
        
        Parameters:
            M_solar: Mass in solar masses
            
        Returns:
            Energy in GeV
        """
        E_J = M_solar * M_sun * c**2
        return E_J / (e * 1e9)  # J to GeV
    
    @staticmethod
    def log_mass_to_energy_GeV(log_M: float) -> float:
        """
        Convert log10(M/M_sun) to log10(E/GeV).
        
        log10(E) = log10(M) + log10(M_sun c² / 1 GeV)
        ≈ log10(M) + 65.85
        
        Parameters:
            log_M: log10(M/M_sun)
            
        Returns:
            log10(E/GeV)
        """
        # M_sun c² in GeV
        M_sun_c2_GeV = M_sun * c**2 / (e * 1e9)
        return log_M + math.log10(M_sun_c2_GeV)
    
    def calculate_UQFF_scaler(self, M_solar: float, Q_res: float = 1e12) -> float:
        """
        Calculate UQFF scaler for BH mass.
        
        S = E_BH / E_LEP × Q_res
        
        Parameters:
            M_solar: Mass in solar masses
            Q_res: Resonance factor (THz from Colman-Gillespie)
            
        Returns:
            UQFF scaler (dimensionless)
        """
        E_BH = self.mass_to_energy_joules(M_solar)
        return (E_BH / E_LEP) * Q_res
    
    def calculate_scaled_Um(self, M_solar: float, Um_lab: float = 1e3) -> float:
        """
        Scale laboratory Um to BH-scale.
        
        Um_BH = Um_lab × S
        
        Parameters:
            M_solar: Mass in solar masses
            Um_lab: Laboratory Um (T·pm³ from Colman-Gillespie)
            
        Returns:
            Scaled Um (T·pm³)
        """
        S = self.calculate_UQFF_scaler(M_solar)
        return Um_lab * S
    
    def calculate_scaled_neutron_rate(self, M_solar: float, eta_lab: float = 1e13) -> float:
        """
        Scale laboratory neutron rate to BH environment.
        
        η_BH = η_lab × S × surface_factor
        
        Parameters:
            M_solar: Mass in solar masses
            eta_lab: Laboratory rate (cm^-2/s)
            
        Returns:
            Scaled rate (very large for SMBHs)
        """
        S = self.calculate_UQFF_scaler(M_solar)
        # Surface area enhancement (relative to lab scale)
        r_BH = 2 * G * M_solar * M_sun / c**2  # Schwarzschild radius
        surface_factor = 4 * math.pi * r_BH**2 / 1e-20  # vs ~fm² lab
        return eta_lab * S * surface_factor
    
    def calculate_buoyancy_force(self, M_solar: float) -> float:
        """
        Calculate UQFF buoyancy force for BH.
        
        F_U_Bi_i = -F_rel × S × g(r)
        
        Negative for coherence stabilization (prevents overgrowth).
        
        Parameters:
            M_solar: Mass in solar masses
            
        Returns:
            Buoyancy force (N)
        """
        S = self.calculate_UQFF_scaler(M_solar)
        
        # Gravity at horizon
        r_s = 2 * G * M_solar * M_sun / c**2
        g_horizon = c**4 / (4 * G * M_solar * M_sun)  # Surface gravity
        
        # Scaled and normalized
        F_UBii = -F_rel * S * g_horizon / 1e50
        
        return F_UBii
    
    def interpret_distribution_peaks(self) -> Dict[str, Any]:
        """
        Interpret mass distribution peaks via UQFF.
        
        Returns:
            Dictionary with peak interpretations
        """
        # Find peaks in distribution
        peaks = []
        for i in range(1, len(self.distribution.values) - 1):
            if (self.distribution.values[i] > self.distribution.values[i-1] and
                self.distribution.values[i] > self.distribution.values[i+1]):
                peaks.append({
                    "log_M": self.distribution.bins[i],
                    "phi": self.distribution.values[i],
                    "M_solar": 10**self.distribution.bins[i],
                })
        
        # Add UQFF interpretation
        for peak in peaks:
            peak["E_scaler"] = self.calculate_UQFF_scaler(peak["M_solar"])
            peak["F_UBii_N"] = self.calculate_buoyancy_force(peak["M_solar"])
            peak["interpretation"] = "coherence resonance" if peak["phi"] > 0.03 else "buoyancy threshold"
        
        return {
            "n_peaks": len(peaks),
            "peaks": peaks,
            "integral": sum(self.distribution.values) * 0.1,  # ~dx
            "source": self.distribution.source,
        }


# =============================================================================
# BH THERMODYNAMICS CALCULATOR
# =============================================================================

class BHThermodynamicsCalculator:
    """
    Calculate BH thermodynamics with UQFF extensions.
    
    Standard Relations:
    - Hawking temperature: T_H = ℏc³ / (8πGMk_B)
    - Bekenstein-Hawking entropy: S = k_B c³ A / (4Gℏ)
    - Evaporation lifetime: τ = 5120 π G² M³ / (ℏc⁴)
    
    UQFF Extensions:
    - Information channels via 26D framework
    - Page curve via buoyancy modulation
    - Coherence terms for information preservation
    """
    
    def __init__(self, bh: BHSystem):
        self.bh = bh
    
    @classmethod
    def from_mass(cls, M_solar: float, name: str = "BH") -> 'BHThermodynamicsCalculator':
        """Create calculator from solar mass."""
        bh = BHSystem(
            name=name,
            mass_solar=M_solar,
            spin=0.0,
            redshift=0.0,
        )
        return cls(bh)
    
    def calculate_hawking_temperature(self) -> float:
        """
        Calculate Hawking temperature.
        
        T_H = ℏc³ / (8πGMk_B)
        
        Returns:
            Temperature (K)
        """
        M = self.bh.mass_solar * M_sun
        return hbar * c**3 / (8 * math.pi * G * M * k_B)
    
    def calculate_bekenstein_hawking_entropy(self) -> float:
        """
        Calculate Bekenstein-Hawking entropy.
        
        S = k_B c³ A / (4Gℏ) = 4π k_B G M² / (ℏc)
        
        Returns:
            Entropy (J/K)
        """
        M = self.bh.mass_solar * M_sun
        return 4 * math.pi * k_B * G * M**2 / (hbar * c)
    
    def calculate_entropy_bits(self) -> float:
        """
        Calculate entropy in bits.
        
        S_bits = S / (k_B ln(2))
        
        Returns:
            Entropy (bits)
        """
        S = self.calculate_bekenstein_hawking_entropy()
        return S / (k_B * math.log(2))
    
    def calculate_evaporation_lifetime(self) -> float:
        """
        Calculate evaporation lifetime.
        
        τ = 5120 π G² M³ / (ℏc⁴)
        
        Returns:
            Lifetime (s)
        """
        M = self.bh.mass_solar * M_sun
        return 5120 * math.pi * G**2 * M**3 / (hbar * c**4)
    
    def calculate_schwarzschild_radius(self) -> float:
        """
        Calculate Schwarzschild radius.
        
        r_s = 2GM/c²
        
        Returns:
            Radius (m)
        """
        M = self.bh.mass_solar * M_sun
        return 2 * G * M / c**2
    
    def calculate_UQFF_information_channels(self) -> int:
        """
        Calculate number of UQFF information channels.
        
        From 26D framework: n_channels = 26 × floor(log10(M/M_sun))
        
        Returns:
            Number of channels
        """
        log_M = math.log10(self.bh.mass_solar)
        return 26 * max(1, int(log_M))
    
    def calculate_page_curve_midpoint(self) -> float:
        """
        Calculate Page curve midpoint (half entropy emitted).
        
        t_Page ≈ τ_evap / 2 (approximately)
        
        Returns:
            Page time (s)
        """
        return self.calculate_evaporation_lifetime() / 2
    
    def calculate_UQFF_coherence_factor(self) -> float:
        """
        Calculate UQFF coherence factor for information preservation.
        
        C = exp(-[SSq]^n / 26) × exp(-κ t)
        
        At t=0: C = exp(-0.57/26) ≈ 0.978 (high coherence)
        
        Returns:
            Coherence factor (0 to 1)
        """
        n = 1  # First order
        return math.exp(-SSq_factor**n / 26)
    
    def compute_full_analysis(self) -> Dict[str, Any]:
        """
        Complete BH thermodynamics analysis.
        
        Returns:
            Dictionary with all BH thermodynamic quantities
        """
        T_H = self.calculate_hawking_temperature()
        S_BH = self.calculate_bekenstein_hawking_entropy()
        tau = self.calculate_evaporation_lifetime()
        r_s = self.calculate_schwarzschild_radius()
        
        return {
            "name": self.bh.name,
            "M_solar": self.bh.mass_solar,
            "spin": self.bh.spin,
            "r_s_m": r_s,
            "r_s_km": r_s / 1000,
            "T_Hawking_K": T_H,
            "S_BH_J_K": S_BH,
            "S_bits": self.calculate_entropy_bits(),
            "tau_evap_s": tau,
            "tau_evap_yr": tau / (365.25 * 86400),
            "n_info_channels": self.calculate_UQFF_information_channels(),
            "t_Page_s": self.calculate_page_curve_midpoint(),
            "coherence_factor": self.calculate_UQFF_coherence_factor(),
            "information_preserved": self.calculate_UQFF_coherence_factor() > 0.9,
        }


# =============================================================================
# MERGER HIERARCHY CALCULATOR
# =============================================================================

class MergerHierarchyCalculator:
    """
    Calculate BH merger hierarchies with UQFF buoyancy.
    
    Explains:
    - Mass gaps (PISN, pair instability)
    - Hierarchical pile-ups at certain masses
    - UQFF buoyancy preventing overgrowth
    """
    
    @staticmethod
    def calculate_ISCO_radius(M_solar: float, spin: float = 0.0) -> float:
        """
        Calculate innermost stable circular orbit.
        
        r_ISCO = 6 GM/c² (Schwarzschild)
        Reduced for spinning BHs.
        
        Parameters:
            M_solar: Mass in solar masses
            spin: Dimensionless spin (0 to 1)
            
        Returns:
            ISCO radius (m)
        """
        r_s = 2 * G * M_solar * M_sun / c**2
        # Kerr correction (prograde)
        Z_1 = 1 + (1 - spin**2)**(1/3) * ((1 + spin)**(1/3) + (1 - spin)**(1/3))
        Z_2 = (3 * spin**2 + Z_1**2)**0.5
        r_ISCO_factor = 3 + Z_2 - ((3 - Z_1) * (3 + Z_1 + 2*Z_2))**0.5
        return r_ISCO_factor * r_s / 2
    
    @staticmethod
    def calculate_merger_timescale(M1_solar: float, M2_solar: float, 
                                   a_AU: float, e: float = 0.0) -> float:
        """
        Calculate gravitational wave merger timescale.
        
        Peters formula for circular orbits.
        
        Parameters:
            M1_solar, M2_solar: Component masses (solar masses)
            a_AU: Semi-major axis (AU)
            e: Eccentricity
            
        Returns:
            Merger time (s)
        """
        a = a_AU * 1.496e11  # AU to m
        M = (M1_solar + M2_solar) * M_sun
        mu = M1_solar * M2_solar / (M1_solar + M2_solar) * M_sun
        
        # Circular approximation
        tau = 5 * c**5 * a**4 / (256 * G**3 * M**2 * mu)
        
        # Eccentricity enhancement
        if e > 0:
            tau *= (1 - e**2)**3.5
        
        return tau
    
    @staticmethod
    def calculate_final_mass(M1_solar: float, M2_solar: float, 
                            spin1: float = 0.0, spin2: float = 0.0) -> Tuple[float, float]:
        """
        Calculate final mass and radiated energy from merger.
        
        Approximate fits to numerical relativity.
        
        Returns:
            (M_final_solar, energy_radiated_solar)
        """
        M_total = M1_solar + M2_solar
        q = min(M1_solar, M2_solar) / max(M1_solar, M2_solar)
        eta = q / (1 + q)**2  # Symmetric mass ratio
        
        # Radiated fraction (fit to NR)
        f_rad = 0.0547 * eta * (1 + 0.5 * (spin1 + spin2))
        
        E_rad = f_rad * M_total
        M_final = M_total - E_rad
        
        return M_final, E_rad
    
    @staticmethod
    def calculate_buoyancy_threshold(M_solar: float) -> Dict[str, float]:
        """
        Calculate UQFF buoyancy threshold for merger.
        
        Negative buoyancy prevents certain mergers (mass gaps).
        
        Returns:
            Dictionary with threshold analysis
        """
        scaler = BHMassEnergyScaler()
        F_UBii = scaler.calculate_buoyancy_force(M_solar)
        
        # Threshold is when buoyancy balances gravitational binding
        E_bind = 0.1 * M_solar * M_sun * c**2  # ~10% rest mass binding
        F_bind = E_bind / (2 * G * M_solar * M_sun / c**2)  # Force scale
        
        return {
            "M_solar": M_solar,
            "F_UBii_N": F_UBii,
            "F_binding_N": F_bind,
            "ratio": abs(F_UBii) / F_bind if F_bind > 0 else float('inf'),
            "merger_allowed": abs(F_UBii) < F_bind,
            "interpretation": "buoyancy blocks merger" if abs(F_UBii) >= F_bind else "merger proceeds",
        }
    
    @classmethod
    def analyze_mass_gap(cls, gap_name: str = "PISN") -> Dict[str, Any]:
        """
        Analyze mass gap via UQFF buoyancy.
        
        Parameters:
            gap_name: "PISN" (50-130 M_sun) or "lower" (2-5 M_sun)
            
        Returns:
            Gap analysis dictionary
        """
        if gap_name == "PISN":
            M_range = [50, 130]  # PISN gap
        elif gap_name == "lower":
            M_range = [2, 5]     # Lower mass gap
        else:
            M_range = [10, 100]  # Generic
        
        results = []
        for M in range(int(M_range[0]), int(M_range[1]) + 1, 10):
            results.append(cls.calculate_buoyancy_threshold(float(M)))
        
        return {
            "gap_name": gap_name,
            "M_range_solar": M_range,
            "n_blocked": sum(1 for r in results if not r["merger_allowed"]),
            "n_allowed": sum(1 for r in results if r["merger_allowed"]),
            "samples": results,
            "interpretation": f"UQFF buoyancy {'explains' if any(not r['merger_allowed'] for r in results) else 'allows'} {gap_name} gap",
        }


# =============================================================================
# REGISTRY
# =============================================================================

BH_MASS_CALCULATORS = {
    "BHMassEnergyScaler": BHMassEnergyScaler,
    "BHThermodynamicsCalculator": BHThermodynamicsCalculator,
    "MergerHierarchyCalculator": MergerHierarchyCalculator,
}

# Pre-defined BH systems
BH_SYSTEMS = {
    "stellar_10": BHSystem("stellar_10M", 10, 0.7, 0.1),
    "IMBH_1000": BHSystem("IMBH_1000M", 1000, 0.5, 0.5),
    "SgrA_star": BHSystem("SgrA*", 4.3e6, 0.5, 0.0),
    "M87_star": BHSystem("M87*", 6.5e9, 0.9, 0.0044),
    "TON_618": BHSystem("TON_618", 6.6e10, 0.8, 1.0),  # Most massive known
}


if __name__ == "__main__":
    # Demo: Mass to energy scaling
    print("=" * 70)
    print("BH MASS → ENERGY SCALING")
    print("=" * 70)
    
    scaler = BHMassEnergyScaler()
    for M in [10, 1e3, 1e6, 1e9]:
        print(f"  M = {M:.0e} M_sun:")
        print(f"    E = {scaler.mass_to_energy_GeV(M):.2e} GeV")
        print(f"    S = {scaler.calculate_UQFF_scaler(M):.2e}")
        print(f"    F_UBii = {scaler.calculate_buoyancy_force(M):.2e} N")
    
    # Demo: Distribution peaks
    print("\n" + "=" * 70)
    print("HARVARD DISTRIBUTION PEAK ANALYSIS")
    print("=" * 70)
    
    peaks = scaler.interpret_distribution_peaks()
    print(f"  Number of peaks: {peaks['n_peaks']}")
    for peak in peaks['peaks'][:3]:
        print(f"  Peak at log M = {peak['log_M']}: φ = {peak['phi']:.3f}, {peak['interpretation']}")
    
    # Demo: BH thermodynamics
    print("\n" + "=" * 70)
    print("BH THERMODYNAMICS: Sgr A*")
    print("=" * 70)
    
    sgra = BHThermodynamicsCalculator(BH_SYSTEMS["SgrA_star"])
    results = sgra.compute_full_analysis()
    
    for key in ["M_solar", "r_s_km", "T_Hawking_K", "tau_evap_yr", "n_info_channels", "coherence_factor"]:
        print(f"  {key}: {results[key]:.3e}" if isinstance(results[key], float) else f"  {key}: {results[key]}")
    
    # Demo: Merger hierarchy
    print("\n" + "=" * 70)
    print("MERGER HIERARCHY: PISN Gap Analysis")
    print("=" * 70)
    
    gap_analysis = MergerHierarchyCalculator.analyze_mass_gap("PISN")
    print(f"  Gap: {gap_analysis['gap_name']}, range: {gap_analysis['M_range_solar']} M_sun")
    print(f"  Blocked: {gap_analysis['n_blocked']}, Allowed: {gap_analysis['n_allowed']}")
    print(f"  {gap_analysis['interpretation']}")
