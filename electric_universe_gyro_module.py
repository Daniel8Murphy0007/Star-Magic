"""
Electric Universe and Gyroscopic Integration Module
=====================================================

Extracted from Grok UQFF conversation (August 2025).
Implements EU theory validation via UQFF and gyroscopic torque nullification.

Key Physics:
- Electric Universe: F_EM / F_g ratio derivations (R ~ 10^71 locally)
- EM dominance via vacuum-modulated Um
- Gyroscopic torque nullification: τ ~ 10^-25 N·m
- Buoyancy as mediator between EU plasma and gravity
- Inertia (Ui) as gyro operator for rotational stability

Theoretical Framework:
- EU claims EM ~10^39 times stronger than gravity
- UQFF proves R ~10^71 at nuclear scales (exceeds EU claim)
- Buoyancy resolves cosmic stability (prevents EM chaos)
- Gyros stabilize rotations via torque nullification

Author: Daniel T. Murphy
Framework: Universal Quantum Field Superconductive Framework (UQFF)
"""

import math
from typing import Dict, Any, Optional, Tuple
from dataclasses import dataclass

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

# Fundamental
G = 6.674e-11        # Gravitational constant (m³/kg/s²)
c = 2.998e8          # Speed of light (m/s)
hbar = 1.055e-34     # Reduced Planck constant (J·s)
e = 1.602e-19        # Elementary charge (C)
m_p = 1.673e-27      # Proton mass (kg)
m_e = 9.109e-31      # Electron mass (kg)
k_B = 1.381e-23      # Boltzmann constant (J/K)
mu_0 = 4 * math.pi * 1e-7  # Permeability of free space (H/m)

# Nuclear/Atomic
fm = 1e-15           # Femtometer (m)
u = 1.661e-27        # Atomic mass unit (kg)
a_0 = 5.29e-11       # Bohr radius (m)

# UQFF Calibrated Constants
F_rel = 4.30e33              # Relativistic coherence force (N) from LEP 1998
E_LEP = 200e9 * e            # LEP baseline energy (J) = 200 GeV
rho_vac_SCm = 7.09e-37       # Superconductive vacuum density (J/m³)
rho_vac_UA = 7.09e-36        # Aether vacuum density (J/m³)
gamma_decay = 5e-5           # Decay constant (day^-1)
omega_c = 1.585e-8           # Critical angular frequency (rad/s)


@dataclass
class ElectricUniverseSystem:
    """Parameters for EU validation calculation."""
    name: str
    mass: float           # kg
    charge: float         # C (can be effective charge)
    radius: float         # m (characteristic scale)
    velocity: float       # m/s (characteristic velocity)
    magnetic_moment: float = 0.0  # A·m² (for Um contribution)


@dataclass
class GyroscopicSystem:
    """Parameters for gyroscopic stabilization calculation."""
    name: str
    moment_of_inertia: float  # kg·m²
    angular_velocity: float   # rad/s
    precession_rate: float    # rad/s²
    radius: float             # m (lever arm for buoyancy)


# =============================================================================
# ELECTRIC UNIVERSE CALCULATORS
# =============================================================================

class ElectricUniverseCalculator:
    """
    Calculate Electric Universe theory validation via UQFF.
    
    EU Theory Claims:
    - EM forces dominate cosmic structure, not gravity
    - F_EM / F_g ~ 10^39 (standard physics comparison)
    - Plasma currents (Birkeland) drive galaxies/stars
    
    UQFF Resolution:
    - R = F_EM / F_g ~ 10^71 at nuclear scales (proves EU locally)
    - Um vacuum amplification enhances EM beyond EU claims
    - Buoyancy (F_U_Bi_i) resolves macro stability
    - Gravity "emerges" at cosmic scales via buoyancy modulation
    """
    
    # Pre-defined test systems
    SYSTEMS = {
        "alpha_particle": ElectricUniverseSystem(
            name="alpha_particle",
            mass=4 * u,            # 4 u
            charge=2 * e,          # +2e
            radius=2 * fm,         # ~2 fm
            velocity=0.1 * c,      # nuclear velocities
        ),
        "hydrogen_atom": ElectricUniverseSystem(
            name="hydrogen_atom",
            mass=m_p,
            charge=e,
            radius=a_0,
            velocity=2.2e6,        # Bohr velocity
        ),
        "magnetar_surface": ElectricUniverseSystem(
            name="magnetar_surface",
            mass=1e30,             # ~0.5 M_sun
            charge=1e20,           # Effective plasma charge
            radius=1e4,            # 10 km
            velocity=0.1 * c,
            magnetic_moment=1e30,  # Extreme magnetic
        ),
        "solar_corona": ElectricUniverseSystem(
            name="solar_corona",
            mass=2e30,             # M_sun
            charge=1e10,           # Coronal plasma
            radius=7e8,            # R_sun
            velocity=1e6,          # Solar wind
        ),
    }
    
    def __init__(self, system: ElectricUniverseSystem):
        self.system = system
    
    @classmethod
    def from_name(cls, name: str) -> 'ElectricUniverseCalculator':
        """Create calculator from predefined system name."""
        if name not in cls.SYSTEMS:
            raise ValueError(f"Unknown system: {name}. Available: {list(cls.SYSTEMS.keys())}")
        return cls(cls.SYSTEMS[name])
    
    def calculate_gravitational_force(self) -> float:
        """
        Calculate gravitational force between two identical masses.
        
        F_g = G M² / r²
        
        Returns:
            Gravitational force (N)
        """
        return G * self.system.mass**2 / self.system.radius**2
    
    def calculate_coulomb_force(self) -> float:
        """
        Calculate Coulomb force between charges.
        
        F_C = k q² / r² (repulsive for same sign)
        
        Returns:
            Coulomb force magnitude (N)
        """
        k_e = 8.99e9  # Coulomb constant (N·m²/C²)
        return k_e * self.system.charge**2 / self.system.radius**2
    
    def calculate_Um_field(self, t: float = 1e-15) -> float:
        """
        Calculate Universal Magnetism field strength.
        
        Um(t,r) = μ_j(t) / r × (1 - exp(-γt cos(πtn))) × E_react × factors
        
        Parameters:
            t: Time scale (s), default nuclear timescale
            
        Returns:
            Um field strength (T·pm³)
        """
        # Magnetic moment oscillation (scaled for time)
        t_days = t / 86400  # Convert to days
        mu_j = (1e3 + 0.4 * math.sin(omega_c * t)) * 3.38e20  # T·pm³
        
        # Exponential term (adjust for short times)
        n = 0
        gamma_t = gamma_decay * t_days
        if gamma_t < 1e-10:
            exp_term = gamma_t  # Taylor expansion for small argument
        else:
            exp_term = 1 - math.exp(-gamma_t * math.cos(math.pi * t_days * n))
        
        # Ensure non-zero for calculations
        exp_term = max(exp_term, 1e-20)
        
        # Energy reactant (scaled)
        E_react = 1e46 * math.exp(-0.0005 * t_days)
        
        # Enhancement factors
        f_Heav = 0.01
        f_quasi = 0.01
        P_SCm = 1.0
        
        Um = (mu_j / self.system.radius) * exp_term * P_SCm * E_react 
        Um *= (1 + 1e13 * f_Heav) * (1 + f_quasi)
        
        return Um
    
    def calculate_Um_electric_field(self, Um: float) -> float:
        """
        Calculate electric field from Um via UQFF.
        
        E = Um × ρ_vac,[UA] / r
        
        Parameters:
            Um: Universal Magnetism field strength
            
        Returns:
            Electric field (V/m)
        """
        return Um * rho_vac_UA / self.system.radius
    
    def calculate_electromagnetic_force(self, Um: Optional[float] = None) -> float:
        """
        Calculate total EM force from UQFF.
        
        F_EM = q × E × v (Lorentz force approximation)
        
        Parameters:
            Um: Precomputed Um field (optional)
            
        Returns:
            EM force (N)
        """
        if Um is None:
            Um = self.calculate_Um_field()
        
        E_field = self.calculate_Um_electric_field(Um)
        
        # Lorentz force
        F_EM = self.system.charge * E_field * self.system.velocity
        
        return abs(F_EM)
    
    def calculate_EU_ratio(self) -> Dict[str, float]:
        """
        Calculate F_EM / F_g ratio (EU validation).
        
        Standard physics: R ~ 10^39 (e²/G m_p²)
        UQFF extension: R ~ 10^71 at nuclear scales
        
        Returns:
            Dictionary with ratios and interpretation
        """
        F_g = self.calculate_gravitational_force()
        F_C = self.calculate_coulomb_force()
        F_EM_UQFF = self.calculate_electromagnetic_force()
        
        # Standard EU ratio
        R_standard = F_C / F_g if F_g > 0 else float('inf')
        
        # UQFF enhanced ratio
        R_UQFF = F_EM_UQFF / F_g if F_g > 0 else float('inf')
        
        return {
            "F_g_N": F_g,
            "F_Coulomb_N": F_C,
            "F_EM_UQFF_N": F_EM_UQFF,
            "R_standard": R_standard,
            "R_UQFF": R_UQFF,
            "EU_claim_10_39": 1e39,
            "exceeds_EU": R_UQFF > 1e39,
        }
    
    def calculate_buoyancy_resolution(self) -> float:
        """
        Calculate buoyancy force that resolves EU at cosmic scales.
        
        F_U_Bi_i = -F_rel × (E_cm / E_LEP) × Q_wave × g(r)
        
        Returns:
            Buoyancy force (N), negative for stabilization
        """
        # Characteristic energy
        E_char = 0.5 * self.system.mass * self.system.velocity**2
        E_ratio = E_char / E_LEP
        
        # Wave factor
        Q_wave = 1e12
        
        # Local gravity
        g_local = G * self.system.mass / self.system.radius**2
        
        # Buoyancy (negative for stability)
        F_UBii = -F_rel * E_ratio * Q_wave * g_local / 1e30
        
        return F_UBii
    
    def compute_full_EU_analysis(self) -> Dict[str, Any]:
        """
        Complete Electric Universe analysis via UQFF.
        
        Returns:
            Dictionary with all EU validation metrics
        """
        ratios = self.calculate_EU_ratio()
        F_buoyancy = self.calculate_buoyancy_resolution()
        
        return {
            "system": self.system.name,
            "mass_kg": self.system.mass,
            "radius_m": self.system.radius,
            **ratios,
            "F_UBii_N": F_buoyancy,
            "EU_local_validity": "PROVEN" if ratios["R_UQFF"] > 1e39 else "PARTIAL",
            "cosmic_resolution": "Buoyancy prevents EM chaos",
            "interpretation": (
                f"EU proven locally (R={ratios['R_UQFF']:.2e} >> 10^39), "
                f"but buoyancy (F={F_buoyancy:.2e} N) stabilizes cosmic scales"
            ),
        }


# =============================================================================
# GYROSCOPIC CALCULATORS
# =============================================================================

class GyroscopicCalculator:
    """
    Calculate gyroscopic torque and UQFF nullification.
    
    Gyroscopic Effects:
    - Precession: τ = I × ω × α
    - Stabilization: Gyros resist torque via angular momentum
    - Nuclear: Spin precession (~10^21 rad/s)
    - Cosmic: Planet precession (26,000-year cycles)
    
    UQFF Integration:
    - Ui (Universal Inertia) as gyro operator
    - Torque nullification via buoyancy: τ + F_U_Bi_i × r = 0
    - Explains stable rotations in EU plasma (Birkeland currents)
    """
    
    # Pre-defined gyroscopic systems
    SYSTEMS = {
        "alpha_cluster": GyroscopicSystem(
            name="alpha_cluster",
            moment_of_inertia=1.07e-56,  # kg·m² (sphere, ~4u, ~2fm)
            angular_velocity=1e21,        # rad/s (nuclear rotation)
            precession_rate=1e10,         # rad/s² (speculative)
            radius=2e-15,                 # m
        ),
        "neutron_star": GyroscopicSystem(
            name="neutron_star",
            moment_of_inertia=1e45,       # kg·m² (~M_sun, 10 km)
            angular_velocity=1e3,          # rad/s (millisecond pulsar)
            precession_rate=1e-8,          # rad/s² (spin-down)
            radius=1e4,                    # m
        ),
        "planet_nine": GyroscopicSystem(
            name="planet_nine",
            moment_of_inertia=1e42,       # kg·m² (Neptune-class)
            angular_velocity=1e-4,         # rad/s (orbital)
            precession_rate=1e-12,         # rad/s² (epoch cycle)
            radius=9e13,                   # m (600 AU)
        ),
        "earth_precession": GyroscopicSystem(
            name="earth_precession",
            moment_of_inertia=8e37,       # kg·m² (Earth)
            angular_velocity=7.3e-5,       # rad/s (daily rotation)
            precession_rate=7.7e-12,       # rad/s² (26,000-year cycle)
            radius=6.4e6,                  # m
        ),
    }
    
    def __init__(self, system: GyroscopicSystem):
        self.system = system
    
    @classmethod
    def from_name(cls, name: str) -> 'GyroscopicCalculator':
        """Create calculator from predefined system name."""
        if name not in cls.SYSTEMS:
            raise ValueError(f"Unknown system: {name}. Available: {list(cls.SYSTEMS.keys())}")
        return cls(cls.SYSTEMS[name])
    
    def calculate_torque(self) -> float:
        """
        Calculate gyroscopic torque.
        
        τ = I × ω × α
        
        Returns:
            Torque (N·m)
        """
        return self.system.moment_of_inertia * self.system.angular_velocity * self.system.precession_rate
    
    def calculate_angular_momentum(self) -> float:
        """
        Calculate angular momentum.
        
        L = I × ω
        
        Returns:
            Angular momentum (kg·m²/s)
        """
        return self.system.moment_of_inertia * self.system.angular_velocity
    
    def calculate_precession_period(self) -> float:
        """
        Calculate precession period.
        
        T_prec = 2π / (α / ω)  [simplified]
        
        Returns:
            Precession period (s)
        """
        if self.system.precession_rate == 0:
            return float('inf')
        return 2 * math.pi * self.system.angular_velocity / self.system.precession_rate
    
    def calculate_buoyancy_nullification(self) -> Dict[str, float]:
        """
        Calculate buoyancy force required to nullify torque.
        
        τ + F_U_Bi_i × r × sin(θ) = 0
        
        For θ = π/2 (maximum coupling):
        F_required = -τ / r
        
        Returns:
            Dictionary with nullification parameters
        """
        tau = self.calculate_torque()
        
        # Required buoyancy force for nullification
        F_required = -tau / self.system.radius
        
        return {
            "tau_N_m": tau,
            "F_required_N": F_required,
            "r_m": self.system.radius,
            "nullified": abs(F_required) < 1e10,  # Reasonable for UQFF
        }
    
    def calculate_Ui_operator(self) -> float:
        """
        Calculate Universal Inertia operator.
        
        Ui = Δρ_vac / ρ_LEP × Q_res
        
        Returns:
            Ui dimensionless operator
        """
        # Vacuum density differential
        delta_rho = rho_vac_UA - rho_vac_SCm
        
        # Resonance factor
        Q_res = 1e12
        
        # Normalize to LEP scale
        rho_LEP = E_LEP / c**2  # effective density
        
        Ui = delta_rho / rho_LEP * Q_res
        
        return Ui
    
    def compute_full_gyro_analysis(self) -> Dict[str, Any]:
        """
        Complete gyroscopic analysis with UQFF integration.
        
        Returns:
            Dictionary with all gyro metrics
        """
        tau = self.calculate_torque()
        L = self.calculate_angular_momentum()
        T_prec = self.calculate_precession_period()
        null_data = self.calculate_buoyancy_nullification()
        Ui = self.calculate_Ui_operator()
        
        return {
            "system": self.system.name,
            "I_kg_m2": self.system.moment_of_inertia,
            "omega_rad_s": self.system.angular_velocity,
            "alpha_rad_s2": self.system.precession_rate,
            "tau_N_m": tau,
            "L_kg_m2_s": L,
            "T_prec_s": T_prec,
            "T_prec_yr": T_prec / (365.25 * 86400) if T_prec < 1e20 else float('inf'),
            "F_nullify_N": null_data["F_required_N"],
            "Ui_operator": Ui,
            "stability": "STABLE" if null_data["nullified"] else "REQUIRES_DAMPING",
            "interpretation": f"Torque τ={tau:.2e} N·m nullified by buoyancy F={null_data['F_required_N']:.2e} N",
        }


# =============================================================================
# EU-GYRO COMBINED CALCULATOR
# =============================================================================

class EUGyroUnifiedCalculator:
    """
    Combined Electric Universe + Gyroscopic calculation.
    
    Unifies:
    - EU plasma dominance (F_EM >> F_g locally)
    - Gyro stabilization (torque nullification)
    - Buoyancy resolution (cosmic coherence)
    
    Applications:
    - Birkeland currents as gyro-stabilized plasma
    - Magnetar fields with precession
    - Planetary systems with orbital gyros
    """
    
    def __init__(self, eu_system: ElectricUniverseSystem, gyro_system: GyroscopicSystem):
        self.eu_calc = ElectricUniverseCalculator(eu_system)
        self.gyro_calc = GyroscopicCalculator(gyro_system)
    
    @classmethod
    def magnetar_analysis(cls) -> 'EUGyroUnifiedCalculator':
        """Create unified calculator for magnetar system."""
        eu = ElectricUniverseCalculator.SYSTEMS["magnetar_surface"]
        gyro = GyroscopicCalculator.SYSTEMS["neutron_star"]
        return cls(eu, gyro)
    
    @classmethod
    def planet_nine_analysis(cls) -> 'EUGyroUnifiedCalculator':
        """Create unified calculator for Planet Nine orbital dynamics."""
        eu = ElectricUniverseCalculator.SYSTEMS["solar_corona"]
        gyro = GyroscopicCalculator.SYSTEMS["planet_nine"]
        return cls(eu, gyro)
    
    def compute_unified_analysis(self) -> Dict[str, Any]:
        """
        Complete EU + Gyro unified analysis.
        
        Returns:
            Dictionary with combined metrics
        """
        eu_results = self.eu_calc.compute_full_EU_analysis()
        gyro_results = self.gyro_calc.compute_full_gyro_analysis()
        
        # Combined interpretation
        eu_local = eu_results["EU_local_validity"] == "PROVEN"
        gyro_stable = gyro_results["stability"] == "STABLE"
        
        return {
            "EU_analysis": eu_results,
            "gyro_analysis": gyro_results,
            "combined_validity": "FULL" if eu_local and gyro_stable else "PARTIAL",
            "unified_interpretation": (
                f"EU proven locally (R={eu_results['R_UQFF']:.2e}), "
                f"gyro {'stabilized' if gyro_stable else 'requires damping'} "
                f"(τ={gyro_results['tau_N_m']:.2e} N·m), "
                f"buoyancy resolves cosmic coherence"
            ),
        }


# =============================================================================
# REGISTRY
# =============================================================================

EU_GYRO_CALCULATORS = {
    "ElectricUniverseCalculator": ElectricUniverseCalculator,
    "GyroscopicCalculator": GyroscopicCalculator,
    "EUGyroUnifiedCalculator": EUGyroUnifiedCalculator,
}

EU_SYSTEMS = ElectricUniverseCalculator.SYSTEMS
GYRO_SYSTEMS = GyroscopicCalculator.SYSTEMS


if __name__ == "__main__":
    # Demo: Electric Universe
    print("=" * 70)
    print("ELECTRIC UNIVERSE ANALYSIS: Alpha Particle")
    print("=" * 70)
    
    eu_calc = ElectricUniverseCalculator.from_name("alpha_particle")
    eu_results = eu_calc.compute_full_EU_analysis()
    
    for key, val in eu_results.items():
        if not isinstance(val, dict):
            print(f"  {key}: {val}")
    
    # Demo: Gyroscopic
    print("\n" + "=" * 70)
    print("GYROSCOPIC ANALYSIS: Alpha Cluster")
    print("=" * 70)
    
    gyro_calc = GyroscopicCalculator.from_name("alpha_cluster")
    gyro_results = gyro_calc.compute_full_gyro_analysis()
    
    for key, val in gyro_results.items():
        print(f"  {key}: {val}")
    
    # Demo: Unified
    print("\n" + "=" * 70)
    print("UNIFIED EU+GYRO ANALYSIS: Magnetar")
    print("=" * 70)
    
    unified = EUGyroUnifiedCalculator.magnetar_analysis()
    unified_results = unified.compute_unified_analysis()
    
    print(f"  combined_validity: {unified_results['combined_validity']}")
    print(f"  interpretation: {unified_results['unified_interpretation']}")
