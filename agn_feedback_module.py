#!/usr/bin/env python3
"""
AGN Feedback Module - Active Galactic Nuclei Outflows and Jet Physics

From Grok Deep Analysis (Feb 2026):
- Equations 23-25: AGN Feedback (terminal momentum, jet power, duty cycle)
- Equations 70-72: Quasar Feedback (wind velocity, ionization parameter, coupling)
- Equations 14-15: Quasar Jets (BZ power, relativistic velocity profile)

Physics domains:
- Blandford-Znajek Mechanism for jet power extraction from BH spin
- Radiation-driven winds and momentum coupling
- AGN duty cycles and cooling flow regulation
- Ionization cones and feedback efficiency

UQFF Integration:
- Um (Universal Magnetism) governs jet launching via magnetic torque
- F_U_Bi_i buoyancy modulates feedback-gravity balance
- Electric Universe validation via EM dominance in jet dynamics
"""

import math
from typing import Dict, Optional, Tuple

# ============== Physical Constants ==============
G = 6.674e-11           # Gravitational constant [m³/(kg·s²)]
c = 2.998e8             # Speed of light [m/s]
M_sun = 1.989e30        # Solar mass [kg]
SIGMA_T = 6.652e-29     # Thomson cross-section [m²]
m_p = 1.673e-27         # Proton mass [kg]
k_B = 1.381e-23         # Boltzmann constant [J/K]
h_bar = 1.055e-34       # Reduced Planck constant [J·s]
eV_to_J = 1.602e-19     # eV to Joules
erg_to_J = 1e-7         # erg to Joules
pc_to_m = 3.086e16      # parsec to meters
year_to_s = 3.154e7     # year to seconds

# UQFF Constants
F_rel = 4.30e33         # Relativistic coherence force [N]
E_LEP = 200e9 * eV_to_J # LEP 1998 baseline energy [J]
rho_vac_SCm = 7.09e-37  # Vacuum density SCm [J/m³]
rho_vac_UA = 7.09e-36   # Vacuum density UA [J/m³]


class BlandfordZnajekCalculator:
    """
    Blandford-Znajek jet power from black hole spin extraction.
    
    Equation 14 (BZ Power):
    P_BZ = (1/32) B² R_H^4 Ω_H² / c
    
    Where:
    - R_H = GM/c² (horizon radius)
    - Ω_H = a c / (2 R_H) (spin angular velocity)
    - B: poloidal magnetic field (~10⁴ G for M87*)
    - a: dimensionless spin parameter (0-1)
    
    Derivation: Force-free MHD in Kerr spacetime, Poynting flux
    from twisted magnetic fields extracting rotational energy.
    """
    
    def compute(self, M_BH: float, a_spin: float, B_field: float) -> Dict:
        """
        Compute BZ jet power.
        
        Args:
            M_BH: Black hole mass [kg]
            a_spin: Dimensionless spin (0 to 1)
            B_field: Poloidal magnetic field [T]
        
        Returns:
            Dict with BZ power and related quantities
        """
        # Horizon radius (Schwarzschild for simplicity, Kerr correction small)
        R_H = G * M_BH / c**2
        
        # Spin angular velocity
        Omega_H = a_spin * c / (2 * R_H)
        
        # BZ power: P = (1/32) B² R_H^4 Ω_H² / c
        P_BZ = (1/32) * B_field**2 * R_H**4 * Omega_H**2 / c
        
        # Eddington luminosity for comparison
        L_Edd = 4 * math.pi * G * M_BH * m_p * c / SIGMA_T
        
        # Efficiency (P/L_Edd)
        efficiency = P_BZ / L_Edd if L_Edd > 0 else 0
        
        # UQFF magnetism contribution
        Um_contribution = self._compute_Um_jet(M_BH, a_spin, B_field, R_H)
        
        return {
            'P_BZ_W': P_BZ,
            'P_BZ_erg_s': P_BZ / erg_to_J,
            'R_horizon_m': R_H,
            'Omega_H_rad_s': Omega_H,
            'L_Eddington_W': L_Edd,
            'efficiency_P_Edd': efficiency,
            'Um_UQFF': Um_contribution,
            'equation': 'P_BZ = (1/32) B² R_H^4 Ω_H² / c'
        }
    
    def _compute_Um_jet(self, M_BH: float, a_spin: float, B_field: float, R_H: float) -> float:
        """UQFF Universal Magnetism for jet launching."""
        # Um = μ(t, ρ_vac) × (1 - e^{-γt}) × a c³ / (2 G M) × B² R_H⁴ / (4π c)
        mu_j = 3.38e20  # T pm³ base
        gamma = 5e-5 / (24 * 3600)  # per second
        t = 1e6 * year_to_s  # Typical AGN age
        
        oscillation = 1 - math.exp(-gamma * t)
        spin_factor = a_spin * c**3 / (2 * G * M_BH)
        power_factor = B_field**2 * R_H**4 / (4 * math.pi * c)
        
        Um = mu_j * oscillation * spin_factor * power_factor * rho_vac_SCm
        return Um


class RelativisticJetVelocityCalculator:
    """
    Relativistic jet velocity profile along propagation axis.
    
    Equation 15:
    Γ(z) ≈ Γ_0 + (z/z_acc)(Γ_∞ - Γ_0)(1 - e^{-z/z_acc})
    
    Where:
    - z: distance along jet axis
    - z_acc: acceleration length (~10 R_s)
    - Γ_0: initial Lorentz factor (~1-5 at base)
    - Γ_∞: terminal Lorentz factor (~10-100)
    
    Derivation: Collimation-acceleration via magnetic nozzle,
    energy conservation Γ ∝ 1/θ² (opening angle).
    """
    
    def compute(self, z: float, z_acc: float, Gamma_0: float = 2.0, 
                Gamma_inf: float = 50.0) -> Dict:
        """
        Compute Lorentz factor at distance z.
        
        Args:
            z: Distance along jet [m]
            z_acc: Acceleration length [m]
            Gamma_0: Initial Lorentz factor
            Gamma_inf: Terminal Lorentz factor
        
        Returns:
            Dict with Lorentz factor and velocity
        """
        # Velocity profile equation
        x = z / z_acc
        exp_term = 1 - math.exp(-x) if x < 100 else 1.0
        
        Gamma = Gamma_0 + x * (Gamma_inf - Gamma_0) * exp_term
        
        # Cap at terminal value
        Gamma = min(Gamma, Gamma_inf)
        
        # Velocity from Lorentz factor: v = c √(1 - 1/Γ²)
        if Gamma > 1:
            beta = math.sqrt(1 - 1/Gamma**2)
            v = beta * c
        else:
            beta = 0
            v = 0
        
        return {
            'Gamma': Gamma,
            'beta': beta,
            'velocity_m_s': v,
            'velocity_c': beta,
            'z_m': z,
            'z_acc_m': z_acc,
            'equation': 'Γ(z) = Γ_0 + (z/z_acc)(Γ_∞ - Γ_0)(1 - e^{-z/z_acc})'
        }


class AGNOutflowMomentumCalculator:
    """
    Terminal momentum from AGN feedback.
    
    Equation 23:
    p_term = Ṁ_out v_out = f(v_out) L_AGN / c
    
    Where:
    - Ṁ_out: outflow mass rate (~10-100 M_☉/yr)
    - v_out: outflow velocity (~1000 km/s)
    - L_AGN: AGN luminosity (~10^{44-46} erg/s)
    - f(v): velocity-dependent factor
    
    Derivation: Balance radiation pressure L/c with ram pressure ρv²,
    with optical depth correction for thick media.
    """
    
    def compute(self, L_AGN: float, v_out: float, tau: float = 1.0,
                M_dot_out: Optional[float] = None) -> Dict:
        """
        Compute AGN outflow momentum.
        
        Args:
            L_AGN: AGN luminosity [W]
            v_out: Outflow velocity [m/s]
            tau: Optical depth (default 1)
            M_dot_out: Optional mass outflow rate [kg/s]
        
        Returns:
            Dict with momentum and coupling efficiency
        """
        # Momentum flux from radiation: p = τ L / c
        p_rad = tau * L_AGN / c
        
        # If mass rate given, compute momentum directly
        if M_dot_out is not None:
            p_term = M_dot_out * v_out
        else:
            # Estimate mass rate from momentum: Ṁ = p / v
            p_term = p_rad
            M_dot_out = p_term / v_out
        
        # Kinetic power
        E_dot_kin = 0.5 * M_dot_out * v_out**2
        
        # Coupling efficiency
        epsilon_f = E_dot_kin / L_AGN if L_AGN > 0 else 0
        
        # UQFF buoyancy contribution
        F_UBii = self._compute_buoyancy_feedback(L_AGN, v_out)
        
        return {
            'p_term_kg_m_s': p_term,
            'M_dot_out_kg_s': M_dot_out,
            'M_dot_out_Msun_yr': M_dot_out / M_sun * year_to_s,
            'E_dot_kin_W': E_dot_kin,
            'epsilon_coupling': epsilon_f,
            'F_UBii_UQFF_N': F_UBii,
            'equation': 'p_term = Ṁ_out·v_out = f(v)·L_AGN/c'
        }
    
    def _compute_buoyancy_feedback(self, L_AGN: float, v_out: float) -> float:
        """UQFF Buoyancy force modulating feedback."""
        # F_UBii = F_rel × (E_AGN / E_LEP) × Q_wave × g(r,t)
        E_AGN = L_AGN * 1e7 * year_to_s  # Typical AGN lifetime energy
        Q_wave = 1e12  # THz resonance
        g_typical = 1e-8  # Cluster gravity ~10^{-8} m/s²
        
        F_UBii = -F_rel * (E_AGN / E_LEP) * Q_wave * g_typical
        return F_UBii


class JetPowerFromSpinCalculator:
    """
    Extended BZ jet power with EHT-consistent parameters.
    
    Equation 24:
    P_jet = (κ / 16π) Φ_BH² Ω_BH² / c
    
    Where:
    - Φ_BH: magnetic flux ~ B M² G² / c⁴
    - Ω_BH = a c³ / (2 G M)
    - κ: efficiency factor (~0.05-1)
    
    Derivation: Force-free electrodynamics in Kerr spacetime,
    normalized magnetic flux from EHT constraints.
    """
    
    def compute(self, M_BH: float, a_spin: float, B_field: float,
                kappa: float = 0.3) -> Dict:
        """
        Compute jet power using EHT-calibrated method.
        
        Args:
            M_BH: Black hole mass [kg]
            a_spin: Dimensionless spin
            B_field: Magnetic field at horizon [T]
            kappa: Efficiency factor
        
        Returns:
            Dict with jet power components
        """
        # BH parameters
        R_H = G * M_BH / c**2
        
        # Angular velocity
        Omega_BH = a_spin * c**3 / (2 * G * M_BH)
        
        # Magnetic flux (normalized)
        Phi_BH = B_field * math.pi * R_H**2
        
        # Jet power
        P_jet = (kappa / (16 * math.pi)) * Phi_BH**2 * Omega_BH**2 / c
        
        return {
            'P_jet_W': P_jet,
            'P_jet_erg_s': P_jet / erg_to_J,
            'Phi_BH_Wb': Phi_BH,
            'Omega_BH_rad_s': Omega_BH,
            'kappa': kappa,
            'R_horizon_m': R_H,
            'equation': 'P_jet = (κ/16π) Φ_BH² Ω_BH² / c'
        }


class FeedbackDutyCycleCalculator:
    """
    AGN feedback duty cycle from simulations.
    
    Equation 25:
    f_duty(t) = (1 - exp(-t/τ_cool)) × (1 + Ṁ_acc/Ṁ_Edd)^{-1}
    
    Where:
    - τ_cool: cooling time (~10^8 yr for clusters)
    - Ṁ_acc: accretion rate
    - Ṁ_Edd: Eddington rate
    
    Derivation: Probabilistic from TNG-Cluster simulations,
    feedback on when gas cooling exceeds time threshold.
    """
    
    def compute(self, t: float, tau_cool: float, M_dot_acc: float,
                M_BH: float) -> Dict:
        """
        Compute feedback duty cycle fraction.
        
        Args:
            t: Time since last feedback event [s]
            tau_cool: Cooling time [s]
            M_dot_acc: Accretion rate [kg/s]
            M_BH: Black hole mass [kg]
        
        Returns:
            Dict with duty cycle parameters
        """
        # Eddington accretion rate
        L_Edd = 4 * math.pi * G * M_BH * m_p * c / SIGMA_T
        epsilon_r = 0.1  # Radiative efficiency
        M_dot_Edd = L_Edd / (epsilon_r * c**2)
        
        # Duty cycle
        cooling_term = 1 - math.exp(-t / tau_cool) if tau_cool > 0 else 1
        accretion_suppression = 1 / (1 + M_dot_acc / M_dot_Edd) if M_dot_Edd > 0 else 0
        
        f_duty = cooling_term * accretion_suppression
        
        return {
            'f_duty': f_duty,
            'cooling_factor': cooling_term,
            'accretion_suppression': accretion_suppression,
            'M_dot_Edd_kg_s': M_dot_Edd,
            'M_dot_ratio': M_dot_acc / M_dot_Edd if M_dot_Edd > 0 else float('inf'),
            'equation': 'f_duty = (1 - e^{-t/τ_cool}) × (1 + Ṁ_acc/Ṁ_Edd)^{-1}'
        }


class WindTerminalVelocityCalculator:
    """
    Radiation-driven wind terminal velocity.
    
    Equation 70:
    v_∞ = √(2GM(1-Γ)/r_launch)
    
    Where:
    - Γ = L/L_Edd (Eddington ratio, ~1 for quasars)
    - r_launch: launch radius (~100 R_s)
    
    Derivation: Momentum balance with radiation pressure L/c
    against ram pressure ρv², energy limit for optically thick.
    """
    
    def compute(self, M_BH: float, L_AGN: float, r_launch: float) -> Dict:
        """
        Compute wind terminal velocity.
        
        Args:
            M_BH: Black hole mass [kg]
            L_AGN: AGN luminosity [W]
            r_launch: Launch radius [m]
        
        Returns:
            Dict with velocity and related parameters
        """
        # Eddington luminosity
        L_Edd = 4 * math.pi * G * M_BH * m_p * c / SIGMA_T
        
        # Eddington ratio
        Gamma_Edd = L_AGN / L_Edd if L_Edd > 0 else 0
        
        # Effective gravity reduction factor
        effective_factor = max(0, 1 - Gamma_Edd)
        
        # Terminal velocity
        if effective_factor > 0:
            v_inf = math.sqrt(2 * G * M_BH * effective_factor / r_launch)
        else:
            # Super-Eddington: use radiation-driven estimate
            v_inf = math.sqrt(2 * L_AGN / (4 * math.pi * r_launch**2 * c))
        
        return {
            'v_inf_m_s': v_inf,
            'v_inf_km_s': v_inf / 1000,
            'v_inf_c': v_inf / c,
            'Gamma_Eddington': Gamma_Edd,
            'L_Eddington_W': L_Edd,
            'r_launch_m': r_launch,
            'equation': 'v_∞ = √(2GM(1-Γ)/r_launch)'
        }


class IonizationParameterCalculator:
    """
    Ionization parameter for AGN feedback efficiency.
    
    Equation 71:
    U = Q_H / (4π r² n_H c)
    
    Where:
    - Q_H: ionizing photon rate (~10^{56} /s)
    - n_H: hydrogen density (~10² cm⁻³)
    - r: distance from AGN
    
    Derivation: Balance photoionization rate with recombination,
    defines ionization state of gas.
    """
    
    def compute(self, Q_H: float, r: float, n_H: float) -> Dict:
        """
        Compute ionization parameter U.
        
        Args:
            Q_H: Ionizing photon rate [photons/s]
            r: Distance from AGN [m]
            n_H: Hydrogen number density [m⁻³]
        
        Returns:
            Dict with ionization parameter
        """
        # Ionization parameter
        U = Q_H / (4 * math.pi * r**2 * n_H * c)
        
        # Log U (commonly used)
        log_U = math.log10(U) if U > 0 else float('-inf')
        
        # Ionization fraction estimate (simplified)
        alpha_B = 2.6e-19  # Case B recombination at 10^4 K [m³/s]
        Gamma_ion = Q_H / (4 * math.pi * r**2 * c)  # Ionization rate
        x_ion = Gamma_ion / (Gamma_ion + alpha_B * n_H)
        
        return {
            'U': U,
            'log_U': log_U,
            'x_ionized': x_ion,
            'Q_H_s': Q_H,
            'n_H_m3': n_H,
            'r_m': r,
            'equation': 'U = Q_H / (4π r² n_H c)'
        }


class FeedbackEnergyCouplingCalculator:
    """
    Feedback energy coupling efficiency to ISM.
    
    Equation 72:
    ε_f = Ė_kin / (Ṁ_acc c²) ≈ 0.05-0.1
    
    Where:
    - Ė_kin = (1/2) Ṁ_w v_w² (kinetic power)
    - Ṁ_acc: accretion rate
    
    Derivation: Couple fraction of rest-mass energy to kinetic,
    regulates M-σ relation.
    """
    
    def compute(self, M_dot_wind: float, v_wind: float, 
                M_dot_acc: float) -> Dict:
        """
        Compute feedback coupling efficiency.
        
        Args:
            M_dot_wind: Wind mass rate [kg/s]
            v_wind: Wind velocity [m/s]
            M_dot_acc: Accretion rate [kg/s]
        
        Returns:
            Dict with coupling efficiency
        """
        # Kinetic power
        E_dot_kin = 0.5 * M_dot_wind * v_wind**2
        
        # Accretion power
        epsilon_rad = 0.1  # Radiative efficiency
        P_acc = M_dot_acc * c**2
        
        # Coupling efficiency
        epsilon_f = E_dot_kin / P_acc if P_acc > 0 else 0
        
        # Mass loading factor
        eta_M = M_dot_wind / M_dot_acc if M_dot_acc > 0 else float('inf')
        
        return {
            'epsilon_coupling': epsilon_f,
            'E_dot_kin_W': E_dot_kin,
            'P_accretion_W': P_acc,
            'mass_loading': eta_M,
            'v_wind_km_s': v_wind / 1000,
            'equation': 'ε_f = Ė_kin / (Ṁ_acc c²)'
        }


class AGNFeedbackCalculator:
    """
    Master calculator combining all AGN feedback physics.
    
    Integrates:
    - Blandford-Znajek jet power
    - Relativistic jet profiles
    - Outflow momentum
    - Duty cycles
    - Wind velocities
    - Ionization and coupling
    
    UQFF Alignment:
    - Um dominates jet launching (EM >> gravity locally)
    - F_UBii buoyancy modulates cooling flow balance
    - Electric Universe validated for jet plasma dynamics
    """
    
    def __init__(self):
        self.bz_calc = BlandfordZnajekCalculator()
        self.jet_velocity_calc = RelativisticJetVelocityCalculator()
        self.outflow_calc = AGNOutflowMomentumCalculator()
        self.jet_power_calc = JetPowerFromSpinCalculator()
        self.duty_calc = FeedbackDutyCycleCalculator()
        self.wind_calc = WindTerminalVelocityCalculator()
        self.ionization_calc = IonizationParameterCalculator()
        self.coupling_calc = FeedbackEnergyCouplingCalculator()
    
    def compute_full_analysis(self, M_BH: float, a_spin: float, B_field: float,
                              L_AGN: float, M_dot_acc: float,
                              r_observation: float = None) -> Dict:
        """
        Full AGN feedback analysis.
        
        Args:
            M_BH: Black hole mass [kg]
            a_spin: Dimensionless spin (0-1)
            B_field: Magnetic field [T]
            L_AGN: AGN luminosity [W]
            M_dot_acc: Accretion rate [kg/s]
            r_observation: Observation radius [m]
        
        Returns:
            Comprehensive feedback analysis
        """
        # Default observation radius
        if r_observation is None:
            R_H = G * M_BH / c**2
            r_observation = 1000 * R_H  # 1000 R_s
        
        # BZ power
        bz_result = self.bz_calc.compute(M_BH, a_spin, B_field)
        
        # Jet power (EHT-calibrated)
        jet_power = self.jet_power_calc.compute(M_BH, a_spin, B_field)
        
        # Jet velocity at observation point
        R_H = G * M_BH / c**2
        z_acc = 10 * R_H
        jet_velocity = self.jet_velocity_calc.compute(r_observation, z_acc)
        
        # Wind properties
        wind = self.wind_calc.compute(M_BH, L_AGN, 100 * R_H)
        
        # Outflow momentum
        outflow = self.outflow_calc.compute(L_AGN, wind['v_inf_m_s'])
        
        # Duty cycle (typical cluster cooling time)
        tau_cool = 1e8 * year_to_s  # 100 Myr
        t_AGN = 1e7 * year_to_s  # 10 Myr AGN age
        duty = self.duty_calc.compute(t_AGN, tau_cool, M_dot_acc, M_BH)
        
        # Ionization parameter (typical NLR)
        Q_H = L_AGN / (20 * eV_to_J)  # ~20 eV per ionizing photon
        n_H = 1e8  # 100 cm⁻³
        ionization = self.ionization_calc.compute(Q_H, r_observation, n_H)
        
        # Coupling efficiency
        M_dot_wind = outflow['M_dot_out_kg_s']
        coupling = self.coupling_calc.compute(M_dot_wind, wind['v_inf_m_s'], M_dot_acc)
        
        return {
            'black_hole': {
                'M_BH_kg': M_BH,
                'M_BH_Msun': M_BH / M_sun,
                'a_spin': a_spin,
                'R_horizon_m': R_H
            },
            'jet_power': {
                'P_BZ_W': bz_result['P_BZ_W'],
                'P_jet_EHT_W': jet_power['P_jet_W'],
                'efficiency': bz_result['efficiency_P_Edd']
            },
            'jet_kinematics': {
                'Gamma': jet_velocity['Gamma'],
                'v_c': jet_velocity['velocity_c'],
                'z_observation_m': r_observation
            },
            'wind': {
                'v_inf_km_s': wind['v_inf_km_s'],
                'Gamma_Eddington': wind['Gamma_Eddington']
            },
            'feedback': {
                'p_term_kg_m_s': outflow['p_term_kg_m_s'],
                'epsilon_coupling': coupling['epsilon_coupling'],
                'mass_loading': coupling['mass_loading'],
                'f_duty': duty['f_duty']
            },
            'ionization': {
                'log_U': ionization['log_U'],
                'x_ionized': ionization['x_ionized']
            },
            'UQFF': {
                'F_UBii_buoyancy_N': outflow['F_UBii_UQFF_N'],
                'Um_magnetism': bz_result['Um_UQFF']
            }
        }


# ============== Pre-defined AGN Systems ==============

# M87* - EHT-observed SMBH
M87_STAR = {
    'name': 'M87*',
    'M_BH': 6.5e9 * M_sun,
    'a_spin': 0.9,
    'B_field': 1.0,  # ~10 Gauss
    'L_AGN': 1e42 / erg_to_J,  # ~10^{42} erg/s
    'M_dot_acc': 1e-4 * M_sun / year_to_s
}

# 3C 273 - Nearby bright quasar
QSO_3C273 = {
    'name': '3C 273',
    'M_BH': 8.9e8 * M_sun,
    'a_spin': 0.7,
    'B_field': 0.1,  # ~1 Gauss
    'L_AGN': 4e46 / erg_to_J,  # ~4×10^{46} erg/s
    'M_dot_acc': 10 * M_sun / year_to_s
}

# Sgr A* - Galactic center
SGR_A_STAR = {
    'name': 'Sgr A*',
    'M_BH': 4.3e6 * M_sun,
    'a_spin': 0.5,
    'B_field': 3.0,  # ~30 Gauss
    'L_AGN': 1e36 / erg_to_J,  # Very low luminosity
    'M_dot_acc': 1e-9 * M_sun / year_to_s
}

# NGC 1275 (Perseus A) - Classic feedback system
NGC_1275 = {
    'name': 'NGC 1275 (Perseus A)',
    'M_BH': 8e8 * M_sun,
    'a_spin': 0.6,
    'B_field': 0.5,
    'L_AGN': 1e45 / erg_to_J,
    'M_dot_acc': 1 * M_sun / year_to_s
}

# Registry of AGN systems
AGN_SYSTEMS = {
    'M87*': M87_STAR,
    '3C_273': QSO_3C273,
    'Sgr_A*': SGR_A_STAR,
    'NGC_1275': NGC_1275
}

# Calculator registry
AGN_FEEDBACK_CALCULATORS = {
    'BlandfordZnajek': BlandfordZnajekCalculator,
    'RelativisticJetVelocity': RelativisticJetVelocityCalculator,
    'AGNOutflowMomentum': AGNOutflowMomentumCalculator,
    'JetPowerFromSpin': JetPowerFromSpinCalculator,
    'FeedbackDutyCycle': FeedbackDutyCycleCalculator,
    'WindTerminalVelocity': WindTerminalVelocityCalculator,
    'IonizationParameter': IonizationParameterCalculator,
    'FeedbackEnergyCoupling': FeedbackEnergyCouplingCalculator,
    'AGNFeedback': AGNFeedbackCalculator
}


def run_demo():
    """Demonstrate AGN feedback calculations."""
    print("=" * 80)
    print("AGN FEEDBACK MODULE - Grok Deep Analysis Equations")
    print("=" * 80)
    
    calc = AGNFeedbackCalculator()
    
    for name, system in AGN_SYSTEMS.items():
        print(f"\n{'='*60}")
        print(f"System: {system['name']}")
        print(f"{'='*60}")
        
        result = calc.compute_full_analysis(
            M_BH=system['M_BH'],
            a_spin=system['a_spin'],
            B_field=system['B_field'],
            L_AGN=system['L_AGN'],
            M_dot_acc=system['M_dot_acc']
        )
        
        print(f"\nBlack Hole: {result['black_hole']['M_BH_Msun']:.2e} M_☉, spin={result['black_hole']['a_spin']}")
        print(f"BZ Jet Power: {result['jet_power']['P_BZ_W']:.2e} W")
        print(f"Jet Lorentz Factor: Γ = {result['jet_kinematics']['Gamma']:.1f}")
        print(f"Wind Velocity: {result['wind']['v_inf_km_s']:.0f} km/s")
        print(f"Eddington Ratio: Γ_Edd = {result['wind']['Gamma_Eddington']:.3f}")
        print(f"Coupling Efficiency: ε = {result['feedback']['epsilon_coupling']:.4f}")
        print(f"Duty Cycle: f = {result['feedback']['f_duty']:.3f}")
        print(f"Ionization: log U = {result['ionization']['log_U']:.2f}")
        print(f"UQFF Buoyancy Force: {result['UQFF']['F_UBii_buoyancy_N']:.2e} N")


if __name__ == '__main__':
    run_demo()
