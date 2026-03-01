#!/usr/bin/env python3
"""
Stellar Evolution Module - Main Sequence to Compact Remnants

From Grok Deep Analysis (Feb 2026):
- Equations 42-44: Main sequence lifetime, mass-luminosity, convective turnover
- Equations 58-60: Type Ia/Core-collapse SN, nucleosynthesis yields
- Equations 61-63: Planetary nebulae, ionization, wind dynamics

Physics domains covered:
- Main sequence stellar structure
- Nuclear burning and luminosity
- Convective mixing and turnover
- Supernova explosions and ejecta
- Planetary nebula evolution
- Chemical enrichment

UQFF Integration:
- Buoyancy forces in stellar interiors
- Vacuum energy contribution to stellar structure
- Validates compressed gravity in stellar cores
"""

import math
from typing import Dict, Optional

# ============== Physical Constants ==============
G = 6.674e-11           # Gravitational constant [m³/(kg·s²)]
c = 2.998e8             # Speed of light [m/s]
M_sun = 1.989e30        # Solar mass [kg]
R_sun = 6.957e8         # Solar radius [m]
L_sun = 3.828e26        # Solar luminosity [W]
T_sun = 5778            # Solar effective temperature [K]
k_B = 1.381e-23         # Boltzmann constant [J/K]
m_p = 1.673e-27         # Proton mass [kg]
m_e = 9.109e-31         # Electron mass [kg]
h_bar = 1.055e-34       # Reduced Planck constant [J·s]
sigma_SB = 5.670e-8     # Stefan-Boltzmann constant [W/(m²·K⁴)]
eV_to_J = 1.602e-19     # eV to Joules
MeV_to_J = 1.602e-13    # MeV to Joules
year_to_s = 3.154e7     # year to seconds

# Nuclear physics
Q_pp = 26.7 * MeV_to_J  # pp-chain energy release [J]
Q_CNO = 25.0 * MeV_to_J # CNO cycle energy release [J]
E_bind_Fe = 8.8 * MeV_to_J  # Binding energy per nucleon for Fe [J]

# UQFF Constants
F_rel = 4.30e33         # Relativistic coherence force [N]
rho_vac_SCm = 7.09e-37  # Vacuum density SCm [J/m³]


class MainSequenceLifetimeCalculator:
    """
    Main sequence lifetime from nuclear burning.
    
    Equation 42:
    τ_MS ≈ 10^{10} (M/M_☉)^{-2.5} yr
    
    Or more precisely:
    τ_MS = ε_nuc M / L = 0.007 × 0.1 × M c² / L
    
    Where:
    - ε_nuc ≈ 0.007 (H → He mass fraction converted)
    - 0.1 is core mass fraction
    
    Derivation: Nuclear fuel reservoir / luminosity,
    with L ∝ M^{3.5} giving τ ∝ M / M^{3.5}.
    """
    
    def compute(self, M: float, L: Optional[float] = None) -> Dict:
        """
        Compute main sequence lifetime.
        
        Args:
            M: Stellar mass [kg]
            L: Luminosity [W] (computed if not provided)
        
        Returns:
            Dict with lifetime parameters
        """
        M_ratio = M / M_sun
        
        # Mass-luminosity if not provided
        if L is None:
            if M_ratio < 0.43:
                L = L_sun * 0.23 * M_ratio**2.3
            elif M_ratio < 2:
                L = L_sun * M_ratio**4
            elif M_ratio < 55:
                L = L_sun * 1.4 * M_ratio**3.5
            else:
                L = L_sun * 32000 * M_ratio  # Eddington limit
        
        # Nuclear efficiency
        epsilon_nuc = 0.007  # H → He mass conversion
        f_core = 0.1        # Core mass fraction burning
        
        # Nuclear energy reservoir
        E_nuc = epsilon_nuc * f_core * M * c**2
        
        # Lifetime
        tau_MS = E_nuc / L
        
        # Simple scaling
        tau_simple = 1e10 * year_to_s * M_ratio**(-2.5)
        
        return {
            'tau_MS_s': tau_MS,
            'tau_MS_yr': tau_MS / year_to_s,
            'tau_MS_Gyr': tau_MS / (1e9 * year_to_s),
            'tau_simple_Gyr': tau_simple / (1e9 * year_to_s),
            'M_Msun': M_ratio,
            'L_Lsun': L / L_sun,
            'E_nuc_J': E_nuc,
            'equation': 'τ_MS ≈ 10^{10} (M/M_☉)^{-2.5} yr'
        }


class MassLuminosityCalculator:
    """
    Mass-luminosity relation for main sequence stars.
    
    Equation 43:
    L/L_☉ = (M/M_☉)^α, α ≈ 3.5 (1-50 M_☉)
    
    With variations:
    - α ≈ 2.3 for M < 0.43 M_☉
    - α ≈ 4 for 0.43-2 M_☉
    - α ≈ 3.5 for 2-55 M_☉
    - α ≈ 1 for M > 55 M_☉ (Eddington)
    
    Derivation: Hydrostatic equilibrium + radiative transport,
    L ∝ M³ μ⁴ / κ where κ is opacity.
    """
    
    def compute(self, M: float) -> Dict:
        """
        Compute luminosity from mass.
        
        Args:
            M: Stellar mass [kg]
        
        Returns:
            Dict with luminosity parameters
        """
        M_ratio = M / M_sun
        
        # Piecewise mass-luminosity
        if M_ratio < 0.43:
            L = L_sun * 0.23 * M_ratio**2.3
            alpha = 2.3
            regime = 'low-mass'
        elif M_ratio < 2:
            L = L_sun * M_ratio**4
            alpha = 4.0
            regime = 'solar-type'
        elif M_ratio < 55:
            L = L_sun * 1.4 * M_ratio**3.5
            alpha = 3.5
            regime = 'intermediate'
        else:
            L = L_sun * 32000 * M_ratio
            alpha = 1.0
            regime = 'Eddington-limited'
        
        # Effective temperature (L = 4πR²σT⁴)
        # Use R ∝ M^{0.8} approximation
        R = R_sun * M_ratio**0.8
        T_eff = (L / (4 * math.pi * R**2 * sigma_SB))**0.25
        
        return {
            'L_W': L,
            'L_Lsun': L / L_sun,
            'M_Msun': M_ratio,
            'alpha': alpha,
            'regime': regime,
            'R_Rsun': R / R_sun,
            'T_eff_K': T_eff,
            'equation': 'L/L_☉ = (M/M_☉)^α'
        }


class ConvectiveTurnoverCalculator:
    """
    Convective turnover time and mixing length.
    
    Equation 44:
    τ_conv = H_p / v_conv ≈ α_MLT H_p / (g T / ∇_ad)^{1/2}
    
    Where:
    - H_p = k_B T / (μ m_p g): pressure scale height
    - v_conv ≈ (F_conv / ρ)^{1/3}: convective velocity
    - α_MLT ≈ 1.5-2: mixing length parameter
    
    Derivation: Mixing length theory (MLT), convective blob
    travels ~H_p before dissipating.
    """
    
    def compute(self, M: float, R: float, L: float,
                alpha_MLT: float = 1.7) -> Dict:
        """
        Compute convective turnover time.
        
        Args:
            M: Stellar mass [kg]
            R: Stellar radius [m]
            L: Luminosity [W]
            alpha_MLT: Mixing length parameter
        
        Returns:
            Dict with convection parameters
        """
        # Surface gravity
        g = G * M / R**2
        
        # Effective temperature
        T_eff = (L / (4 * math.pi * R**2 * sigma_SB))**0.25
        
        # Pressure scale height (assuming μ ≈ 0.6)
        mu = 0.6
        H_p = k_B * T_eff / (mu * m_p * g)
        
        # Adiabatic gradient
        nabla_ad = 0.4  # For ideal gas
        
        # Convective velocity estimate
        F_conv = L / (4 * math.pi * R**2)  # Flux
        rho_phot = 3.5e-4  # kg/m³ (photospheric)
        v_conv = (F_conv / rho_phot)**(1/3)
        
        # Turnover time
        tau_conv = alpha_MLT * H_p / v_conv
        
        # Rossby number for activity (P_rot / τ_conv)
        # Assume P_rot scales with τ_MS
        P_rot_sun = 25 * 24 * 3600  # 25 days
        M_ratio = M / M_sun
        P_rot = P_rot_sun * M_ratio**0.5  # Rough scaling
        Ro = P_rot / tau_conv
        
        return {
            'tau_conv_s': tau_conv,
            'tau_conv_days': tau_conv / (24 * 3600),
            'H_p_m': H_p,
            'H_p_Rsun': H_p / R_sun,
            'v_conv_m_s': v_conv,
            'v_conv_km_s': v_conv / 1000,
            'g_m_s2': g,
            'T_eff_K': T_eff,
            'Rossby': Ro,
            'alpha_MLT': alpha_MLT,
            'equation': 'τ_conv = α_MLT H_p / v_conv'
        }


class TypeIaSupernovaCalculator:
    """
    Type Ia supernova energetics and light curve.
    
    Equation 58:
    E_Ia ≈ 10^{51} erg, M_Ni ≈ 0.6 M_☉
    L_peak ∝ M_Ni^{0.8} (Phillips relation)
    
    Where:
    - M_Ni: ⁵⁶Ni mass (powers light curve)
    - Decay: ⁵⁶Ni → ⁵⁶Co → ⁵⁶Fe
    
    Derivation: Thermonuclear runaway at M_Ch,
    C+O → ⁵⁶Ni releases ~10⁵¹ erg.
    """
    
    def compute(self, M_Ni: float = 0.6 * M_sun,
                M_ej: float = 1.4 * M_sun,
                t: float = 0) -> Dict:
        """
        Compute Type Ia parameters.
        
        Args:
            M_Ni: ⁵⁶Ni mass [kg]
            M_ej: Ejecta mass [kg]
            t: Time since explosion [s]
        
        Returns:
            Dict with SN Ia parameters
        """
        # Energy release (~10⁵¹ erg = 10⁴⁴ J)
        E_Ia = 1e44  # J
        
        # Kinetic energy to ejecta
        v_ej = math.sqrt(2 * E_Ia / M_ej)
        
        # ⁵⁶Ni decay (τ = 8.8 days)
        tau_Ni = 8.8 * 24 * 3600  # s
        # ⁵⁶Co decay (τ = 111.3 days)
        tau_Co = 111.3 * 24 * 3600  # s
        
        # Heating rate
        Q_Ni = 3.9e10  # W/kg
        Q_Co = 6.8e9   # W/kg
        
        # Luminosity (Arnett rule at peak)
        L_peak = M_Ni * Q_Ni  # At t ~ τ_Ni
        
        # Light curve
        if t > 0:
            x_Ni = math.exp(-t / tau_Ni)
            x_Co = tau_Co / (tau_Co - tau_Ni) * (math.exp(-t/tau_Ni) - math.exp(-t/tau_Co))
            L_t = M_Ni * (Q_Ni * x_Ni + Q_Co * x_Co)
        else:
            L_t = L_peak
        
        # Peak absolute magnitude
        M_B_peak = -19.3 - 2.5 * math.log10(M_Ni / (0.6 * M_sun))
        
        return {
            'E_Ia_J': E_Ia,
            'E_Ia_erg': E_Ia * 1e7,
            'M_Ni_Msun': M_Ni / M_sun,
            'M_ej_Msun': M_ej / M_sun,
            'v_ej_km_s': v_ej / 1000,
            'L_peak_W': L_peak,
            'L_peak_Lsun': L_peak / L_sun,
            'L_t_W': L_t,
            'M_B_peak': M_B_peak,
            't_days': t / (24 * 3600),
            'equation': 'E_Ia ≈ 10^{51} erg, L ∝ M_Ni'
        }


class CoreCollapseSupernovaCalculator:
    """
    Core-collapse supernova (Type II/Ib/Ic).
    
    Equation 59:
    E_cc ≈ 3 × 10^{53} erg (gravitational)
    E_kin ≈ 10^{51} erg (kinetic)
    E_ν ≈ 3 × 10^{53} erg (99% in neutrinos)
    
    Where:
    - E_bind = 3/5 G M_NS² / R_NS
    - Most energy escapes as neutrinos
    
    Derivation: Core collapse to NS, binding energy
    ~ 10⁵³ erg, 1% couples to shock → ejecta.
    """
    
    def compute(self, M_prog: float, M_NS: float = 1.4 * M_sun,
                R_NS: float = 12e3) -> Dict:
        """
        Compute core-collapse SN parameters.
        
        Args:
            M_prog: Progenitor mass [kg]
            M_NS: Neutron star mass [kg]
            R_NS: Neutron star radius [m]
        
        Returns:
            Dict with CC-SN parameters
        """
        # Binding energy of NS
        E_bind = 0.6 * G * M_NS**2 / R_NS
        
        # Neutrino energy (99%)
        E_nu = 0.99 * E_bind
        
        # Kinetic energy (~1%)
        E_kin = 0.01 * E_bind
        
        # Ejecta mass
        M_ej = M_prog - M_NS
        
        # Ejecta velocity
        v_ej = math.sqrt(2 * E_kin / M_ej) if M_ej > 0 else 0
        
        # Nickel mass estimate (scales with progenitor)
        M_Ni = 0.07 * M_sun * (M_prog / (15 * M_sun))**1.5
        
        # Peak luminosity
        L_peak = 1e9 * L_sun * (M_Ni / M_sun)
        
        return {
            'E_bind_J': E_bind,
            'E_bind_erg': E_bind * 1e7,
            'E_nu_J': E_nu,
            'E_kin_J': E_kin,
            'E_kin_erg': E_kin * 1e7,
            'M_prog_Msun': M_prog / M_sun,
            'M_NS_Msun': M_NS / M_sun,
            'M_ej_Msun': M_ej / M_sun,
            'M_Ni_Msun': M_Ni / M_sun,
            'v_ej_km_s': v_ej / 1000,
            'L_peak_Lsun': L_peak / L_sun,
            'equation': 'E_bind = 3/5 G M_NS²/R_NS ≈ 3×10^{53} erg'
        }


class NucleosynthesisYieldCalculator:
    """
    Nucleosynthesis yields from stellar evolution.
    
    Equation 60:
    M_i = ∫ X_i(m) dm (ejected mass of element i)
    
    Key yields:
    - C, O from He-burning
    - Si, S from O-burning (massive stars)
    - Fe-peak from explosive burning
    - r-process from NS mergers
    
    Derivation: Stellar models + nuclear networks,
    integrated over IMF for chemical evolution.
    """
    
    def compute(self, M_prog: float, metallicity: float = 0.02) -> Dict:
        """
        Compute nucleosynthesis yields.
        
        Args:
            M_prog: Progenitor mass [kg]
            metallicity: Initial Z
        
        Returns:
            Dict with yield parameters
        """
        M_ratio = M_prog / M_sun
        
        # Simplified yields (Nomoto et al. fitting)
        # Massive stars (M > 8 M_☉)
        if M_ratio > 8:
            # Carbon
            M_C = 0.1 * M_sun * (M_ratio / 20)**1.5
            # Oxygen
            M_O = 0.4 * M_sun * (M_ratio / 20)**2
            # Silicon
            M_Si = 0.05 * M_sun * (M_ratio / 20)**1.5
            # Iron
            M_Fe = 0.07 * M_sun * (M_ratio / 20)**0.5
            channel = 'CCSN'
        else:
            # AGB stars (1-8 M_☉)
            M_C = 0.01 * M_sun * M_ratio
            M_O = 0.02 * M_sun * M_ratio
            M_Si = 0
            M_Fe = 0
            channel = 'AGB'
        
        # Total metals ejected
        M_metals = M_C + M_O + M_Si + M_Fe
        
        # Enrichment factor
        delta_Z = M_metals / M_prog
        
        return {
            'M_C_Msun': M_C / M_sun,
            'M_O_Msun': M_O / M_sun,
            'M_Si_Msun': M_Si / M_sun,
            'M_Fe_Msun': M_Fe / M_sun,
            'M_metals_Msun': M_metals / M_sun,
            'delta_Z': delta_Z,
            'M_prog_Msun': M_ratio,
            'channel': channel,
            'initial_Z': metallicity,
            'equation': 'M_i = ∫ X_i(m) dm'
        }


class PlanetaryNebulaCalculator:
    """
    Planetary nebula expansion and ionization.
    
    Equation 61:
    R_PN(t) = v_exp × t, v_exp ≈ 20-30 km/s
    
    Equation 62:
    n_e² V = Q_H / α_B (ionization equilibrium)
    
    Equation 63:
    Ṁ_wind = L / (v_∞ c) (radiation-driven wind)
    
    Derivation: AGB mass loss → central star UV ionizes ejected envelope.
    """
    
    def compute(self, M_cs: float, L_cs: float, t: float,
                M_PN: float = 0.3 * M_sun,
                v_exp: float = 25e3) -> Dict:
        """
        Compute PN parameters.
        
        Args:
            M_cs: Central star mass [kg]
            L_cs: Central star luminosity [W]
            t: Age since ejection [s]
            M_PN: Nebula mass [kg]
            v_exp: Expansion velocity [m/s]
        
        Returns:
            Dict with PN parameters
        """
        # Nebula radius
        R_PN = v_exp * t
        
        # Nebula volume (thin shell approximation)
        dR = 0.2 * R_PN  # Shell thickness
        V_PN = 4/3 * math.pi * (R_PN**3 - (R_PN - dR)**3)
        
        # Average density
        rho_PN = M_PN / V_PN if V_PN > 0 else 0
        n_e = rho_PN / m_p  # Number density
        
        # Ionizing photon rate (central star T ~ 100,000 K)
        T_cs = 1e5  # K
        Q_H = L_cs * 0.3 / (13.6 * eV_to_J)  # Simple estimate
        
        # Recombination coefficient (case B)
        alpha_B = 2.6e-19  # m³/s at 10⁴ K
        
        # Stromgren radius
        R_S = (3 * Q_H / (4 * math.pi * alpha_B * n_e**2))**(1/3) if n_e > 0 else 0
        
        # Is PN fully ionized?
        ionization_bounded = R_PN < R_S
        
        return {
            'R_PN_m': R_PN,
            'R_PN_pc': R_PN / 3.086e16,
            't_yr': t / year_to_s,
            'v_exp_km_s': v_exp / 1000,
            'n_e_m3': n_e,
            'n_e_cm3': n_e / 1e6,
            'Q_H_s': Q_H,
            'R_S_pc': R_S / 3.086e16,
            'ionization_bounded': ionization_bounded,
            'L_cs_Lsun': L_cs / L_sun,
            'M_PN_Msun': M_PN / M_sun,
            'equation': 'R_PN = v_exp t, n_e² V = Q_H/α_B'
        }


class StellarWindCalculator:
    """
    Radiation-driven stellar wind.
    
    Equation (AGB/massive star winds):
    Ṁ = L / (v_∞ c) × Γ / (1 - Γ)
    
    Where:
    - Γ = L/L_Edd (Eddington ratio)
    - v_∞ ≈ v_esc × (Γ/(1-Γ))^{1/2}
    
    Derivation: Momentum transfer from radiation to dust/ions,
    CAK theory for line-driven winds.
    """
    
    def compute(self, M: float, L: float, R: float) -> Dict:
        """
        Compute stellar wind parameters.
        
        Args:
            M: Stellar mass [kg]
            L: Luminosity [W]
            R: Radius [m]
        
        Returns:
            Dict with wind parameters
        """
        # Eddington luminosity
        kappa_es = 0.034  # m²/kg (electron scattering)
        L_Edd = 4 * math.pi * G * M * c / kappa_es
        
        # Eddington ratio
        Gamma = L / L_Edd
        
        # Escape velocity
        v_esc = math.sqrt(2 * G * M / R)
        
        # Terminal velocity (CAK scaling)
        if Gamma < 1:
            v_inf = v_esc * math.sqrt(Gamma / (1 - Gamma)) * 2.6
        else:
            v_inf = 3 * v_esc  # Super-Eddington
        
        # Mass-loss rate
        if Gamma < 1:
            M_dot = L / (v_inf * c) * Gamma / (1 - Gamma)
        else:
            M_dot = L / (v_inf * c)  # Photon-momentum limited
        
        # Wind momentum
        p_dot = M_dot * v_inf
        
        return {
            'M_dot_kg_s': M_dot,
            'M_dot_Msun_yr': M_dot * year_to_s / M_sun,
            'v_inf_km_s': v_inf / 1000,
            'v_esc_km_s': v_esc / 1000,
            'Gamma': Gamma,
            'L_Edd_Lsun': L_Edd / L_sun,
            'p_dot_N': p_dot,
            'wind_efficiency': p_dot * c / L,
            'equation': 'Ṁ = L/(v_∞ c) × Γ/(1-Γ)'
        }


class StellarEvolutionCalculator:
    """
    Master calculator for stellar evolution physics.
    
    Integrates:
    - Main sequence structure
    - Nuclear burning and lifetime
    - Convection and mixing
    - Supernova explosions
    - Nucleosynthesis
    - Planetary nebulae and winds
    
    UQFF Integration:
    - Buoyancy effects in stellar cores
    - Vacuum energy contributions to luminosity
    """
    
    def __init__(self):
        self.ms_lifetime = MainSequenceLifetimeCalculator()
        self.ml_relation = MassLuminosityCalculator()
        self.convection = ConvectiveTurnoverCalculator()
        self.sn_ia = TypeIaSupernovaCalculator()
        self.sn_cc = CoreCollapseSupernovaCalculator()
        self.yields = NucleosynthesisYieldCalculator()
        self.pn = PlanetaryNebulaCalculator()
        self.wind = StellarWindCalculator()
    
    def compute_stellar_fate(self, M: float) -> Dict:
        """
        Compute complete stellar evolution fate.
        
        Args:
            M: Initial stellar mass [kg]
        
        Returns:
            Comprehensive evolution analysis
        """
        M_ratio = M / M_sun
        
        # Mass-luminosity
        ml = self.ml_relation.compute(M)
        L = ml['L_W']
        R = ml['R_Rsun'] * R_sun
        
        # Main sequence lifetime
        ms = self.ms_lifetime.compute(M, L)
        
        # Convection
        conv = self.convection.compute(M, R, L)
        
        # Determine fate
        if M_ratio < 0.08:
            fate = 'brown_dwarf'
            remnant_mass = M
        elif M_ratio < 8:
            fate = 'AGB_to_WD'
            remnant_mass = 0.5 * M_sun + 0.1 * (M - M_sun) if M > M_sun else 0.5 * M_sun
        elif M_ratio < 25:
            fate = 'CCSN_to_NS'
            remnant_mass = 1.4 * M_sun
        else:
            fate = 'BH_formation'
            remnant_mass = 0.3 * M
        
        # Yields
        yields = self.yields.compute(M)
        
        # Wind (for massive stars)
        wind = None
        if M_ratio > 10:
            wind = self.wind.compute(M, L, R)
        
        return {
            'M_initial_Msun': M_ratio,
            'fate': fate,
            'remnant_Msun': remnant_mass / M_sun,
            'mass_luminosity': ml,
            'main_sequence': ms,
            'convection': conv,
            'nucleosynthesis': yields,
            'wind': wind,
            'UQFF': {
                'F_core_N': F_rel * (M / M_sun)**2 * (rho_vac_SCm / 1e-37),
                'note': 'Buoyancy in degenerate cores'
            }
        }


# ============== Pre-defined Systems ==============

SUN = {
    'name': 'Sun',
    'M': 1.0 * M_sun,
    'R': 1.0 * R_sun,
    'L': 1.0 * L_sun,
    'T_eff': 5778
}

BETELGEUSE = {
    'name': 'Betelgeuse',
    'M': 15 * M_sun,
    'R': 900 * R_sun,
    'L': 1e5 * L_sun,
    'fate': 'CCSN_imminent'
}

SIRIUS_A = {
    'name': 'Sirius A',
    'M': 2.1 * M_sun,
    'R': 1.7 * R_sun,
    'L': 25 * L_sun
}

ETA_CARINAE = {
    'name': 'Eta Carinae',
    'M': 100 * M_sun,
    'L': 5e6 * L_sun,
    'fate': 'pair_instability_SN'
}

STELLAR_SYSTEMS = {
    'Sun': SUN,
    'Betelgeuse': BETELGEUSE,
    'Sirius_A': SIRIUS_A,
    'Eta_Carinae': ETA_CARINAE
}

STELLAR_EVOLUTION_CALCULATORS = {
    'MainSequenceLifetime': MainSequenceLifetimeCalculator,
    'MassLuminosity': MassLuminosityCalculator,
    'ConvectiveTurnover': ConvectiveTurnoverCalculator,
    'TypeIaSupernova': TypeIaSupernovaCalculator,
    'CoreCollapseSupernova': CoreCollapseSupernovaCalculator,
    'NucleosynthesisYield': NucleosynthesisYieldCalculator,
    'PlanetaryNebula': PlanetaryNebulaCalculator,
    'StellarWind': StellarWindCalculator,
    'StellarEvolution': StellarEvolutionCalculator
}


def run_demo():
    """Demonstrate stellar evolution calculations."""
    print("=" * 80)
    print("STELLAR EVOLUTION MODULE - Grok Deep Analysis")
    print("=" * 80)
    
    calc = StellarEvolutionCalculator()
    
    # Mass range analysis
    print("\n--- Main Sequence Lifetimes ---")
    for M_ratio in [0.5, 1.0, 3.0, 10.0, 30.0]:
        M = M_ratio * M_sun
        result = calc.ms_lifetime.compute(M)
        print(f"M = {M_ratio:.1f} M_☉: τ_MS = {result['tau_MS_Gyr']:.2f} Gyr, "
              f"L = {result['L_Lsun']:.1f} L_☉")
    
    # Stellar fates
    print("\n--- Stellar Fates ---")
    for M_ratio in [0.5, 1.0, 5.0, 15.0, 40.0]:
        M = M_ratio * M_sun
        fate = calc.compute_stellar_fate(M)
        print(f"M = {M_ratio:.0f} M_☉: {fate['fate']} → {fate['remnant_Msun']:.2f} M_☉")
    
    # Type Ia light curve
    print("\n--- Type Ia Supernova ---")
    sn = calc.sn_ia
    for t_days in [0, 15, 30, 60, 100]:
        t = t_days * 24 * 3600
        result = sn.compute(t=t)
        L_ratio = result['L_t_W'] / result['L_peak_W']
        print(f"t = {t_days} days: L/L_peak = {L_ratio:.3f}")
    
    # Core-collapse SN
    print("\n--- Core-Collapse SN (Betelgeuse proxy) ---")
    cc = calc.sn_cc.compute(15 * M_sun)
    print(f"E_bind = {cc['E_bind_erg']:.2e} erg")
    print(f"E_ν (99%) = {cc['E_nu_J']:.2e} J")
    print(f"v_ej = {cc['v_ej_km_s']:.0f} km/s")
    
    # Nucleosynthesis
    print("\n--- Nucleosynthesis Yields ---")
    yields = calc.yields.compute(20 * M_sun)
    print(f"M_O = {yields['M_O_Msun']:.3f} M_☉")
    print(f"M_Fe = {yields['M_Fe_Msun']:.3f} M_☉")


if __name__ == '__main__':
    run_demo()
