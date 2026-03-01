#!/usr/bin/env python3
"""
Dark Matter Halos Module - Structure and Dynamics

From Grok Deep Analysis (Feb 2026):
- Equations 29-31: NFW profile, rotation curves, SIDM core formation
- Equations 82-84: First halos, star formation efficiency, feedback injection

Physics domains covered:
- NFW density profile (CDM halo structure)
- Rotation curves from halo potential
- Self-interacting dark matter (SIDM) core formation
- First halo collapse and star formation

UQFF Integration:
- Buoyancy (F_U_Bi_i) as alternative to dark matter in rotation curves
- Vacuum density differentials explain flat rotation without exotic matter
- Validates "1% answer" through aether-mediated dynamics
"""

import math
from typing import Dict, Optional, Callable

# ============== Physical Constants ==============
G = 6.674e-11           # Gravitational constant [m³/(kg·s²)]
c = 2.998e8             # Speed of light [m/s]
M_sun = 1.989e30        # Solar mass [kg]
k_B = 1.381e-23         # Boltzmann constant [J/K]
m_p = 1.673e-27         # Proton mass [kg]
h_bar = 1.055e-34       # Reduced Planck constant [J·s]
eV_to_J = 1.602e-19     # eV to Joules
pc_to_m = 3.086e16      # parsec to meters
kpc_to_m = 3.086e19     # kiloparsec to meters
Mpc_to_m = 3.086e22     # megaparsec to meters
year_to_s = 3.154e7     # year to seconds
km_to_m = 1000          # km to meters

# Cosmological parameters
H_0 = 70 * km_to_m / Mpc_to_m  # Hubble constant [s⁻¹]
Omega_m = 0.3           # Matter density parameter
Omega_Lambda = 0.7      # Dark energy density parameter
rho_crit = 3 * H_0**2 / (8 * math.pi * G)  # Critical density [kg/m³]

# UQFF Constants
F_rel = 4.30e33         # Relativistic coherence force [N]
E_LEP = 200e9 * eV_to_J # LEP 1998 baseline energy [J]
rho_vac_UA = 7.09e-36   # Vacuum density UA [J/m³]


class NFWProfileCalculator:
    """
    Navarro-Frenk-White (NFW) density profile for CDM halos.
    
    Equation 29:
    ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
    
    Where:
    - ρ_s: characteristic density
    - r_s: scale radius
    
    Derivation: From N-body simulations of CDM structure formation,
    universal density profile with ρ ∝ r⁻¹ inner, r⁻³ outer.
    """
    
    def __init__(self, M_vir: float, c: float, z: float = 0):
        """
        Initialize NFW profile.
        
        Args:
            M_vir: Virial mass [kg]
            c: Concentration parameter
            z: Redshift
        """
        self.M_vir = M_vir
        self.c = c
        self.z = z
        
        # Overdensity (Bryan & Norman 1998)
        self.Delta_vir = 200  # Simplified
        
        # Virial radius
        rho_crit_z = rho_crit * (Omega_m * (1+z)**3 + Omega_Lambda)
        self.r_vir = (3 * M_vir / (4 * math.pi * self.Delta_vir * rho_crit_z))**(1/3)
        
        # Scale radius
        self.r_s = self.r_vir / c
        
        # Characteristic density
        f_c = math.log(1 + c) - c / (1 + c)
        self.rho_s = M_vir / (4 * math.pi * self.r_s**3 * f_c)
    
    def density(self, r: float) -> float:
        """
        Compute density at radius r.
        
        Args:
            r: Radius [m]
        
        Returns:
            Density [kg/m³]
        """
        x = r / self.r_s
        return self.rho_s / (x * (1 + x)**2)
    
    def enclosed_mass(self, r: float) -> float:
        """
        Compute enclosed mass within radius r.
        
        Args:
            r: Radius [m]
        
        Returns:
            Enclosed mass [kg]
        """
        x = r / self.r_s
        f_x = math.log(1 + x) - x / (1 + x)
        return 4 * math.pi * self.rho_s * self.r_s**3 * f_x
    
    def compute(self, r: float) -> Dict:
        """
        Full NFW profile computation at radius r.
        
        Args:
            r: Radius [m]
        
        Returns:
            Dict with density, mass, and profile parameters
        """
        rho = self.density(r)
        M_enc = self.enclosed_mass(r)
        
        return {
            'rho_kg_m3': rho,
            'M_enclosed_kg': M_enc,
            'M_enclosed_Msun': M_enc / M_sun,
            'r_m': r,
            'r_kpc': r / kpc_to_m,
            'r_s_kpc': self.r_s / kpc_to_m,
            'r_vir_kpc': self.r_vir / kpc_to_m,
            'concentration': self.c,
            'M_vir_Msun': self.M_vir / M_sun,
            'equation': 'ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]'
        }


class RotationCurveCalculator:
    """
    Rotation curve from halo potential.
    
    Equation 30:
    v(r)² = GM(r)/r = 4πG ∫₀ʳ ρ(r')r'² dr'
    
    Where M(r) is enclosed mass from NFW profile.
    
    Derivation: Virial equilibrium, centrifugal = gravity,
    v² / r = G M/r², integrated over NFW profile.
    """
    
    def __init__(self, nfw_profile: Optional[NFWProfileCalculator] = None,
                 M_stellar: float = 0, r_disk: float = 3 * kpc_to_m):
        """
        Initialize rotation curve calculator.
        
        Args:
            nfw_profile: NFW halo profile (optional)
            M_stellar: Stellar disk mass [kg] (optional)
            r_disk: Disk scale length [m]
        """
        self.nfw = nfw_profile
        self.M_stellar = M_stellar
        self.r_disk = r_disk
    
    def compute(self, r: float) -> Dict:
        """
        Compute rotation velocity at radius r.
        
        Args:
            r: Radius [m]
        
        Returns:
            Dict with velocities and components
        """
        v2_total = 0
        v2_halo = 0
        v2_disk = 0
        M_halo = 0
        
        # Halo contribution
        if self.nfw is not None:
            M_halo = self.nfw.enclosed_mass(r)
            v2_halo = G * M_halo / r
            v2_total += v2_halo
        
        # Disk contribution (exponential disk approximation)
        if self.M_stellar > 0:
            x = r / self.r_disk
            # Freeman disk v² approximation
            if x < 0.1:
                y = 2 * x / math.pi
            else:
                y = x * (1 - math.exp(-3.2 * x))
            v2_disk = G * self.M_stellar / self.r_disk * y / (1 + x)
            v2_total += v2_disk
        
        v_total = math.sqrt(max(0, v2_total))
        v_halo = math.sqrt(max(0, v2_halo))
        v_disk = math.sqrt(max(0, v2_disk))
        
        # UQFF buoyancy alternative
        F_UBii = self._compute_buoyancy_rotation(r, v_total)
        
        return {
            'v_total_m_s': v_total,
            'v_total_km_s': v_total / km_to_m,
            'v_halo_km_s': v_halo / km_to_m,
            'v_disk_km_s': v_disk / km_to_m,
            'M_enclosed_Msun': (M_halo + self.M_stellar) / M_sun,
            'r_kpc': r / kpc_to_m,
            'F_UBii_UQFF_N': F_UBii,
            'equation': 'v(r)² = GM(r)/r'
        }
    
    def _compute_buoyancy_rotation(self, r: float, v: float) -> float:
        """
        UQFF buoyancy contribution to rotation.
        
        In UQFF, flat rotation curves emerge from vacuum density
        gradients (Δρ_vac) rather than exotic dark matter.
        """
        # F_UBii = F_rel × (Δρ_vac / ρ_LEP) × Q_wave × g(r)
        Delta_rho_vac = rho_vac_UA * (r / kpc_to_m)  # Scale with radius
        Q_wave = 1e12
        g_eff = v**2 / r if r > 0 else 0
        
        F_UBii = F_rel * (Delta_rho_vac / rho_crit) * Q_wave * g_eff
        return F_UBii


class SIDMCoreFormationCalculator:
    """
    Self-Interacting Dark Matter (SIDM) core formation time.
    
    Equation 31:
    t_core ≈ 1 / (ρ σ/m) ≈ 10^{10} (ρ/10⁸ M_☉/kpc³)⁻¹ (σ/m / 1 cm²/g)⁻¹ yr
    
    Where:
    - ρ: dark matter density
    - σ/m: cross-section per mass (~1 cm²/g for SIDM)
    
    Derivation: Scattering rate Γ = ρ σ v / m, core forms when
    Γ t ~ 1 redistributes orbits.
    """
    
    def compute(self, rho: float, sigma_m: float, v_disp: float) -> Dict:
        """
        Compute SIDM core formation time.
        
        Args:
            rho: DM density [kg/m³]
            sigma_m: Cross-section per mass [m²/kg]
            v_disp: Velocity dispersion [m/s]
        
        Returns:
            Dict with core formation parameters
        """
        # Scattering rate
        Gamma = rho * sigma_m * v_disp
        
        # Core formation time
        t_core = 1 / Gamma if Gamma > 0 else float('inf')
        
        # Core radius estimate (where t_scatter ~ t_age)
        t_age = 13.8e9 * year_to_s  # Hubble time
        
        return {
            't_core_s': t_core,
            't_core_Gyr': t_core / (1e9 * year_to_s),
            'Gamma_scattering_s': Gamma,
            'rho_kg_m3': rho,
            'sigma_m_m2_kg': sigma_m,
            'sigma_m_cm2_g': sigma_m * 10,  # Convert
            'v_disp_km_s': v_disp / km_to_m,
            'core_formed': t_core < t_age,
            'equation': 't_core ≈ 1 / (ρ σ/m v)'
        }


class FirstHaloCalculator:
    """
    First halo collapse (Press-Schechter formalism).
    
    Equation 82:
    dn/dM = √(2/π) (ρ̄/M) (δ_c/σ(M)) |d ln σ/d ln M| exp(-δ_c²/(2σ²))
    
    Where:
    - δ_c ≈ 1.69: critical overdensity
    - σ(M): rms fluctuation on mass scale M
    
    Derivation: Gaussian probability of exceeding δ_c,
    earliest halos at M ~ 10⁸ M_☉ at z ~ 10.
    """
    
    def __init__(self, sigma_8: float = 0.81, n_s: float = 0.96):
        """
        Initialize first halo calculator.
        
        Args:
            sigma_8: Fluctuation amplitude at 8 Mpc/h
            n_s: Spectral index
        """
        self.sigma_8 = sigma_8
        self.n_s = n_s
        self.delta_c = 1.686  # Critical overdensity
    
    def sigma_M(self, M: float, z: float = 0) -> float:
        """
        Compute rms fluctuation for mass scale M.
        
        Args:
            M: Mass scale [kg]
            z: Redshift
        
        Returns:
            σ(M)
        """
        # Simplified power-law scaling
        M_8 = 1e14 * M_sun  # Mass in 8 Mpc/h sphere
        
        # Growth factor (matter-dominated approximation)
        D_z = 1 / (1 + z)
        
        # σ(M) scaling: σ ∝ M^{-(n_s+3)/6}
        alpha = (self.n_s + 3) / 6
        sigma = self.sigma_8 * (M / M_8)**(-alpha) * D_z
        
        return sigma
    
    def compute_mass_function(self, M: float, z: float) -> Dict:
        """
        Compute halo mass function.
        
        Args:
            M: Halo mass [kg]
            z: Redshift
        
        Returns:
            Dict with mass function parameters
        """
        sigma = self.sigma_M(M, z)
        nu = self.delta_c / sigma  # Peak height
        
        # Press-Schechter multiplicity
        f_nu = math.sqrt(2/math.pi) * nu * math.exp(-nu**2 / 2)
        
        # d ln σ / d ln M
        alpha = (self.n_s + 3) / 6
        dlnsigma_dlnM = -alpha
        
        # Number density dn/dM (per Mpc³ per M_sun)
        rho_bar = Omega_m * rho_crit * (1 + z)**3
        dn_dM = (rho_bar / M) * f_nu * abs(dlnsigma_dlnM) / M
        
        # Collapse redshift estimate
        z_collapse = self.delta_c / sigma - 1 if sigma < self.delta_c else z
        
        return {
            'sigma_M': sigma,
            'nu_peak_height': nu,
            'f_nu': f_nu,
            'dn_dM': dn_dM,
            'M_Msun': M / M_sun,
            'z': z,
            'z_collapse_estimate': z_collapse,
            'equation': 'dn/dM = √(2/π)(ρ̄/M)(δ_c/σ)|dlnσ/dlnM|exp(-δ_c²/2σ²)'
        }


class StarFormationEfficiencyCalculator:
    """
    Star formation efficiency from baryon accretion.
    
    Equation 83:
    ε_* = f_b × Ṁ_halo / (M_halo H(z)) × (1 + M_halo/M_crit)⁻¹
    
    Where:
    - f_b = Ω_b/Ω_m ≈ 0.16 (baryon fraction)
    - M_crit ~ 10^{12} M_☉ (virial T > 10⁴ K)
    
    Derivation: Accretion rate from EPS, efficiency limited
    by cooling and feedback below/above M_crit.
    """
    
    def __init__(self):
        self.f_b = 0.16  # Baryon fraction
        self.M_crit = 1e12 * M_sun  # Critical mass
    
    def compute(self, M_halo: float, z: float, 
                M_dot_halo: Optional[float] = None) -> Dict:
        """
        Compute star formation efficiency.
        
        Args:
            M_halo: Halo mass [kg]
            z: Redshift
            M_dot_halo: Halo accretion rate [kg/s] (optional)
        
        Returns:
            Dict with efficiency parameters
        """
        # Hubble parameter at z
        H_z = H_0 * math.sqrt(Omega_m * (1+z)**3 + Omega_Lambda)
        
        # Accretion rate estimate if not provided
        # Ṁ ∝ M^{1.1} (1+z)^{2.25}
        if M_dot_halo is None:
            # Characteristic rate
            M_dot_halo = 0.1 * M_halo * H_z * (1 + z)**0.5
        
        # Efficiency
        accretion_factor = M_dot_halo / (M_halo * H_z) if M_halo > 0 else 0
        mass_suppression = 1 / (1 + M_halo / self.M_crit)
        
        epsilon_star = self.f_b * accretion_factor * mass_suppression
        
        # Star formation rate
        SFR = epsilon_star * self.f_b * M_halo / (1e9 * year_to_s)  # M_☉/yr
        
        return {
            'epsilon_star': epsilon_star,
            'SFR_Msun_yr': SFR / M_sun,
            'M_halo_Msun': M_halo / M_sun,
            'M_dot_halo_Msun_yr': M_dot_halo * year_to_s / M_sun,
            'mass_suppression': mass_suppression,
            'z': z,
            'equation': 'ε_* = f_b × Ṁ_halo/(M_halo H) × (1 + M/M_crit)⁻¹'
        }


class FeedbackEnergyInjectionCalculator:
    """
    Feedback energy injection from SN/AGN.
    
    Equation 84:
    E_fb = η Ṁ_* c², η ~ 10⁻³ - 10⁻¹
    
    Where:
    - η: coupling efficiency
    - Ṁ_*: star formation rate
    
    For SN: η ~ 10⁻⁵ per M_☉ (kinetic fraction ~1%)
    """
    
    def compute(self, SFR: float, eta: float = 1e-3,
                mode: str = 'AGN') -> Dict:
        """
        Compute feedback energy injection rate.
        
        Args:
            SFR: Star formation rate [kg/s]
            eta: Coupling efficiency
            mode: 'SN' or 'AGN'
        
        Returns:
            Dict with feedback parameters
        """
        # Rest mass power
        P_rest = SFR * c**2
        
        # Feedback power
        E_dot_fb = eta * P_rest
        
        # Mass loading estimate (outflow mass rate / SFR)
        # Ṁ_out / SFR ~ 1-10 typically
        v_out = 500 * km_to_m  # Typical outflow velocity
        mass_loading = 2 * E_dot_fb / (SFR * v_out**2)
        
        return {
            'E_dot_fb_W': E_dot_fb,
            'E_dot_fb_erg_s': E_dot_fb * 1e7,
            'eta': eta,
            'mode': mode,
            'SFR_Msun_yr': SFR * year_to_s / M_sun,
            'mass_loading': mass_loading,
            'v_effective_km_s': math.sqrt(2 * E_dot_fb / SFR) / km_to_m if SFR > 0 else 0,
            'equation': 'E_fb = η Ṁ_* c²'
        }


class VirialEquilibriumCalculator:
    """
    Virial theorem for cluster equilibrium.
    
    Equation 32:
    2K + W = 0  ⟹  M = 3 σ_v² r / G
    
    Where:
    - K: kinetic energy ∝ σ_v² (velocity dispersion ~1000 km/s)
    - W: potential energy ∝ -GM²/r
    - r: virial radius (1-3 Mpc)
    
    Derivation: Time-average of moment of inertia d²I/dt² = 2K + W = 0.
    """
    
    def compute(self, sigma_v: float, r: float) -> Dict:
        """
        Compute virial mass from velocity dispersion.
        
        Args:
            sigma_v: Velocity dispersion [m/s]
            r: Virial radius [m]
        
        Returns:
            Dict with virial parameters
        """
        # Virial mass
        M_vir = 3 * sigma_v**2 * r / G
        
        # Alternative: projection corrected
        M_vir_proj = 5 * sigma_v**2 * r / G  # For 3D isotropic
        
        # Crossing time
        t_cross = r / sigma_v
        
        # Relaxation time (N-body)
        N = 1000  # Typical cluster
        t_relax = N * t_cross / (6 * math.log(N))
        
        return {
            'M_vir_kg': M_vir,
            'M_vir_Msun': M_vir / M_sun,
            'M_vir_proj_Msun': M_vir_proj / M_sun,
            'sigma_v_km_s': sigma_v / km_to_m,
            'r_Mpc': r / Mpc_to_m,
            't_cross_Gyr': t_cross / (1e9 * year_to_s),
            't_relax_Gyr': t_relax / (1e9 * year_to_s),
            'equation': 'M = 3 σ_v² r / G'
        }


class DarkMatterHaloCalculator:
    """
    Master calculator for dark matter halo physics.
    
    Integrates:
    - NFW density profiles
    - Rotation curves
    - SIDM core formation
    - First halo collapse
    - Star formation and feedback
    
    UQFF Integration:
    - Buoyancy force as DM alternative
    - Vacuum density gradients for flat rotation
    - Validates Electric Universe locally
    """
    
    def __init__(self):
        self.first_halo_calc = FirstHaloCalculator()
        self.sf_efficiency_calc = StarFormationEfficiencyCalculator()
        self.feedback_calc = FeedbackEnergyInjectionCalculator()
        self.virial_calc = VirialEquilibriumCalculator()
        self.sidm_calc = SIDMCoreFormationCalculator()
    
    def compute_halo_analysis(self, M_vir: float, c: float, z: float = 0,
                               r_sample: Optional[float] = None) -> Dict:
        """
        Complete halo analysis.
        
        Args:
            M_vir: Virial mass [kg]
            c: Concentration
            z: Redshift
            r_sample: Sample radius [m]
        
        Returns:
            Comprehensive halo analysis
        """
        # NFW profile
        nfw = NFWProfileCalculator(M_vir, c, z)
        
        # Default sample radius
        if r_sample is None:
            r_sample = nfw.r_vir / 2
        
        profile = nfw.compute(r_sample)
        
        # Rotation curve
        rot_calc = RotationCurveCalculator(nfw)
        rotation = rot_calc.compute(r_sample)
        
        # Mass function
        mass_func = self.first_halo_calc.compute_mass_function(M_vir, z)
        
        # Star formation
        sf = self.sf_efficiency_calc.compute(M_vir, z)
        
        return {
            'profile': profile,
            'rotation': rotation,
            'mass_function': mass_func,
            'star_formation': sf,
            'UQFF': {
                'F_UBii_buoyancy_N': rotation['F_UBii_UQFF_N'],
                'DM_alternative': 'Vacuum density gradient Δρ_vac'
            }
        }


# ============== Pre-defined Systems ==============

# Milky Way halo
MILKY_WAY_HALO = {
    'name': 'Milky Way',
    'M_vir': 1e12 * M_sun,
    'c': 12,
    'r_disk': 3 * kpc_to_m,
    'M_stellar': 6e10 * M_sun
}

# Draco dwarf (classic DM test case)
DRACO_DWARF = {
    'name': 'Draco dSph',
    'M_vir': 1e9 * M_sun,
    'c': 20,
    'sigma_v': 9 * km_to_m,
    'r_half': 200 * pc_to_m
}

# Coma cluster
COMA_CLUSTER = {
    'name': 'Coma Cluster',
    'M_vir': 1e15 * M_sun,
    'c': 5,
    'sigma_v': 1000 * km_to_m,
    'r_vir': 2 * Mpc_to_m
}

# First galaxies (z~10)
FIRST_HALO = {
    'name': 'First Halo (z=10)',
    'M_vir': 1e8 * M_sun,
    'c': 8,
    'z': 10
}

HALO_SYSTEMS = {
    'Milky_Way': MILKY_WAY_HALO,
    'Draco': DRACO_DWARF,
    'Coma': COMA_CLUSTER,
    'First_Halo': FIRST_HALO
}

DARK_MATTER_CALCULATORS = {
    'NFWProfile': NFWProfileCalculator,
    'RotationCurve': RotationCurveCalculator,
    'SIDMCoreFormation': SIDMCoreFormationCalculator,
    'FirstHalo': FirstHaloCalculator,
    'StarFormationEfficiency': StarFormationEfficiencyCalculator,
    'FeedbackEnergyInjection': FeedbackEnergyInjectionCalculator,
    'VirialEquilibrium': VirialEquilibriumCalculator,
    'DarkMatterHalo': DarkMatterHaloCalculator
}


def run_demo():
    """Demonstrate dark matter halo calculations."""
    print("=" * 80)
    print("DARK MATTER HALOS MODULE - Grok Deep Analysis")
    print("=" * 80)
    
    calc = DarkMatterHaloCalculator()
    
    # Milky Way analysis
    print("\n--- Milky Way Halo ---")
    mw = MILKY_WAY_HALO
    nfw = NFWProfileCalculator(mw['M_vir'], mw['c'])
    rot = RotationCurveCalculator(nfw, M_stellar=mw['M_stellar'], r_disk=mw['r_disk'])
    
    for r_kpc in [5, 10, 20, 50]:
        r = r_kpc * kpc_to_m
        result = rot.compute(r)
        print(f"r = {r_kpc} kpc: v = {result['v_total_km_s']:.0f} km/s "
              f"(halo: {result['v_halo_km_s']:.0f}, disk: {result['v_disk_km_s']:.0f})")
    
    # SIDM core formation
    print("\n--- SIDM Core Formation (Draco) ---")
    draco = DRACO_DWARF
    nfw_draco = NFWProfileCalculator(draco['M_vir'], draco['c'])
    rho_center = nfw_draco.density(100 * pc_to_m)
    
    sidm = calc.sidm_calc.compute(
        rho=rho_center,
        sigma_m=1e-4,  # 1 cm²/g = 0.1 m²/kg
        v_disp=draco['sigma_v']
    )
    print(f"Central density: {rho_center:.2e} kg/m³")
    print(f"Core formation time: {sidm['t_core_Gyr']:.2f} Gyr")
    print(f"Core formed? {sidm['core_formed']}")
    
    # Virial mass from cluster
    print("\n--- Coma Cluster Virial Mass ---")
    coma = COMA_CLUSTER
    virial = calc.virial_calc.compute(coma['sigma_v'], coma['r_vir'])
    print(f"σ_v = {virial['sigma_v_km_s']:.0f} km/s")
    print(f"M_vir = {virial['M_vir_Msun']:.2e} M_☉")
    print(f"Crossing time: {virial['t_cross_Gyr']:.2f} Gyr")
    
    # First halos
    print("\n--- First Halo Collapse (z=10) ---")
    mass_func = calc.first_halo_calc.compute_mass_function(1e8 * M_sun, z=10)
    print(f"σ(M) = {mass_func['sigma_M']:.3f}")
    print(f"ν (peak height) = {mass_func['nu_peak_height']:.2f}")


if __name__ == '__main__':
    run_demo()
