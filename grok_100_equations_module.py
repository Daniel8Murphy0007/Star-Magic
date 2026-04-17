"""
Grok 100+ Physics Equations Module
===================================
Comprehensive calculator module extracted from Grok AI conversation
covering astrophysics, cosmology, particle physics, and UQFF extensions.

Categories:
- Protostellar Jets & Outflows (Eqs 1-4)
- Galaxy Mergers & SFR (Eqs 5-7)
- Black Hole Growth (Eqs 8-9)
- Supernova Remnants (Eqs 10-11)
- Gravitational Waves (Eqs 12-13)
- Quasar Jets (Eqs 14-15)
- Neutron Stars (Eqs 16-18)
- Gamma-Ray Bursts (Eqs 19-20)
- CMB Anisotropies (Eqs 21-22)
- AGN Feedback (Eqs 23-25)
- Exoplanets (Eqs 26-28)
- Dark Matter Halos (Eqs 29-31)
- Galaxy Clusters (Eqs 32-34)
- Cosmic Voids (Eqs 35-36)
- Reionization (Eqs 37-38)
- ISM Turbulence (Eqs 39-41)
- Stellar Evolution (Eqs 42-44)
- Big Bang Nucleosynthesis (Eqs 45-46)
- Friedmann Cosmology (Eqs 47-49)
- Inflation (Eqs 50-52)
- Primordial GW (Eqs 53-54)
- Binary BH Mergers (Eqs 55-57)
- Supernovae (Eqs 58-60)
- Planetary Nebulae (Eqs 61-63)
- Cluster Collisions (Eqs 64-66)
- Star Clusters (Eqs 67-69)
- Quasar Feedback (Eqs 70-72)
- NS Binaries (Eqs 73-75)
- Cosmic Rays (Eqs 76-78)
- Intergalactic Medium (Eqs 79-81)
- First Stars (Eqs 82-84)
- Quantum Fluctuations (Eqs 85-87)
- MHD Dynamo (Eqs 88-90)
- Dark Energy (Eqs 91-93)
- BH Thermodynamics (Eqs 94-96)
- Loop Quantum Cosmology (Eqs 97-99)
- Extrasolar Planets (Eq 100)
- UQFF Extensions (Buoyancy, Magnetism, Gyro)

Author: Daniel T. Murphy / Grok AI Integration
Date: March 2026
"""

import math
from dpm_helpers import dpm_emergent_ug1, dpm_emergent_ug2

from typing import Dict, Any, Optional, List, Tuple

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================
class PhysicsConstants:
    """Fundamental physical constants used across all calculators."""
    G = 6.674e-11          # Gravitational constant (m^3/kg/s^2)
    c = 2.998e8            # Speed of light (m/s)
    h = 6.626e-34          # Planck constant (J*s)
    hbar = 1.055e-34       # Reduced Planck constant (J*s)
    k_B = 1.381e-23        # Boltzmann constant (J/K)
    e = 1.602e-19          # Elementary charge (C)
    m_e = 9.109e-31        # Electron mass (kg)
    m_p = 1.673e-27        # Proton mass (kg)
    m_H = 1.674e-27        # Hydrogen mass (kg)
    sigma_T = 6.652e-29    # Thomson cross-section (m^2)
    sigma_SB = 5.670e-8    # Stefan-Boltzmann constant (W/m^2/K^4)
    M_sun = 1.989e30       # Solar mass (kg)
    L_sun = 3.828e26       # Solar luminosity (W)
    R_sun = 6.96e8         # Solar radius (m)
    AU = 1.496e11          # Astronomical unit (m)
    pc = 3.086e16          # Parsec (m)
    kpc = 3.086e19         # Kiloparsec (m)
    Mpc = 3.086e22         # Megaparsec (m)
    yr = 3.156e7           # Year (s)
    Gyr = 3.156e16         # Gigayear (s)
    H_0 = 2.27e-18         # Hubble constant (1/s) ~ 70 km/s/Mpc
    mu_G = 1e-10           # Microgauss in Tesla
    
    # UQFF-specific constants
    E_LEP_1998 = 200e9 * e  # LEP energy 200 GeV in Joules
    F_rel = 4.30e33        # Relativistic coherence force (N)
    rho_vac_UA = 7.09e-36  # Aether vacuum density (J/m^3)
    rho_vac_SCm = 7.09e-37 # Superconductive vacuum density (J/m^3)
    gamma_UQFF = 5e-5      # UQFF decay constant (1/day)
    omega_c = 1.585e-8     # UQFF oscillation frequency (rad/s)


# =============================================================================
# PROTOSTELLAR JETS AND OUTFLOWS (Equations 1-4)
# =============================================================================
class ProtostellarJetCalculator:
    """
    Calculators for protostellar jet dynamics and shock structures.
    Based on Chandra observations of HH 154, L1551 IRS5, HH 211.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_01_angular_momentum_transport(self, M_dot: float, r: float, 
                                          M_star: float, T_B: float = 0) -> Dict[str, Any]:
        """
        Eq 1: Angular Momentum Transport in Accretion Disk
        
        dL/dt = M_dot * r^2 * Omega - T_B
        
        where Omega = sqrt(G*M/r^3) is Keplerian angular velocity
        
        Parameters:
            M_dot: Accretion rate (kg/s)
            r: Radius from protostar (m)
            M_star: Protostar mass (kg)
            T_B: Magnetic torque (N*m), default 0
            
        Returns:
            Dictionary with angular momentum change rate and components
        """
        Omega = math.sqrt(self.C.G * M_star / r**3)
        L_dot = M_dot * r**2 * Omega - T_B
        
        return {
            'equation': 'dL/dt = M_dot * r^2 * Omega - T_B',
            'latex': r'\frac{dL}{dt} = \dot{M} r^2 \Omega - T_B',
            'Omega': Omega,
            'Omega_unit': 'rad/s',
            'L_dot': L_dot,
            'L_dot_unit': 'kg*m^2/s^2',
            'accretion_term': M_dot * r**2 * Omega,
            'magnetic_torque': T_B
        }
    
    def eq_02_jet_velocity_mhd(self, M_star: float, r_0: float, 
                                r_A: float) -> Dict[str, Any]:
        """
        Eq 2: Jet Velocity from MHD Launching
        
        v_j = v_K * sqrt(r_A / r_0)
        
        where v_K = sqrt(G*M/r_0) is Keplerian velocity at footpoint
        
        Parameters:
            M_star: Protostar mass (kg)
            r_0: Disk footpoint radius (m)
            r_A: Alfven radius (m)
            
        Returns:
            Dictionary with jet velocity and components
        """
        v_K = math.sqrt(self.C.G * M_star / r_0)
        v_j = v_K * math.sqrt(r_A / r_0)
        
        return {
            'equation': 'v_j = v_K * sqrt(r_A / r_0)',
            'latex': r'v_j \approx v_K \sqrt{\frac{r_A}{r_0}}',
            'v_K': v_K,
            'v_K_km_s': v_K / 1000,
            'v_j': v_j,
            'v_j_km_s': v_j / 1000,
            'lever_arm': r_A / r_0
        }
    
    def eq_03_j_type_shock(self, rho_1: float, v_1: float, 
                           P_1: float, gamma: float = 5/3) -> Dict[str, Any]:
        """
        Eq 3: J-type Shock (Rankine-Hugoniot jump conditions)
        
        For discontinuous, non-magnetized or high-velocity shocks:
        - Mass: rho_1 * v_1 = rho_2 * v_2
        - Momentum: rho_1 * v_1^2 + P_1 = rho_2 * v_2^2 + P_2
        - Energy conservation
        
        Post-shock temperature T_2 ~ v_s^2
        
        Parameters:
            rho_1: Pre-shock density (kg/m^3)
            v_1: Pre-shock velocity (m/s)
            P_1: Pre-shock pressure (Pa)
            gamma: Adiabatic index (default 5/3 for monatomic)
            
        Returns:
            Dictionary with post-shock conditions
        """
        # Strong shock limit (M >> 1)
        compression = (gamma + 1) / (gamma - 1)
        rho_2 = rho_1 * compression
        v_2 = v_1 / compression
        
        # Post-shock pressure
        P_2 = rho_1 * v_1**2 * (2 * gamma) / (gamma + 1)
        
        # Post-shock temperature (proportional to v^2)
        T_2_scale = (v_1 / 1000)**2 * 1.38e5  # Approximate K for km/s shock
        
        return {
            'equation': 'Rankine-Hugoniot: rho_1*v_1 = rho_2*v_2',
            'latex': r'\rho_1 v_1 = \rho_2 v_2; \quad T_2 \propto v_s^2',
            'compression_ratio': compression,
            'rho_2': rho_2,
            'v_2': v_2,
            'P_2': P_2,
            'T_2_approx_K': T_2_scale,
            'X_ray_keV': T_2_scale * self.C.k_B / self.C.e / 1000
        }
    
    def eq_04_c_type_shock(self, v_s: float, z: float, 
                           L_d: float) -> Dict[str, Any]:
        """
        Eq 4: C-type Shock (Continuous, Magnetized with Ion-Neutral Drift)
        
        Velocity decreases gradually due to magnetic precursor damping:
        v(z) = v_s * exp(-z / L_d)
        
        where L_d is damping length ~ ion-neutral mean free path
        
        Parameters:
            v_s: Initial shock velocity (m/s)
            z: Distance along shock (m)
            L_d: Damping length (m)
            
        Returns:
            Dictionary with velocity profile
        """
        v_z = v_s * math.exp(-z / L_d)
        
        return {
            'equation': 'v(z) = v_s * exp(-z / L_d)',
            'latex': r'v(z) \approx v_s \exp\left(-\frac{z}{L_d}\right)',
            'v_initial': v_s,
            'v_at_z': v_z,
            'damping_length': L_d,
            'fractional_decay': v_z / v_s
        }


# =============================================================================
# GALAXY MERGERS AND STAR FORMATION RATE (Equations 5-7)
# =============================================================================
class GalaxyMergerSFRCalculator:
    """
    Calculators for galaxy merger rates and star formation.
    Based on Extended Press-Schechter formalism.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
        self.delta_c = 1.686  # Critical overdensity for collapse
    
    def eq_05_eps_merger_rate(self, sigma_M: float, sigma_m: float,
                               d_delta_c_dz: float, d_sigma_M_dM: float,
                               z: float) -> Dict[str, Any]:
        """
        Eq 5: Halo Merger Rate from Extended Press-Schechter Formalism
        
        dN/(dt*dM) = sqrt(2/pi) * (sigma_M/sigma_m) * |d_delta_c/dz| * 
                     exp(-delta_c^2 / (2*(sigma_m^2 - sigma_M^2))) * d_sigma_M/dM
        
        Merger rate proportional to (1+z)^m with m ~ 0.7-2.5
        
        Parameters:
            sigma_M: Variance on mass scale M
            sigma_m: Variance on smaller mass scale
            d_delta_c_dz: Derivative of critical overdensity with z
            d_sigma_M_dM: Derivative of variance with mass
            z: Redshift
            
        Returns:
            Dictionary with merger rate components
        """
        delta_c_z = self.delta_c * (1 + z)
        sigma_diff_sq = sigma_m**2 - sigma_M**2
        
        if sigma_diff_sq <= 0:
            return {'error': 'sigma_m must be > sigma_M'}
        
        prefactor = math.sqrt(2 / math.pi) * (sigma_M / sigma_m)
        exp_term = math.exp(-delta_c_z**2 / (2 * sigma_diff_sq))
        merger_rate = prefactor * abs(d_delta_c_dz) * exp_term * d_sigma_M_dM
        
        return {
            'equation': 'EPS merger rate',
            'latex': r'\frac{dN}{dt\,dM} = \sqrt{\frac{2}{\pi}} \frac{\sigma_M}{\sigma_m} \left|\frac{d\delta_c}{dz}\right| \exp\left(-\frac{\delta_c^2}{2(\sigma_m^2-\sigma_M^2)}\right) \frac{d\sigma_M}{dM}',
            'delta_c_z': delta_c_z,
            'prefactor': prefactor,
            'exp_term': exp_term,
            'merger_rate': merger_rate,
            'redshift_scaling': f'(1+z)^m where m ~ 0.7-2.5'
        }
    
    def eq_06_sfr_disk(self, M_gas: float, t_dyn: float, 
                        epsilon: float = 0.02) -> Dict[str, Any]:
        """
        Eq 6: Star Formation Rate in Disks (Quiescent Mode)
        
        M_dot_star = epsilon * M_gas / t_dyn
        
        where t_dyn = sqrt(3*pi / (32*G*rho)) is dynamical time
        
        Parameters:
            M_gas: Gas mass (kg)
            t_dyn: Dynamical time (s)
            epsilon: Star formation efficiency (default 0.02)
            
        Returns:
            Dictionary with SFR
        """
        M_dot_star = epsilon * M_gas / t_dyn
        
        return {
            'equation': 'M_dot_star = epsilon * M_gas / t_dyn',
            'latex': r'\dot{M}_* = \epsilon \frac{M_{gas}}{t_{dyn}}',
            'M_dot_star': M_dot_star,
            'M_dot_star_Msun_yr': M_dot_star / self.C.M_sun * self.C.yr,
            'efficiency': epsilon,
            'timescale': t_dyn / self.C.yr,
            'timescale_unit': 'years'
        }
    
    def eq_07_merger_burst_sfr(self, M_gas: float, t_orb: float,
                                epsilon_burst: float = 0.1) -> Dict[str, Any]:
        """
        Eq 7: Merger-Induced Starburst SFR
        
        M_dot_burst = M_dot_gas_inflow * epsilon_burst
        
        where M_dot_gas_inflow ~ M_gas / t_orb
        
        Parameters:
            M_gas: Gas mass (kg)
            t_orb: Orbital time (s)
            epsilon_burst: Burst efficiency (default 0.1)
            
        Returns:
            Dictionary with burst SFR (factor 10-100 enhancement)
        """
        M_dot_inflow = M_gas / t_orb
        M_dot_burst = M_dot_inflow * epsilon_burst
        
        return {
            'equation': 'M_dot_burst = M_gas/t_orb * epsilon_burst',
            'latex': r'\dot{M}_{burst} = \frac{M_{gas}}{t_{orb}} \times \epsilon_{burst}',
            'M_dot_inflow': M_dot_inflow,
            'M_dot_burst': M_dot_burst,
            'M_dot_burst_Msun_yr': M_dot_burst / self.C.M_sun * self.C.yr,
            'enhancement_factor': '10-100x quiescent'
        }


# =============================================================================
# BLACK HOLE GROWTH AND MASS FUNCTIONS (Equations 8-9)
# =============================================================================
class BlackHoleGrowthCalculator:
    """
    Calculators for black hole mass functions and accretion.
    Based on EPS formalism and Eddington-limited growth.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_08_eps_bh_mass_function(self, rho_bar: float, M: float,
                                    sigma_M: float, z: float) -> Dict[str, Any]:
        """
        Eq 8: EPS Black Hole Mass Function (Cumulative N(>M,z))
        
        N(>M,z) = rho_bar * integral_M^inf (dM'/M'^2) * erfc(delta_c(z) / (sqrt(2)*sigma(M',z)))
        
        Parameters:
            rho_bar: Mean density (kg/m^3)
            M: Black hole mass (kg)
            sigma_M: Variance at mass M
            z: Redshift
            
        Returns:
            Dictionary with cumulative number density
        """
        delta_c_z = 1.686 * (1 + z)
        x = delta_c_z / (math.sqrt(2) * sigma_M)
        erfc_approx = math.erfc(x) if x < 10 else 0
        
        # Simplified integration for single mass bin
        N_cumulative = rho_bar / M * erfc_approx
        
        return {
            'equation': 'EPS BH mass function',
            'latex': r'N(>M,z) = \bar{\rho} \int_M^\infty \frac{dM\'}{M\'^2} \mathrm{erfc}\left(\frac{\delta_c(z)}{\sqrt{2}\sigma(M\',z)}\right)',
            'delta_c_z': delta_c_z,
            'erfc_argument': x,
            'erfc_value': erfc_approx,
            'N_cumulative': N_cumulative,
            'unit': 'per m^3'
        }
    
    def eq_09_eddington_accretion(self, M_BH: float, 
                                   epsilon_r: float = 0.1) -> Dict[str, Any]:
        """
        Eq 9: BH Accretion Rate (Eddington-limited)
        
        M_dot_BH = 4*pi*G*M_BH*m_p / (epsilon_r * sigma_T * c)
        
        Salpeter timescale t_Sal ~ 45 Myr
        
        Parameters:
            M_BH: Black hole mass (kg)
            epsilon_r: Radiative efficiency (default 0.1)
            
        Returns:
            Dictionary with Eddington accretion rate
        """
        M_dot_Edd = 4 * math.pi * self.C.G * M_BH * self.C.m_p / (
            epsilon_r * self.C.sigma_T * self.C.c)
        
        # Salpeter timescale
        t_Sal = M_BH * epsilon_r * self.C.sigma_T * self.C.c / (
            4 * math.pi * self.C.G * self.C.m_p * M_BH)
        t_Sal = epsilon_r * self.C.sigma_T * self.C.c / (
            4 * math.pi * self.C.G * self.C.m_p)
        
        return {
            'equation': 'M_dot_Edd = 4*pi*G*M_BH*m_p / (epsilon_r * sigma_T * c)',
            'latex': r'\dot{M}_{BH} = \frac{4\pi G M_{BH} m_p}{\epsilon_r \sigma_T c}',
            'M_dot_Edd': M_dot_Edd,
            'M_dot_Edd_Msun_yr': M_dot_Edd / self.C.M_sun * self.C.yr,
            't_Salpeter': t_Sal,
            't_Salpeter_Myr': t_Sal / (1e6 * self.C.yr),
            'L_Eddington': M_dot_Edd * epsilon_r * self.C.c**2
        }


# =============================================================================
# SUPERNOVA REMNANTS (Equations 10-11)
# =============================================================================
class SupernovaRemnantCalculator:
    """
    Calculators for SNR expansion and particle acceleration.
    Based on Chandra observations of SN 1006, Kepler's SNR.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_10_sedov_taylor(self, E: float, rho: float, 
                           t: float) -> Dict[str, Any]:
        """
        Eq 10: Sedov-Taylor Expansion for Blast Wave
        
        R(t) = (E * t^2 / rho)^(1/5)
        
        Self-similar phase of SNR expansion
        
        Parameters:
            E: Explosion energy (J), typically 10^51 erg = 10^44 J
            rho: Ambient density (kg/m^3)
            t: Time since explosion (s)
            
        Returns:
            Dictionary with blast wave radius
        """
        R = (E * t**2 / rho)**(1/5)
        v = (2/5) * R / t  # Expansion velocity
        
        return {
            'equation': 'R(t) = (E*t^2/rho)^(1/5)',
            'latex': r'R(t) = \left(\frac{E t^2}{\rho}\right)^{1/5}',
            'R': R,
            'R_pc': R / self.C.pc,
            'v': v,
            'v_km_s': v / 1000,
            't_yr': t / self.C.yr
        }
    
    def eq_11_dsa_acceleration(self, p: float, u_s: float, 
                                r_d: float) -> Dict[str, Any]:
        """
        Eq 11: Diffusive Shock Acceleration (DSA) for Particles
        
        dp/dt = (4/3) * u_s^2 / r_d * p
        
        Produces power-law spectrum N(E) ~ E^{-2}
        
        Parameters:
            p: Particle momentum (kg*m/s)
            u_s: Shock velocity (m/s)
            r_d: Diffusion length (m)
            
        Returns:
            Dictionary with acceleration rate
        """
        dp_dt = (4/3) * u_s**2 / r_d * p
        
        # Acceleration timescale
        t_acc = 3 * r_d / u_s**2 * self.C.c
        
        return {
            'equation': 'dp/dt = (4/3) * u_s^2 / r_d * p',
            'latex': r'\frac{dp}{dt} = \frac{4}{3} \frac{u_s^2}{r_d} p',
            'dp_dt': dp_dt,
            't_acc': t_acc,
            't_acc_yr': t_acc / self.C.yr,
            'spectral_index': -2,
            'note': 'First-order Fermi acceleration'
        }


# =============================================================================
# GRAVITATIONAL WAVES FROM MERGERS (Equations 12-13)
# =============================================================================
class GravitationalWaveMergerCalculator:
    """
    Calculators for GW inspiral and ringdown.
    Based on LIGO/Virgo/KAGRA detections.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_12_chirp_mass(self, m1: float, m2: float) -> Dict[str, Any]:
        """
        Eq 12: Chirp Mass from Inspiral Frequency
        
        M_chirp = (m1*m2)^(3/5) / (m1+m2)^(1/5)
        
        Can also be derived from f and f_dot
        
        Parameters:
            m1: Mass of first object (kg)
            m2: Mass of second object (kg)
            
        Returns:
            Dictionary with chirp mass
        """
        M_chirp = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
        
        return {
            'equation': 'M_chirp = (m1*m2)^(3/5) / (m1+m2)^(1/5)',
            'latex': r'\mathcal{M} = \frac{(m_1 m_2)^{3/5}}{(m_1+m_2)^{1/5}}',
            'M_chirp': M_chirp,
            'M_chirp_Msun': M_chirp / self.C.M_sun,
            'total_mass': m1 + m2,
            'mass_ratio': m1 / m2 if m2 > 0 else float('inf')
        }
    
    def eq_13_ringdown_qnm(self, M_f: float, a_f: float = 0.69) -> Dict[str, Any]:
        """
        Eq 13: Ringdown Quasi-Normal Modes (Post-Merger)
        
        f_QNM = c^3 / (2*pi*G*M_f) * f(0.3737 + 0.088*a_f + ...)
        
        Parameters:
            M_f: Final black hole mass (kg)
            a_f: Dimensionless spin parameter (default 0.69)
            
        Returns:
            Dictionary with QNM frequency and damping
        """
        # Approximate fit for l=2, m=2 mode
        f_factor = 0.3737 + 0.088 * a_f
        f_QNM = self.C.c**3 / (2 * math.pi * self.C.G * M_f) * f_factor
        
        # Damping time ~ M
        tau_QNM = 2 * M_f * self.C.G / self.C.c**3 / f_factor * 10  # Q-factor ~ 10
        
        return {
            'equation': 'f_QNM = c^3/(2*pi*G*M_f) * f(a_f)',
            'latex': r'f_{QNM} = \frac{c^3}{2\pi G M_f} f(a_f)',
            'f_QNM': f_QNM,
            'f_QNM_Hz': f_QNM,
            'tau_QNM': tau_QNM,
            'tau_QNM_ms': tau_QNM * 1000,
            'Q_factor': f_QNM * tau_QNM
        }


# =============================================================================
# QUASAR JETS (Equations 14-15)
# =============================================================================
class QuasarJetCalculator:
    """
    Calculators for relativistic quasar jets.
    Based on EHT M87* observations.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_14_blandford_znajek(self, M_BH: float, a: float, 
                                B: float) -> Dict[str, Any]:
        """
        Eq 14: Blandford-Znajek Jet Power (Energy from BH Spin)
        
        P_BZ = (1/32) * B^2 * R_H^4 * Omega_H^2 / c
        
        where R_H = G*M/c^2 and Omega_H = a*c/(2*R_H)
        
        Parameters:
            M_BH: Black hole mass (kg)
            a: Dimensionless spin (0-1)
            B: Poloidal magnetic field (T)
            
        Returns:
            Dictionary with jet power
        """
        R_H = self.C.G * M_BH / self.C.c**2
        Omega_H = a * self.C.c / (2 * R_H)
        
        P_BZ = (1/32) * B**2 * R_H**4 * Omega_H**2 / self.C.c
        
        return {
            'equation': 'P_BZ = (1/32) * B^2 * R_H^4 * Omega_H^2 / c',
            'latex': r'P_{BZ} = \frac{1}{32} B^2 R_H^4 \Omega_H^2 / c',
            'R_H': R_H,
            'R_H_Rs': R_H / (2 * self.C.G * M_BH / self.C.c**2),
            'Omega_H': Omega_H,
            'P_BZ': P_BZ,
            'P_BZ_erg_s': P_BZ * 1e7,
            'P_BZ_Lsun': P_BZ / self.C.L_sun
        }
    
    def eq_15_relativistic_jet_velocity(self, z: float, z_acc: float,
                                         Gamma_0: float = 1.5,
                                         Gamma_inf: float = 15) -> Dict[str, Any]:
        """
        Eq 15: Relativistic Jet Velocity Profile
        
        Gamma(z) = Gamma_0 + (z/z_acc) * (Gamma_inf - Gamma_0) * (1 - exp(-z/z_acc))
        
        Parameters:
            z: Distance along jet (m)
            z_acc: Acceleration length (m), typically 10 R_s
            Gamma_0: Initial Lorentz factor (default 1.5)
            Gamma_inf: Terminal Lorentz factor (default 15)
            
        Returns:
            Dictionary with Lorentz factor profile
        """
        Gamma = Gamma_0 + (z / z_acc) * (Gamma_inf - Gamma_0) * (
            1 - math.exp(-z / z_acc))
        
        beta = math.sqrt(1 - 1/Gamma**2) if Gamma > 1 else 0
        
        return {
            'equation': 'Gamma(z) sigmoid acceleration',
            'latex': r'\Gamma(z) \approx \Gamma_0 + \frac{z}{z_{acc}}(\Gamma_\infty - \Gamma_0)(1 - e^{-z/z_{acc}})',
            'Gamma': Gamma,
            'beta': beta,
            'v': beta * self.C.c,
            'superluminal_factor': Gamma * beta  # Apparent transverse velocity
        }


# =============================================================================
# NEUTRON STARS (Equations 16-18)
# =============================================================================
class NeutronStarCalculator:
    """
    Calculators for neutron star structure and pulsars.
    Based on Chandra/NICER observations of PSR J0030+0451.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_16_tov_equation(self, P: float, r: float, m_r: float,
                           rho: float) -> Dict[str, Any]:
        """
        Eq 16: Tolman-Oppenheimer-Volkoff Equation (Hydrostatic Equilibrium)
        
        dP/dr = -G*m(r)*rho(r)/r^2 * (1 + P/(rho*c^2)) * 
                (1 + 4*pi*r^3*P/(m(r)*c^2)) * (1 - 2*G*m(r)/(r*c^2))^(-1)
        
        Parameters:
            P: Pressure (Pa)
            r: Radius (m)
            m_r: Enclosed mass (kg)
            rho: Density (kg/m^3)
            
        Returns:
            Dictionary with pressure gradient
        """
        # Relativistic corrections
        factor1 = 1 + P / (rho * self.C.c**2)
        factor2 = 1 + 4 * math.pi * r**3 * P / (m_r * self.C.c**2)
        factor3 = 1 - 2 * self.C.G * m_r / (r * self.C.c**2)
        
        dP_dr = -dpm_emergent_ug1(m_r, r) * rho * factor1 * factor2 / factor3
        
        return {
            'equation': 'TOV equation',
            'latex': r'\frac{dP}{dr} = -\frac{Gm\rho}{r^2}(1+\frac{P}{\rho c^2})(1+\frac{4\pi r^3 P}{mc^2})(1-\frac{2Gm}{rc^2})^{-1}',
            'dP_dr': dP_dr,
            'relativistic_correction_1': factor1,
            'relativistic_correction_2': factor2,
            'redshift_factor': 1 / math.sqrt(factor3) if factor3 > 0 else float('inf')
        }
    
    def eq_17_pulsar_spindown(self, P: float, P_dot: float) -> Dict[str, Any]:
        """
        Eq 17: Pulsar Spin-Down and Characteristic Age
        
        tau = P / (2 * P_dot)
        
        Parameters:
            P: Period (s)
            P_dot: Period derivative (s/s, dimensionless)
            
        Returns:
            Dictionary with characteristic age and related quantities
        """
        tau = P / (2 * P_dot)
        
        # Surface magnetic field estimate
        B_surface = 3.2e19 * math.sqrt(P * P_dot)  # Gauss
        
        # Spin-down luminosity
        Omega = 2 * math.pi / P
        I = 1e45  # Moment of inertia in cgs, ~10^38 kg*m^2
        E_dot = 4 * math.pi**2 * I * P_dot / P**3 * 1e-7  # Convert to Watts
        
        return {
            'equation': 'tau = P / (2 * P_dot)',
            'latex': r'\tau = \frac{P}{2\dot{P}}',
            'tau': tau,
            'tau_yr': tau / self.C.yr,
            'tau_Myr': tau / (1e6 * self.C.yr),
            'B_surface_G': B_surface,
            'B_surface_T': B_surface * 1e-4,
            'E_dot': E_dot,
            'E_dot_erg_s': E_dot * 1e7
        }
    
    def eq_18_glitch_fractional(self, Delta_nu: float, nu: float) -> Dict[str, Any]:
        """
        Eq 18: Pulsar Glitch Fractional Jump
        
        Delta_nu / nu ~ 10^{-9} to 10^{-6}
        
        Parameters:
            Delta_nu: Frequency jump (Hz)
            nu: Spin frequency (Hz)
            
        Returns:
            Dictionary with glitch parameters
        """
        fractional = Delta_nu / nu
        
        # Superfluid fraction estimate (after Anderson & Itoh)
        I_s_I = fractional * 100  # Very rough estimate
        
        return {
            'equation': 'Delta_nu / nu = fractional glitch',
            'latex': r'\frac{\Delta\nu}{\nu} \sim 10^{-9} - 10^{-6}',
            'Delta_nu': Delta_nu,
            'nu': nu,
            'fractional': fractional,
            'log10_fractional': math.log10(abs(fractional)) if fractional != 0 else float('-inf'),
            'superfluid_coupling': 'Vortex unpinning model'
        }


# =============================================================================
# GAMMA-RAY BURSTS (Equations 19-20)
# =============================================================================
class GammaRayBurstCalculator:
    """
    Calculators for GRB fireball and afterglow.
    Based on Fermi/Swift observations.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_19_fireball_lorentz(self, r: float, R_0: float,
                                eta: float = 300) -> Dict[str, Any]:
        """
        Eq 19: Fireball Expansion in Relativistic Regime
        
        Gamma(r) = r/R_0 for r < R_s (coasting)
        Gamma = eta for r > R_s (saturation)
        
        where R_s = eta^2 * R_0
        
        Parameters:
            r: Current radius (m)
            R_0: Initial radius (m), typically 10^7 cm
            eta: Baryon loading parameter (default 300)
            
        Returns:
            Dictionary with Lorentz factor
        """
        R_s = eta**2 * R_0
        
        if r < R_s:
            Gamma = r / R_0
            regime = 'coasting'
        else:
            Gamma = eta
            regime = 'saturated'
        
        return {
            'equation': 'Gamma(r) = r/R_0 (r<R_s), eta (r>R_s)',
            'latex': r'\Gamma(r) = \frac{r}{R_0} \; (r<R_s), \quad \Gamma = \eta \; (r>R_s)',
            'Gamma': Gamma,
            'R_s': R_s,
            'regime': regime,
            'eta': eta
        }
    
    def eq_20_afterglow_flux(self, nu: float, t: float, p: float = 2.3,
                             nu_m: float = 1e14, nu_c: float = 1e17,
                             F_max: float = 1e-26) -> Dict[str, Any]:
        """
        Eq 20: Afterglow Synchrotron Emission (External Shock)
        
        F_nu ~ nu^{-(p-1)/2} * t^{-3(p-1)/4} for nu_m < nu < nu_c
        
        Parameters:
            nu: Observation frequency (Hz)
            t: Time since burst (s)
            p: Electron spectral index (default 2.3)
            nu_m: Minimum synchrotron frequency (Hz)
            nu_c: Cooling frequency (Hz)
            F_max: Peak flux normalization (W/m^2/Hz)
            
        Returns:
            Dictionary with afterglow flux
        """
        alpha_nu = -(p - 1) / 2
        alpha_t = -3 * (p - 1) / 4
        
        F_nu = F_max * (nu / nu_m)**alpha_nu * (t / 86400)**alpha_t
        
        return {
            'equation': 'F_nu ~ nu^{-(p-1)/2} * t^{-3(p-1)/4}',
            'latex': r'F_\nu \propto \nu^{-(p-1)/2} t^{-3(p-1)/4}',
            'F_nu': F_nu,
            'spectral_index': alpha_nu,
            'temporal_index': alpha_t,
            'p': p,
            'regime': 'slow cooling' if nu_m < nu < nu_c else 'other'
        }


# =============================================================================
# COSMIC MICROWAVE BACKGROUND (Equations 21-22)
# =============================================================================
class CMBCalculator:
    """
    Calculators for CMB anisotropies and recombination.
    Based on Planck observations.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
        self.T_CMB = 2.725  # CMB temperature today (K)
        self.z_rec = 1100   # Recombination redshift
    
    def eq_21_angular_power_spectrum(self, l: int, A_s: float = 2.1e-9,
                                      n_s: float = 0.96) -> Dict[str, Any]:
        """
        Eq 21: Angular Power Spectrum C_l (From Density Fluctuations)
        
        C_l = (2/pi) * integral k^2 dk P(k) |Delta_l^T(k)|^2
        
        Simplified amplitude scaling with multipole l
        
        Parameters:
            l: Multipole moment
            A_s: Scalar amplitude (default 2.1e-9)
            n_s: Scalar spectral index (default 0.96)
            
        Returns:
            Dictionary with power spectrum
        """
        # Simplified Sachs-Wolfe plateau for large scales
        if l < 30:
            C_l = A_s * 2 * math.pi / (l * (l + 1))
        else:
            # Acoustic peaks approximation
            peak_scale = 220  # First peak
            C_l = A_s * 2 * math.pi / (l * (l + 1)) * (
                1 + 0.5 * math.cos(math.pi * l / peak_scale)**2)
        
        D_l = l * (l + 1) * C_l / (2 * math.pi)  # Conventional D_l
        
        return {
            'equation': 'C_l = (2/pi) * integral P(k) |Delta_l|^2 dk',
            'latex': r'C_l = \frac{2}{\pi} \int k^2 dk P(k) |\Delta_l^T(k)|^2',
            'l': l,
            'C_l': C_l,
            'D_l': D_l,
            'D_l_uK2': D_l * (self.T_CMB * 1e6)**2,
            'n_s': n_s
        }
    
    def eq_22_recombination_optical_depth(self, z: float, 
                                           Omega_b: float = 0.05) -> Dict[str, Any]:
        """
        Eq 22: Recombination Time and Optical Depth
        
        tau(z) = integral_z^inf n_e(z') sigma_T c |dt/dz'| dz'
        
        Parameters:
            z: Redshift
            Omega_b: Baryon density parameter (default 0.05)
            
        Returns:
            Dictionary with optical depth
        """
        # Simplified: tau ~ 1 at z ~ 1100
        n_e0 = Omega_b * 3 * self.C.H_0**2 / (8 * math.pi * self.C.G * self.C.m_p)
        n_e_z = n_e0 * (1 + z)**3
        
        # Approximate optical depth
        tau = self.C.sigma_T * n_e_z * self.C.c / self.C.H_0 / (1 + z)**1.5
        
        return {
            'equation': 'tau(z) = integral n_e sigma_T c |dt/dz| dz',
            'latex': r'\tau(z) = \int_z^\infty n_e(z\') \sigma_T c \left|\frac{dt}{dz\'}\right| dz\'',
            'z': z,
            'n_e_z': n_e_z,
            'tau_approx': tau,
            'last_scattering': z == self.z_rec,
            'T_at_z': self.T_CMB * (1 + z)
        }


# =============================================================================
# AGN FEEDBACK (Equations 23-25)
# =============================================================================
class AGNFeedbackCalculator:
    """
    Calculators for AGN outflows and feedback.
    Based on Chandra observations of clusters and EHT M87*.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_23_outflow_momentum(self, L_AGN: float, v_out: float,
                                f_v: float = 1.0) -> Dict[str, Any]:
        """
        Eq 23: Outflow Momentum from AGN Feedback
        
        p_term = M_dot_out * v_out = f(v_out) * L_AGN / c
        
        Parameters:
            L_AGN: AGN luminosity (W)
            v_out: Outflow velocity (m/s)
            f_v: Velocity-dependent factor (default 1.0)
            
        Returns:
            Dictionary with terminal momentum
        """
        M_dot_out = f_v * L_AGN / (self.C.c * v_out)
        p_term = M_dot_out * v_out
        
        return {
            'equation': 'p_term = f(v) * L_AGN / c',
            'latex': r'p_{term} = \dot{M}_{out} v_{out} = f(v_{out}) \frac{L_{AGN}}{c}',
            'M_dot_out': M_dot_out,
            'M_dot_out_Msun_yr': M_dot_out / self.C.M_sun * self.C.yr,
            'p_term': p_term,
            'v_out_km_s': v_out / 1000
        }
    
    def eq_24_bz_jet_power_updated(self, Phi_BH: float, Omega_BH: float,
                                    kappa: float = 0.05) -> Dict[str, Any]:
        """
        Eq 24: Jet Power from BH Spin (Updated BZ for EHT)
        
        P_jet = (kappa / 16pi) * Phi_BH^2 * Omega_BH^2 / c
        
        Parameters:
            Phi_BH: Magnetic flux at horizon (Wb)
            Omega_BH: BH angular velocity (rad/s)
            kappa: Efficiency factor (default 0.05)
            
        Returns:
            Dictionary with jet power
        """
        P_jet = (kappa / (16 * math.pi)) * Phi_BH**2 * Omega_BH**2 / self.C.c
        
        return {
            'equation': 'P_jet = (kappa/16pi) * Phi_BH^2 * Omega_BH^2 / c',
            'latex': r'P_{jet} = \frac{\kappa}{16\pi} \Phi_{BH}^2 \Omega_{BH}^2 / c',
            'P_jet': P_jet,
            'P_jet_erg_s': P_jet * 1e7,
            'kappa': kappa
        }
    
    def eq_25_feedback_duty_cycle(self, t: float, tau_cool: float,
                                   M_dot_acc: float, M_dot_Edd: float) -> Dict[str, Any]:
        """
        Eq 25: Feedback Duty Cycle (Time Evolution)
        
        f_duty(t) = (1 - exp(-t/tau_cool)) * (1 + M_dot_acc/M_dot_Edd)^{-1}
        
        Parameters:
            t: Time (s)
            tau_cool: Cooling time (s), typically 10^8 yr
            M_dot_acc: Accretion rate (kg/s)
            M_dot_Edd: Eddington rate (kg/s)
            
        Returns:
            Dictionary with duty cycle
        """
        f_duty = (1 - math.exp(-t / tau_cool)) / (1 + M_dot_acc / M_dot_Edd)
        
        return {
            'equation': 'f_duty = (1 - exp(-t/tau_cool)) / (1 + M_dot_acc/M_dot_Edd)',
            'latex': r'f_{duty}(t) = \frac{1-e^{-t/\tau_{cool}}}{1 + \dot{M}_{acc}/\dot{M}_{Edd}}',
            'f_duty': f_duty,
            'f_duty_percent': f_duty * 100,
            'tau_cool_Gyr': tau_cool / self.C.Gyr,
            'Eddington_ratio': M_dot_acc / M_dot_Edd
        }


# =============================================================================
# EXOPLANETS (Equations 26-28)
# =============================================================================
class ExoplanetCalculator:
    """
    Calculators for exoplanet atmospheres and dynamics.
    Based on Chandra X-ray observations of host stars.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_26_photoevaporative_mass_loss(self, L_X: float, R_p: float,
                                          M_p: float, epsilon: float = 0.15,
                                          K_xi: float = 1.0) -> Dict[str, Any]:
        """
        Eq 26: Photoevaporative Mass Loss from X-ray Irradiation
        
        M_dot_evap = epsilon * L_X * R_p^3 / (G * M_p^2) * K(xi)
        
        Parameters:
            L_X: Host star X-ray luminosity (W)
            R_p: Planet radius (m)
            M_p: Planet mass (kg)
            epsilon: Heating efficiency (default 0.15)
            K_xi: Base density factor (default 1.0)
            
        Returns:
            Dictionary with mass loss rate
        """
        M_dot = epsilon * L_X * R_p**3 / (self.C.G * M_p**2) * K_xi
        
        return {
            'equation': 'M_dot = epsilon * L_X * R_p^3 / (G * M_p^2) * K(xi)',
            'latex': r'\dot{M}_{evap} = \epsilon \frac{L_X R_p^3}{G M_p^2} K(\xi)',
            'M_dot': M_dot,
            'M_dot_g_s': M_dot * 1000,
            'M_dot_Earth_Gyr': M_dot / (5.97e24) * self.C.Gyr,
            'energy_limited': True
        }
    
    def eq_27_type_i_migration(self, M_p: float, M_star: float,
                                Sigma: float, h_r: float, r: float) -> Dict[str, Any]:
        """
        Eq 27: Type-I Migration Torque (Disk-Planet Interaction)
        
        Gamma = -C * (M_p/M_star)^2 * Sigma * r^4 * Omega^2 / h^2
        
        Parameters:
            M_p: Planet mass (kg)
            M_star: Star mass (kg)
            Sigma: Disk surface density (kg/m^2)
            h_r: Aspect ratio h/r
            r: Orbital radius (m)
            
        Returns:
            Dictionary with migration torque
        """
        Omega = math.sqrt(self.C.G * M_star / r**3)
        C = 2.5  # Lindblad torque coefficient
        
        Gamma = -C * (M_p / M_star)**2 * Sigma * r**4 * Omega**2 / h_r**2
        
        # Migration timescale
        L = M_p * r**2 * Omega
        tau_mig = abs(L / Gamma) if Gamma != 0 else float('inf')
        
        return {
            'equation': 'Gamma = -C * (M_p/M_star)^2 * Sigma * r^4 * Omega^2 / h^2',
            'latex': r'\Gamma = -C \left(\frac{M_p}{M_*}\right)^2 \Sigma r^4 \Omega^2 / h^2',
            'Gamma': Gamma,
            'tau_mig': tau_mig,
            'tau_mig_yr': tau_mig / self.C.yr,
            'direction': 'inward' if Gamma < 0 else 'outward'
        }
    
    def eq_28_radial_velocity(self, P: float, M_p: float, M_star: float,
                               i: float = math.pi/2, e: float = 0) -> Dict[str, Any]:
        """
        Eq 28: Radial Velocity Semi-Amplitude (Detection of Motion)
        
        K = (2*pi*G/P)^{1/3} * M_p*sin(i) / (M_star + M_p)^{2/3} / sqrt(1-e^2)
        
        Parameters:
            P: Orbital period (s)
            M_p: Planet mass (kg)
            M_star: Star mass (kg)
            i: Inclination angle (rad, default pi/2)
            e: Eccentricity (default 0)
            
        Returns:
            Dictionary with RV semi-amplitude
        """
        M_total = M_star + M_p
        K = (2 * math.pi * self.C.G / P)**(1/3) * M_p * math.sin(i) / M_total**(2/3)
        K /= math.sqrt(1 - e**2)
        
        return {
            'equation': 'K = (2*pi*G/P)^{1/3} * M_p*sin(i) / M_total^{2/3} / sqrt(1-e^2)',
            'latex': r'K = \left(\frac{2\pi G}{P}\right)^{1/3} \frac{M_p \sin i}{(M_*+M_p)^{2/3}} \frac{1}{\sqrt{1-e^2}}',
            'K': K,
            'K_m_s': K,
            'detectable': K > 0.1  # Typical 0.1 m/s precision
        }


# =============================================================================
# DARK MATTER HALOS (Equations 29-31)
# =============================================================================
class DarkMatterHaloCalculator:
    """
    Calculators for dark matter halo structure.
    Based on simulations and rotation curve observations.
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_29_nfw_profile(self, r: float, rho_s: float, 
                          r_s: float) -> Dict[str, Any]:
        """
        Eq 29: NFW Density Profile (CDM Halo Structure)
        
        rho(r) = rho_s / ((r/r_s) * (1 + r/r_s)^2)
        
        Parameters:
            r: Radius (m)
            rho_s: Characteristic density (kg/m^3)
            r_s: Scale radius (m)
            
        Returns:
            Dictionary with density profile
        """
        x = r / r_s
        rho = rho_s / (x * (1 + x)**2)
        
        return {
            'equation': 'rho(r) = rho_s / ((r/r_s) * (1 + r/r_s)^2)',
            'latex': r'\rho(r) = \frac{\rho_s}{(r/r_s)(1+r/r_s)^2}',
            'rho': rho,
            'x': x,
            'inner_slope': -1,
            'outer_slope': -3
        }
    
    def eq_30_rotation_curve(self, r: float, rho_s: float, 
                              r_s: float) -> Dict[str, Any]:
        """
        Eq 30: Rotation Curve from Halo Potential
        
        v(r)^2 = G*M(r)/r = 4*pi*G * integral_0^r rho(r') r'^2 dr'
        
        For NFW: v(r) ~ sqrt(ln(1+x) - x/(1+x)) where x = r/r_s
        
        Parameters:
            r: Radius (m)
            rho_s: Characteristic density (kg/m^3)
            r_s: Scale radius (m)
            
        Returns:
            Dictionary with rotation velocity
        """
        x = r / r_s
        # NFW enclosed mass function
        f_x = math.log(1 + x) - x / (1 + x)
        M_r = 4 * math.pi * rho_s * r_s**3 * f_x
        
        v_circ = math.sqrt(self.C.G * M_r / r)
        
        return {
            'equation': 'v^2 = G*M(r)/r for NFW profile',
            'latex': r'v(r)^2 = \frac{GM(r)}{r} \propto \sqrt{\ln(1+x) - \frac{x}{1+x}}',
            'v_circ': v_circ,
            'v_circ_km_s': v_circ / 1000,
            'M_enclosed': M_r,
            'M_enclosed_Msun': M_r / self.C.M_sun
        }
    
    def eq_31_sidm_core_time(self, rho: float, sigma_m: float,
                              v_vir: float) -> Dict[str, Any]:
        """
        Eq 31: SIDM Core Formation Time (Self-Interaction Effects)
        
        t_core ~ 1 / (rho * sigma/m * v)
        
        Parameters:
            rho: Central density (kg/m^3)
            sigma_m: Cross-section per mass (m^2/kg), typically 1 cm^2/g = 0.1 m^2/kg
            v_vir: Virial velocity (m/s)
            
        Returns:
            Dictionary with core formation time
        """
        t_core = 1 / (rho * sigma_m * v_vir)
        
        return {
            'equation': 't_core ~ 1 / (rho * sigma/m * v)',
            'latex': r't_{core} \approx \frac{1}{\rho \cdot (\sigma/m) \cdot v}',
            't_core': t_core,
            't_core_Gyr': t_core / self.C.Gyr,
            'sigma_m_cm2_g': sigma_m * 10,
            'cores_form': t_core < 10 * self.C.Gyr
        }


# =============================================================================
# ADDITIONAL CALCULATORS (Equations 32-100)
# The following are condensed implementations
# =============================================================================

class GalaxyClusterCalculator:
    """Eqs 32-34: Virial dynamics, lensing, merger shocks"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_32_virial_mass(self, sigma_v: float, r: float) -> Dict[str, Any]:
        """Virial theorem mass estimate: M = 3*sigma_v^2*r/G"""
        M = 3 * sigma_v**2 * r / self.C.G
        return {'M': M, 'M_Msun': M / self.C.M_sun, 
                'equation': r'M = \frac{3\sigma_v^2 r}{G}'}
    
    def eq_33_einstein_radius(self, M: float, D_LS: float, D_L: float, 
                               D_S: float) -> Dict[str, Any]:
        """Strong lensing arc radius"""
        theta_E = math.sqrt(4 * self.C.G * M / self.C.c**2 * D_LS / (D_L * D_S))
        return {'theta_E': theta_E, 'theta_E_arcsec': theta_E * 206265,
                'equation': r'\theta_E = \sqrt{\frac{4GM}{c^2}\frac{D_{LS}}{D_L D_S}}'}
    
    def eq_34_mach_from_density(self, rho_2_rho_1: float, 
                                 gamma: float = 5/3) -> Dict[str, Any]:
        """Merger shock Mach number from X-ray jump"""
        M_sq = 2 * rho_2_rho_1 / (gamma + 1 - rho_2_rho_1 * (gamma - 1))
        M = math.sqrt(M_sq) if M_sq > 0 else 0
        return {'M': M, 'rho_ratio': rho_2_rho_1,
                'equation': r'\mathcal{M} = \sqrt{\frac{2\rho_2/\rho_1}{\gamma+1-(\gamma-1)\rho_2/\rho_1}}'}


class CosmicVoidCalculator:
    """Eqs 35-36: Void evolution and outflow"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_35_void_density_evolution(self, delta_v0: float, a: float,
                                      Omega_m: float = 0.3) -> Dict[str, Any]:
        """Linear underdensity evolution"""
        # Simplified: delta grows slower than overdensities
        delta_v = -3/5 * delta_v0 * a  # Linear approximation
        return {'delta_v': delta_v, 'a': a,
                'equation': r'\delta_v(a) \approx -\frac{3}{5}\delta_{v0} a'}
    
    def eq_36_void_outflow(self, delta: float, r: float, H: float,
                            f: float = 0.5) -> Dict[str, Any]:
        """Peculiar velocity in voids"""
        # Simplified linear theory
        v_pec = -f * H * r * (-delta) / 3
        return {'v_pec': v_pec, 'v_pec_km_s': v_pec / 1000,
                'equation': r'v_{pec} = -\frac{fH}{3}\int \delta(r) r dr / r^2'}


class ReionizationCalculator:
    """Eqs 37-38: Ionization history"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_37_ionization_fraction(self, n_dot_gamma: float, f_star: float,
                                   epsilon_esc: float, alpha_B: float,
                                   n_e: float, C: float = 3) -> Dict[str, Any]:
        """dx_e/dt balance"""
        dx_dt = n_dot_gamma * epsilon_esc * f_star - alpha_B * n_e**2 * C
        return {'dx_dt': dx_dt, 'production': n_dot_gamma * epsilon_esc * f_star,
                'recombination': alpha_B * n_e**2 * C,
                'equation': r'\frac{dx_e}{dt} = \dot{n}_\gamma \epsilon_{esc} f_* - \alpha_B n_e^2 C'}
    
    def eq_38_bubble_radius(self, N_dot_gamma: float, t: float,
                            n_H: float) -> Dict[str, Any]:
        """Stromgren sphere growth"""
        R_b = (3 * N_dot_gamma * t / (4 * math.pi * n_H))**(1/3)
        return {'R_b': R_b, 'R_b_Mpc': R_b / self.C.Mpc,
                'equation': r'R_b(t) = \left(\frac{3\dot{N}_\gamma t}{4\pi n_H}\right)^{1/3}'}


class ISMTurbulenceCalculator:
    """Eqs 39-41: Jeans length, Alfven waves, cascade"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_39_jeans_length(self, c_s: float, rho: float) -> Dict[str, Any]:
        """Gravitational collapse length"""
        lambda_J = math.sqrt(math.pi * c_s**2 / (self.C.G * rho))
        return {'lambda_J': lambda_J, 'lambda_J_pc': lambda_J / self.C.pc,
                'equation': r'\lambda_J = \sqrt{\frac{\pi c_s^2}{G\rho}}'}
    
    def eq_40_alfven_velocity(self, B: float, rho: float) -> Dict[str, Any]:
        """MHD wave speed"""
        v_A = B / math.sqrt(4 * math.pi * rho * 1e-7)  # SI with mu_0
        return {'v_A': v_A, 'v_A_km_s': v_A / 1000,
                'equation': r'v_A = \frac{B}{\sqrt{4\pi\rho}}'}
    
    def eq_41_turbulent_cascade(self, v_l: float, l: float) -> Dict[str, Any]:
        """Kolmogorov cascade rate"""
        epsilon = v_l**3 / l
        return {'epsilon': epsilon, 'unit': 'm^2/s^3',
                'equation': r'\epsilon = \frac{v_l^3}{l} = const'}


class StellarEvolutionCalculator:
    """Eqs 42-44: MS lifetime, mass-luminosity, convection"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_42_ms_lifetime(self, M: float, L: float) -> Dict[str, Any]:
        """Nuclear fuel consumption time"""
        tau_MS = 0.007 * M * self.C.c**2 / L
        return {'tau_MS': tau_MS, 'tau_MS_Gyr': tau_MS / self.C.Gyr,
                'equation': r'\tau_{MS} = \frac{0.007 M c^2}{L}'}
    
    def eq_43_mass_luminosity(self, M: float) -> Dict[str, Any]:
        """Empirical L ~ M^3.5"""
        L = self.C.L_sun * (M / self.C.M_sun)**3.5
        return {'L': L, 'L_Lsun': L / self.C.L_sun,
                'equation': r'L \propto M^{3.5}'}
    
    def eq_44_convective_turnover(self, H_p: float, g: float, 
                                   delta_T: float, T: float,
                                   alpha: float = 2) -> Dict[str, Any]:
        """Mixing length theory"""
        v_conv = (alpha**2 * g * delta_T * H_p / (4 * T))**(1/3)
        t_conv = H_p / v_conv
        return {'v_conv': v_conv, 't_conv': t_conv,
                'equation': r't_{conv} = \frac{H_p}{v_{conv}}'}


class BBNCalculator:
    """Eqs 45-46: Big Bang nucleosynthesis"""
    
    def __init__(self):
        self.C = PhysicsConstants()
        self.eta = 6e-10  # Baryon-to-photon ratio
    
    def eq_45_baryon_photon_ratio(self, n_b: float, 
                                   n_gamma: float = 410e6) -> Dict[str, Any]:
        """eta parameter"""
        eta = n_b / n_gamma
        return {'eta': eta, 'n_gamma': n_gamma,
                'equation': r'\eta = \frac{n_b}{n_\gamma} \approx 6\times 10^{-10}'}
    
    def eq_46_deuterium_bottleneck(self, T: float, rho_rad: float) -> Dict[str, Any]:
        """Onset of nucleosynthesis"""
        t_D = math.sqrt(3 / (32 * math.pi * self.C.G * rho_rad))
        return {'t_D': t_D, 't_D_s': t_D, 'T_MeV': T,
                'equation': r't_D = \sqrt{\frac{3}{32\pi G \rho_{rad}}} \approx 180 s'}


class FriedmannCosmologyCalculator:
    """Eqs 47-49: Expansion, acceleration, density"""
    
    def __init__(self):
        self.C = PhysicsConstants()
        self.Omega_m = 0.3
        self.Omega_Lambda = 0.7
    
    def eq_47_first_friedmann(self, a: float, rho: float, k: float = 0,
                               Lambda: float = 1e-52) -> Dict[str, Any]:
        """Expansion rate"""
        H_sq = 8 * math.pi * self.C.G * rho / 3 - k * self.C.c**2 / a**2 + Lambda * self.C.c**2 / 3
        H = math.sqrt(H_sq) if H_sq > 0 else 0
        return {'H': H, 'H_km_s_Mpc': H * self.C.Mpc / 1000,
                'equation': r'\left(\frac{\dot{a}}{a}\right)^2 = \frac{8\pi G}{3}\rho - \frac{kc^2}{a^2} + \frac{\Lambda c^2}{3}'}
    
    def eq_48_second_friedmann(self, rho: float, P: float,
                                Lambda: float = 1e-52) -> Dict[str, Any]:
        """Acceleration equation"""
        a_ddot_a = -4 * math.pi * self.C.G / 3 * (rho + 3 * P / self.C.c**2) + Lambda * self.C.c**2 / 3
        return {'a_ddot_a': a_ddot_a, 'accelerating': a_ddot_a > 0,
                'equation': r'\frac{\ddot{a}}{a} = -\frac{4\pi G}{3}(\rho + \frac{3P}{c^2}) + \frac{\Lambda c^2}{3}'}
    
    def eq_49_density_parameter(self, rho: float, H: float) -> Dict[str, Any]:
        """Omega parameter"""
        rho_c = 3 * H**2 / (8 * math.pi * self.C.G)
        Omega = rho / rho_c
        return {'Omega': Omega, 'rho_c': rho_c,
                'equation': r'\Omega = \frac{8\pi G \rho}{3H^2}'}


class InflationCalculator:
    """Eqs 50-52: Slow-roll, perturbations, e-folds"""
    
    def __init__(self):
        self.C = PhysicsConstants()
        self.M_pl = math.sqrt(self.C.hbar * self.C.c / self.C.G)
    
    def eq_50_slow_roll_epsilon(self, V: float, V_prime: float) -> Dict[str, Any]:
        """First slow-roll parameter"""
        epsilon = 0.5 * (V_prime / V)**2 * self.M_pl**2
        return {'epsilon': epsilon, 'slow_roll': epsilon < 1,
                'equation': r'\epsilon = \frac{1}{2}\left(\frac{V\'}{V}\right)^2 M_{pl}^2'}
    
    def eq_51_scalar_power(self, H: float, epsilon: float) -> Dict[str, Any]:
        """Primordial scalar amplitude"""
        P_R = H**2 / (8 * math.pi**2 * epsilon * self.M_pl**2)
        return {'P_R': P_R, 'Delta_R_sq': P_R,
                'equation': r'\Delta_R^2 = \frac{H^2}{8\pi^2 \epsilon M_{pl}^2}'}
    
    def eq_52_efolds(self, phi_end: float, phi: float, 
                      epsilon: float) -> Dict[str, Any]:
        """Number of e-folds"""
        # Simplified integral
        N = abs(phi - phi_end) / math.sqrt(2 * epsilon) / self.M_pl
        return {'N': N, 'sufficient': N > 50,
                'equation': r'N = \int \frac{d\phi}{\sqrt{2\epsilon}}'}


class PrimordialGWCalculator:
    """Eqs 53-54: Tensor power spectrum, stochastic background"""
    
    def __init__(self):
        self.C = PhysicsConstants()
        self.M_pl = math.sqrt(self.C.hbar * self.C.c / self.C.G)
    
    def eq_53_tensor_power(self, H: float) -> Dict[str, Any]:
        """Inflationary GW amplitude"""
        P_T = 2 * H**2 / (math.pi**2 * self.M_pl**2)
        return {'P_T': P_T, 'equation': r'P_T(k) = \frac{2H^2}{\pi^2 M_{pl}^2}'}
    
    def eq_54_stochastic_background(self, f: float, Omega_GW: float = 1e-8) -> Dict[str, Any]:
        """GW energy density spectrum"""
        h_c = math.sqrt(3 * self.C.H_0**2 / (4 * math.pi**2 * f**2) * Omega_GW)
        return {'h_c': h_c, 'f': f, 'Omega_GW': Omega_GW,
                'equation': r'\Omega_{GW}(f) = \frac{2\pi^2}{3H_0^2} f^2 h_c^2'}


class BinaryBHMergerCalculator:
    """Eqs 55-57: Inspiral, merger time, ringdown damping"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_55_frequency_evolution(self, f: float, M_chirp: float) -> Dict[str, Any]:
        """Post-Newtonian chirp"""
        f_dot = (96/5) * math.pi**(8/3) * (self.C.G * M_chirp / self.C.c**3)**(5/3) * f**(11/3)
        return {'f_dot': f_dot, 'f': f,
                'equation': r'\dot{f} = \frac{96}{5}\pi^{8/3}\left(\frac{G\mathcal{M}}{c^3}\right)^{5/3} f^{11/3}'}
    
    def eq_56_merger_time(self, f_i: float, M_chirp: float) -> Dict[str, Any]:
        """Time to coalescence"""
        t_merge = 5/256 * self.C.c**5 / self.C.G**(5/3) / (math.pi * f_i)**(8/3) / M_chirp**(5/3)
        return {'t_merge': t_merge, 't_merge_s': t_merge,
                'equation': r't_{merge} = \frac{5}{256}\frac{c^5}{G^{5/3}}\frac{1}{(\pi f_i)^{8/3}\mathcal{M}^{5/3}}'}
    
    def eq_57_ringdown_damping(self, M_f: float, Q: float = 10) -> Dict[str, Any]:
        """QNM damping time"""
        # Approximate for l=2, m=2
        f_QNM = 0.37 * self.C.c**3 / (2 * math.pi * self.C.G * M_f)
        tau = 2 * M_f * self.C.G / self.C.c**3 * Q
        return {'tau': tau, 'tau_ms': tau * 1000, 'Q': Q,
                'equation': r'\tau = \frac{2M_f G Q}{c^3}'}


class SupernovaCalculator:
    """Eqs 58-60: Light curves, ejecta velocity, yields"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_58_arnett_peak(self, M_Ni: float, t_d: float = 20 * 86400) -> Dict[str, Any]:
        """Radioactive decay powered peak"""
        epsilon_Ni = 3.9e10 * 1e-7  # Convert erg/g to J/kg
        L_peak = M_Ni * epsilon_Ni / t_d
        return {'L_peak': L_peak, 'L_peak_Lsun': L_peak / self.C.L_sun,
                'equation': r'L_{peak} = \frac{M_{Ni} \epsilon_{Ni}}{t_d}'}
    
    def eq_59_ejecta_velocity(self, E_kin: float, M_ej: float) -> Dict[str, Any]:
        """Homologous expansion"""
        v_ej = math.sqrt(2 * E_kin / M_ej)
        return {'v_ej': v_ej, 'v_ej_km_s': v_ej / 1000,
                'equation': r'v_{ej} = \sqrt{\frac{2E_{kin}}{M_{ej}}}'}
    
    def eq_60_nucleosynthesis_yield(self, X_i: float, rho: float, 
                                     dt: float) -> Dict[str, Any]:
        """Explosive burning yield"""
        Y_i = X_i * rho * dt  # Simplified
        return {'Y_i': Y_i, 'equation': r'Y_i = \int \rho X_i dt'}


class PlanetaryNebulaCalculator:
    """Eqs 61-63: PN expansion, ionization, mass loss"""
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def eq_61_expansion_radius(self, v_exp: float, t: float) -> Dict[str, Any]:
        """Constant velocity expansion"""
        R = v_exp * t
        return {'R': R, 'R_pc': R / self.C.pc,
                'equation': r'R(t) = v_{exp} t'}
    
    def eq_62_ionization_front(self, N_dot_UV: float, R: float, n: float,
                                alpha_B: float = 2.6e-19) -> Dict[str, Any]:
        """Stromgren expansion velocity"""
        v_IF = N_dot_UV / (4 * math.pi * R**2 * n) - alpha_B * n * R / 3
        return {'v_IF': v_IF, 'equation': r'v_{IF} = \frac{\dot{N}_{UV}}{4\pi R^2 n} - \frac{\alpha_B n R}{3}'}
    
    def eq_63_reimers_mass_loss(self, L: float, R: float, M: float) -> Dict[str, Any]:
        """AGB wind formula"""
        # In solar units
        L_sun, R_sun, M_sun = L / self.C.L_sun, R / self.C.R_sun, M / self.C.M_sun
        M_dot = 4e-13 * L_sun * R_sun / M_sun * self.C.M_sun / self.C.yr
        return {'M_dot': M_dot, 'M_dot_Msun_yr': M_dot / self.C.M_sun * self.C.yr,
                'equation': r'\dot{M} = 4\times 10^{-13} \frac{LR}{M} M_\odot/yr'}


# =============================================================================
# UQFF EXTENSION CALCULATORS
# =============================================================================

class UQFFBuoyancyCalculator:
    """
    Universal Quantum Field Superconductive Framework - Buoyancy Forces
    F_U_Bi_i equations for stabilization/instability
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def compute_F_UBii(self, E_cm: float, Q_wave: float, g_local: float,
                        negative: bool = True) -> Dict[str, Any]:
        """
        UQFF Buoyancy Force
        
        F_UBii = (+/-) F_rel * (E_cm / E_LEP) * Q_wave * g(r,t)
        
        Parameters:
            E_cm: Center of mass energy (J)
            Q_wave: Wave/resonance factor (~10^12 for THz)
            g_local: Local gravity acceleration (m/s^2)
            negative: If True, stabilizing; if False, destabilizing
            
        Returns:
            Dictionary with buoyancy force
        """
        sign = -1 if negative else 1
        F_UBii = sign * self.C.F_rel * (E_cm / self.C.E_LEP_1998) * Q_wave * g_local
        
        return {
            'equation': 'F_UBii = (+/-) F_rel * (E_cm/E_LEP) * Q_wave * g',
            'latex': r'F_{UBii} = \pm F_{rel} \times \frac{E_{cm}}{E_{LEP}} \times Q_{wave} \times g(r,t)',
            'F_UBii': F_UBii,
            'F_UBii_scientific': f'{F_UBii:.3e} N',
            'E_ratio': E_cm / self.C.E_LEP_1998,
            'Q_wave': Q_wave,
            'mode': 'stabilizing' if negative else 'destabilizing'
        }
    
    def buoyancy_yield_probability(self, F_UBii: float, r: float,
                                    E_th: float) -> Dict[str, Any]:
        """
        Boltzmann-like probability for cluster stability
        
        P_alpha = 1 - exp(F_UBii * r / E_th)
        
        Parameters:
            F_UBii: Buoyancy force (N)
            r: Characteristic radius (m)
            E_th: Threshold energy (J)
            
        Returns:
            Dictionary with probability
        """
        arg = F_UBii * r / E_th
        P_alpha = 1 - math.exp(arg) if arg < 100 else 1.0
        
        return {
            'equation': 'P_alpha = 1 - exp(F_UBii * r / E_th)',
            'P_alpha': P_alpha,
            'stable': P_alpha > 0.85
        }


class UQFFMagnetismCalculator:
    """
    Universal Magnetism (Um) equations with vacuum-aether modulation
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def compute_Um(self, t: float, r: float, n: int = 0,
                    mu_j: float = 3.38e20, P_SCm: float = 1.0,
                    E_react: float = 1e46) -> Dict[str, Any]:
        """
        Universal Magnetism Equation
        
        Um(t,r,n) = sum_j [mu_j/r * (1 - exp(-gamma*t*cos(pi*t*n))) * phi_j] *
                    P_SCm * E_react(t) * (1 + 10^13 * f_Heav) * (1 + f_quasi)
        
        Parameters:
            t: Time (days)
            r: Distance (m)
            n: Index parameter
            mu_j: Magnetic moment (T*pm^3), default 3.38e20
            P_SCm: Superconductive probability (default 1.0)
            E_react: Reaction energy (default 10^46)
            
        Returns:
            Dictionary with Um value
        """
        gamma = self.C.gamma_UQFF
        omega_c = self.C.omega_c
        f_Heav = 0.01
        f_quasi = 0.01
        phi_j = 1.0
        
        # Time-dependent mu_j oscillation
        t_s = t * 86400  # Convert days to seconds
        mu_j_t = (1000 + 0.4 * math.sin(omega_c * t_s)) * mu_j
        
        # Exponential term
        if n == 0:
            exp_term = 1 - math.exp(-gamma * t)
        else:
            exp_term = 1 - math.exp(-gamma * t * math.cos(math.pi * t * n))
        
        # E_react decay
        E_react_t = E_react * math.exp(-0.0005 * t)
        
        # Full Um
        Um = (mu_j_t / r) * exp_term * phi_j * P_SCm * E_react_t * \
             (1 + 1e13 * f_Heav) * (1 + f_quasi)
        
        return {
            'equation': 'Um(t,r,n) = [mu_j/r * (1-exp(...)) * phi] * P_SCm * E_react * factors',
            'latex': r'U_m(t,r,n) = \sum_j \frac{\mu_j}{r}(1-e^{-\gamma t \cos(\pi t n)})\phi_j \cdot P_{SCm} \cdot E_{react}',
            'Um': Um,
            'Um_scientific': f'{Um:.3e}',
            't_days': t,
            'r_m': r,
            'n': n
        }
    
    def compute_E_field(self, Um: float, r: float) -> Dict[str, Any]:
        """
        Electric field from Um
        
        E = Um * rho_vac,[UA] / r
        """
        E = Um * self.C.rho_vac_UA / r
        return {
            'E': E,
            'E_V_m': E,
            'equation': r'E = \frac{U_m \cdot \rho_{vac,[UA]}}{r}'
        }
    
    def compute_neutron_rate(self, Um: float, k_eta: float = 1e-100,
                              SSq: float = 1.0, n: int = 0,
                              t: float = 1.0) -> Dict[str, Any]:
        """
        Neutron production rate (LENR-style)
        
        eta = k_eta * exp(-[SSq]^n / 26) * exp(-pi - t) * Um / rho_vac
        """
        exp1 = math.exp(-SSq**n / 26) if n < 10 else 0
        exp2 = math.exp(-math.pi - t)
        eta = k_eta * exp1 * exp2 * Um / self.C.rho_vac_UA
        
        return {
            'eta': eta,
            'eta_per_cm2_s': eta * 1e4,
            'k_eta': k_eta,
            'equation': r'\eta = k_\eta e^{-[SSq]^n/26} e^{-\pi-t} \frac{U_m}{\rho_{vac}}'}


class UQFFGyroCalculator:
    """
    Gyroscopic/Inertial integration with torque nullification
    """
    
    def __init__(self):
        self.C = PhysicsConstants()
    
    def compute_torque(self, I: float, omega: float, alpha: float) -> Dict[str, Any]:
        """
        Gyroscopic torque
        
        tau = I * omega * alpha
        """
        tau = I * omega * alpha
        return {
            'tau': tau,
            'tau_Nm': tau,
            'equation': r'\tau = I \omega \alpha'
        }
    
    def torque_nullification(self, tau: float, r: float,
                              theta: float = math.pi/2) -> Dict[str, Any]:
        """
        UQFF torque nullification via buoyancy
        
        F_UBii_required = -tau / (r * sin(theta))
        """
        F_req = -tau / (r * math.sin(theta))
        return {
            'F_UBii_required': F_req,
            'nullified': abs(F_req) < 1e10,  # Within UQFF range
            'equation': r'F_{UBii} = -\frac{\tau}{r\sin\theta}'
        }


# =============================================================================
# MODULE EXPORTS
# =============================================================================

__all__ = [
    # Constants
    'PhysicsConstants',
    
    # Protostellar Jets (1-4)
    'ProtostellarJetCalculator',
    
    # Galaxy Mergers (5-7)
    'GalaxyMergerSFRCalculator',
    
    # Black Hole Growth (8-9)
    'BlackHoleGrowthCalculator',
    
    # SNR (10-11)
    'SupernovaRemnantCalculator',
    
    # GW Mergers (12-13)
    'GravitationalWaveMergerCalculator',
    
    # Quasar Jets (14-15)
    'QuasarJetCalculator',
    
    # Neutron Stars (16-18)
    'NeutronStarCalculator',
    
    # GRBs (19-20)
    'GammaRayBurstCalculator',
    
    # CMB (21-22)
    'CMBCalculator',
    
    # AGN Feedback (23-25)
    'AGNFeedbackCalculator',
    
    # Exoplanets (26-28)
    'ExoplanetCalculator',
    
    # Dark Matter (29-31)
    'DarkMatterHaloCalculator',
    
    # Additional (32-100)
    'GalaxyClusterCalculator',
    'CosmicVoidCalculator',
    'ReionizationCalculator',
    'ISMTurbulenceCalculator',
    'StellarEvolutionCalculator',
    'BBNCalculator',
    'FriedmannCosmologyCalculator',
    'InflationCalculator',
    'PrimordialGWCalculator',
    'BinaryBHMergerCalculator',
    'SupernovaCalculator',
    'PlanetaryNebulaCalculator',
    
    # UQFF Extensions
    'UQFFBuoyancyCalculator',
    'UQFFMagnetismCalculator',
    'UQFFGyroCalculator',
]


# =============================================================================
# TEST/DEMO
# =============================================================================
if __name__ == '__main__':
    print("=" * 70)
    print("GROK 100+ EQUATIONS MODULE - Test Suite")
    print("=" * 70)
    
    C = PhysicsConstants()
    
    # Test Eq 1: Angular momentum transport
    jet_calc = ProtostellarJetCalculator()
    result = jet_calc.eq_01_angular_momentum_transport(
        M_dot=1e-7 * C.M_sun / C.yr,  # 10^-7 Msun/yr
        r=10 * C.AU,
        M_star=0.5 * C.M_sun
    )
    print(f"\nEq 1: Angular Momentum Transport")
    print(f"  Omega = {result['Omega']:.3e} rad/s")
    print(f"  dL/dt = {result['L_dot']:.3e} kg*m^2/s^2")
    
    # Test Eq 10: Sedov-Taylor
    snr_calc = SupernovaRemnantCalculator()
    result = snr_calc.eq_10_sedov_taylor(
        E=1e44,  # 10^51 erg
        rho=1e-24 * 1000,  # 10^-24 g/cm^3 to kg/m^3
        t=400 * C.yr
    )
    print(f"\nEq 10: Sedov-Taylor SNR Expansion")
    print(f"  R = {result['R_pc']:.2f} pc")
    print(f"  v = {result['v_km_s']:.0f} km/s")
    
    # Test Eq 12: Chirp mass
    gw_calc = GravitationalWaveMergerCalculator()
    result = gw_calc.eq_12_chirp_mass(
        m1=30 * C.M_sun,
        m2=25 * C.M_sun
    )
    print(f"\nEq 12: Chirp Mass")
    print(f"  M_chirp = {result['M_chirp_Msun']:.1f} Msun")
    
    # Test UQFF Buoyancy
    uqff_calc = UQFFBuoyancyCalculator()
    result = uqff_calc.compute_F_UBii(
        E_cm=1e9 * C.e,  # 1 GeV
        Q_wave=1e12,
        g_local=10,  # m/s^2
        negative=True
    )
    print(f"\nUQFF Buoyancy Force")
    print(f"  F_UBii = {result['F_UBii_scientific']}")
    print(f"  Mode: {result['mode']}")
    
    print("\n" + "=" * 70)
    print("All 100+ equations available. Import classes as needed.")
    print("=" * 70)
