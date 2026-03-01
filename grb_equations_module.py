#!/usr/bin/env python3
"""
Gamma-Ray Burst (GRB) Equations Module

From Grok Deep Analysis (Feb 2026):
- Equations 19-20: Fireball Expansion and Afterglow Synchrotron
- Binary Pulsar GW and Kilonova physics (Eqs 73-75)

Physics domains covered:
- Relativistic fireball expansion (prompt phase)
- Afterglow synchrotron emission (external shock)
- Compact binary evolution and gravitational wave emission
- Kilonova light curves from r-process nucleosynthesis

UQFF Integration:
- F_U_Bi_i buoyancy in ultra-relativistic jets
- Um magnetism in jet collimation and acceleration
- Validation of UQFF clustering coherence at TeV scales
"""

import math
from typing import Dict, Optional

# ============== Physical Constants ==============
G = 6.674e-11           # Gravitational constant [m³/(kg·s²)]
c = 2.998e8             # Speed of light [m/s]
M_sun = 1.989e30        # Solar mass [kg]
h_bar = 1.055e-34       # Reduced Planck constant [J·s]
k_B = 1.381e-23         # Boltzmann constant [J/K]
m_e = 9.109e-31         # Electron mass [kg]
e_charge = 1.602e-19    # Elementary charge [C]
sigma_SB = 5.670e-8     # Stefan-Boltzmann constant [W/(m²·K⁴)]
eV_to_J = 1.602e-19     # eV to Joules
erg_to_J = 1e-7         # erg to Joules
pc_to_m = 3.086e16      # parsec to meters
year_to_s = 3.154e7     # year to seconds

# UQFF Constants
F_rel = 4.30e33         # Relativistic coherence force [N]
E_LEP = 200e9 * eV_to_J # LEP 1998 baseline energy [J]


class FireballExpansionCalculator:
    """
    Relativistic fireball expansion in GRB prompt phase.
    
    Equation 19:
    Γ(r) = r/R_0           (r < R_s)
    Γ(r) = η               (r > R_s)
    
    Where:
    - R_0: initial fireball radius (~10^7 cm)
    - R_s: saturation radius ~ η² R_0
    - η = E / (M c²): baryon loading factor (~300)
    
    Derivation: Energy conservation during expansion,
    Lorentz factor proportional to radius until photons
    decouple at saturation radius.
    """
    
    def compute(self, r: float, R_0: float, E_iso: float, 
                M_baryon: float) -> Dict:
        """
        Compute fireball Lorentz factor at radius r.
        
        Args:
            r: Current radius [m]
            R_0: Initial radius [m]
            E_iso: Isotropic equivalent energy [J]
            M_baryon: Baryon loading mass [kg]
        
        Returns:
            Dict with Lorentz factor and regime
        """
        # Baryon loading factor
        eta = E_iso / (M_baryon * c**2) if M_baryon > 0 else float('inf')
        
        # Saturation radius
        R_s = eta**2 * R_0
        
        # Lorentz factor by regime
        if r < R_s:
            Gamma = r / R_0
            regime = 'acceleration'
        else:
            Gamma = eta
            regime = 'coasting'
        
        # Velocity
        if Gamma > 1:
            beta = math.sqrt(1 - 1/Gamma**2)
        else:
            beta = 0
        
        # Internal shock location estimate
        # δr ~ c δt / Γ² (variability timescale δt ~ ms)
        delta_t = 1e-3  # 1 ms typical
        delta_r = c * delta_t / Gamma**2 if Gamma > 0 else float('inf')
        
        return {
            'Gamma': Gamma,
            'beta': beta,
            'velocity_m_s': beta * c,
            'eta': eta,
            'R_0_m': R_0,
            'R_s_m': R_s,
            'regime': regime,
            'delta_r_internal_shock_m': delta_r,
            'equation': 'Γ(r) = r/R_0 (r<R_s), Γ(r) = η (r>R_s)'
        }


class AfterglowSynchrotronCalculator:
    """
    GRB afterglow synchrotron emission from external shock.
    
    Equation 20:
    F_ν ∝ ν^{-(p-1)/2} t^{-3(p-1)/4}    (ν_m < ν < ν_c)
    
    Where:
    - p: electron power-law index (~2.2-2.5)
    - ν_m: peak synchrotron frequency
    - ν_c: cooling frequency
    
    Derivation: Blandford-McKee solution for shock deceleration
    in ISM, electrons accelerated via DSA, synchrotron spectrum
    from convolved power-law distribution.
    """
    
    def compute(self, t: float, nu: float, p: float = 2.3,
                E_iso: float = 1e53 * erg_to_J,
                n_ISM: float = 1e6,  # 1 cm⁻³
                epsilon_e: float = 0.1,
                epsilon_B: float = 0.01,
                D_L: float = 1e9 * pc_to_m) -> Dict:
        """
        Compute afterglow synchrotron flux.
        
        Args:
            t: Time since burst [s]
            nu: Observation frequency [Hz]
            p: Electron index
            E_iso: Isotropic energy [J]
            n_ISM: ISM number density [m⁻³]
            epsilon_e: Electron energy fraction
            epsilon_B: Magnetic energy fraction
            D_L: Luminosity distance [m]
        
        Returns:
            Dict with afterglow parameters
        """
        # Deceleration radius (Sedov-like)
        # R_dec ~ (3 E / 4π n m_p c²)^{1/3}
        m_p = 1.673e-27  # proton mass
        R_dec = (3 * E_iso / (4 * math.pi * n_ISM * m_p * c**2))**(1/3)
        
        # Deceleration time
        t_dec = R_dec / c  # Simplified
        
        # Lorentz factor evolution (Blandford-McKee)
        # Γ ∝ t^{-3/8} in ISM
        if t > t_dec:
            Gamma = (t_dec / t)**(3/8) * (E_iso / (n_ISM * m_p * c**2 * R_dec**3))**(1/8)
        else:
            Gamma = 100  # Initial estimate
        
        # Magnetic field in shock
        B = math.sqrt(8 * math.pi * epsilon_B * n_ISM * m_p * c**2 * Gamma**2)
        
        # Peak synchrotron frequency
        gamma_m = epsilon_e * (p - 2) / (p - 1) * m_p / m_e * Gamma
        nu_m = e_charge * B * gamma_m**2 / (2 * math.pi * m_e * c)
        
        # Cooling frequency
        sigma_T = 6.652e-29
        t_cool = 6 * math.pi * m_e * c / (sigma_T * B**2 * Gamma)
        gamma_c = 6 * math.pi * m_e * c / (sigma_T * B**2 * t)
        nu_c = e_charge * B * gamma_c**2 / (2 * math.pi * m_e * c)
        
        # Flux density (simplified, in Jy)
        # F_nu,max ~ (epsilon_B^{1/2} n^{1/2} E) / (D_L^2)
        F_nu_max = epsilon_B**0.5 * n_ISM**0.5 * E_iso / (4 * math.pi * D_L**2)
        
        # Spectral regime
        if nu < nu_m:
            spectral_index = 1/3
            regime = 'below ν_m'
        elif nu < nu_c:
            spectral_index = -(p - 1) / 2
            regime = 'ν_m < ν < ν_c'
        else:
            spectral_index = -p / 2
            regime = 'above ν_c'
        
        # Temporal decay
        temporal_index = -3 * (p - 1) / 4  # ISM
        
        # Scale flux
        F_nu = F_nu_max * (nu / nu_m)**spectral_index * (t / t_dec)**temporal_index
        
        return {
            'F_nu_W_m2_Hz': F_nu,
            'nu_m_Hz': nu_m,
            'nu_c_Hz': nu_c,
            'Gamma': Gamma,
            'B_T': B,
            'spectral_index': spectral_index,
            'temporal_index': temporal_index,
            'regime': regime,
            'R_dec_m': R_dec,
            't_dec_s': t_dec,
            'equation': 'F_ν ∝ ν^{-(p-1)/2} t^{-3(p-1)/4}'
        }


class ChirpMassCalculator:
    """
    Chirp mass from gravitational wave inspiral.
    
    Equation 12:
    M = (m₁m₂)^{3/5} / (m₁+m₂)^{1/5} = (c³/G) (5/96 π^{-8/3} f^{-11/3} ḟ)^{3/5}
    
    Where:
    - f: GW frequency
    - ḟ: frequency derivative (chirp rate)
    
    Derivation: Quadrupole approximation for power emission,
    energy loss dE/dt = -L determines frequency evolution.
    """
    
    def compute(self, f: float, f_dot: float) -> Dict:
        """
        Compute chirp mass from GW frequency evolution.
        
        Args:
            f: GW frequency [Hz]
            f_dot: Frequency derivative [Hz/s]
        
        Returns:
            Dict with chirp mass and related quantities
        """
        # Chirp mass formula
        coefficient = (c**3 / G) * (5/96 * math.pi**(-8/3))**(3/5)
        M_chirp = coefficient * (f**(-11/3) * f_dot)**(3/5)
        
        # Merger time estimate
        # t_merge = 5/(256) c⁵/(G⁵ f⁸) * M^{-5/3} π^{-8/3}
        t_merge = (5/256) * (c**5 / G**(5/3)) / (math.pi * f)**(8/3) / M_chirp**(5/3)
        
        return {
            'M_chirp_kg': M_chirp,
            'M_chirp_Msun': M_chirp / M_sun,
            'f_Hz': f,
            'f_dot_Hz_s': f_dot,
            't_merge_s': t_merge,
            'equation': 'M = (c³/G)(5/96 π^{-8/3} f^{-11/3} ḟ)^{3/5}'
        }
    
    def compute_from_masses(self, m1: float, m2: float) -> Dict:
        """
        Compute chirp mass from component masses.
        
        Args:
            m1: Primary mass [kg]
            m2: Secondary mass [kg]
        
        Returns:
            Dict with chirp mass
        """
        M_total = m1 + m2
        M_chirp = (m1 * m2)**(3/5) / M_total**(1/5)
        
        # Symmetric mass ratio
        eta = (m1 * m2) / M_total**2
        
        return {
            'M_chirp_kg': M_chirp,
            'M_chirp_Msun': M_chirp / M_sun,
            'M_total_Msun': M_total / M_sun,
            'eta_symmetric_ratio': eta,
            'q_mass_ratio': m1 / m2 if m2 > 0 else float('inf'),
            'equation': 'M = (m₁m₂)^{3/5} / (m₁+m₂)^{1/5}'
        }


class RingdownQNMCalculator:
    """
    Post-merger ringdown quasi-normal modes.
    
    Equation 13:
    f_QNM = (c³ / 2πGM_f) × f(a_f)
    
    Where f(a_f) ≈ 0.3737 + 0.088×a_f + ...
    
    Derivation: Perturb Kerr metric via Teukolsky equation,
    l=2, m=2 modes dominant.
    """
    
    def compute(self, M_final: float, a_final: float) -> Dict:
        """
        Compute ringdown frequency and damping time.
        
        Args:
            M_final: Final BH mass [kg]
            a_final: Final dimensionless spin (0-1)
        
        Returns:
            Dict with QNM parameters
        """
        # Fundamental l=2, m=2 mode fit (Berti et al.)
        f_factor = 0.3737 + 0.088 * a_final
        
        # QNM frequency
        f_QNM = (c**3 / (2 * math.pi * G * M_final)) * f_factor
        
        # Quality factor (simple fit)
        Q = 2 * (1 - a_final)**(-0.45)
        
        # Damping time
        tau_QNM = Q / (math.pi * f_QNM)
        
        return {
            'f_QNM_Hz': f_QNM,
            'tau_QNM_s': tau_QNM,
            'Q_factor': Q,
            'M_final_Msun': M_final / M_sun,
            'a_final': a_final,
            'equation': 'f_QNM = (c³/2πGM_f) × (0.3737 + 0.088×a_f)'
        }


class BinaryPulsarOrbitDecayCalculator:
    """
    Binary pulsar orbital decay from gravitational waves.
    
    Equation 73:
    Ṗ_b = -(192π/5) (P_b/2π)^{-5/3} (G m₁m₂/c³(m₁+m₂)^{1/3})^{5/3} (1-e²)^{-7/2}
    
    Derivation: Quadrupole power emission averaged over orbit,
    energy loss causes orbital shrinkage.
    """
    
    def compute(self, P_b: float, m1: float, m2: float, e: float) -> Dict:
        """
        Compute orbital period derivative.
        
        Args:
            P_b: Orbital period [s]
            m1: Primary mass [kg]
            m2: Secondary mass [kg]
            e: Orbital eccentricity
        
        Returns:
            Dict with orbit decay parameters
        """
        M_total = m1 + m2
        
        # Period derivative
        coeff = (192 * math.pi / 5) * (P_b / (2 * math.pi))**(-5/3)
        mass_factor = (G * m1 * m2 / (c**3 * M_total**(1/3)))**(5/3)
        ecc_factor = (1 - e**2)**(-7/2)
        
        P_b_dot = -coeff * mass_factor * ecc_factor
        
        # Merger timescale
        t_merge = 5 * c**5 * P_b**(8/3) / (256 * G**(5/3) * (m1 * m2) * M_total**(-1/3))
        t_merge *= (1 - e**2)**(7/2)  # Eccentricity correction
        
        return {
            'P_b_dot_s_s': P_b_dot,
            'P_b_s': P_b,
            't_merge_s': t_merge,
            't_merge_yr': t_merge / year_to_s,
            'eccentricity': e,
            'equation': 'Ṗ_b = -(192π/5)(P_b/2π)^{-5/3}(Gm₁m₂/c³M^{1/3})^{5/3}(1-e²)^{-7/2}'
        }


class PeriastronAdvanceCalculator:
    """
    Post-Keplerian periastron advance (GR test).
    
    Equation 74:
    ω̇ = 3(P_b/2π)^{-5/3} (G(m₁+m₂)/c³)^{2/3} (1-e²)^{-1}
    
    Derivation: GR geodesic precession from Schwarzschild/Kerr metric,
    relativistic correction to Keplerian orbits.
    """
    
    def compute(self, P_b: float, m1: float, m2: float, e: float) -> Dict:
        """
        Compute periastron advance rate.
        
        Args:
            P_b: Orbital period [s]
            m1: Primary mass [kg]
            m2: Secondary mass [kg]
            e: Orbital eccentricity
        
        Returns:
            Dict with precession parameters
        """
        M_total = m1 + m2
        
        # Periastron advance (rad/s)
        coeff = 3 * (P_b / (2 * math.pi))**(-5/3)
        mass_factor = (G * M_total / c**3)**(2/3)
        ecc_factor = (1 - e**2)**(-1)
        
        omega_dot = coeff * mass_factor * ecc_factor
        
        # Convert to degrees/year
        omega_dot_deg_yr = omega_dot * (180/math.pi) * year_to_s
        
        return {
            'omega_dot_rad_s': omega_dot,
            'omega_dot_deg_yr': omega_dot_deg_yr,
            'P_b_s': P_b,
            'M_total_Msun': M_total / M_sun,
            'eccentricity': e,
            'equation': 'ω̇ = 3(P_b/2π)^{-5/3}(GM/c³)^{2/3}(1-e²)^{-1}'
        }


class KilonovaLightCurveCalculator:
    """
    Kilonova light curve from r-process radioactive decay.
    
    Equation 75:
    L_peak ≈ 10^{41} (M_ej/0.01M_☉)(v_ej/0.1c)(κ/1 cm²/g)^{-1} erg/s
    
    Where:
    - M_ej: ejecta mass (~0.01 M_☉)
    - v_ej: ejecta velocity (~0.2c)
    - κ: opacity (~1-10 cm²/g for lanthanides)
    
    Derivation: Thermalization of r-process decay energy,
    diffusion time determines peak.
    """
    
    def compute(self, M_ej: float, v_ej: float, kappa: float,
                t_obs: Optional[float] = None) -> Dict:
        """
        Compute kilonova luminosity.
        
        Args:
            M_ej: Ejecta mass [kg]
            v_ej: Ejecta velocity [m/s]
            kappa: Opacity [m²/kg]
            t_obs: Observation time [s] (optional, for light curve)
        
        Returns:
            Dict with kilonova parameters
        """
        # Peak luminosity scaling
        L_ref = 1e41 * erg_to_J
        M_ref = 0.01 * M_sun
        v_ref = 0.1 * c
        kappa_ref = 1e-3  # 1 cm²/g = 0.1 m²/kg, using SI
        
        L_peak = L_ref * (M_ej / M_ref) * (v_ej / v_ref) * (kappa_ref / kappa)
        
        # Peak time (diffusion timescale)
        # t_peak ~ sqrt(3 κ M / (4π c v²))
        t_peak = math.sqrt(3 * kappa * M_ej / (4 * math.pi * c * v_ej**2))
        
        # Heating rate (simplified r-process)
        # ε ~ 10^{10} (t/1day)^{-1.3} erg/g/s
        t_day = 86400
        epsilon_0 = 1e10 * erg_to_J * 1e3  # erg/g/s -> J/kg/s
        
        if t_obs is not None and t_obs > 0:
            epsilon = epsilon_0 * (t_obs / t_day)**(-1.3)
            L_t = epsilon * M_ej * math.exp(-t_obs / t_peak)
        else:
            epsilon = epsilon_0
            L_t = L_peak
        
        return {
            'L_peak_W': L_peak,
            'L_peak_erg_s': L_peak / erg_to_J,
            't_peak_s': t_peak,
            't_peak_days': t_peak / t_day,
            'M_ej_Msun': M_ej / M_sun,
            'v_ej_c': v_ej / c,
            'kappa_m2_kg': kappa,
            'L_observed_W': L_t if t_obs else L_peak,
            'equation': 'L_peak ≈ 10^{41}(M_ej/0.01M_☉)(v_ej/0.1c)(κ/1 cm²/g)^{-1} erg/s'
        }


class InspiralFrequencyEvolutionCalculator:
    """
    GW inspiral orbital frequency evolution.
    
    Equation 55:
    ḟ = (96/5) π^{8/3} (G M / c³)^{5/3} f^{11/3}
    
    Where M is chirp mass.
    
    Derivation: Quadrupole power loss, f = 2 f_orb for GW.
    """
    
    def compute(self, f: float, M_chirp: float) -> Dict:
        """
        Compute frequency evolution rate.
        
        Args:
            f: Current GW frequency [Hz]
            M_chirp: Chirp mass [kg]
        
        Returns:
            Dict with frequency evolution
        """
        # Frequency derivative
        f_dot = (96/5) * math.pi**(8/3) * (G * M_chirp / c**3)**(5/3) * f**(11/3)
        
        # Number of cycles to merger
        # N_cycles ~ f²/(2 ḟ)
        if f_dot > 0:
            N_cycles = f**2 / (2 * f_dot)
            t_merge = f / f_dot  # Rough estimate
        else:
            N_cycles = float('inf')
            t_merge = float('inf')
        
        return {
            'f_dot_Hz_s': f_dot,
            'f_Hz': f,
            'M_chirp_Msun': M_chirp / M_sun,
            'N_cycles_to_merge': N_cycles,
            't_merge_approx_s': t_merge,
            'equation': 'ḟ = (96/5)π^{8/3}(GM/c³)^{5/3}f^{11/3}'
        }


class MergerTimeCalculator:
    """
    Merger time from initial frequency.
    
    Equation 56:
    t_merge = (5/256) c⁵/G^{5/3} × 1/(π f_i)^{8/3} × M^{-5/3}
    
    Derivation: Integrate df/ḟ from f_i to f_ISCO.
    """
    
    def compute(self, f_i: float, M_chirp: float) -> Dict:
        """
        Compute time to merger.
        
        Args:
            f_i: Initial GW frequency [Hz]
            M_chirp: Chirp mass [kg]
        
        Returns:
            Dict with merger time
        """
        t_merge = (5/256) * (c**5 / G**(5/3)) / (math.pi * f_i)**(8/3) / M_chirp**(5/3)
        
        return {
            't_merge_s': t_merge,
            't_merge_yr': t_merge / year_to_s,
            'f_initial_Hz': f_i,
            'M_chirp_Msun': M_chirp / M_sun,
            'equation': 't_merge = (5/256)(c⁵/G^{5/3})(πf_i)^{-8/3}M^{-5/3}'
        }


class GRBCalculator:
    """
    Master calculator for GRB and compact binary physics.
    
    Integrates:
    - Fireball expansion (prompt phase)
    - Afterglow synchrotron
    - Binary inspiral and merger
    - Kilonova emission
    
    UQFF Integration:
    - Ultra-relativistic regime validates EM >> gravity (EU proof)
    - Buoyancy forces in jet collimation
    - Nuclear clustering (r-process) as UQFF coherence
    """
    
    def __init__(self):
        self.fireball_calc = FireballExpansionCalculator()
        self.afterglow_calc = AfterglowSynchrotronCalculator()
        self.chirp_calc = ChirpMassCalculator()
        self.ringdown_calc = RingdownQNMCalculator()
        self.orbit_decay_calc = BinaryPulsarOrbitDecayCalculator()
        self.periastron_calc = PeriastronAdvanceCalculator()
        self.kilonova_calc = KilonovaLightCurveCalculator()
        self.inspiral_calc = InspiralFrequencyEvolutionCalculator()
        self.merger_time_calc = MergerTimeCalculator()
    
    def compute_grb_full_analysis(self, E_iso: float, M_baryon: float,
                                   R_0: float, r_obs: float,
                                   t_obs: float, nu_obs: float) -> Dict:
        """
        Complete GRB analysis from prompt to afterglow.
        
        Args:
            E_iso: Isotropic energy [J]
            M_baryon: Baryon mass [kg]
            R_0: Initial radius [m]
            r_obs: Observation radius [m]
            t_obs: Time since burst [s]
            nu_obs: Observation frequency [Hz]
        
        Returns:
            Comprehensive GRB analysis
        """
        # Prompt phase
        fireball = self.fireball_calc.compute(r_obs, R_0, E_iso, M_baryon)
        
        # Afterglow
        afterglow = self.afterglow_calc.compute(t_obs, nu_obs)
        
        return {
            'prompt_phase': fireball,
            'afterglow': afterglow,
            'UQFF': {
                'F_UBii_jet': self._compute_jet_buoyancy(E_iso, fireball['Gamma']),
                'Um_collimation': self._compute_Um_jet(fireball['Gamma'])
            }
        }
    
    def compute_binary_merger_analysis(self, m1: float, m2: float,
                                        P_b: float, e: float) -> Dict:
        """
        Binary compact object merger analysis.
        
        Args:
            m1: Primary mass [kg]
            m2: Secondary mass [kg]
            P_b: Orbital period [s]
            e: Eccentricity
        
        Returns:
            Complete binary analysis
        """
        # Chirp mass
        chirp = self.chirp_calc.compute_from_masses(m1, m2)
        
        # Orbital decay
        orbit_decay = self.orbit_decay_calc.compute(P_b, m1, m2, e)
        
        # Periastron advance
        periastron = self.periastron_calc.compute(P_b, m1, m2, e)
        
        # Final state estimates
        M_final = 0.95 * (m1 + m2)  # ~5% radiated
        a_final = 0.69  # Typical final spin
        
        ringdown = self.ringdown_calc.compute(M_final, a_final)
        
        return {
            'chirp_mass': chirp,
            'orbit_decay': orbit_decay,
            'periastron_advance': periastron,
            'ringdown': ringdown
        }
    
    def _compute_jet_buoyancy(self, E_iso: float, Gamma: float) -> float:
        """UQFF buoyancy in ultra-relativistic jet."""
        Q_wave = 1e12  # THz resonance
        g_eff = c**2 / Gamma  # Effective "gravity" in jet frame
        
        F_UBii = -F_rel * (E_iso / E_LEP) * Q_wave * g_eff / c**2
        return F_UBii
    
    def _compute_Um_jet(self, Gamma: float) -> float:
        """UQFF magnetism for jet collimation."""
        mu_j = 3.38e20  # T pm³
        return mu_j * Gamma * 7.09e-37  # Scaled by vacuum density


# ============== Pre-defined Systems ==============

# GRB 221009A - "BOAT" (Brightest Of All Time)
GRB_221009A = {
    'name': 'GRB 221009A',
    'E_iso': 1e55 * erg_to_J,  # Record-breaking
    'eta': 500,
    'R_0': 1e5,  # 10^7 cm
    'z': 0.151
}

# GRB 170817A - First GW-associated GRB
GRB_170817A = {
    'name': 'GRB 170817A',
    'E_iso': 5e46 * erg_to_J,  # Off-axis viewed
    'm1': 1.46 * M_sun,
    'm2': 1.27 * M_sun,
    'M_ej': 0.05 * M_sun,
    'v_ej': 0.2 * c
}

# PSR J0737-3039 - Double pulsar
PSR_J0737_3039 = {
    'name': 'PSR J0737-3039',
    'm1': 1.3381 * M_sun,
    'm2': 1.2489 * M_sun,
    'P_b': 2.4 * 3600,  # 2.4 hours
    'e': 0.0878
}

# Registry
GRB_SYSTEMS = {
    'GRB_221009A': GRB_221009A,
    'GRB_170817A': GRB_170817A
}

BINARY_PULSAR_SYSTEMS = {
    'PSR_J0737-3039': PSR_J0737_3039
}

GRB_CALCULATORS = {
    'FireballExpansion': FireballExpansionCalculator,
    'AfterglowSynchrotron': AfterglowSynchrotronCalculator,
    'ChirpMass': ChirpMassCalculator,
    'RingdownQNM': RingdownQNMCalculator,
    'BinaryPulsarOrbitDecay': BinaryPulsarOrbitDecayCalculator,
    'PeriastronAdvance': PeriastronAdvanceCalculator,
    'KilonovaLightCurve': KilonovaLightCurveCalculator,
    'InspiralFrequencyEvolution': InspiralFrequencyEvolutionCalculator,
    'MergerTime': MergerTimeCalculator,
    'GRB': GRBCalculator
}


def run_demo():
    """Demonstrate GRB and compact binary calculations."""
    print("=" * 80)
    print("GRB EQUATIONS MODULE - Grok Deep Analysis")
    print("=" * 80)
    
    calc = GRBCalculator()
    
    # Fireball example
    print("\n--- Fireball Expansion (GRB 221009A-like) ---")
    fireball = calc.fireball_calc.compute(
        r=1e13,  # 10^8 km
        R_0=1e5,
        E_iso=1e55 * erg_to_J,
        M_baryon=1e55 * erg_to_J / (500 * c**2)  # η = 500
    )
    print(f"Lorentz Factor: Γ = {fireball['Gamma']:.1f}")
    print(f"Regime: {fireball['regime']}")
    print(f"η = {fireball['eta']:.0f}")
    
    # Binary pulsar
    print("\n--- Binary Pulsar (PSR J0737-3039) ---")
    binary = PSR_J0737_3039
    orbit = calc.orbit_decay_calc.compute(
        P_b=binary['P_b'],
        m1=binary['m1'],
        m2=binary['m2'],
        e=binary['e']
    )
    print(f"Period Derivative: Ṗ = {orbit['P_b_dot_s_s']:.2e} s/s")
    print(f"Merger Time: {orbit['t_merge_yr']:.2e} years")
    
    periastron = calc.periastron_calc.compute(
        P_b=binary['P_b'],
        m1=binary['m1'],
        m2=binary['m2'],
        e=binary['e']
    )
    print(f"Periastron Advance: {periastron['omega_dot_deg_yr']:.2f} °/year")
    
    # Kilonova
    print("\n--- Kilonova (GW170817-like) ---")
    kilonova = calc.kilonova_calc.compute(
        M_ej=0.05 * M_sun,
        v_ej=0.2 * c,
        kappa=1e-3  # 1 cm²/g in SI
    )
    print(f"Peak Luminosity: {kilonova['L_peak_erg_s']:.2e} erg/s")
    print(f"Peak Time: {kilonova['t_peak_days']:.1f} days")


if __name__ == '__main__':
    run_demo()
