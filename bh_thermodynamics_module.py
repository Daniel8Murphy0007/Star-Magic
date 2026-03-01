#!/usr/bin/env python3
"""
Black Hole Thermodynamics Module - Hawking Radiation and LQC Bounce

From Grok Deep Analysis (Feb 2026):
- Equations 94-96: Hawking temperature, Bekenstein-Hawking entropy, evaporation time
- Equations 97-99: Loop Quantum Cosmology bounce, Planck density, quantum corrections

Physics domains covered:
- Black hole thermodynamics
- Hawking radiation and temperature
- Bekenstein-Hawking entropy
- Black hole evaporation
- Loop quantum gravity corrections
- LQC bounce cosmology
- Information paradox aspects

UQFF Integration:
- Vacuum density rho_vac near event horizon
- Quantum corrections to gravitational coupling
- Validates holographic principles via [SSq]
"""

import math
from typing import Dict, Optional

# ============== Physical Constants ==============
G = 6.674e-11           # Gravitational constant [m³/(kg·s²)]
c = 2.998e8             # Speed of light [m/s]
h_bar = 1.055e-34       # Reduced Planck constant [J·s]
k_B = 1.381e-23         # Boltzmann constant [J/K]
sigma_SB = 5.670e-8     # Stefan-Boltzmann [W/(m²·K⁴)]
M_sun = 1.989e30        # Solar mass [kg]
eV_to_J = 1.602e-19     # eV to Joules
year_to_s = 3.154e7     # year to seconds

# Planck units
l_P = math.sqrt(h_bar * G / c**3)  # Planck length [m]
t_P = l_P / c                       # Planck time [s]
m_P = math.sqrt(h_bar * c / G)     # Planck mass [kg]
T_P = m_P * c**2 / k_B             # Planck temperature [K]
rho_P = m_P / l_P**3               # Planck density [kg/m³]

# UQFF Constants
F_rel = 4.30e33         # Relativistic coherence force [N]
SSq = 0.57              # Holographic entropy factor
rho_vac_SCm = 7.09e-37  # Vacuum density SCm [J/m³]


class HawkingTemperatureCalculator:
    """
    Hawking temperature of black holes.
    
    Equation 94:
    T_H = ℏc³ / (8πGMk_B)
    
    Simplified:
    T_H = 6.17 × 10⁻⁸ (M_☉/M) K
    
    Derivation: Quantum field theory in curved spacetime,
    pair creation at event horizon with thermal spectrum.
    """
    
    def compute(self, M: float) -> Dict:
        """
        Compute Hawking temperature.
        
        Args:
            M: Black hole mass [kg]
        
        Returns:
            Dict with temperature parameters
        """
        # Hawking temperature
        T_H = h_bar * c**3 / (8 * math.pi * G * M * k_B)
        
        # Schwarzschild radius
        r_s = 2 * G * M / c**2
        
        # Surface gravity
        kappa = c**4 / (4 * G * M)  # = c²/2r_s
        
        # Peak wavelength of radiation
        lambda_peak = 2.898e-3 / T_H if T_H > 0 else float('inf')  # Wien's law
        
        # Luminosity (Stefan-Boltzmann)
        A_horizon = 4 * math.pi * r_s**2
        L_H = sigma_SB * A_horizon * T_H**4
        
        # Energy of typical Hawking photon
        E_photon = k_B * T_H
        
        return {
            'T_H_K': T_H,
            'T_H_T_P': T_H / T_P,
            'r_s_m': r_s,
            'r_s_l_P': r_s / l_P,
            'kappa_m_s2': kappa,
            'lambda_peak_m': lambda_peak,
            'L_H_W': L_H,
            'E_photon_eV': E_photon / eV_to_J,
            'M_kg': M,
            'M_Msun': M / M_sun,
            'equation': 'T_H = ℏc³/(8πGMk_B)'
        }


class BekensteinHawkingEntropyCalculator:
    """
    Bekenstein-Hawking entropy of black holes.
    
    Equation 95:
    S_BH = A / (4 l_P²) = 4π G M² / (ℏc)
    
    Where:
    - A = 4πr_s² = 16π(GM/c²)²: horizon area
    - l_P = √(ℏG/c³): Planck length
    
    Derivation: Bekenstein's area theorem + thermodynamics,
    Hawking's exact coefficient from QFT.
    """
    
    def compute(self, M: float) -> Dict:
        """
        Compute Bekenstein-Hawking entropy.
        
        Args:
            M: Black hole mass [kg]
        
        Returns:
            Dict with entropy parameters
        """
        # Schwarzschild radius
        r_s = 2 * G * M / c**2
        
        # Horizon area
        A = 4 * math.pi * r_s**2
        
        # Bekenstein-Hawking entropy
        S_BH = A / (4 * l_P**2)
        
        # In physical units
        S_BH_JK = S_BH * k_B  # In J/K
        
        # Information content
        N_bits = S_BH / math.log(2)
        
        # Compare to thermal entropy of star
        # S_star ~ N_particles k_B ~ M/m_p k_B
        S_star = (M / (M_sun * 1e-30)) * k_B / k_B  # dimensionless comparison
        S_ratio = S_BH / S_star if S_star > 0 else float('inf')
        
        # UQFF: Holographic contribution
        S_UQFF = SSq * S_BH  # Effective quantum entropy
        
        return {
            'S_BH': S_BH,
            'S_BH_kB': S_BH,
            'S_BH_JK': S_BH_JK,
            'A_m2': A,
            'A_l_P2': A / l_P**2,
            'N_bits': N_bits,
            'M_Msun': M / M_sun,
            'S_UQFF': S_UQFF,
            'SSq_factor': SSq,
            'equation': 'S_BH = A/(4l_P²) = πr_s²/l_P²'
        }


class BlackHoleEvaporationCalculator:
    """
    Black hole evaporation time.
    
    Equation 96:
    τ_evap = 5120 π G² M³ / (ℏ c⁴)
    
    Simplified:
    τ_evap ≈ 2.1 × 10⁶⁷ (M/M_☉)³ years
    
    Derivation: dM/dt = -L_H/c² with L_H ∝ T⁴ ∝ M⁻⁴,
    integrate M³ dM = const × dt.
    """
    
    def compute(self, M: float, t: float = 0) -> Dict:
        """
        Compute evaporation time and mass evolution.
        
        Args:
            M: Initial black hole mass [kg]
            t: Time elapsed [s] (for mass at time t)
        
        Returns:
            Dict with evaporation parameters
        """
        # Evaporation time (initial M)
        tau_evap = 5120 * math.pi * G**2 * M**3 / (h_bar * c**4)
        
        # Mass at time t
        # M(t)³ = M_0³ - 3 × const × t
        const_evap = h_bar * c**4 / (5120 * math.pi * G**2)
        M_cubed_t = M**3 - 3 * const_evap * t
        
        if M_cubed_t > 0:
            M_t = M_cubed_t**(1/3)
            evaporated = False
        else:
            M_t = 0
            evaporated = True
        
        # Evaporation rate (current)
        if M_t > 0:
            T_H = h_bar * c**3 / (8 * math.pi * G * M_t * k_B)
            r_s = 2 * G * M_t / c**2
            A = 4 * math.pi * r_s**2
            L_H = sigma_SB * A * T_H**4
            dMdt = -L_H / c**2
        else:
            dMdt = 0
            T_H = float('inf')
        
        # Final burst energy
        E_burst = M * c**2  # Total rest mass energy
        
        return {
            'tau_evap_s': tau_evap,
            'tau_evap_yr': tau_evap / year_to_s,
            'M_0_kg': M,
            'M_0_Msun': M / M_sun,
            'M_t_kg': M_t,
            't_s': t,
            'dMdt_kg_s': dMdt,
            'T_H_K': T_H,
            'E_burst_J': E_burst,
            'evaporated': evaporated,
            'equation': 'τ = 5120πG²M³/(ℏc⁴)'
        }


class LQCBounceCalculator:
    """
    Loop Quantum Cosmology (LQC) bounce dynamics.
    
    Equation 97:
    ρ_bounce ≈ ρ_Planck ≈ c⁵/(ℏG²) ≈ 5.16 × 10⁹⁶ kg/m³
    
    Equation 98:
    H² = 8πGρ/3 × (1 - ρ/ρ_c) (modified Friedmann)
    
    Equation 99:
    a_min = a_0 (ρ_0/ρ_c)^{1/6} (bounce scale factor)
    
    Derivation: Holonomy corrections from loop quantization,
    prevents singularity via quantum geometry effects.
    """
    
    def __init__(self, gamma: float = 0.2375):
        """
        Initialize LQC calculator.
        
        Args:
            gamma: Barbero-Immirzi parameter
        """
        self.gamma = gamma
        
        # Critical density (LQC)
        # ρ_c = √3 / (32π² γ³ G² ℏ)
        self.rho_c = math.sqrt(3) / (32 * math.pi**2 * gamma**3 * G**2 * h_bar) * c**5
    
    def compute_bounce(self, rho: float, a: float) -> Dict:
        """
        Compute LQC bounce parameters.
        
        Args:
            rho: Energy density [kg/m³]
            a: Scale factor
        
        Returns:
            Dict with LQC bounce parameters
        """
        # Modified Friedmann equation
        # H² = 8πGρ/3 × (1 - ρ/ρ_c)
        correction = 1 - rho / self.rho_c
        H_squared_classical = 8 * math.pi * G * rho / 3
        H_squared_LQC = H_squared_classical * correction
        
        # At bounce: H = 0, ρ = ρ_c
        at_bounce = abs(correction) < 1e-10
        
        H = math.sqrt(max(0, H_squared_LQC))
        
        # Minimum scale factor
        # a_min = a × (ρ/ρ_c)^{1/6}
        a_min = a * (rho / self.rho_c)**(1/6) if rho < self.rho_c else a
        
        # Quantum correction to expansion
        quantum_correction = correction
        
        # Time to bounce (if contracting)
        if H_squared_LQC > 0:
            # Approximate: τ ~ 1/H × (1-ρ/ρ_c)^{-1/2}
            tau_bounce = 1 / H if correction > 0 else 0
        else:
            tau_bounce = 0
        
        return {
            'rho_kg_m3': rho,
            'rho_rho_c': rho / self.rho_c,
            'rho_c_kg_m3': self.rho_c,
            'rho_P_kg_m3': rho_P,
            'H_Hz': H,
            'H_squared_classical': H_squared_classical,
            'H_squared_LQC': H_squared_LQC,
            'quantum_correction': quantum_correction,
            'a': a,
            'a_min': a_min,
            'at_bounce': at_bounce,
            'gamma': self.gamma,
            'equation': 'H² = 8πGρ/3 × (1 - ρ/ρ_c)'
        }
    
    def compute_planck_regime(self) -> Dict:
        """
        Compute Planck regime parameters.
        
        Returns:
            Dict with Planck scale parameters
        """
        return {
            'l_P_m': l_P,
            't_P_s': t_P,
            'm_P_kg': m_P,
            'T_P_K': T_P,
            'rho_P_kg_m3': rho_P,
            'rho_c_kg_m3': self.rho_c,
            'rho_c_rho_P': self.rho_c / rho_P,
            'E_P_J': m_P * c**2,
            'E_P_GeV': m_P * c**2 / (1e9 * eV_to_J)
        }


class PrimordialBlackHoleCalculator:
    """
    Primordial black hole formation and constraints.
    
    PBH mass from horizon at formation:
    M_PBH ~ c³ t / G ~ 10⁵ M_☉ (t/s)
    
    At matter-radiation equality: M ~ 10¹⁵ g (10⁻¹⁸ M_☉)
    evaporates now (τ ~ 10¹⁰ years).
    """
    
    def compute(self, t_form: float) -> Dict:
        """
        Compute PBH parameters.
        
        Args:
            t_form: Formation time after Big Bang [s]
        
        Returns:
            Dict with PBH parameters
        """
        # Horizon mass at formation
        M_horizon = c**3 * t_form / G
        
        # PBH mass (order of horizon mass)
        M_PBH = M_horizon
        
        # Evaporation time
        tau_evap = 5120 * math.pi * G**2 * M_PBH**3 / (h_bar * c**4)
        
        # Constraint masses
        M_evap_now = (h_bar * c**4 * 13.8e9 * year_to_s / 
                      (5120 * math.pi * G**2))**(1/3)
        
        # Current status
        age_universe = 13.8e9 * year_to_s
        if tau_evap < age_universe:
            status = 'evaporated'
        elif tau_evap < 10 * age_universe:
            status = 'evaporating_now'
        else:
            status = 'stable'
        
        return {
            'M_PBH_kg': M_PBH,
            'M_PBH_Msun': M_PBH / M_sun,
            'M_PBH_g': M_PBH * 1000,
            't_form_s': t_form,
            'tau_evap_s': tau_evap,
            'tau_evap_yr': tau_evap / year_to_s,
            'M_evap_now_kg': M_evap_now,
            'M_evap_now_g': M_evap_now * 1000,
            'status': status,
            'equation': 'M_PBH ~ c³t/G'
        }


class HolographicPrincipleCalculator:
    """
    Holographic principle bounds.
    
    Bekenstein bound:
    S ≤ 2πRE/(ℏc) = 2πMc R/(ℏ)
    
    Covariant entropy bound:
    S ≤ A/(4l_P²)
    
    Holographic principle: maximum entropy ~ area, not volume.
    """
    
    def compute(self, M: float, R: float) -> Dict:
        """
        Compute holographic bounds.
        
        Args:
            M: Mass [kg]
            R: Radius [m]
        
        Returns:
            Dict with holographic parameters
        """
        # Bekenstein bound
        E = M * c**2
        S_Bekenstein = 2 * math.pi * R * E / (h_bar * c)
        
        # Area bound (covariant)
        A = 4 * math.pi * R**2
        S_area = A / (4 * l_P**2)
        
        # Which bound is tighter?
        S_max = min(S_Bekenstein, S_area)
        tighter = 'Bekenstein' if S_Bekenstein < S_area else 'area'
        
        # Schwarzschild radius
        r_s = 2 * G * M / c**2
        
        # Is it a black hole?
        is_bh = R <= r_s
        
        # Information in qubits
        N_qubits = S_max / math.log(2)
        
        # UQFF [SSq] factor
        S_UQFF = SSq * S_max
        
        return {
            'S_Bekenstein': S_Bekenstein,
            'S_area': S_area,
            'S_max': S_max,
            'tighter_bound': tighter,
            'N_qubits': N_qubits,
            'A_m2': A,
            'R_m': R,
            'r_s_m': r_s,
            'is_black_hole': is_bh,
            'S_UQFF': S_UQFF,
            'equation': 'S ≤ min(2πRE/ℏc, A/4l_P²)'
        }


class BlackHoleThermodynamicsCalculator:
    """
    Master calculator for black hole thermodynamics.
    
    Integrates:
    - Hawking temperature
    - Bekenstein-Hawking entropy
    - Evaporation dynamics
    - LQC bounce
    - Primordial black holes
    - Holographic bounds
    
    UQFF Integration:
    - [SSq] holographic entropy factor
    - Vacuum density near horizons
    - Quantum gravity corrections
    """
    
    def __init__(self):
        self.hawking = HawkingTemperatureCalculator()
        self.entropy = BekensteinHawkingEntropyCalculator()
        self.evaporation = BlackHoleEvaporationCalculator()
        self.lqc = LQCBounceCalculator()
        self.pbh = PrimordialBlackHoleCalculator()
        self.holographic = HolographicPrincipleCalculator()
    
    def compute_full_thermodynamics(self, M: float) -> Dict:
        """
        Complete black hole thermodynamics.
        
        Args:
            M: Black hole mass [kg]
        
        Returns:
            Comprehensive analysis
        """
        # Hawking temperature
        hawking = self.hawking.compute(M)
        
        # Entropy
        entropy = self.entropy.compute(M)
        
        # Evaporation
        evap = self.evaporation.compute(M)
        
        # Holographic bounds
        r_s = 2 * G * M / c**2
        holo = self.holographic.compute(M, r_s)
        
        # UQFF vacuum energy at horizon
        # ρ_vac ~ T_H⁴ / (ℏc)³
        T_H = hawking['T_H_K']
        rho_vac_horizon = k_B**4 * T_H**4 / (h_bar**3 * c**3)
        
        return {
            'hawking': hawking,
            'entropy': entropy,
            'evaporation': evap,
            'holographic': holo,
            'UQFF': {
                'rho_vac_horizon_kg_m3': rho_vac_horizon,
                'SSq': SSq,
                'S_quantum': entropy['S_UQFF'],
                'note': 'Holographic entropy with [SSq] correction'
            }
        }


# ============== Pre-defined Systems ==============

SGR_A_STAR = {
    'name': 'Sagittarius A*',
    'M': 4e6 * M_sun,
    'type': 'SMBH'
}

M87_STAR = {
    'name': 'M87*',
    'M': 6.5e9 * M_sun,
    'type': 'SMBH'
}

STELLAR_BH_10 = {
    'name': 'Typical Stellar BH',
    'M': 10 * M_sun,
    'type': 'stellar'
}

PBH_ASTEROID = {
    'name': 'Asteroid-mass PBH',
    'M': 1e15,  # kg, ~10¹² g
    'type': 'primordial',
    'note': 'Dark matter candidate'
}

PLANCK_MASS_BH = {
    'name': 'Planck Mass BH',
    'M': m_P,
    'type': 'quantum',
    'note': 'Minimum BH mass'
}

BH_SYSTEMS = {
    'Sgr_A': SGR_A_STAR,
    'M87': M87_STAR,
    'Stellar_10': STELLAR_BH_10,
    'PBH_asteroid': PBH_ASTEROID,
    'Planck': PLANCK_MASS_BH
}

BH_THERMODYNAMICS_CALCULATORS = {
    'HawkingTemperature': HawkingTemperatureCalculator,
    'BekensteinHawkingEntropy': BekensteinHawkingEntropyCalculator,
    'BlackHoleEvaporation': BlackHoleEvaporationCalculator,
    'LQCBounce': LQCBounceCalculator,
    'PrimordialBlackHole': PrimordialBlackHoleCalculator,
    'HolographicPrinciple': HolographicPrincipleCalculator,
    'BlackHoleThermodynamics': BlackHoleThermodynamicsCalculator
}


def run_demo():
    """Demonstrate black hole thermodynamics calculations."""
    print("=" * 80)
    print("BLACK HOLE THERMODYNAMICS MODULE - Grok Deep Analysis")
    print("=" * 80)
    
    calc = BlackHoleThermodynamicsCalculator()
    
    # Hawking temperature across mass range
    print("\n--- Hawking Temperature ---")
    for M_Msun in [1, 10, 1e6, 1e9]:
        M = M_Msun * M_sun
        result = calc.hawking.compute(M)
        print(f"M = {M_Msun:.0e} M_☉: T_H = {result['T_H_K']:.2e} K")
    
    # Bekenstein-Hawking entropy
    print("\n--- Bekenstein-Hawking Entropy ---")
    for name, system in [('Sgr A*', SGR_A_STAR), ('M87*', M87_STAR)]:
        M = system['M']
        result = calc.entropy.compute(M)
        print(f"{name}: S = {result['S_BH']:.2e} k_B, N_bits = {result['N_bits']:.2e}")
    
    # Evaporation time
    print("\n--- Evaporation Time ---")
    for M_kg in [1e12, 1e15, M_sun, 10*M_sun]:
        result = calc.evaporation.compute(M_kg)
        if M_kg < M_sun:
            print(f"M = {M_kg:.0e} kg: τ = {result['tau_evap_yr']:.2e} years")
        else:
            print(f"M = {M_kg/M_sun:.0f} M_☉: τ = {result['tau_evap_yr']:.2e} years")
    
    # LQC bounce
    print("\n--- LQC Bounce (Big Bounce) ---")
    lqc = calc.lqc
    planck = lqc.compute_planck_regime()
    print(f"Planck density: {planck['rho_P_kg_m3']:.2e} kg/m³")
    print(f"Critical density (LQC): {planck['rho_c_kg_m3']:.2e} kg/m³")
    
    # Test bounce condition
    for rho_ratio in [0.1, 0.5, 0.9, 0.99, 1.0]:
        rho = rho_ratio * lqc.rho_c
        result = lqc.compute_bounce(rho, a=1)
        print(f"ρ/ρ_c = {rho_ratio}: H² correction = {result['quantum_correction']:.3f}")
    
    # Primordial black holes
    print("\n--- Primordial Black Holes ---")
    for t_form in [1e-23, 1e-5, 1, 1e10]:
        result = calc.pbh.compute(t_form)
        print(f"t_form = {t_form:.0e} s: M = {result['M_PBH_g']:.2e} g, "
              f"status: {result['status']}")
    
    # Holographic principle
    print("\n--- Holographic Bounds (M87*) ---")
    M87 = M87_STAR['M']
    r_s = 2 * G * M87 / c**2
    holo = calc.holographic.compute(M87, r_s)
    print(f"S_Bekenstein = {holo['S_Bekenstein']:.2e}")
    print(f"S_area = {holo['S_area']:.2e}")
    print(f"Qubits = {holo['N_qubits']:.2e}")


if __name__ == '__main__':
    run_demo()
