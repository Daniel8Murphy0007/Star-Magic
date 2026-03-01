#!/usr/bin/env python3
"""
MHD Dynamo Module - Magnetic Field Generation and Evolution

From Grok Deep Analysis (Feb 2026):
- Equations 88-90: Kazantsev dynamo, Alfvén Mach number, field reversals
- Equations 39-41: ISM turbulence cascade, Larson relations

Physics domains covered:
- Small-scale turbulent dynamo (Kazantsev)
- Large-scale mean-field dynamo (α-Ω)
- Alfvénic turbulence and MHD cascades
- Magnetic field amplification
- Geomagnetic/solar reversals

UQFF Integration:
- Aether field as magnetic carrier medium
- Vacuum permittivity gradients in dynamo action
- Validates Electric Universe magnetism locally
"""

import math
from typing import Dict, Optional

# ============== Physical Constants ==============
G = 6.674e-11           # Gravitational constant [m³/(kg·s²)]
c = 2.998e8             # Speed of light [m/s]
M_sun = 1.989e30        # Solar mass [kg]
mu_0 = 1.257e-6         # Vacuum permeability [H/m]
epsilon_0 = 8.854e-12   # Vacuum permittivity [F/m]
k_B = 1.381e-23         # Boltzmann constant [J/K]
m_p = 1.673e-27         # Proton mass [kg]
e = 1.602e-19           # Elementary charge [C]
pc_to_m = 3.086e16      # parsec to meters
kpc_to_m = 3.086e19     # kiloparsec to meters
year_to_s = 3.154e7     # year to seconds
km_to_m = 1000          # km to meters

# UQFF Constants
F_rel = 4.30e33         # Relativistic coherence force [N]
rho_vac_UA = 7.09e-36   # Vacuum density UA [J/m³]


class KazantsevDynamoCalculator:
    """
    Kazantsev small-scale turbulent dynamo.
    
    Equation 88:
    B_sat ~ √(4πρ) v_turb, γ_K ~ v_turb / l_turb × Re_m^{1/2}
    
    Where:
    - Re_m = v l / η (magnetic Reynolds number)
    - η: magnetic diffusivity
    - Growth rate γ_K amplifies seed field to equipartition
    
    Derivation: Random stretching of field lines by turbulent
    motions, kinematic dynamo in high Re_m limit.
    """
    
    def compute(self, rho: float, v_turb: float, l_turb: float,
                eta: float = 1e4, B_seed: float = 1e-15) -> Dict:
        """
        Compute Kazantsev dynamo parameters.
        
        Args:
            rho: Density [kg/m³]
            v_turb: Turbulent velocity [m/s]
            l_turb: Turbulent scale [m]
            eta: Magnetic diffusivity [m²/s]
            B_seed: Seed magnetic field [T]
        
        Returns:
            Dict with dynamo parameters
        """
        # Magnetic Reynolds number
        Re_m = v_turb * l_turb / eta
        
        # Eddy turnover time
        tau_eddy = l_turb / v_turb
        
        # Kazantsev growth rate (kinematic phase)
        gamma_K = v_turb / l_turb * math.sqrt(Re_m) if Re_m > 1 else 0
        
        # Saturation field (equipartition)
        # B² / (2μ₀) = ρ v²/2
        B_sat = math.sqrt(mu_0 * rho) * v_turb
        
        # e-folding time
        t_efold = 1 / gamma_K if gamma_K > 0 else float('inf')
        
        # Time to saturation
        if B_seed > 0 and gamma_K > 0:
            n_efold = math.log(B_sat / B_seed)
            t_sat = n_efold / gamma_K
        else:
            t_sat = float('inf')
        
        # Resistive scale (where diffusion dominates)
        l_resist = l_turb / Re_m**0.5 if Re_m > 1 else l_turb
        
        return {
            'B_sat_T': B_sat,
            'B_sat_G': B_sat * 1e4,  # Gauss
            'gamma_K_Hz': gamma_K,
            't_efold_s': t_efold,
            't_efold_yr': t_efold / year_to_s,
            't_sat_yr': t_sat / year_to_s,
            'Re_m': Re_m,
            'tau_eddy_s': tau_eddy,
            'l_resist_m': l_resist,
            'B_seed_T': B_seed,
            'equation': 'B_sat ~ √(4πρ) v_turb, γ ~ v/l × Re_m^{1/2}'
        }


class AlfvenMachNumberCalculator:
    """
    Alfvén Mach number and MHD turbulence.
    
    Equation 89:
    M_A = v / v_A = v √(μ₀ρ) / B
    
    Where:
    - v_A = B / √(μ₀ρ): Alfvén velocity
    - M_A < 1: sub-Alfvénic (magnetically dominated)
    - M_A > 1: super-Alfvénic (kinetic dominated)
    
    Derivation: Ratio of kinetic to magnetic energy flux,
    determines MHD cascade anisotropy.
    """
    
    def compute(self, B: float, rho: float, v: float) -> Dict:
        """
        Compute Alfvén Mach number.
        
        Args:
            B: Magnetic field [T]
            rho: Density [kg/m³]
            v: Flow velocity [m/s]
        
        Returns:
            Dict with MHD parameters
        """
        # Alfvén velocity
        v_A = B / math.sqrt(mu_0 * rho) if B > 0 else float('inf')
        
        # Alfvén Mach number
        M_A = v / v_A if v_A > 0 and v_A != float('inf') else float('inf')
        
        # Magnetic pressure
        P_B = B**2 / (2 * mu_0)
        
        # Plasma beta (thermal to magnetic pressure)
        # Assuming T ~ 10⁴ K for ISM
        T = 1e4
        P_th = rho * k_B * T / m_p
        beta = P_th / P_B if P_B > 0 else float('inf')
        
        # Regime classification
        if M_A < 1:
            regime = 'sub-Alfvénic (magnetically dominated)'
        elif M_A < 10:
            regime = 'trans-Alfvénic'
        else:
            regime = 'super-Alfvénic (kinetically dominated)'
        
        # Goldreich-Sridhar anisotropy
        # l_perp / l_para ~ M_A for strong MHD turbulence
        anisotropy = 1 / M_A if M_A > 0 else 0
        
        return {
            'M_A': M_A,
            'v_A_m_s': v_A,
            'v_A_km_s': v_A / km_to_m,
            'v_m_s': v,
            'B_T': B,
            'B_uG': B * 1e10,  # microgauss
            'rho_kg_m3': rho,
            'P_B_Pa': P_B,
            'beta': beta,
            'regime': regime,
            'GS_anisotropy': anisotropy,
            'equation': 'M_A = v/v_A = v√(μ₀ρ)/B'
        }


class FieldReversalCalculator:
    """
    Magnetic field polarity reversals (geodynamo/solar).
    
    Equation 90:
    τ_rev ~ τ_ohm × (E_m / E_k)^α ∝ Rm^{-β}
    
    Where:
    - τ_ohm = l² / η (ohmic diffusion time)
    - E_m / E_k ~ 1 at reversal threshold
    - α, β ~ 1-2 from simulations
    
    For Earth: τ_rev ~ 10⁵-10⁶ years
    For Sun: τ_rev ~ 11 years (activity cycle)
    
    Derivation: Chaotic dynamo behavior near marginal
    stability, stochastic reversals from turbulent fluctuations.
    """
    
    def compute(self, l: float, eta: float, v: float,
                object_type: str = 'planet') -> Dict:
        """
        Compute reversal timescale.
        
        Args:
            l: Characteristic length scale [m]
            eta: Magnetic diffusivity [m²/s]
            v: Flow velocity [m/s]
            object_type: 'planet', 'star', or 'disk'
        
        Returns:
            Dict with reversal parameters
        """
        # Ohmic diffusion time
        tau_ohm = l**2 / eta
        
        # Magnetic Reynolds number
        Rm = v * l / eta
        
        # Turnover time
        tau_conv = l / v
        
        # Reversal timescale (empirical scaling)
        # Stochastic process with mean interval
        if object_type == 'planet':
            # Earth-like: ~10⁵-10⁶ years
            alpha = 1.5
            tau_rev = tau_ohm * Rm**(-alpha/2)
        elif object_type == 'star':
            # Solar-like: ~11 years (but modulated)
            tau_rev = tau_conv * 10  # Approx period
        else:
            # Galactic disk
            tau_rev = tau_ohm * 0.1
        
        # Reversal rate
        rate = 1 / tau_rev if tau_rev > 0 else 0
        
        return {
            'tau_rev_s': tau_rev,
            'tau_rev_yr': tau_rev / year_to_s,
            'tau_ohm_s': tau_ohm,
            'tau_ohm_yr': tau_ohm / year_to_s,
            'tau_conv_s': tau_conv,
            'Rm': Rm,
            'reversal_rate_per_Myr': rate * 1e6 * year_to_s,
            'object_type': object_type,
            'equation': 'τ_rev ~ τ_ohm × Rm^{-α}'
        }


class MeanFieldDynamoCalculator:
    """
    Mean-field (α-Ω) dynamo theory.
    
    ∂B/∂t = ∇×(v×B + αB) - η∇²B
    
    Key numbers:
    - α-effect: α ~ v·ω / k (helical turbulence)
    - Ω-effect: differential rotation shearing poloidal → toroidal
    - Dynamo number: D = αΔΩ L³/η² (must exceed critical D_c)
    
    Derivation: Mean-field electrodynamics, separation of
    large-scale and fluctuating components.
    """
    
    def compute(self, alpha: float, Delta_Omega: float, L: float,
                eta: float) -> Dict:
        """
        Compute mean-field dynamo parameters.
        
        Args:
            alpha: α-effect coefficient [m/s]
            Delta_Omega: Differential rotation [rad/s]
            L: Characteristic scale [m]
            eta: Magnetic diffusivity [m²/s]
        
        Returns:
            Dict with α-Ω dynamo parameters
        """
        # Dynamo number
        D = alpha * Delta_Omega * L**3 / eta**2
        
        # Critical dynamo number (order ~10)
        D_c = 10
        
        # Is dynamo active?
        active = abs(D) > D_c
        
        # Oscillation period (if α-Ω type)
        # P ~ η / (α L)
        P_cyc = eta / (abs(alpha) * L) if alpha != 0 else float('inf')
        
        # Growth rate
        gamma = math.sqrt(abs(D) - D_c) * eta / L**2 if active else 0
        
        return {
            'D': D,
            'D_c': D_c,
            'active': active,
            'alpha_m_s': alpha,
            'Delta_Omega_rad_s': Delta_Omega,
            'P_cyc_s': P_cyc,
            'P_cyc_yr': P_cyc / year_to_s,
            'gamma_Hz': gamma,
            'L_m': L,
            'eta_m2_s': eta,
            'equation': 'D = α ΔΩ L³/η² > D_c'
        }


class ISMTurbulenceCascadeCalculator:
    """
    ISM turbulence energy cascade (Kolmogorov/Goldreich-Sridhar).
    
    Equation 39:
    E(k) ∝ k^{-5/3} (Kolmogorov, hydrodynamic)
    E(k) ∝ k_perp^{-5/3} (GS95, MHD)
    
    Equation 40:
    v(l) ~ v_L (l/L)^{1/3} (velocity-size relation)
    
    Equation 41:
    σ_v ~ 1 km/s (L/pc)^{0.5} (Larson relation)
    
    Derivation: Energy injected at L, cascades to dissipation
    scale with constant energy flux ε ~ v³/l.
    """
    
    def compute(self, L_inj: float, v_inj: float, rho: float,
                B: Optional[float] = None) -> Dict:
        """
        Compute turbulence cascade parameters.
        
        Args:
            L_inj: Injection scale [m]
            v_inj: Injection velocity [m/s]
            rho: Density [kg/m³]
            B: Magnetic field [T] (optional, for MHD)
        
        Returns:
            Dict with cascade parameters
        """
        # Energy flux (energy cascade rate)
        epsilon = v_inj**3 / L_inj
        
        # Kolmogorov dissipation scale
        # l_diss ~ (ν³/ε)^{1/4}
        nu_kinematic = 1e-5  # m²/s (estimate for warm ISM)
        l_diss = (nu_kinematic**3 / epsilon)**(1/4)
        
        # Turbulent energy
        E_turb = 0.5 * rho * v_inj**2
        
        # Dissipation time
        tau_diss = L_inj / v_inj
        
        # MHD extensions if B provided
        if B is not None and B > 0:
            v_A = B / math.sqrt(mu_0 * rho)
            M_A = v_inj / v_A
            
            # GS critical balance: k_para ~ k_perp^{2/3}
            # Anisotropy at scale l
            anisotropy = M_A**(4/3)
            
            # MHD cascade rate
            epsilon_mhd = v_inj**3 / L_inj * min(1, M_A**2)
        else:
            v_A = None
            M_A = None
            anisotropy = 1
            epsilon_mhd = epsilon
        
        # Larson relation check
        L_pc = L_inj / pc_to_m
        sigma_Larson = 1 * km_to_m * math.sqrt(L_pc)
        
        return {
            'epsilon_W_kg': epsilon,
            'E_turb_J_m3': E_turb,
            'l_diss_m': l_diss,
            'tau_diss_s': tau_diss,
            'tau_diss_Myr': tau_diss / (1e6 * year_to_s),
            'L_inj_pc': L_pc,
            'v_inj_km_s': v_inj / km_to_m,
            'sigma_Larson_km_s': sigma_Larson / km_to_m,
            'M_A': M_A,
            'GS_anisotropy': anisotropy,
            'equation': 'E(k) ∝ k^{-5/3}, v(l) ~ v_L(l/L)^{1/3}'
        }


class MagneticFluxFreezeInCalculator:
    """
    Magnetic flux freezing (ideal MHD).
    
    d/dt (Φ_B) = d/dt ∫ B·dA = 0 (ideal limit)
    
    ⟹ B ∝ ρ^{2/3} for isotropic compression
    ⟹ B ∝ ρ for flux tubes (1D compression)
    
    Derivation: Faraday's law with E = -v×B,
    magnetic field advected with fluid.
    """
    
    def compute(self, B_0: float, rho_0: float, rho: float,
                compression_mode: str = 'isotropic') -> Dict:
        """
        Compute flux frozen field evolution.
        
        Args:
            B_0: Initial field [T]
            rho_0: Initial density [kg/m³]
            rho: Final density [kg/m³]
            compression_mode: 'isotropic' or 'flux_tube'
        
        Returns:
            Dict with field evolution
        """
        compression_ratio = rho / rho_0
        
        if compression_mode == 'isotropic':
            # 3D compression: B ∝ ρ^{2/3}
            B = B_0 * compression_ratio**(2/3)
            exponent = 2/3
        else:
            # 1D flux tube: B ∝ ρ
            B = B_0 * compression_ratio
            exponent = 1.0
        
        # Magnetic pressure
        P_B_0 = B_0**2 / (2 * mu_0)
        P_B = B**2 / (2 * mu_0)
        
        # Magnetic energy density
        u_B_0 = B_0**2 / (2 * mu_0)
        u_B = B**2 / (2 * mu_0)
        
        return {
            'B_T': B,
            'B_G': B * 1e4,
            'B_0_T': B_0,
            'compression_ratio': compression_ratio,
            'B_exponent': exponent,
            'P_B_Pa': P_B,
            'P_B_0_Pa': P_B_0,
            'u_B_J_m3': u_B,
            'compression_mode': compression_mode,
            'equation': 'B ∝ ρ^{2/3} (isotropic) or B ∝ ρ (1D)'
        }


class MHDDynamoCalculator:
    """
    Master calculator for MHD dynamo physics.
    
    Integrates:
    - Kazantsev small-scale dynamo
    - Mean-field α-Ω dynamo
    - Alfvénic turbulence
    - Field reversals
    - ISM cascade
    - Flux freezing
    
    UQFF Integration:
    - Aether field as magnetic carrier
    - Vacuum permittivity gradients
    - Validates Ug3 rotational magnetism
    """
    
    def __init__(self):
        self.kazantsev = KazantsevDynamoCalculator()
        self.alfven_mach = AlfvenMachNumberCalculator()
        self.reversal = FieldReversalCalculator()
        self.mean_field = MeanFieldDynamoCalculator()
        self.cascade = ISMTurbulenceCascadeCalculator()
        self.flux_freeze = MagneticFluxFreezeInCalculator()
    
    def compute_full_mhd_analysis(self, B: float, rho: float,
                                   v: float, L: float,
                                   eta: float = 1e4) -> Dict:
        """
        Complete MHD analysis.
        
        Args:
            B: Magnetic field [T]
            rho: Density [kg/m³]
            v: Velocity [m/s]
            L: Scale [m]
            eta: Diffusivity [m²/s]
        
        Returns:
            Comprehensive MHD analysis
        """
        # Kazantsev dynamo
        kaz = self.kazantsev.compute(rho, v, L, eta)
        
        # Alfvén Mach
        mach = self.alfven_mach.compute(B, rho, v)
        
        # Cascade
        cascade = self.cascade.compute(L, v, rho, B)
        
        # UQFF integration
        # Ug3 rotational magnetism: Um = Ug3 = |m| × Ω × r²
        Omega = v / L  # Angular velocity estimate
        m_aether = 1e-30  # Effective aether mass [kg]
        Um = abs(m_aether) * Omega * L**2
        
        return {
            'kazantsev': kaz,
            'alfven_mach': mach,
            'cascade': cascade,
            'UQFF': {
                'Um_magnetism_J': Um,
                'Omega_Hz': Omega,
                'note': 'Ug3 rotational magnetism via aether'
            }
        }


# ============== Pre-defined Systems ==============

SOLAR_DYNAMO = {
    'name': 'Solar Dynamo',
    'L': 7e8,           # Solar radius [m]
    'v': 50,            # Convective velocity [m/s]
    'B': 1e-4,          # Poloidal field [T]
    'eta': 1e3,         # Magnetic diffusivity [m²/s]
    'Delta_Omega': 3e-6 # Differential rotation [rad/s]
}

EARTH_GEODYNAMO = {
    'name': 'Earth Geodynamo',
    'L': 3.5e6,         # Outer core radius [m]
    'v': 1e-3,          # Core convection [m/s]
    'B': 3e-5,          # Dipole field [T]
    'eta': 2,           # Magnetic diffusivity [m²/s]
    'tau_rev': 3e5 * year_to_s  # Reversal interval [s]
}

ISM_TURBULENCE = {
    'name': 'ISM Warm Phase',
    'L': 100 * pc_to_m, # Injection scale
    'v': 10 * km_to_m,  # Turbulent velocity
    'rho': 1e-21,       # Density [kg/m³]
    'B': 3e-10          # Field [T] = 3 μG
}

GALAXY_DISK_DYNAMO = {
    'name': 'Milky Way Disk Dynamo',
    'L': 1 * kpc_to_m,  # Scale height
    'v': 50 * km_to_m,  # Turbulent velocity
    'B': 6e-10,         # Field [T] = 6 μG
    'eta': 1e26         # Turbulent diffusivity [m²/s]
}

MHD_SYSTEMS = {
    'Solar': SOLAR_DYNAMO,
    'Earth': EARTH_GEODYNAMO,
    'ISM': ISM_TURBULENCE,
    'Galaxy': GALAXY_DISK_DYNAMO
}

MHD_DYNAMO_CALCULATORS = {
    'KazantsevDynamo': KazantsevDynamoCalculator,
    'AlfvenMachNumber': AlfvenMachNumberCalculator,
    'FieldReversal': FieldReversalCalculator,
    'MeanFieldDynamo': MeanFieldDynamoCalculator,
    'ISMTurbulenceCascade': ISMTurbulenceCascadeCalculator,
    'MagneticFluxFreezeIn': MagneticFluxFreezeInCalculator,
    'MHDDynamo': MHDDynamoCalculator
}


def run_demo():
    """Demonstrate MHD dynamo calculations."""
    print("=" * 80)
    print("MHD DYNAMO MODULE - Grok Deep Analysis")
    print("=" * 80)
    
    calc = MHDDynamoCalculator()
    
    # Solar dynamo
    print("\n--- Solar Dynamo (α-Ω) ---")
    solar = SOLAR_DYNAMO
    mf = calc.mean_field.compute(
        alpha=10,  # m/s
        Delta_Omega=solar['Delta_Omega'],
        L=solar['L'],
        eta=solar['eta']
    )
    print(f"Dynamo number D = {mf['D']:.1f}")
    print(f"Active: {mf['active']}")
    print(f"Cycle period: {mf['P_cyc_yr']:.1f} years")
    
    # Earth geodynamo reversals
    print("\n--- Earth Geodynamo Reversals ---")
    earth = EARTH_GEODYNAMO
    rev = calc.reversal.compute(earth['L'], earth['eta'], earth['v'], 'planet')
    print(f"Ohmic diffusion time: {rev['tau_ohm_yr']:.2e} years")
    print(f"Reversal interval: {rev['tau_rev_yr']:.2e} years")
    print(f"Reversal rate: {rev['reversal_rate_per_Myr']:.1f} per Myr")
    
    # ISM turbulence
    print("\n--- ISM Turbulence Cascade ---")
    ism = ISM_TURBULENCE
    cascade = calc.cascade.compute(ism['L'], ism['v'], ism['rho'], ism['B'])
    print(f"Injection scale: {cascade['L_inj_pc']:.0f} pc")
    print(f"Injection velocity: {cascade['v_inj_km_s']:.1f} km/s")
    print(f"Larson prediction: {cascade['sigma_Larson_km_s']:.1f} km/s")
    print(f"Alfvén Mach: {cascade['M_A']:.2f}")
    
    # Kazantsev dynamo
    print("\n--- Kazantsev Dynamo (Proto-galaxy) ---")
    proto = calc.kazantsev.compute(
        rho=1e-23,        # Early universe
        v_turb=30 * km_to_m,
        l_turb=10 * kpc_to_m,
        eta=1e20,
        B_seed=1e-20      # Primordial seed
    )
    print(f"Saturation field: {proto['B_sat_G']:.2e} G")
    print(f"Re_m = {proto['Re_m']:.2e}")
    print(f"e-folding time: {proto['t_efold_yr']:.2e} years")
    
    # Flux freezing
    print("\n--- Flux Freezing (Star Formation) ---")
    freeze = calc.flux_freeze.compute(
        B_0=3e-10,     # ISM field
        rho_0=1e-21,   # ISM density
        rho=1e-13,     # Protostellar core
        compression_mode='isotropic'
    )
    print(f"Initial B: {freeze['B_0_T']*1e10:.1f} μG")
    print(f"Compression: {freeze['compression_ratio']:.2e}×")
    print(f"Final B: {freeze['B_G']:.2e} G")


if __name__ == '__main__':
    run_demo()
