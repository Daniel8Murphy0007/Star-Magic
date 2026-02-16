"""
SOURCE86: Extended Fields MUGE (Master Universal Gravity Equation) Module

Extraction Date: 2026-02-14
Source: source86.cpp (MUGEModule, 715 lines, 29 KB)
"""
import math
from typing import Dict, Any, Optional
from enum import Enum

class SystemType(Enum):
    """Supported astrophysical systems."""
    MAGNETAR_SGR_1745_2900 = 0
    SAGITTARIUS_A = 1
    TAPESTRY_BLAZING_STARBIRTH = 2
    WESTERLUND_2 = 3
    PILLARS_CREATION = 4
    RINGS_RELATIVITY = 5
    STUDENTS_GUIDE_UNIVERSE = 6

class Source86_Extended:
    """
    Master Universal Gravity Equation (MUGE) calculator supporting 7 astrophysical systems.
    
    Implements 12 physics functions from source86.cpp (MUGEModule):
      1. calculate_hubble_expansion() - H(t,z) cosmological expansion
      2. calculate_ug_sum() - Σ Ugi (UQFF subterms)
      3. calculate_lambda_term() - Λc²/3 cosmological constant
      4. calculate_quantum_term() - Quantum entanglement integral
      5. calculate_fluid_term() - Fluid dynamics
      6. calculate_dm_term() - Dark matter perturbation + curvature
      7. calculate_system_specific_term() - 7 system-specific physics
      8. calculate_adpm() - Base resonance term (a_DPM)
      9. calculate_athz() - THz frequency resonance (a_THz)
      10. calculate_osc_term() - Oscillation term 2A cos(ωt)
      11. calculate_muge_compressed() - Master compressed MUGE
      12. calculate_muge_resonance() - Master resonance MUGE
    
    Physics Model:
    MUGE unifies compressed UQFF (base gravity + corrections) and resonance UQFF
    (frequency-domain terms) for multiple systems including magnetars, SMBHs,
    star-forming regions, and gravitational lenses.
    
    Compressed MUGE:
        g(r,t) = [G M / r²] × [1+H(t,z)] × [1-B/B_crit] × [1+F_env]
                 + Σ Ugi + Λc²/3 + quantum_term + EM_term + fluid_term
                 + resonant_term + DM_term + system_specific_term
    
    Resonance MUGE:
        g(r,t) = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res
                 + Ug4i + a_quantum_freq + a_Aether_freq + a_fluid_freq
                 + Osc_term + a_exp_freq + f_TRZ
    
    Supported Systems (7):
        1. Magnetar SGR 1745-2900: Ultra-strong B-field, spin-down
        2. Sagittarius A*: SMBH, frame-dragging, spin dynamics
        3. Tapestry of Blazing Starbirth: Stellar winds, high SFR
        4. Westerlund 2: Young cluster, wind ram pressure
        5. Pillars of Creation: Photoevaporation, erosion (1-E_t factor)
        6. Rings of Relativity: Gravitational lensing (1+L_t factor)
        7. Student Guide Universe: Cosmological, simplified
    
    UQFF Innovations:
        1. **Dual Computation Methods**: Compressed (direct) vs Resonance (frequency)
        2. **System-Specific Physics**: 7 distinct physical regimes
        3. **Frequency Resonances**: THz, quantum, aether, fluid, expansion modes
        4. **a_DPM Base Term**: c × V_sys × F_DPM × f_DPM × E_vac,neb (foundation)
    
    Validation:
        - SGR 1745-2900: g_comp ~ 1.79e12 m/s² (base dominated)
        - SGR 1745-2900: g_res ~ 1e-10 m/s² (resonant scaled)
        - Nebulae: Effective g ~ 1e-11 m/s² (fluid/resonant dominated)
        - Compressed suitable for strong-field regimes
        - Resonance suitable for weak-field/nebular regimes
    
    References:
        - source86.cpp (MUGEModule, 715 lines, 29 KB)
        - Supports: Magnetars, SMBHs, star-forming regions, lensing
        - Copyright: Daniel T. Murphy, analyzed Oct 10, 2025
    """
    
    DEFAULT_PARAMS = {
        # System selection
        'system': SystemType.MAGNETAR_SGR_1745_2900,
        
        # Universal constants
        'G': 6.6743e-11,           # m³ kg⁻¹ s⁻²
        'c': 3e8,                  # m/s
        'hbar': 1.0546e-34,        # J·s
        'Lambda': 1.1e-52,         # m⁻² (cosmological constant)
        'q': 1.602e-19,            # C
        'pi': 3.141592653589793,
        'M_sun': 1.989e30,         # kg
        
        # Cosmology
        'H0': 2.269e-18,           # s⁻¹ (70 km/s/Mpc)
        'Omega_m': 0.3,
        'Omega_Lambda': 0.7,
        't_Hubble': 4.35e17,       # s (13.8 Gyr)
        
        # System parameters (Magnetar SGR 1745-2900 defaults)
        'M': 3 * 1.989e30,         # kg (3 M_sun)
        'M_visible': 3 * 1.989e30,
        'M_DM': 0.0,               # No DM for magnetar
        'r': 1e4,                  # m (10 km radius)
        'z': 0.0,                  # Local (Galactic center)
        'B': 1e11,                 # T (magnetar field)
        'B_crit': 4.4e13,          # T (critical field)
        'F_env': 0.0,              # No envelope
        
        # UQFF terms
        'Ug1': 0.0,                # Negligible in compressed
        'Ug2': 0.0,
        'Ug3_prime': 0.0,          # External Ug3
        'Ug4': 0.0,
        
        # Quantum parameters
        'Delta_x': 1e-10,          # m (position uncertainty)
        'Delta_p': 1.0546e-24,     # kg·m/s (momentum uncertainty)
        'integral_psi': 2.176e-18, # J (normalized wavefunction integral)
        
        # Fluid & DM
        'rho_fluid': 1e-20,        # kg/m³
        'V': 1e3,                  # m³ (volume)
        'g_local': 9.8,            # m/s² (reference)
        'delta_rho_over_rho': 1e-5, # DM perturbation
        
        # Resonance parameters
        'Evac_neb': 7.09e-36,      # J/m³ (nebula vacuum energy)
        'Evac_ISM': 7.09e-37,      # J/m³ (ISM vacuum energy)
        'Delta_Evac': 6.381e-36,   # J/m³ (vacuum difference)
        'v_exp': 1e3,              # m/s (expansion velocity)
        'f_THz': 1e12,             # Hz (THz frequency)
        'f_DPM': 1e9,              # Hz (DPM frequency)
        'FDPM': 6.284e29,          # A·m² (DPM flux)
        'F_super': 6.287e-19,      # Superconductor factor
        'UA_SCm': 10.0,            # UA/SCm scaling
        'omega_i': 1e-8,           # rad/s (inertia frequency)
        'k4': 1.0,                 # Reactor coupling
        'f_react': 1e10,           # Hz (reactor frequency)
        'E_react': 1e-20,          # J (reactor energy)
        'f_quantum': 1.445e-17,    # Hz
        'f_Aether': 1.576e-35,     # Hz
        'f_fluid': 1.269e-14,      # Hz
        'f_osc': 4.57e14,          # Hz (oscillation)
        'f_exp': 1e-18,            # Hz (expansion)
        'f_TRZ': 0.1,              # TRZ factor
        
        # Oscillation parameters
        'A': 1e-10,                # Amplitude
        'k': 1e20,                 # Wave number
        'omega': 1e15,             # rad/s
        
        # System-specific parameters
        'v_wind': 1e6,             # m/s (stellar wind)
        'rho': 1e-20,              # kg/m³ (wind density)
        'E_t': 0.5,                # Erosion factor (Pillars)
        'L_t': 0.1,                # Lensing magnification (Rings)
        'dOmega_dt': 1e-10,        # rad/s² (SgrA* spin)
        'spin_adjust': 0.5,        # sin(30°)
        'scale_macro': 1e-12,      # EM scaling
    }
    
    @staticmethod
    def calculate_hubble_expansion(
        t: float,
        z: float = DEFAULT_PARAMS['z'],
        H0: float = DEFAULT_PARAMS['H0'],
        Omega_m: float = DEFAULT_PARAMS['Omega_m'],
        Omega_Lambda: float = DEFAULT_PARAMS['Omega_Lambda']
    ) -> float:
        """Calculate Hubble expansion factor H(t,z)."""
        Hz = H0 * math.sqrt(Omega_m * ((1 + z) ** 3) + Omega_Lambda)
        return Hz * t
    
    @staticmethod
    def calculate_ug_sum(
        Ug1: float = DEFAULT_PARAMS['Ug1'],
        Ug2: float = DEFAULT_PARAMS['Ug2'],
        Ug3_prime: float = DEFAULT_PARAMS['Ug3_prime'],
        Ug4: float = DEFAULT_PARAMS['Ug4']
    ) -> float:
        """Calculate Ug sum (Ug3' for external systems, others ~0 in compressed)."""
        return Ug1 + Ug2 + Ug3_prime + Ug4
    
    @staticmethod
    def calculate_lambda_term(
        Lambda: float = DEFAULT_PARAMS['Lambda'],
        c: float = DEFAULT_PARAMS['c']
    ) -> float:
        """Calculate cosmological constant term Λc²/3."""
        return (Lambda * c * c) / 3.0
    
    @staticmethod
    def calculate_quantum_term(
        hbar: float = DEFAULT_PARAMS['hbar'],
        Delta_x: float = DEFAULT_PARAMS['Delta_x'],
        Delta_p: float = DEFAULT_PARAMS['Delta_p'],
        integral_psi: float = DEFAULT_PARAMS['integral_psi'],
        pi: float = DEFAULT_PARAMS['pi'],
        t_Hubble: float = DEFAULT_PARAMS['t_Hubble']
    ) -> float:
        """Calculate quantum entanglement term (ℏ/√(ΔxΔp)) × ∫|ψ|² × (2π/t_H)."""
        unc = math.sqrt(Delta_x * Delta_p)
        return (hbar / unc) * integral_psi * (2 * pi / t_Hubble)
    
    @staticmethod
    def calculate_fluid_term(
        g_base: float,
        rho_fluid: float = DEFAULT_PARAMS['rho_fluid'],
        V: float = DEFAULT_PARAMS['V']
    ) -> float:
        """Calculate fluid dynamics term ρ_fluid × V × g_base."""
        return rho_fluid * V * g_base
    
    @staticmethod
    def calculate_dm_term(
        M: float = DEFAULT_PARAMS['M'],
        M_visible: float = DEFAULT_PARAMS['M_visible'],
        M_DM: float = DEFAULT_PARAMS['M_DM'],
        r: float = DEFAULT_PARAMS['r'],
        G: float = DEFAULT_PARAMS['G'],
        delta_rho_over_rho: float = DEFAULT_PARAMS['delta_rho_over_rho']
    ) -> float:
        """Calculate dark matter term (M_vis+M_DM) × (δρ/ρ + 3GM/r³)."""
        pert = delta_rho_over_rho
        curv = 3 * G * M / (r ** 3)
        return (M_visible + M_DM) * (pert + curv)
    
    @staticmethod
    def calculate_system_specific_term(
        system: SystemType,
        t: float,
        **params
    ) -> float:
        """
        Calculate system-specific physics term.
        
        Systems:
        --------
        MAGNETAR_SGR_1745_2900: ρ × v_wind² (stellar wind ram pressure)
        SAGITTARIUS_A: (GM²/c⁴r) × (dΩ/dt)² × sin(30°) (frame-dragging)
        TAPESTRY_BLAZING_STARBIRTH: ρ × v_wind² (stellar winds)
        WESTERLUND_2: ρ × v_wind² (wind ram pressure)
        PILLARS_CREATION: ρ × v_wind² × (1-E_t) (photoevaporation erosion)
        RINGS_RELATIVITY: ρ_fluid × V × g_local × (1+L_t) (lensing magnification)
        STUDENTS_GUIDE_UNIVERSE: 0.0 (simplified cosmological)
        """
        G = params.get('G', Source86_Extended.DEFAULT_PARAMS['G'])
        M = params.get('M', Source86_Extended.DEFAULT_PARAMS['M'])
        r = params.get('r', Source86_Extended.DEFAULT_PARAMS['r'])
        c = params.get('c', Source86_Extended.DEFAULT_PARAMS['c'])
        rho = params.get('rho', Source86_Extended.DEFAULT_PARAMS['rho'])
        rho_fluid = params.get('rho_fluid', Source86_Extended.DEFAULT_PARAMS['rho_fluid'])
        v_wind = params.get('v_wind', Source86_Extended.DEFAULT_PARAMS['v_wind'])
        V = params.get('V', Source86_Extended.DEFAULT_PARAMS['V'])
        g_local = params.get('g_local', Source86_Extended.DEFAULT_PARAMS['g_local'])
        E_t = params.get('E_t', Source86_Extended.DEFAULT_PARAMS['E_t'])
        L_t = params.get('L_t', Source86_Extended.DEFAULT_PARAMS['L_t'])
        dOmega_dt = params.get('dOmega_dt', Source86_Extended.DEFAULT_PARAMS['dOmega_dt'])
        spin_adjust = params.get('spin_adjust', Source86_Extended.DEFAULT_PARAMS['spin_adjust'])
        
        term = 0.0
        
        if system == SystemType.SAGITTARIUS_A:
            # Frame-dragging term
            term = (G * M * M / (c ** 4 * r)) * (dOmega_dt ** 2) * spin_adjust
        elif system in [SystemType.TAPESTRY_BLAZING_STARBIRTH, SystemType.WESTERLUND_2]:
            # Stellar wind ram pressure
            term = rho * (v_wind ** 2)
        elif system == SystemType.PILLARS_CREATION:
            # Photoevaporation with erosion factor
            term = rho * (v_wind ** 2) * (1 - E_t)
        elif system == SystemType.RINGS_RELATIVITY:
            # Gravitational lensing magnification
            term = rho_fluid * V * g_local * (1 + L_t)
        elif system == SystemType.STUDENTS_GUIDE_UNIVERSE:
            # Simplified cosmological
            term = 0.0
        else:
            # Default: wind ram pressure
            term = rho_fluid * (v_wind ** 2)
        
        return term
    
    @staticmethod
    def calculate_adpm(
        c: float = DEFAULT_PARAMS['c'],
        V: float = DEFAULT_PARAMS['V'],
        FDPM: float = DEFAULT_PARAMS['FDPM'],
        f_DPM: float = DEFAULT_PARAMS['f_DPM'],
        Evac_neb: float = DEFAULT_PARAMS['Evac_neb']
    ) -> float:
        """
        Calculate base resonance term a_DPM (Dipole Moment).
        
        Foundation for all resonance components:
        a_DPM = c × V_sys × F_DPM × f_DPM × E_vac,neb
        """
        return c * V * FDPM * f_DPM * Evac_neb
    
    @staticmethod
    def calculate_athz(
        adpm: float,
        Evac_ISM: float = DEFAULT_PARAMS['Evac_ISM'],
        c: float = DEFAULT_PARAMS['c'],
        f_THz: float = DEFAULT_PARAMS['f_THz'],
        Evac_neb: float = DEFAULT_PARAMS['Evac_neb'],
        v_exp: float = DEFAULT_PARAMS['v_exp']
    ) -> float:
        """
        Calculate THz frequency resonance a_THz.
        
        a_THz = (E_vac,ISM / c) × f_THz × E_vac,neb × v_exp × a_DPM
        """
        return (Evac_ISM / c) * f_THz * Evac_neb * v_exp * adpm
    
    @staticmethod
    def calculate_osc_term(
        t: float,
        A: float = DEFAULT_PARAMS['A'],
        f_osc: float = DEFAULT_PARAMS['f_osc'],
        pi: float = DEFAULT_PARAMS['pi']
    ) -> float:
        """
        Calculate oscillation term 2A cos(ωt).
        
        Simplified harmonic oscillation with frequency f_osc.
        """
        omega = f_osc * 2 * pi
        return 2 * A * math.cos(omega * t)
    
    @staticmethod
    def calculate_muge_compressed(
        t: float,
        params: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Master function: Calculate compressed MUGE g(r,t).
        
        Complete UQFF compressed model:
        g(r,t) = [G M / r²] × [1+H(t,z)] × [1-B/B_crit] × [1+F_env]
                 + Σ Ugi + Λc²/3 + quantum_term + EM_term + fluid_term
                 + resonant_term + DM_term + system_specific_term
        
        Suitable for strong-field regimes (magnetars, SMBHs, compact objects).
        
        Returns all components for validation and analysis.
        """
        if params is None:
            params = {}
        
        # Extract parameters with defaults
        G = params.get('G', Source86_Extended.DEFAULT_PARAMS['G'])
        M = params.get('M', Source86_Extended.DEFAULT_PARAMS['M'])
        r = params.get('r', Source86_Extended.DEFAULT_PARAMS['r'])
        z = params.get('z', Source86_Extended.DEFAULT_PARAMS['z'])
        B = params.get('B', Source86_Extended.DEFAULT_PARAMS['B'])
        B_crit = params.get('B_crit', Source86_Extended.DEFAULT_PARAMS['B_crit'])
        F_env = params.get('F_env', Source86_Extended.DEFAULT_PARAMS['F_env'])
        c = params.get('c', Source86_Extended.DEFAULT_PARAMS['c'])
        q = params.get('q', Source86_Extended.DEFAULT_PARAMS['q'])
        v_wind = params.get('v_wind', Source86_Extended.DEFAULT_PARAMS['v_wind'])
        scale_macro = params.get('scale_macro', Source86_Extended.DEFAULT_PARAMS['scale_macro'])
        system = params.get('system', Source86_Extended.DEFAULT_PARAMS['system'])
        
        # 1. Hubble expansion
        Hz_t = Source86_Extended.calculate_hubble_expansion(t, z, **params)
        expansion = 1.0 + Hz_t
        
        # 2. Superconductor correction
        sc_correction = 1.0 - (B / B_crit)
        
        # 3. Envelope factor
        env_factor = 1.0 + F_env
        
        # 4. Base gravity
        g_base = (G * M / (r ** 2)) * expansion * sc_correction * env_factor
        
        # 5. Ug sum
        ug_sum = Source86_Extended.calculate_ug_sum(**params)
        
        # 6. Lambda term
        lambda_term = Source86_Extended.calculate_lambda_term(**params)
        
        # 7. Quantum term
        quantum_term = Source86_Extended.calculate_quantum_term(**params)
        
        # 8. EM term q(v × B)
        em_term = (q * v_wind * B) / 1.673e-27 * scale_macro
        
        # 9. Fluid term
        fluid_term = Source86_Extended.calculate_fluid_term(g_base, **params)
        
        # 10. Resonant term (simplified cos + complex exp)
        A = params.get('A', Source86_Extended.DEFAULT_PARAMS['A'])
        k = params.get('k', Source86_Extended.DEFAULT_PARAMS['k'])
        omega = params.get('omega', Source86_Extended.DEFAULT_PARAMS['omega'])
        pi = params.get('pi', Source86_Extended.DEFAULT_PARAMS['pi'])
        cos_term = 2 * A * math.cos(k * 0.0) * math.cos(omega * t)
        exp_factor = (2 * pi / 13.8)
        real_exp = A * math.cos(k * 0.0 - omega * t)
        resonant_term = cos_term + exp_factor * real_exp
        
        # 11. DM term
        dm_term = Source86_Extended.calculate_dm_term(**params)
        
        # 12. System-specific term
        sys_term = Source86_Extended.calculate_system_specific_term(system, t, **params)
        
        # Total gravity
        g_total = g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + sys_term
        
        return {
            'g_total': g_total,
            'g_base': g_base,
            'expansion_factor': expansion,
            'sc_correction': sc_correction,
            'env_factor': env_factor,
            'Hz_t': Hz_t,
            'ug_sum': ug_sum,
            'lambda_term': lambda_term,
            'quantum_term': quantum_term,
            'em_term': em_term,
            'fluid_term': fluid_term,
            'resonant_term': resonant_term,
            'dm_term': dm_term,
            'sys_term': sys_term,
        }
    
    @staticmethod
    def calculate_muge_resonance(
        t: float,
        params: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Master function: Calculate resonance MUGE g(r,t).
        
        Complete UQFF resonance model (frequency-domain):
        g(r,t) = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res
                 + Ug4i + a_quantum_freq + a_Aether_freq + a_fluid_freq
                 + Osc_term + a_exp_freq + f_TRZ
        
        Suitable for weak-field regimes (nebulae, star-forming regions, diffuse ISM).
        
        Returns all resonance components for validation and analysis.
        """
        if params is None:
            params = {}
        
        # Extract resonance-specific parameters
        c = params.get('c', Source86_Extended.DEFAULT_PARAMS['c'])
        Evac_ISM = params.get('Evac_ISM', Source86_Extended.DEFAULT_PARAMS['Evac_ISM'])
        Evac_neb = params.get('Evac_neb', Source86_Extended.DEFAULT_PARAMS['Evac_neb'])
        Delta_Evac = params.get('Delta_Evac', Source86_Extended.DEFAULT_PARAMS['Delta_Evac'])
        V = params.get('V', Source86_Extended.DEFAULT_PARAMS['V'])
        v_exp = params.get('v_exp', Source86_Extended.DEFAULT_PARAMS['v_exp'])
        f_THz = params.get('f_THz', Source86_Extended.DEFAULT_PARAMS['f_THz'])
        F_super = params.get('F_super', Source86_Extended.DEFAULT_PARAMS['F_super'])
        UA_SCm = params.get('UA_SCm', Source86_Extended.DEFAULT_PARAMS['UA_SCm'])
        omega_i = params.get('omega_i', Source86_Extended.DEFAULT_PARAMS['omega_i'])
        k4 = params.get('k4', Source86_Extended.DEFAULT_PARAMS['k4'])
        E_react = params.get('E_react', Source86_Extended.DEFAULT_PARAMS['E_react'])
        f_react = params.get('f_react', Source86_Extended.DEFAULT_PARAMS['f_react'])
        f_quantum = params.get('f_quantum', Source86_Extended.DEFAULT_PARAMS['f_quantum'])
        f_Aether = params.get('f_Aether', Source86_Extended.DEFAULT_PARAMS['f_Aether'])
        f_fluid = params.get('f_fluid', Source86_Extended.DEFAULT_PARAMS['f_fluid'])
        f_exp = params.get('f_exp', Source86_Extended.DEFAULT_PARAMS['f_exp'])
        f_TRZ = params.get('f_TRZ', Source86_Extended.DEFAULT_PARAMS['f_TRZ'])
        
        # 1. Base resonance term
        adpm = Source86_Extended.calculate_adpm(**params)
        
        # 2. THz frequency resonance
        athz = Source86_Extended.calculate_athz(adpm, **params)
        
        # 3. Vacuum difference term
        avac_diff = (Evac_neb / (c * c)) * Delta_Evac * (v_exp ** 2) * adpm
        
        # 4. Super frequency term
        asuper_freq = (Evac_neb / c) * F_super * f_THz * adpm
        
        # 5. Aether resonance
        aaether_res = UA_SCm * omega_i * f_THz * adpm * (1 + f_TRZ)
        
        # 6. Ug4i reactor term
        ug4i = k4 * E_react * f_react * adpm / (Evac_neb * c)
        
        # 7. Quantum frequency
        aquantum_freq = (Evac_ISM / c) * f_quantum * Evac_neb * adpm
        
        # 8. Aether frequency
        aaether_freq = (Evac_ISM / c) * f_Aether * Evac_neb * adpm
        
        # 9. Fluid frequency
        afluid_freq = (Evac_ISM / c) * f_fluid * Evac_neb * V
        
        # 10. Oscillation term
        osc_term = Source86_Extended.calculate_osc_term(t, **params)
        
        # 11. Expansion frequency
        aexp_freq = (Evac_ISM / c) * f_exp * Evac_neb * adpm
        
        # 12. f_TRZ factor
        ftrz = f_TRZ
        
        # Total resonance gravity
        g_total = adpm + athz + avac_diff + asuper_freq + aaether_res + ug4i + aquantum_freq + aaether_freq + afluid_freq + osc_term + aexp_freq + ftrz
        
        return {
            'g_total': g_total,
            'adpm': adpm,
            'athz': athz,
            'avac_diff': avac_diff,
            'asuper_freq': asuper_freq,
            'aaether_res': aaether_res,
            'ug4i': ug4i,
            'aquantum_freq': aquantum_freq,
            'aaether_freq': aaether_freq,
            'afluid_freq': afluid_freq,
            'osc_term': osc_term,
            'aexp_freq': aexp_freq,
            'ftrz': ftrz,
        }
