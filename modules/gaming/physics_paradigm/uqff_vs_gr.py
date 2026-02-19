#!/usr/bin/env python3
"""
Physics Paradigm - UQFF vs GR Comparison
========================================

Interactive comparison between UQFF triadic gravity framework and
standard General Relativity predictions.

Key Comparisons:
- Galaxy rotation curves (dark matter interpretation)
- Gravitational lensing predictions
- Gravitational wave propagation
- Extreme field regimes (magnetars, SMBHs)
- Cosmological expansion

Physics Paradigms:
1. General Relativity (Einstein 1915)
2. UQFF Triadic Gravity (Murphy 2024-2026)
3. MUGE Compressed (Newtonian + corrections)
4. MUGE Resonance (frequency-domain)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic Plug/Play Architecture v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import sys
import math
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field
from enum import Enum
import json

# Add modules to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from modules.module_interface import (
    GamingModule, ModuleType, ModuleFormat, GameFrame
)


# Physical constants
CONSTANTS = {
    'G': 6.67430e-11,      # Gravitational constant (m³/kg/s²)
    'c': 2.998e8,          # Speed of light (m/s)
    'h': 6.626e-34,        # Planck constant (J·s)
    'hbar': 1.055e-34,     # Reduced Planck constant
    'k_B': 1.381e-23,      # Boltzmann constant
    'mu_0': 1.257e-6,      # Vacuum permeability
    'epsilon_0': 8.854e-12, # Vacuum permittivity
    'M_sun': 1.989e30,     # Solar mass (kg)
    'R_sun': 6.96e8,       # Solar radius (m)
    'pc': 3.086e16,        # Parsec (m)
    'kpc': 3.086e19,       # Kiloparsec (m)
    'ly': 9.461e15,        # Light-year (m)
    
    # UQFF calibration constants
    'kappa': 0.0005,       # κ per day
    'SSq': 0.57,           # [SSq] calibration
    'H_SCm': 0.99,         # SCm horizon factor
    'U_UA': 0.0001,        # UA factor
    'k_eta': 1e-113,       # Coupling constant
    'beta_i': 0.603,       # β_i factor
}


class Paradigm(Enum):
    """Physics paradigm selection."""
    GR = "General Relativity"
    UQFF = "UQFF Triadic"
    MUGE_COMPRESSED = "MUGE Compressed"
    MUGE_RESONANCE = "MUGE Resonance"


@dataclass
class ComparisonResult:
    """Result of paradigm comparison."""
    system_name: str
    observable: str
    gr_prediction: float
    uqff_prediction: float
    muge_prediction: float
    observed_value: Optional[float]
    gr_residual: Optional[float]
    uqff_residual: Optional[float]
    muge_residual: Optional[float]
    winner: str
    units: str
    notes: str = ""


@dataclass
class SystemParams:
    """Physical system parameters."""
    name: str
    M: float  # Mass (kg)
    r: float  # Distance/radius (m)
    z: float = 0.0  # Redshift
    B: float = 0.0  # Magnetic field (T)
    omega: float = 0.0  # Angular velocity (rad/s)
    T: float = 0.0  # Temperature (K)


class PhysicsParadigmGame(GamingModule):
    """
    Interactive physics paradigm comparison game.
    
    Players explore different physical systems and compare predictions
    from GR, UQFF, MUGE Compressed, and MUGE Resonance.
    """
    
    def __init__(self):
        super().__init__()
        
        # Set metadata
        self.metadata.name = "PhysicsParadigmGame"
        self.metadata.description = "UQFF vs GR comparison game"
        self.metadata.version = "1.0.0"
        self.metadata.module_type = ModuleType.GAMING_PHYSICS_PARADIGM
        self.metadata.format = ModuleFormat.PYTHON
        
        # Capabilities
        self.capabilities.can_hot_reload = True
        self.capabilities.requires_sandbox = False
        
        # Game state
        self._systems: List[SystemParams] = []
        self._comparisons: List[ComparisonResult] = []
        self._score = {'GR': 0, 'UQFF': 0, 'MUGE': 0}
        self._current_system_idx = 0
    
    def load(self) -> bool:
        """Initialize with preset systems."""
        self._load_preset_systems()
        self.state.is_loaded = True
        return True
    
    def unload(self) -> bool:
        """Cleanup."""
        self._systems.clear()
        self._comparisons.clear()
        self.state.is_loaded = False
        return True
    
    def verify(self) -> bool:
        """Verify module."""
        return len(self._systems) > 0
    
    def _load_preset_systems(self):
        """Load preset astrophysical systems."""
        self._systems = [
            SystemParams(
                name="Sun",
                M=CONSTANTS['M_sun'],
                r=CONSTANTS['R_sun'],
            ),
            SystemParams(
                name="Earth (surface)",
                M=5.972e24,
                r=6.371e6,
            ),
            SystemParams(
                name="SGR 1745-2900 (magnetar)",
                M=1.4 * CONSTANTS['M_sun'],
                r=1e4,
                B=2.3e10,  # Tesla
                omega=2 * math.pi / 3.764,
            ),
            SystemParams(
                name="Sagittarius A*",
                M=4e6 * CONSTANTS['M_sun'],
                r=6e9,  # 6 billion m
            ),
            SystemParams(
                name="M87* (SMBH)",
                M=6.5e9 * CONSTANTS['M_sun'],
                r=1e13,
            ),
            SystemParams(
                name="NGC 1365 (galaxy)",
                M=1e42,
                r=10 * CONSTANTS['kpc'],
                z=0.0055,
            ),
            SystemParams(
                name="NGC 3596 (galaxy)",
                M=3e41,
                r=8 * CONSTANTS['kpc'],
                z=0.0038,
            ),
            SystemParams(
                name="Neutron Star Binary",
                M=2.8 * CONSTANTS['M_sun'],
                r=1e5,
                omega=100,
            ),
        ]
    
    def update(self, dt: float, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Update game state."""
        if 'select_system' in inputs:
            self._current_system_idx = inputs['select_system']
        
        if 'compare' in inputs:
            observable = inputs.get('observable', 'gravity')
            result = self._compare_paradigms(
                self._systems[self._current_system_idx],
                observable
            )
            self._comparisons.append(result)
            self._update_score(result)
        
        return {
            'current_system': self._systems[self._current_system_idx].name,
            'score': self._score,
            'comparisons': len(self._comparisons),
        }
    
    def render(self) -> GameFrame:
        """Render current frame."""
        frame = GameFrame()
        
        system = self._systems[self._current_system_idx]
        
        frame.entities = {
            'systems': [s.name for s in self._systems],
            'current': self._current_system_idx,
            'current_params': {
                'name': system.name,
                'M': system.M,
                'r': system.r,
                'B': system.B,
                'omega': system.omega,
            },
        }
        
        frame.ui_state = {
            'score': self._score,
            'comparisons': len(self._comparisons),
            'last_winner': self._comparisons[-1].winner if self._comparisons else "None",
        }
        
        return frame
    
    def execute(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute comparison operations.
        
        Args:
            params: Operation parameters with 'operation' key.
        """
        operation = params.get('operation', 'info')
        
        if operation == 'compare':
            system_idx = params.get('system_idx', 0)
            observable = params.get('observable', 'gravity')
            system = self._systems[system_idx]
            result = self._compare_paradigms(system, observable)
            return {
                'comparison': {
                    'system': result.system_name,
                    'observable': result.observable,
                    'GR': result.gr_prediction,
                    'UQFF': result.uqff_prediction,
                    'MUGE': result.muge_prediction,
                    'observed': result.observed_value,
                    'winner': result.winner,
                    'units': result.units,
                }
            }
        
        elif operation == 'compare_all':
            observable = params.get('observable', 'gravity')
            results = []
            for i, system in enumerate(self._systems):
                result = self._compare_paradigms(system, observable)
                results.append({
                    'system': result.system_name,
                    'GR': result.gr_prediction,
                    'UQFF': result.uqff_prediction,
                    'MUGE': result.muge_prediction,
                    'winner': result.winner,
                })
            return {'comparisons': results}
        
        elif operation == 'systems':
            return {'systems': [
                {'idx': i, 'name': s.name, 'M': s.M, 'r': s.r}
                for i, s in enumerate(self._systems)
            ]}
        
        elif operation == 'score':
            return {'score': self._score}
        
        else:
            return {'error': f'Unknown operation: {operation}'}
    
    # ═══════════════════════════════════════════════════════════════════════════
    # PARADIGM CALCULATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compare_paradigms(self, system: SystemParams, 
                          observable: str) -> ComparisonResult:
        """
        Compare paradigm predictions for a given observable.
        """
        if observable == 'gravity':
            return self._compare_gravity(system)
        elif observable == 'rotation':
            return self._compare_rotation_curve(system)
        elif observable == 'lensing':
            return self._compare_lensing(system)
        elif observable == 'redshift':
            return self._compare_gravitational_redshift(system)
        elif observable == 'frame_drag':
            return self._compare_frame_dragging(system)
        else:
            return self._compare_gravity(system)
    
    def _compare_gravity(self, system: SystemParams) -> ComparisonResult:
        """Compare gravitational acceleration predictions."""
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        M = system.M
        r = system.r
        
        # GR: Schwarzschild surface gravity
        # g = GM/r² * (1 - r_s/r)^(-1/2)
        r_s = 2 * G * M / c**2
        if r > r_s:
            gr_factor = 1 / math.sqrt(1 - r_s / r)
        else:
            gr_factor = float('inf')
        g_gr = G * M / r**2 * gr_factor
        
        # UQFF Triadic: F_U = Ug1 + Ug2 + Ug3 + Ug4
        g_uqff = self._compute_uqff_gravity(system)
        
        # MUGE Compressed
        g_muge = self._compute_muge_compressed(system)
        
        # Winner (closest to Newtonian for everyday systems)
        g_newton = G * M / r**2
        
        return ComparisonResult(
            system_name=system.name,
            observable="Surface gravity",
            gr_prediction=g_gr,
            uqff_prediction=g_uqff,
            muge_prediction=g_muge,
            observed_value=g_newton,
            gr_residual=abs(g_gr - g_newton) / g_newton * 100,
            uqff_residual=abs(g_uqff - g_newton) / g_newton * 100,
            muge_residual=abs(g_muge - g_newton) / g_newton * 100,
            winner=self._determine_winner(g_gr, g_uqff, g_muge, g_newton),
            units="m/s²",
        )
    
    def _compute_uqff_gravity(self, system: SystemParams) -> float:
        """
        Compute UQFF triadic gravity.
        
        F_U = Ug1 + Ug2 + Ug3 + Ug4
        
        Ug1: Magnetic dipole contribution
        Ug2: Charge-reactivity 
        Ug3: String rotation
        Ug4: Vacuum concentration
        """
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        M = system.M
        r = system.r
        B = system.B
        omega = system.omega
        mu_0 = CONSTANTS['mu_0']
        
        # Base Newtonian
        g_base = G * M / r**2
        
        # Ug1: Magnetic contribution
        if B > 0:
            mu = B * r**3  # Magnetic moment
            Ug1 = (mu_0 * mu**2) / (4 * math.pi * r**5) / M
        else:
            Ug1 = 0
        
        # Ug2: Electrostatic/charge-reactivity (~0 for neutral matter)
        Ug2 = 0
        
        # Ug3: Rotational contribution
        if omega > 0:
            Ug3 = omega**2 * r * 0.001  # Small correction factor
        else:
            Ug3 = 0
        
        # Ug4: Vacuum polarization (significant near extreme objects)
        # B_crit = m_e² c³ / (e ℏ) ≈ 4.4e9 T
        B_crit = 4.4e9
        if B > 0:
            Ug4 = (B / B_crit)**2 * g_base * 0.01
        else:
            # Dark energy contribution
            Lambda = 1.1e-52  # m^-2 (cosmological constant)
            Ug4 = -Lambda * c**2 * r / 3  # Anti-gravity at large scales
        
        Ubi = g_base * CONSTANTS['beta_i']  # Buoyancy
        
        F_U = g_base + Ug1 + Ug2 + Ug3 + Ug4
        
        return F_U
    
    def _compute_muge_compressed(self, system: SystemParams) -> float:
        """
        Compute MUGE Compressed gravity.
        
        g = g_Newton * (1 + Σ corrections)
        
        Corrections:
        - Hubble expansion
        - Magnetic suppression
        - Envelope
        - Cosmological Λ
        - Quantum ℏ
        """
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        M = system.M
        r = system.r
        
        g_newton = G * M / r**2
        
        # Hubble correction (expansion at large scales)
        H_0 = 70e3 / 3.086e22  # Hubble constant in s^-1
        delta_H = H_0**2 * r / g_newton if g_newton > 0 else 0
        
        # Magnetic suppression
        if system.B > 0:
            delta_B = -0.001 * (system.B / 4.4e9)**2
        else:
            delta_B = 0
        
        # Quantum correction (negligible except at Planck scale)
        l_P = 1.616e-35  # Planck length
        delta_Q = (l_P / r)**2 if r > l_P else 1
        
        # Cosmological constant
        Lambda = 1.1e-52
        delta_Lambda = -Lambda * c**2 * r**2 / (3 * g_newton) if g_newton > 0 else 0
        
        # Total
        total_correction = 1 + delta_H + delta_B + delta_Q + delta_Lambda
        total_correction = max(0.9, min(1.1, total_correction))  # Clamp
        
        return g_newton * total_correction
    
    def _compare_rotation_curve(self, system: SystemParams) -> ComparisonResult:
        """Compare galaxy rotation curve predictions."""
        G = CONSTANTS['G']
        M = system.M
        r = system.r
        
        # GR/Keplerian: v = sqrt(GM/r)
        v_gr = math.sqrt(G * M / r)
        
        # UQFF: Includes dark matter via UQFF triadic
        # v_UQFF ~ v_Keplerian * (1 + f(r/r_s))
        r_s = r / 10  # Scale radius
        v_uqff = v_gr * math.sqrt(1 + math.log(1 + r / r_s) - r / (r + r_s))
        
        # MUGE: Similar with NFW profile
        c = 10  # Concentration
        x = r / r_s
        f_x = math.log(1 + x) - x / (1 + x)
        v_muge = v_gr * math.sqrt(1 + 5 * f_x / x)  # 5x dark matter factor
        
        # Observed flat rotation (~200 km/s for Milky Way-like)
        v_obs = 200e3 if "galaxy" in system.name.lower() else v_gr
        
        return ComparisonResult(
            system_name=system.name,
            observable="Rotation velocity",
            gr_prediction=v_gr,
            uqff_prediction=v_uqff,
            muge_prediction=v_muge,
            observed_value=v_obs,
            gr_residual=abs(v_gr - v_obs) / v_obs * 100 if v_obs > 0 else 0,
            uqff_residual=abs(v_uqff - v_obs) / v_obs * 100 if v_obs > 0 else 0,
            muge_residual=abs(v_muge - v_obs) / v_obs * 100 if v_obs > 0 else 0,
            winner=self._determine_winner(v_gr, v_uqff, v_muge, v_obs),
            units="m/s",
            notes="Flat rotation curves suggest dark matter or modified gravity",
        )
    
    def _compare_lensing(self, system: SystemParams) -> ComparisonResult:
        """Compare gravitational lensing predictions."""
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        M = system.M
        b = system.r  # Impact parameter
        
        # GR Einstein deflection: α = 4GM/(c²b)
        alpha_gr = 4 * G * M / (c**2 * b)
        
        # UQFF: Slightly modified deflection
        # Additional contribution from Ug4 (vacuum polarization)
        alpha_uqff = alpha_gr * (1 + 0.001 * CONSTANTS['SSq'])
        
        # MUGE: Same as GR for lensing
        alpha_muge = alpha_gr
        
        # Observed (use GR as reference)
        alpha_obs = alpha_gr
        
        return ComparisonResult(
            system_name=system.name,
            observable="Deflection angle",
            gr_prediction=alpha_gr,
            uqff_prediction=alpha_uqff,
            muge_prediction=alpha_muge,
            observed_value=alpha_obs,
            gr_residual=0,
            uqff_residual=abs(alpha_uqff - alpha_obs) / alpha_obs * 100,
            muge_residual=0,
            winner="GR",
            units="radians",
        )
    
    def _compare_gravitational_redshift(self, system: SystemParams) -> ComparisonResult:
        """Compare gravitational redshift predictions."""
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        M = system.M
        r = system.r
        
        # GR: z = 1/sqrt(1 - r_s/r) - 1
        r_s = 2 * G * M / c**2
        if r > r_s:
            z_gr = 1 / math.sqrt(1 - r_s / r) - 1
        else:
            z_gr = float('inf')
        
        # UQFF: Small additional contribution
        z_uqff = z_gr * (1 + CONSTANTS['kappa'] * 1e-6)
        
        # MUGE:
        z_muge = z_gr
        
        return ComparisonResult(
            system_name=system.name,
            observable="Gravitational redshift",
            gr_prediction=z_gr,
            uqff_prediction=z_uqff,
            muge_prediction=z_muge,
            observed_value=z_gr,
            gr_residual=0,
            uqff_residual=abs(z_uqff - z_gr) / max(z_gr, 1e-10) * 100,
            muge_residual=0,
            winner="GR",
            units="Δλ/λ",
        )
    
    def _compare_frame_dragging(self, system: SystemParams) -> ComparisonResult:
        """Compare Lense-Thirring precession predictions."""
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        M = system.M
        r = system.r
        omega = system.omega
        
        # GR: Ω_LT = 2GJ/(c²r³) where J = Iω ~ Mr²ω
        I = 0.4 * M * r**2  # Moment of inertia (sphere)
        J = I * omega
        Omega_gr = 2 * G * J / (c**2 * r**3) if omega > 0 else 0
        
        # UQFF: Enhanced by Ug3 rotational term
        Omega_uqff = Omega_gr * (1 + CONSTANTS['beta_i'] * 0.01)
        
        # MUGE:
        Omega_muge = Omega_gr
        
        return ComparisonResult(
            system_name=system.name,
            observable="Frame dragging rate",
            gr_prediction=Omega_gr,
            uqff_prediction=Omega_uqff,
            muge_prediction=Omega_muge,
            observed_value=Omega_gr,
            gr_residual=0,
            uqff_residual=abs(Omega_uqff - Omega_gr) / max(Omega_gr, 1e-20) * 100,
            muge_residual=0,
            winner="GR",
            units="rad/s",
        )
    
    def _determine_winner(self, gr: float, uqff: float, muge: float, 
                         observed: float) -> str:
        """Determine which paradigm best matches observation."""
        if observed is None or observed == 0:
            return "N/A"
        
        residuals = {
            'GR': abs(gr - observed) / abs(observed),
            'UQFF': abs(uqff - observed) / abs(observed),
            'MUGE': abs(muge - observed) / abs(observed),
        }
        
        return min(residuals, key=residuals.get)
    
    def _update_score(self, result: ComparisonResult):
        """Update score based on comparison result."""
        if result.winner in self._score:
            self._score[result.winner] += 1


# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Physics Paradigm - UQFF vs GR Comparison")
    print("=" * 60)
    
    game = PhysicsParadigmGame()
    game.load()
    
    # List systems
    systems = game.execute({'operation': 'systems'})
    print("\nAvailable Systems:")
    for s in systems['systems']:
        print(f"  {s['idx']}: {s['name']}")
    
    # Compare gravity for all systems
    print("\n" + "=" * 60)
    print("Gravity Comparison")
    print("=" * 60)
    
    for i in range(len(systems['systems'])):
        result = game.execute({
            'operation': 'compare',
            'system_idx': i,
            'observable': 'gravity',
        })
        
        c = result['comparison']
        print(f"\n{c['system']}:")
        print(f"  GR:   {c['GR']:.4e} {c['units']}")
        print(f"  UQFF: {c['UQFF']:.4e} {c['units']}")
        print(f"  MUGE: {c['MUGE']:.4e} {c['units']}")
        print(f"  Winner: {c['winner']}")
    
    # Rotation curves
    print("\n" + "=" * 60)
    print("Rotation Curve Comparison (galaxies)")
    print("=" * 60)
    
    for i, s in enumerate(systems['systems']):
        if 'galaxy' in s['name'].lower():
            result = game.execute({
                'operation': 'compare',
                'system_idx': i,
                'observable': 'rotation',
            })
            
            c = result['comparison']
            print(f"\n{c['system']}:")
            print(f"  GR:   {c['GR']/1e3:.1f} km/s")
            print(f"  UQFF: {c['UQFF']/1e3:.1f} km/s")
            print(f"  MUGE: {c['MUGE']/1e3:.1f} km/s")
            print(f"  Obs:  {c['observed']/1e3:.1f} km/s")
            print(f"  Winner: {c['winner']}")
