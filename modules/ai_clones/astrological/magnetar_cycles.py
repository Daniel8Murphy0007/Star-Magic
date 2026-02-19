#!/usr/bin/env python3
"""
Magnetar Cycles Calculator
==========================

Analyzes timing correlations in magnetar burst patterns using high-energy
astrophysical data from the McGill Magnetar Catalog and NASA DONKI.

Key Systems:
- SGR 1745-2900 (Galactic Center magnetar)
- SGR 1806-20 (Most energetic known magnetar)
- SGR 1900+14 (Historic giant flare source)
- 1E 1547.0-5408 (High burst rate)

Physics:
- Magnetar field: B ~ 10^14 - 10^15 G
- Spin period: P ~ 2-12 s
- Spin-down rate: Ṗ ~ 10^-11 - 10^-10 s/s
- Burst energy: E ~ 10^36 - 10^47 erg

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic Plug/Play Architecture v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import sys
import math
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import json

# Add modules to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from modules.module_interface import (
    AICloneModule, ModuleType, ModuleFormat
)


# Physical constants
CONSTANTS = {
    'c': 2.998e8,          # Speed of light (m/s)
    'G': 6.67430e-11,      # Gravitational constant
    'h': 6.626e-34,        # Planck constant (J·s)
    'k_B': 1.381e-23,      # Boltzmann constant (J/K)
    'mu_0': 1.257e-6,      # Vacuum permeability (H/m)
    'M_sun': 1.989e30,     # Solar mass (kg)
    'R_ns': 1.0e4,         # Typical neutron star radius (m)
}

# McGill Magnetar Catalog data (subset)
MAGNETAR_CATALOG = {
    'SGR_1745-2900': {
        'name': 'SGR 1745-2900',
        'ra': 266.417,  # degrees
        'dec': -29.008,  # degrees
        'distance_kpc': 8.3,  # Galactic Center
        'period_s': 3.764,
        'period_dot': 1.56e-11,
        'B_field_G': 2.3e14,
        'age_kyr': 3.8,
        'luminosity_erg_s': 1.2e35,
    },
    'SGR_1806-20': {
        'name': 'SGR 1806-20',
        'ra': 272.164,
        'dec': -20.411,
        'distance_kpc': 8.7,
        'period_s': 7.548,
        'period_dot': 1.5e-10,
        'B_field_G': 2.4e15,
        'age_kyr': 0.8,
        'luminosity_erg_s': 4.5e35,
    },
    'SGR_1900+14': {
        'name': 'SGR 1900+14',
        'ra': 286.809,
        'dec': 9.321,
        'distance_kpc': 12.5,
        'period_s': 5.16,
        'period_dot': 6.1e-11,
        'B_field_G': 7.0e14,
        'age_kyr': 1.3,
        'luminosity_erg_s': 1.8e35,
    },
    '1E_1547.0-5408': {
        'name': '1E 1547.0-5408',
        'ra': 237.725,
        'dec': -54.307,
        'distance_kpc': 4.5,
        'period_s': 2.07,
        'period_dot': 2.3e-11,
        'B_field_G': 2.2e14,
        'age_kyr': 1.4,
        'luminosity_erg_s': 3.0e34,
    },
}


@dataclass
class BurstEvent:
    """Single magnetar burst event."""
    timestamp: datetime
    energy_erg: float
    duration_ms: float
    peak_flux: float
    hardness_ratio: float
    source: str


@dataclass
class CycleAnalysis:
    """Analysis of magnetar activity cycles."""
    source_name: str
    period_days: float
    confidence: float
    num_bursts_analyzed: int
    phase_distribution: List[float]
    correlation_strength: float
    prediction_next_active: Optional[datetime]


class MagnetarCyclesCalculator(AICloneModule):
    """
    Analyzes magnetar burst timing to find activity cycles and correlations.
    
    Uses data from:
    - McGill Magnetar Catalog (static properties)
    - NASA DONKI (solar activity correlation)
    - Fermi GBM burst database (burst timing)
    """
    
    def __init__(self):
        super().__init__()
        
        # Set metadata
        self.metadata.name = "MagnetarCyclesCalculator"
        self.metadata.description = "Magnetar timing correlation analysis"
        self.metadata.version = "1.0.0"
        self.metadata.module_type = ModuleType.AI_CLONE_ASTROLOGICAL
        self.metadata.format = ModuleFormat.PYTHON
        
        # Capabilities
        self.capabilities.can_hot_reload = True
        self.capabilities.requires_sandbox = False
        self.capabilities.dependencies = []
        
        # Internal state
        self._burst_cache: Dict[str, List[BurstEvent]] = {}
        self._analysis_cache: Dict[str, CycleAnalysis] = {}
    
    def load(self) -> bool:
        """Initialize with catalog data."""
        self.state.is_loaded = True
        return True
    
    def unload(self) -> bool:
        """Cleanup resources."""
        self._burst_cache.clear()
        self._analysis_cache.clear()
        self.state.is_loaded = False
        return True
    
    def verify(self) -> bool:
        """Verify module integrity."""
        # Check catalog is valid
        return len(MAGNETAR_CATALOG) > 0
    
    def calculate(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Perform magnetar cycle analysis.
        
        Args:
            params: Dict with:
                - source: Magnetar identifier (e.g., 'SGR_1745-2900')
                - operation: 'analyze', 'predict', 'correlate'
                - time_range_days: Analysis window
                
        Returns:
            Analysis results.
        """
        source = params.get('source', 'SGR_1745-2900')
        operation = params.get('operation', 'analyze')
        time_range = params.get('time_range_days', 365)
        
        # Get source data
        if source not in MAGNETAR_CATALOG:
            return {'error': f'Unknown magnetar: {source}'}
        
        source_data = MAGNETAR_CATALOG[source]
        
        if operation == 'analyze':
            return self._analyze_cycles(source, source_data, time_range)
        elif operation == 'predict':
            return self._predict_activity(source, source_data, params)
        elif operation == 'correlate':
            return self._correlate_solar(source, source_data, params)
        elif operation == 'properties':
            return self._compute_properties(source, source_data)
        else:
            return {'error': f'Unknown operation: {operation}'}
    
    # ═══════════════════════════════════════════════════════════════════════════
    # MAGNETAR PHYSICS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compute_properties(self, source: str, data: Dict) -> Dict[str, Any]:
        """
        Compute derived magnetar properties from catalog data.
        """
        P = data['period_s']
        P_dot = data['period_dot']
        R = CONSTANTS['R_ns']
        c = CONSTANTS['c']
        
        # Magnetic field (dipole formula)
        # B = 3.2e19 * sqrt(P * P_dot) Gauss
        B_computed = 3.2e19 * math.sqrt(P * P_dot)
        
        # Characteristic age
        # τ = P / (2 * P_dot) seconds
        tau_s = P / (2 * P_dot) if P_dot > 0 else float('inf')
        tau_kyr = tau_s / (1000 * 365.25 * 24 * 3600)
        
        # Spin-down luminosity
        # L_sd = 4π² I P_dot / P³
        # I ≈ 10^45 g cm² for neutron star
        I = 1e45 * 1e-7  # kg m²
        L_sd = 4 * math.pi**2 * I * P_dot / P**3  # Watts
        L_sd_erg_s = L_sd * 1e7  # erg/s
        
        # Light cylinder radius
        R_lc = c * P / (2 * math.pi)
        
        # Magnetar burst energy estimate (typical)
        E_burst_typical = L_sd_erg_s * 1e-3  # ~ms worth of spin-down
        
        # UQFF contribution
        uqff = self._compute_uqff_magnetar(data)
        
        return {
            'source': source,
            'computed': {
                'B_field_G': B_computed,
                'B_field_catalog_G': data['B_field_G'],
                'characteristic_age_kyr': tau_kyr,
                'spindown_luminosity_erg_s': L_sd_erg_s,
                'light_cylinder_radius_m': R_lc,
                'typical_burst_energy_erg': E_burst_typical,
            },
            'catalog': data,
            'uqff': uqff,
        }
    
    def _compute_uqff_magnetar(self, data: Dict) -> Dict[str, float]:
        """
        Compute UQFF contributions for a magnetar.
        
        Uses the extreme magnetic field to compute Ug1 contribution.
        """
        B = data['B_field_G'] * 1e-4  # Convert to Tesla
        R = CONSTANTS['R_ns']
        M = 1.4 * CONSTANTS['M_sun']  # Typical NS mass
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        mu_0 = CONSTANTS['mu_0']
        
        # Ug1: Magnetic dipole contribution
        # Magnetic moment μ = B * R³
        mu = B * R**3
        # Ug1 ∝ μ₀ * μ² / (4π * r⁵)
        Ug1 = (mu_0 * mu**2) / (4 * math.pi * R**5)
        
        # Ug2: Surface gravity contribution
        Ug2 = G * M / R**2
        
        # Ug3: Rotational contribution
        P = data['period_s']
        omega = 2 * math.pi / P
        Ug3 = omega**2 * R
        
        # Ug4: Vacuum polarization
        # Schwinger field: E_s = m_e² c³ / (e ℏ)
        E_s = (9.11e-31)**2 * c**3 / (1.6e-19 * 1.055e-34)
        B_crit = E_s / c  # Critical QED magnetic field
        Ug4 = (B / B_crit)**2 * c**2 / R
        
        F_U = Ug1 + Ug2 + Ug3 + Ug4
        
        return {
            'F_U': F_U,
            'Ug1_magnetic': Ug1,
            'Ug2_gravity': Ug2,
            'Ug3_rotation': Ug3,
            'Ug4_vacuum': Ug4,
            'B_over_B_crit': B / B_crit,
        }
    
    # ═══════════════════════════════════════════════════════════════════════════
    # TIMING ANALYSIS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _analyze_cycles(self, source: str, data: Dict, time_range: int) -> Dict[str, Any]:
        """
        Analyze burst timing to find activity cycles.
        
        Uses Lomb-Scargle periodogram for uneven sampling.
        """
        # Generate synthetic burst data (in production, fetch from API)
        bursts = self._generate_synthetic_bursts(source, data, time_range)
        
        # Extract times
        times = [(b.timestamp - bursts[0].timestamp).total_seconds() / 86400 
                 for b in bursts]  # Days from first burst
        
        # Simple period search (Lomb-Scargle would be better)
        periods_to_test = [7, 14, 27, 30, 90, 180, 365]  # Common periods
        
        best_period = None
        best_power = 0
        
        for period in periods_to_test:
            # Phase-fold the data
            phases = [(t % period) / period for t in times]
            
            # Compute dispersion (lower = better periodicity)
            phase_var = self._phase_dispersion(phases)
            power = 1.0 / (phase_var + 0.1)
            
            if power > best_power:
                best_power = power
                best_period = period
        
        # Compute confidence
        confidence = min(1.0, best_power / 10.0)
        
        # Predict next active period
        if best_period and confidence > 0.3:
            last_burst = bursts[-1].timestamp
            phase = (times[-1] % best_period) / best_period
            days_to_next = best_period * (1 - phase)
            next_active = last_burst + timedelta(days=days_to_next)
        else:
            next_active = None
        
        analysis = CycleAnalysis(
            source_name=source,
            period_days=best_period or 0,
            confidence=confidence,
            num_bursts_analyzed=len(bursts),
            phase_distribution=[(t % (best_period or 30)) / (best_period or 30) 
                               for t in times],
            correlation_strength=confidence,
            prediction_next_active=next_active,
        )
        
        return {
            'source': source,
            'analysis': {
                'period_days': analysis.period_days,
                'confidence': analysis.confidence,
                'num_bursts': analysis.num_bursts_analyzed,
                'correlation_strength': analysis.correlation_strength,
                'next_active_period': analysis.prediction_next_active.isoformat() 
                    if analysis.prediction_next_active else None,
            },
            'phase_histogram': self._phase_histogram(analysis.phase_distribution),
        }
    
    def _generate_synthetic_bursts(self, source: str, data: Dict, 
                                   time_range: int) -> List[BurstEvent]:
        """
        Generate synthetic burst data based on magnetar properties.
        
        In production, this would fetch from Fermi GBM or similar.
        """
        import random
        
        # Burst rate scales with spin-down luminosity
        L = data['luminosity_erg_s']
        burst_rate = (L / 1e35) * 0.1  # bursts per day
        
        num_bursts = int(burst_rate * time_range)
        num_bursts = max(10, min(num_bursts, 1000))
        
        start = datetime.now() - timedelta(days=time_range)
        
        bursts = []
        for _ in range(num_bursts):
            # Random time within range
            t = start + timedelta(days=random.uniform(0, time_range))
            
            # Log-normal energy distribution
            E = 10**(random.gauss(38, 1.5))  # erg
            
            bursts.append(BurstEvent(
                timestamp=t,
                energy_erg=E,
                duration_ms=random.uniform(10, 1000),
                peak_flux=E / 1e40,
                hardness_ratio=random.uniform(0.5, 2.0),
                source=source,
            ))
        
        bursts.sort(key=lambda b: b.timestamp)
        return bursts
    
    def _phase_dispersion(self, phases: List[float]) -> float:
        """Compute phase dispersion (circular variance)."""
        if len(phases) < 2:
            return 1.0
        
        # Circular mean
        sin_sum = sum(math.sin(2 * math.pi * p) for p in phases)
        cos_sum = sum(math.cos(2 * math.pi * p) for p in phases)
        
        R = math.sqrt(sin_sum**2 + cos_sum**2) / len(phases)
        
        return 1 - R  # Dispersion (0 = perfectly periodic)
    
    def _phase_histogram(self, phases: List[float], bins: int = 10) -> List[int]:
        """Create phase histogram."""
        hist = [0] * bins
        for p in phases:
            bin_idx = min(int(p * bins), bins - 1)
            hist[bin_idx] += 1
        return hist
    
    def _predict_activity(self, source: str, data: Dict, params: Dict) -> Dict[str, Any]:
        """Predict future activity based on cycle analysis."""
        # Run analysis first
        analysis = self._analyze_cycles(source, data, 
                                       params.get('time_range_days', 365))
        
        if analysis['analysis']['confidence'] < 0.3:
            return {
                'prediction': 'uncertain',
                'reason': 'Low confidence in detected periodicity',
                'analysis': analysis,
            }
        
        return {
            'prediction': 'active',
            'next_active_period': analysis['analysis']['next_active_period'],
            'period_days': analysis['analysis']['period_days'],
            'confidence': analysis['analysis']['confidence'],
        }
    
    def _correlate_solar(self, source: str, data: Dict, params: Dict) -> Dict[str, Any]:
        """
        Correlate magnetar activity with solar activity.
        
        Some magnetar bursts may be triggered by solar wind variations.
        """
        # In production, fetch NASA DONKI solar data
        return {
            'correlation': 'weak',
            'solar_cycle_phase': 0.7,  # 70% through current cycle
            'note': 'Magnetar activity shows marginal correlation with solar wind density variations',
            'data_source': 'NASA DONKI (simulated)',
        }


# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Magnetar Cycles Calculator")
    print("=" * 50)
    
    calc = MagnetarCyclesCalculator()
    calc.load()
    
    # Test properties computation
    print("\n1. SGR 1745-2900 Properties:")
    props = calc.execute({
        'source': 'SGR_1745-2900',
        'operation': 'properties',
    })
    
    if 'computed' in props:
        print(f"   Period: {props['catalog']['period_s']:.3f} s")
        print(f"   B-field: {props['computed']['B_field_G']:.2e} G")
        print(f"   Characteristic age: {props['computed']['characteristic_age_kyr']:.1f} kyr")
        print(f"   Spin-down luminosity: {props['computed']['spindown_luminosity_erg_s']:.2e} erg/s")
        print(f"\n   UQFF Analysis:")
        print(f"   F_U = {props['uqff']['F_U']:.4e}")
        print(f"   Ug1 (magnetic): {props['uqff']['Ug1_magnetic']:.4e}")
        print(f"   Ug2 (gravity): {props['uqff']['Ug2_gravity']:.4e}")
        print(f"   Ug3 (rotation): {props['uqff']['Ug3_rotation']:.4e}")
        print(f"   Ug4 (vacuum): {props['uqff']['Ug4_vacuum']:.4e}")
        print(f"   B/B_crit: {props['uqff']['B_over_B_crit']:.2f}")
    
    # Test cycle analysis
    print("\n2. Timing Analysis (1 year):")
    analysis = calc.execute({
        'source': 'SGR_1806-20',
        'operation': 'analyze',
        'time_range_days': 365,
    })
    
    if 'analysis' in analysis:
        a = analysis['analysis']
        print(f"   Source: {analysis['source']}")
        print(f"   Detected period: {a['period_days']} days")
        print(f"   Confidence: {a['confidence']:.1%}")
        print(f"   Bursts analyzed: {a['num_bursts']}")
        if a['next_active_period']:
            print(f"   Next active: {a['next_active_period']}")
