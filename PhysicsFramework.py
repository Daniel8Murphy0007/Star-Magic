"""
PhysicsFramework.py - Self-Expanding Physics Term System

Python implementation of the C++ PhysicsTerm framework from source70.cpp.
Provides dynamic parameter updates, runtime term registration, metadata tracking,
adaptive learning, and validation for all UQFF physics modules.

Architecture mirrors C++ implementation:
- PhysicsTerm abstract base class (ABC)
- Concrete term implementations (DynamicVacuumTerm, QuantumCouplingTerm)
- Runtime extensibility via dynamicTerms list
- Parameter mutation via dynamicParameters dict
- Provenance tracking via metadata dict

Usage:
    from PhysicsFramework import PhysicsTerm, DynamicVacuumTerm
    
    # Create physics term
    term = DynamicVacuumTerm(amplitude=1e-10, frequency=1e-15)
    
    # Runtime parameter updates
    term.set_dynamic_parameter('amplitude', 2e-10)
    
    # Add nested terms
    term.register_dynamic_term(QuantumCouplingTerm())
    
    # Compute with metadata
    result = term.compute(t=1e6, params={'M': 1e30, 'r': 1e4})
    
    # Export state
    term.export_state('term_state.json')

Author: Daniel T. Murphy
Date: February 14, 2026
Version: 2.0-Enhanced (Python port from C++)
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any
import json
import math
from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell

import logging
from datetime import datetime

# Physical constants (matches QCalc.py)
G = 6.6743e-11          # m^3 kg^-1 s^-2
c = 2.998e8             # m/s
hbar = 1.0546e-34       # J·s
k_B = 1.381e-23         # J/K
M_sun = 1.989e30        # kg


class PhysicsTerm(ABC):
    """
    Abstract base class for all physics terms in the self-expanding framework.
    
    Mirrors C++ PhysicsTerm class from source70.cpp lines 40-58.
    All derived classes must implement compute(), getName(), getDescription().
    
    Features:
    - Dynamic parameter updates at runtime
    - Nested term registration (composition)
    - Metadata tracking for provenance
    - Optional logging for debugging
    - Validation before computation
    - Adaptive learning rate for optimization
    
    Attributes:
        dynamicParameters (dict): Runtime-mutable parameters
        dynamicTerms (list): Nested PhysicsTerm objects
        metadata (dict): Provenance tracking (version, author, timestamp)
        enableDynamicTerms (bool): Toggle nested term contribution
        enableLogging (bool): Debug logging flag
        learningRate (float): Adaptive parameter optimization rate
    """
    
    def __init__(self):
        """Initialize self-expanding framework members."""
        self.dynamicParameters: Dict[str, float] = {}
        self.dynamicTerms: List[PhysicsTerm] = []
        self.metadata: Dict[str, str] = {
            'version': '2.0-Enhanced',
            'created': datetime.now().isoformat(),
            'framework': 'PhysicsFramework.py'
        }
        self.enableDynamicTerms: bool = False  # Disabled by default (additive to core)
        self.enableLogging: bool = False
        self.learningRate: float = 0.001
        
        # Setup logging
        self.logger = logging.getLogger(self.__class__.__name__)
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            ))
            self.logger.addHandler(handler)
            self.logger.setLevel(logging.INFO)
    
    @abstractmethod
    def compute(self, t: float, params: Dict[str, float]) -> float:
        """
        Compute physics term value at time t with given parameters.
        
        Args:
            t: Time coordinate (seconds or dimensionless)
            params: Dictionary of physical parameters
                   Example: {'M': 1e30, 'r': 1e4, 'z': 0.002}
        
        Returns:
            Computed term value (units depend on term type)
        """
        pass
    
    @abstractmethod
    def getName(self) -> str:
        """Return unique term identifier."""
        pass
    
    @abstractmethod
    def getDescription(self) -> str:
        """Return human-readable physics description."""
        pass
    
    def validate(self, params: Dict[str, float]) -> bool:
        """
        Validate parameters before computation.
        
        Override in derived classes for custom validation logic.
        Default: always returns True (no validation).
        
        Args:
            params: Parameter dictionary to validate
        
        Returns:
            True if valid, False otherwise
        """
        return True
    
    def set_dynamic_parameter(self, name: str, value: float) -> None:
        """
        Set or update a dynamic parameter at runtime.
        
        Args:
            name: Parameter identifier
            value: New parameter value
        """
        if self.enableLogging:
            self.logger.info(f"{self.getName()}: Set {name} = {value}")
        self.dynamicParameters[name] = value
    
    def get_dynamic_parameter(self, name: str, default: float = 0.0) -> float:
        """
        Get dynamic parameter value.
        
        Args:
            name: Parameter identifier
            default: Return value if parameter not found
        
        Returns:
            Parameter value or default
        """
        return self.dynamicParameters.get(name, default)
    
    def register_dynamic_term(self, term: 'PhysicsTerm') -> None:
        """
        Register nested physics term for additive contribution.
        
        Nested terms are computed and added to base term if enableDynamicTerms=True.
        
        Args:
            term: PhysicsTerm instance to register
        """
        if self.enableLogging:
            self.logger.info(f"{self.getName()}: Registered {term.getName()}")
        self.dynamicTerms.append(term)
    
    def compute_with_dynamic_terms(self, t: float, params: Dict[str, float]) -> float:
        """
        Compute term value including nested dynamic terms.
        
        Base term value + sum of all registered dynamic terms (if enabled).
        
        Args:
            t: Time coordinate
            params: Parameter dictionary
        
        Returns:
            Total value (base + dynamic contributions)
        """
        # Base term
        base_value = self.compute(t, params)
        
        # Add dynamic terms if enabled
        if self.enableDynamicTerms:
            for term in self.dynamicTerms:
                if term.validate(params):
                    base_value += term.compute(t, params)
                    if self.enableLogging:
                        self.logger.debug(f"  + {term.getName()}: {term.compute(t, params):.6e}")
        
        return base_value
    
    def set_metadata(self, key: str, value: str) -> None:
        """Set metadata field for provenance tracking."""
        self.metadata[key] = value
    
    def get_metadata(self, key: str) -> Optional[str]:
        """Get metadata field value."""
        return self.metadata.get(key)
    
    def export_state(self, filepath: str) -> None:
        """
        Export term state to JSON file.
        
        Includes: dynamicParameters, metadata, configuration flags.
        Does NOT include dynamicTerms (circular reference issues).
        
        Args:
            filepath: Output JSON file path
        """
        state = {
            'class': self.__class__.__name__,
            'name': self.getName(),
            'description': self.getDescription(),
            'dynamicParameters': self.dynamicParameters,
            'metadata': self.metadata,
            'enableDynamicTerms': self.enableDynamicTerms,
            'enableLogging': self.enableLogging,
            'learningRate': self.learningRate,
            'num_dynamic_terms': len(self.dynamicTerms)
        }
        
        with open(filepath, 'w') as f:
            json.dump(state, f, indent=2)
        
        if self.enableLogging:
            self.logger.info(f"State exported to {filepath}")
    
    def import_state(self, filepath: str) -> None:
        """
        Import term state from JSON file.
        
        Restores: dynamicParameters, metadata, configuration flags.
        Does NOT restore dynamicTerms (must be re-registered).
        
        Args:
            filepath: Input JSON file path
        """
        with open(filepath, 'r') as f:
            state = json.load(f)
        
        self.dynamicParameters = state.get('dynamicParameters', {})
        self.metadata.update(state.get('metadata', {}))
        self.enableDynamicTerms = state.get('enableDynamicTerms', False)
        self.enableLogging = state.get('enableLogging', False)
        self.learningRate = state.get('learningRate', 0.001)
        
        if self.enableLogging:
            self.logger.info(f"State imported from {filepath}")


class DynamicVacuumTerm(PhysicsTerm):
    """
    Time-varying vacuum energy density modulation.
    
    Mirrors C++ DynamicVacuumTerm from source70.cpp lines 60-78.
    
    Physics: Vacuum energy density oscillates with amplitude and frequency.
    Formula: ρ_vac(t) = amplitude * ρ_vac_UA * sin(frequency * t)
    
    Parameters:
        amplitude: Oscillation amplitude (dimensionless)
        frequency: Oscillation frequency (Hz)
        ρ_vac_UA: Baseline vacuum density (kg/m³, default 7.09e-36)
    
    Usage:
        term = DynamicVacuumTerm(amplitude=1e-10, frequency=1e-15)
        rho = term.compute(t=1e6, params={'rho_vac_UA': 7.09e-36})
    """
    
    def __init__(self, amplitude: float = 1e-10, frequency: float = 1e-15):
        """
        Initialize dynamic vacuum term.
        
        Args:
            amplitude: Oscillation amplitude (default 1e-10)
            frequency: Oscillation frequency in Hz (default 1e-15)
        """
        super().__init__()
        self.amplitude = amplitude
        self.frequency = frequency
        self.set_metadata('physics_type', 'vacuum_energy')
        self.set_metadata('author', 'Daniel T. Murphy')
    
    def compute(self, t: float, params: Dict[str, float]) -> float:
        """
        Compute vacuum energy contribution at time t.
        
        Args:
            t: Time in seconds
            params: Must contain 'rho_vac_UA' (baseline vacuum density)
        
        Returns:
            Vacuum energy contribution (kg/m³)
        """
        rho_vac = params.get('rho_vac_UA', 7.09e-36)
        
        # Override with dynamic parameters if set
        amp = self.get_dynamic_parameter('amplitude', self.amplitude)
        freq = self.get_dynamic_parameter('frequency', self.frequency)
        
        result = amp * rho_vac * math.sin(freq * t)
        
        if self.enableLogging:
            self.logger.debug(f"DynamicVacuum: t={t:.3e}, amp={amp:.3e}, "
                            f"freq={freq:.3e}, result={result:.6e}")
        
        return result
    
    def getName(self) -> str:
        return "DynamicVacuumTerm"
    
    def getDescription(self) -> str:
        return f"Time-varying vacuum energy (A={self.amplitude:.2e}, f={self.frequency:.2e} Hz)"
    
    def validate(self, params: Dict[str, float]) -> bool:
        """Validate that vacuum density is positive."""
        rho_vac = params.get('rho_vac_UA', 7.09e-36)
        return rho_vac > 0


class QuantumCouplingTerm(PhysicsTerm):
    """
    Non-local quantum coupling effects.
    
    Mirrors C++ QuantumCouplingTerm from source70.cpp lines 92-110.
    
    Physics: Quantum uncertainty introduces gravity corrections.
    Formula: g_quantum = coupling_strength * (ℏ²) / (M * r²) * cos(t / τ)
    
    Parameters:
        coupling_strength: Dimensionless coupling constant (default 1e-40)
        timescale: Oscillation timescale τ (default 1e6 s)
    
    Usage:
        term = QuantumCouplingTerm(coupling_strength=1e-40)
        g_q = term.compute(t=1e6, params={'M': 1e30, 'r': 1e4, 'hbar': 1.0546e-34})
    """
    
    def __init__(self, coupling_strength: float = 1e-40, timescale: float = 1e6):
        """
        Initialize quantum coupling term.
        
        Args:
            coupling_strength: Coupling constant (default 1e-40)
            timescale: Oscillation timescale in seconds (default 1e6)
        """
        super().__init__()
        self.coupling_strength = coupling_strength
        self.timescale = timescale
        self.set_metadata('physics_type', 'quantum_gravity')
        self.set_metadata('author', 'Daniel T. Murphy')
    
    def compute(self, t: float, params: Dict[str, float]) -> float:
        """
        Compute quantum coupling contribution.
        
        Args:
            t: Time in seconds
            params: Must contain 'M' (mass in kg), 'r' (radius in m), 'hbar' (optional)
        
        Returns:
            Quantum gravity correction (m/s²)
        """
        hbar_val = params.get('hbar', hbar)
        M = params.get('M', 1.989e30)  # Default solar mass
        r = params.get('r', 1e4)       # Default 10 km
        
        # Override with dynamic parameters
        coupling = self.get_dynamic_parameter('coupling_strength', self.coupling_strength)
        tau = self.get_dynamic_parameter('timescale', self.timescale)
        
        if r == 0:
            return 0.0  # Avoid division by zero
        
        result = coupling * (hbar_val ** 2) / (M * r * r) * math.cos(t / tau)
        
        if self.enableLogging:
            self.logger.debug(f"QuantumCoupling: t={t:.3e}, M={M:.3e}, "
                            f"r={r:.3e}, result={result:.6e}")
        
        return result
    
    def getName(self) -> str:
        return "QuantumCouplingTerm"
    
    def getDescription(self) -> str:
        return f"Non-local quantum effects (κ={self.coupling_strength:.2e}, τ={self.timescale:.2e} s)"
    
    def validate(self, params: Dict[str, float]) -> bool:
        """Validate that mass and radius are positive."""
        M = params.get('M', 1.989e30)
        r = params.get('r', 1e4)
        return M > 0 and r > 0


class DarkMatterHaloTerm(PhysicsTerm):
    """
    Dark matter halo gravitational contribution.
    
    Physics: NFW profile with scale radius and density.
    Formula: g_DM = 4πG ρ_s r_s³ / r² * [ln((r_s + r)/r_s) - r/(r_s + r)]
    
    Parameters:
        rho_s: Characteristic density (kg/m³)
        r_s: Scale radius (m)
    
    Usage:
        term = DarkMatterHaloTerm(rho_s=1e-20, r_s=1e4)
        g_DM = term.compute(t=0, params={'r': 1e5})
    """
    
    def __init__(self, rho_s: float = 1e-20, r_s: float = 1e4):
        """
        Initialize dark matter halo term.
        
        Args:
            rho_s: Characteristic density in kg/m³
            r_s: Scale radius in m
        """
        super().__init__()
        self.rho_s = rho_s
        self.r_s = r_s
        self.set_metadata('physics_type', 'dark_matter')
        self.set_metadata('profile', 'NFW')
    
    def compute(self, t: float, params: Dict[str, float]) -> float:
        """
        Compute dark matter halo contribution.
        
        Args:
            t: Time (unused for static halo)
            params: Must contain 'r' (radius in m)
        
        Returns:
            Dark matter gravity (m/s²)
        """
        r = params.get('r', 1e4)
        
        # Override with dynamic parameters
        rho = self.get_dynamic_parameter('rho_s', self.rho_s)
        rs = self.get_dynamic_parameter('r_s', self.r_s)
        
        if r == 0 or rs == 0:
            return 0.0
        
        # NFW profile
        x = r / rs
        f_x = math.log((rs + r) / rs) - r / (rs + r)
        result = 4 * math.pi * dpm_ug1_seed(rho * (rs ** 3), r) * f_x
        
        if self.enableLogging:
            self.logger.debug(f"DarkMatterHalo: r={r:.3e}, x={x:.3f}, result={result:.6e}")
        
        return result
    
    def getName(self) -> str:
        return "DarkMatterHaloTerm"
    
    def getDescription(self) -> str:
        return f"NFW dark matter halo (ρ_s={self.rho_s:.2e} kg/m³, r_s={self.r_s:.2e} m)"
    
    def validate(self, params: Dict[str, float]) -> bool:
        """Validate positive radius."""
        r = params.get('r', 1e4)
        return r > 0


# ===========================================================================================
# MODULE METADATA
# ===========================================================================================

FRAMEWORK_VERSION = "2.0-Enhanced"
FRAMEWORK_AUTHOR = "Daniel T. Murphy"
FRAMEWORK_DATE = "2026-02-14"
FRAMEWORK_DESCRIPTION = "Python port of C++ self-expanding physics framework from source70.cpp"

# Available physics terms
AVAILABLE_TERMS = {
    'DynamicVacuumTerm': DynamicVacuumTerm,
    'QuantumCouplingTerm': QuantumCouplingTerm,
    'DarkMatterHaloTerm': DarkMatterHaloTerm,
}


def create_term(term_type: str, **kwargs) -> PhysicsTerm:
    """
    Factory function to create physics terms.
    
    Args:
        term_type: Term class name (e.g., 'DynamicVacuumTerm')
        **kwargs: Term-specific initialization parameters
    
    Returns:
        PhysicsTerm instance
    
    Example:
        >>> term = create_term('DynamicVacuumTerm', amplitude=1e-10, frequency=1e-15)
        >>> result = term.compute(t=1e6, params={'rho_vac_UA': 7.09e-36})
    """
    if term_type not in AVAILABLE_TERMS:
        raise ValueError(f"Unknown term type: {term_type}. Available: {list(AVAILABLE_TERMS.keys())}")
    
    return AVAILABLE_TERMS[term_type](**kwargs)


if __name__ == "__main__":
    # Demo usage
    print(f"PhysicsFramework.py - Version {FRAMEWORK_VERSION}")
    print(f"Available terms: {list(AVAILABLE_TERMS.keys())}")
    print()
    
    # Create and test DynamicVacuumTerm
    vacuum_term = DynamicVacuumTerm(amplitude=1e-10, frequency=1e-15)
    vacuum_term.enableLogging = True
    result = vacuum_term.compute(t=1e6, params={'rho_vac_UA': 7.09e-36})
    print(f"{vacuum_term.getName()}: {result:.6e}")
    print(f"Description: {vacuum_term.getDescription()}")
    print()
    
    # Test quantum coupling
    quantum_term = QuantumCouplingTerm(coupling_strength=1e-40)
    quantum_term.set_dynamic_parameter('coupling_strength', 2e-40)
    result = quantum_term.compute(t=1e6, params={'M': 1e30, 'r': 1e4})
    print(f"{quantum_term.getName()}: {result:.6e}")
    print()
    
    # Test nested terms
    print("Testing nested term registration:")
    vacuum_term.register_dynamic_term(quantum_term)
    vacuum_term.enableDynamicTerms = True
    total = vacuum_term.compute_with_dynamic_terms(t=1e6, params={
        'rho_vac_UA': 7.09e-36,
        'M': 1e30,
        'r': 1e4
    })
    print(f"Total with dynamic terms: {total:.6e}")
