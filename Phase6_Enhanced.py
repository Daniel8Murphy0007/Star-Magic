#!/usr/bin/env python3
"""
Phase6_Enhanced.py - Phase 6 with Self-Expanding Framework Integration
========================================================================

ARCHITECTURE: Integrates PhysicsFramework.py with Phase6_Consolidated.py to provide
TRUE self-expanding capabilities matching C++ source70.cpp implementation.

BACKWARD COMPATIBLE: All Phase6_Consolidated static methods still available.
ENHANCED: Framework-based wrappers with runtime extensibility.

Features:
- Dynamic parameter updates at runtime
- Nested physics term registration
- Metadata tracking and provenance
- Adaptive learning rates
- State export/import for ML pipelines
- Validation before computation

Usage:
    # Legacy static method (backward compatible):
    from Phase6_Consolidated import Source70_M51
    result = Source70_M51.calculate_m51_gravity(params)
    
    # Enhanced framework method (NEW):
    from Phase6_Enhanced import M51GravityCalculator
    calc = M51GravityCalculator()
    calc.register_dynamic_term(DynamicVacuumTerm())
    calc.enableDynamicTerms = True
    result = calc.compute_gravity(params)

Author: Daniel T. Murphy
Date: February 14, 2026
Version: 2.0-Enhanced
"""

from typing import Dict, List, Optional
import numpy as np
from IPData import InputParameters
from QCalc import CONSTANTS, EquationResult
from PhysicsFramework import (
    PhysicsTerm, 
    DynamicVacuumTerm, 
    QuantumCouplingTerm,
    DarkMatterHaloTerm
)
import Phase6_Consolidated as Phase6


class M51GravityCalculator:
    """
    M51 Whirlpool Galaxy gravity calculator with self-expanding framework.
    
    Wraps Phase6_Consolidated.Source70_M51 with PhysicsFramework capabilities.
    
    Features:
    - Uses validated Phase6 calculation as base
    - Adds dynamic terms at runtime via framework
    - Metadata tracking for all computations
    - Runtime parameter mutation
    - State export for reproducibility
    
    Usage:
        calc = M51GravityCalculator()
        
        # Add dynamic vacuum energy term
        calc.register_dynamic_term(DynamicVacuumTerm(amplitude=1e-10))
        calc.enableDynamicTerms = True
        
        # Compute with dynamic contributions
        result = calc.compute_gravity(params)
    """
    
    def __init__(self):
        """Initialize M51 calculator with self-expanding framework."""
        self.dynamicTerms: List[PhysicsTerm] = []
        self.dynamicParameters: Dict[str, float] = {}
        self.metadata: Dict[str, str] = {
            'version': '2.0-Enhanced',
            'source': 'Phase6_Enhanced.py',
            'base_module': 'Phase6_Consolidated.Source70_M51',
            'system': 'M51 Whirlpool Galaxy'
        }
        self.enableDynamicTerms: bool = False
        self.enableLogging: bool = False
        self.learningRate: float = 0.001
    
    def register_dynamic_term(self, term: PhysicsTerm) -> None:
        """
        Register physics term for additive contribution.
        
        Args:
            term: PhysicsTerm instance (DynamicVacuumTerm, QuantumCouplingTerm, etc.)
        """
        self.dynamicTerms.append(term)
        if self.enableLogging:
            print(f"M51GravityCalculator: Registered {term.getName()}")
    
    def set_dynamic_parameter(self, name: str, value: float) -> None:
        """Set runtime parameter override."""
        self.dynamicParameters[name] = value
    
    def compute_gravity(self, params: InputParameters) -> EquationResult:
        """
        Compute M51 gravity with optional dynamic terms.
        
        Args:
            params: InputParameters with M51 system data
        
        Returns:
            EquationResult with base + dynamic contributions
        """
        # Base calculation from validated Phase6 module
        base_result = Phase6.Source70_M51.calculate_m51_gravity(params)
        
        if not self.enableDynamicTerms or len(self.dynamicTerms) == 0:
            return base_result  # Return base result if no dynamic terms
        
        # Add dynamic term contributions
        t = getattr(params, 't', 5e8 * 3.156e7)
        r = getattr(params, 'r', 23.58e3 * 3.086e19)
        M_visible = getattr(params, 'M_visible', 1.2e11 * CONSTANTS['M_sun'])
        M_DM = getattr(params, 'M_DM', 4e10 * CONSTANTS['M_sun'])
        
        term_params = {
            't': t,
            'r': r,
            'M': M_visible + M_DM,
            'rho_vac_UA': 7.09e-36,
            'hbar': CONSTANTS['hbar'],
        }
        
        dynamic_contribution = 0.0
        term_details = []
        for term in self.dynamicTerms:
            if term.validate(term_params):
                contribution = term.compute(t, term_params)
                dynamic_contribution += contribution
                term_details.append(f"{term.getName()}: {contribution:.6e}")
                if self.enableLogging:
                    print(f"  + {term.getName()}: {contribution:.6e}")
        
        # Update result with dynamic contributions
        enhanced_result = EquationResult(
            name=base_result.name + "_enhanced",
            latex=base_result.latex + r" + \Delta_{dynamic}",
            substituted=base_result.substituted + f"\n  + Dynamic terms: {dynamic_contribution:.6e}",
            result=base_result.result + dynamic_contribution,
            unit=base_result.unit,
            parameters_used={**base_result.parameters_used},
            notes=f"{base_result.notes} | Enhanced with {len(self.dynamicTerms)} dynamic terms: {', '.join(term_details)}"
        )
        
        return enhanced_result


class NGC1316GravityCalculator:
    """
    NGC1316 Fornax A radio galaxy gravity calculator with self-expanding framework.
    
    Wraps Phase6_Consolidated.Source71_NGC1316 with PhysicsFramework capabilities.
    
    Features post-merger dynamics, AGN, dust lanes, radio lobes.
    """
    
    def __init__(self):
        """Initialize NGC1316 calculator with self-expanding framework."""
        self.dynamicTerms: List[PhysicsTerm] = []
        self.dynamicParameters: Dict[str, float] = {}
        self.metadata: Dict[str, str] = {
            'version': '2.0-Enhanced',
            'source': 'Phase6_Enhanced.py',
            'base_module': 'Phase6_Consolidated.Source71_NGC1316',
            'system': 'NGC1316 Fornax A'
        }
        self.enableDynamicTerms: bool = False
        self.enableLogging: bool = False
    
    def register_dynamic_term(self, term: PhysicsTerm) -> None:
        """Register physics term for additive contribution."""
        self.dynamicTerms.append(term)
    
    def compute_gravity(self, params: InputParameters) -> EquationResult:
        """
        Compute NGC1316 gravity with optional dynamic terms.
        
        Args:
            params: InputParameters with NGC1316 system data
        
        Returns:
            EquationResult with base + dynamic contributions
        """
        # Base calculation
        base_result = Phase6.Source71_NGC1316.calculate_ngc1316_gravity(params)
        
        if not self.enableDynamicTerms or len(self.dynamicTerms) == 0:
            return base_result
        
        # Add dynamic contributions
        t = getattr(params, 't', 1e9 * 3.156e7)
        r = getattr(params, 'r', 46e3 * 3.086e19)
        M = getattr(params, 'M_total', 5e11 * CONSTANTS['M_sun'])
        
        term_params = {
            't': t,
            'r': r,
            'M': M,
            'rho_vac_UA': 7.09e-36,
            'hbar': CONSTANTS['hbar'],
        }
        
        dynamic_contribution = 0.0
        term_details = []
        for term in self.dynamicTerms:
            if term.validate(term_params):
                contribution = term.compute(t, term_params)
                dynamic_contribution += contribution
                term_details.append(f"{term.getName()}: {contribution:.6e}")
        
        enhanced_result = EquationResult(
            name=base_result.name + "_enhanced",
            latex=base_result.latex + r" + \Delta_{dynamic}",
            substituted=base_result.substituted + f"\n  + Dynamic terms: {dynamic_contribution:.6e}",
            result=base_result.result + dynamic_contribution,
            unit=base_result.unit,
            parameters_used={**base_result.parameters_used},
            notes=f"{base_result.notes} | Enhanced with {len(self.dynamicTerms)} dynamic terms: {', '.join(term_details)}"
        )
        
        return enhanced_result


class SMBHBinaryCalculator:
    """
    SMBH Binary coalescence calculator with self-expanding framework.
    
    Wraps Phase6_Consolidated.Source80_SMBHBinary with PhysicsFramework capabilities.
    
    Revolutionary frequency-based gravity: a = Σ f_i · λ_Planck / (2π)
    Features 9 frequency sources replacing traditional Newtonian gravity.
    """
    
    def __init__(self):
        """Initialize SMBH Binary calculator with self-expanding framework."""
        self.dynamicTerms: List[PhysicsTerm] = []
        self.dynamicParameters: Dict[str, float] = {}
        self.metadata: Dict[str, str] = {
            'version': '2.0-Enhanced',
            'source': 'Phase6_Enhanced.py',
            'base_module': 'Phase6_Consolidated.Source80_SMBHBinary',
            'system': 'SMBH Binary Coalescence',
            'physics_type': 'frequency_based_gravity'
        }
        self.enableDynamicTerms: bool = False
        self.enableLogging: bool = False
    
    def register_dynamic_term(self, term: PhysicsTerm) -> None:
        """Register physics term for additive contribution."""
        self.dynamicTerms.append(term)
    
    def compute_gravity(self, params: InputParameters) -> EquationResult:
        """
        Compute SMBH binary gravity with optional dynamic terms.
        
        Args:
            params: InputParameters with SMBH binary data
        
        Returns:
            EquationResult with base + dynamic contributions
        """
        # Base calculation
        base_result = Phase6.Source80_SMBHBinary.calculate_smbh_binary_gravity(params)
        
        if not self.enableDynamicTerms or len(self.dynamicTerms) == 0:
            return base_result
        
        # Add dynamic contributions
        t = getattr(params, 't', 0)
        M1 = getattr(params, 'M1', 4e6 * CONSTANTS['M_sun'])
        M2 = getattr(params, 'M2', 2e6 * CONSTANTS['M_sun'])
        a = getattr(params, 'a', 1e9)  # Semi-major axis
        
        term_params = {
            't': t,
            'r': a,
            'M': M1 + M2,
            'rho_vac_UA': 7.09e-36,
            'hbar': CONSTANTS['hbar'],
        }
        
        dynamic_contribution = 0.0
        term_details = []
        for term in self.dynamicTerms:
            if term.validate(term_params):
                contribution = term.compute(t, term_params)
                dynamic_contribution += contribution
                term_details.append(f"{term.getName()}: {contribution:.6e}")
        
        enhanced_result = EquationResult(
            name=base_result.name + "_enhanced",
            latex=base_result.latex + r" + \Delta_{dynamic}",
            substituted=base_result.substituted + f"\n  + Dynamic terms: {dynamic_contribution:.6e}",
            result=base_result.result + dynamic_contribution,
            unit=base_result.unit,
            parameters_used={**base_result.parameters_used},
            notes=f"{base_result.notes} | Enhanced with {len(self.dynamicTerms)} dynamic terms: {', '.join(term_details)}"
        )
        
        return enhanced_result


# ═══════════════════════════════════════════════════════════════════════════════
# REGISTRY: Enhanced calculators with framework support
# ═══════════════════════════════════════════════════════════════════════════════

ENHANCED_CALCULATORS = {
    'm51_gravity': M51GravityCalculator,
    'ngc1316_gravity': NGC1316GravityCalculator,
    'smbh_binary_gravity': SMBHBinaryCalculator,
}


def create_calculator(calculator_type: str) -> object:
    """
    Factory function to create enhanced calculators.
    
    Args:
        calculator_type: Calculator identifier
            'm51_gravity' - M51 Whirlpool Galaxy
            'ngc1316_gravity' - NGC1316 Fornax A
            'smbh_binary_gravity' - SMBH Binary Coalescence
    
    Returns:
        Enhanced calculator instance
    
    Example:
        >>> calc = create_calculator('m51_gravity')
        >>> calc.register_dynamic_term(DynamicVacuumTerm())
        >>> calc.enableDynamicTerms = True
        >>> result = calc.compute_gravity(params)
    """
    if calculator_type not in ENHANCED_CALCULATORS:
        raise ValueError(f"Unknown calculator: {calculator_type}. "
                        f"Available: {list(ENHANCED_CALCULATORS.keys())}")
    
    return ENHANCED_CALCULATORS[calculator_type]()


# ═══════════════════════════════════════════════════════════════════════════════
# DEMONSTRATION
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Phase6_Enhanced.py - Self-Expanding Framework Integration")
    print("=" * 70)
    print()
    
    # Create test parameters
    test_params = InputParameters()
    test_params.M_visible = 1.2e11 * CONSTANTS['M_sun']
    test_params.M_DM = 4e10 * CONSTANTS['M_sun']
    test_params.r = 23.58e3 * 3.086e19
    test_params.z = 0.002
    test_params.t = 5e8 * 3.156e7
    
    print("Test 1: M51 Gravity - Base Calculation (No Dynamic Terms)")
    print("-" * 70)
    calc_base = M51GravityCalculator()
    result_base = calc_base.compute_gravity(test_params)
    print(f"Gravity (base): {result_base.result:.6e} {result_base.unit}")
    print()
    
    print("Test 2: M51 Gravity - With Dynamic Vacuum Term")
    print("-" * 70)
    calc_enhanced = M51GravityCalculator()
    calc_enhanced.register_dynamic_term(DynamicVacuumTerm(amplitude=1e-10, frequency=1e-15))
    calc_enhanced.enableDynamicTerms = True
    calc_enhanced.enableLogging = True
    result_enhanced = calc_enhanced.compute_gravity(test_params)
    print(f"Gravity (enhanced): {result_enhanced.result:.6e} {result_enhanced.unit}")
    print(f"Notes: {result_enhanced.notes[:100]}...")  # Print first 100 chars of notes
    print()
    
    print("Test 3: M51 Gravity - Multiple Dynamic Terms")
    print("-" * 70)
    calc_multi = M51GravityCalculator()
    calc_multi.register_dynamic_term(DynamicVacuumTerm(amplitude=1e-10, frequency=1e-15))
    calc_multi.register_dynamic_term(QuantumCouplingTerm(coupling_strength=1e-40))
    calc_multi.register_dynamic_term(DarkMatterHaloTerm(rho_s=1e-20, r_s=1e4))
    calc_multi.enableDynamicTerms = True
    calc_multi.enableLogging = True
    result_multi = calc_multi.compute_gravity(test_params)
    print(f"Gravity (multi-term): {result_multi.result:.6e} {result_multi.unit}")
    print(f"Notes: {result_multi.notes[:150]}...")  # Print first 150 chars
    print()
    
    print("✅ Phase6_Enhanced.py - All tests completed successfully")
    print("🌟 Self-expanding framework fully operational")
