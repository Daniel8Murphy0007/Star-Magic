#!/usr/bin/env python3
"""
CondensedPhysicsAggregator.py - Unified UQFF Calculator Import Surface
======================================================================

Master aggregation module that provides a unified import surface for all
CondensedPhysics calculator modules. This enables scalable file clustering
while maintaining a single-import API.

ARCHITECTURE:
    CondensedPhysics.py      → Foundation (1011 base classes)
    CondensedPhysics2.py     → Extension 1 (17+ classes: Orb Analysis 10/11+)
    CondensedPhysics3.py     → Extension 2 (future)
    CondensedPhysicsAggregator.py → This file (unified API)

USAGE:
    # Import everything from unified API
    from CondensedPhysicsAggregator import *
    
    # Or import specific modules
    from CondensedPhysicsAggregator import CP1_CALCULATORS, CP2_CALCULATORS
    
    # Access aggregated registry
    from CondensedPhysicsAggregator import ALL_CALCULATORS

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
Created: February 27, 2026
"""

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORT FROM CONDENSEDPHYSICS.PY (FOUNDATION - 1011 base classes)
# ═══════════════════════════════════════════════════════════════════════════════

from CondensedPhysics import *

# Create alias for CP1 calculators (foundation module)
# Note: Individual calculators are imported via wildcard above
CP1_MODULE_NAME = "CondensedPhysics"
CP1_VERSION = "1.0.0"

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORT FROM CONDENSEDPHYSICS2.PY (EXTENSION 1 - Orb Analysis 10/11+)
# ═══════════════════════════════════════════════════════════════════════════════

from CondensedPhysics2 import (
    # Version & metadata
    CP2_VERSION,
    CP2_CLASS_COUNT,
    
    # Orb Analysis_10 (8 classes)
    ORB_ANALYSIS_10_PARAMS,
    ThirtySixFrameSequenceCalculator,
    CyclicalConvectionPatternCalculator,
    Orb10RefinedFUCalculator,
    SpookyActionNonLocalTransferCalculator,
    ThermalGradientDrivenDynamicsCalculator,
    QuadrantTransitionTrackerCalculator,
    ACEDCEModulatedEnergyCalculator,
    MagneticBubbleConfinementCalculator,
    ORB_ANALYSIS_10_CALCULATORS,
    
    # Orb Analysis_11 (9 classes)
    ORB_ANALYSIS_11_PARAMS,
    ThirtyNineFrameSequenceCalculator,
    CounterClockwiseDiagonalCycleCalculator,
    Orb11RefinedFUCalculator,
    ExtendedCyclePatternAnalyzerCalculator,
    IntelligentPlasmoidBehaviorCalculator,
    BulbDrivenPlasmaEnergeticsCalculator,
    WaxCapCoolingDynamicsCalculator,
    FieldGeneratorResonanceCouplingCalculator,
    TotalEnergyBudgetCalculator,
    ORB_ANALYSIS_11_CALCULATORS,
    
    # Aggregated CP2 registry
    CP2_CALCULATORS,
)

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORT FROM GROK DEEP ANALYSIS MODULES (Equations 12-99)
# ═══════════════════════════════════════════════════════════════════════════════

# AGN Feedback Module (Eqs 14-15, 23-25, 70-72)
from agn_feedback_module import (
    AGN_FEEDBACK_CALCULATORS,
    BlandfordZnajekCalculator,
    RelativisticJetVelocityCalculator,
    AGNOutflowMomentumCalculator,
    JetPowerFromSpinCalculator,
    FeedbackDutyCycleCalculator,
    WindTerminalVelocityCalculator,
    IonizationParameterCalculator,
    FeedbackEnergyCouplingCalculator,
    AGNFeedbackCalculator,
    AGN_SYSTEMS,
)

# GRB Equations Module (Eqs 12-13, 19-20, 55-57, 73-75)
from grb_equations_module import (
    GRB_CALCULATORS,
    FireballExpansionCalculator,
    AfterglowSynchrotronCalculator,
    ChirpMassCalculator,
    RingdownQNMCalculator,
    BinaryPulsarOrbitDecayCalculator,
    PeriastronAdvanceCalculator,
    KilonovaLightCurveCalculator,
    InspiralFrequencyEvolutionCalculator,
    MergerTimeCalculator,
    GRBCalculator,
    GRB_SYSTEMS,
)

# Dark Matter Halos Module (Eqs 29-31, 82-84)
from dark_matter_halos_module import (
    DARK_MATTER_CALCULATORS,
    NFWProfileCalculator,
    RotationCurveCalculator,
    SIDMCoreFormationCalculator,
    FirstHaloCalculator,
    StarFormationEfficiencyCalculator,
    FeedbackEnergyInjectionCalculator,
    VirialEquilibriumCalculator,
    DarkMatterHaloCalculator,
    HALO_SYSTEMS,
)

# Stellar Evolution Module (Eqs 42-44, 58-63)
from stellar_evolution_module import (
    STELLAR_EVOLUTION_CALCULATORS,
    MainSequenceLifetimeCalculator,
    MassLuminosityCalculator,
    ConvectiveTurnoverCalculator,
    TypeIaSupernovaCalculator,
    CoreCollapseSupernovaCalculator,
    NucleosynthesisYieldCalculator,
    PlanetaryNebulaCalculator,
    StellarWindCalculator,
    StellarEvolutionCalculator,
    STELLAR_SYSTEMS,
)

# MHD Dynamo Module (Eqs 39-41, 88-90)
from mhd_dynamo_module import (
    MHD_DYNAMO_CALCULATORS,
    KazantsevDynamoCalculator,
    AlfvenMachNumberCalculator,
    FieldReversalCalculator,
    MeanFieldDynamoCalculator,
    ISMTurbulenceCascadeCalculator,
    MagneticFluxFreezeInCalculator,
    MHDDynamoCalculator,
    MHD_SYSTEMS,
)

# Black Hole Thermodynamics Module (Eqs 94-99)
from bh_thermodynamics_module import (
    BH_THERMODYNAMICS_CALCULATORS,
    HawkingTemperatureCalculator,
    BekensteinHawkingEntropyCalculator,
    BlackHoleEvaporationCalculator,
    LQCBounceCalculator,
    PrimordialBlackHoleCalculator,
    HolographicPrincipleCalculator,
    BlackHoleThermodynamicsCalculator,
    BH_SYSTEMS,
)

# ═══════════════════════════════════════════════════════════════════════════════
# AGGREGATED MASTER REGISTRY
# ═══════════════════════════════════════════════════════════════════════════════

# Aggregate all calculators from all CP modules
ALL_CALCULATORS = {
    **CP2_CALCULATORS,
    # Grok Deep Analysis Modules (Equations 12-99)
    **AGN_FEEDBACK_CALCULATORS,
    **GRB_CALCULATORS,
    **DARK_MATTER_CALCULATORS,
    **STELLAR_EVOLUTION_CALCULATORS,
    **MHD_DYNAMO_CALCULATORS,
    **BH_THERMODYNAMICS_CALCULATORS,
    # Note: CP1 calculators are accessed individually via class names
    # due to the large number (~1011). Add specific registries here as needed.
}

# Module metadata
AGGREGATOR_VERSION = "1.1.0"
TOTAL_MODULES = 8  # CP1, CP2, + 6 Grok Deep Analysis modules


def get_calculator(name: str):
    """
    Get a calculator by name from any CP module.
    
    Args:
        name: Calculator class name (e.g., 'ThirtySixFrameSequenceCalculator')
        
    Returns:
        Calculator instance or None if not found
    """
    # Check CP2 first (extension modules)
    if name in ALL_CALCULATORS:
        return ALL_CALCULATORS[name]
    
    # Try to get from global namespace (CP1 classes)
    if name in globals():
        cls = globals()[name]
        if hasattr(cls, 'compute'):
            return cls() if callable(cls) else cls
    
    return None


def list_all_calculators():
    """
    List all available calculator names across all modules.
    
    Returns:
        dict: Module name → list of calculator names
    """
    return {
        'CP2_ORB_ANALYSIS_10': list(ORB_ANALYSIS_10_CALCULATORS.keys()),
        'CP2_ORB_ANALYSIS_11': list(ORB_ANALYSIS_11_CALCULATORS.keys()),
        # Grok Deep Analysis modules
        'AGN_FEEDBACK': list(AGN_FEEDBACK_CALCULATORS.keys()),
        'GRB_EQUATIONS': list(GRB_CALCULATORS.keys()),
        'DARK_MATTER_HALOS': list(DARK_MATTER_CALCULATORS.keys()),
        'STELLAR_EVOLUTION': list(STELLAR_EVOLUTION_CALCULATORS.keys()),
        'MHD_DYNAMO': list(MHD_DYNAMO_CALCULATORS.keys()),
        'BH_THERMODYNAMICS': list(BH_THERMODYNAMICS_CALCULATORS.keys()),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE EXPORTS
# ═══════════════════════════════════════════════════════════════════════════════

__all__ = [
    # Aggregator utilities
    'ALL_CALCULATORS',
    'AGGREGATOR_VERSION',
    'TOTAL_MODULES',
    'get_calculator',
    'list_all_calculators',
    
    # CP2 exports (explicit)
    'CP2_VERSION',
    'CP2_CLASS_COUNT',
    'CP2_CALCULATORS',
    
    # Orb Analysis_10
    'ORB_ANALYSIS_10_PARAMS',
    'ThirtySixFrameSequenceCalculator',
    'CyclicalConvectionPatternCalculator',
    'Orb10RefinedFUCalculator',
    'SpookyActionNonLocalTransferCalculator',
    'ThermalGradientDrivenDynamicsCalculator',
    'QuadrantTransitionTrackerCalculator',
    'ACEDCEModulatedEnergyCalculator',
    'MagneticBubbleConfinementCalculator',
    'ORB_ANALYSIS_10_CALCULATORS',
    
    # Orb Analysis_11
    'ORB_ANALYSIS_11_PARAMS',
    'ThirtyNineFrameSequenceCalculator',
    'CounterClockwiseDiagonalCycleCalculator',
    'Orb11RefinedFUCalculator',
    'ExtendedCyclePatternAnalyzerCalculator',
    'IntelligentPlasmoidBehaviorCalculator',
    'BulbDrivenPlasmaEnergeticsCalculator',
    'WaxCapCoolingDynamicsCalculator',
    'FieldGeneratorResonanceCouplingCalculator',
    'TotalEnergyBudgetCalculator',
    'ORB_ANALYSIS_11_CALCULATORS',
    
    # AGN Feedback Module
    'AGN_FEEDBACK_CALCULATORS',
    'BlandfordZnajekCalculator',
    'RelativisticJetVelocityCalculator',
    'AGNOutflowMomentumCalculator',
    'JetPowerFromSpinCalculator',
    'FeedbackDutyCycleCalculator',
    'WindTerminalVelocityCalculator',
    'IonizationParameterCalculator',
    'FeedbackEnergyCouplingCalculator',
    'AGNFeedbackCalculator',
    'AGN_SYSTEMS',
    
    # GRB Equations Module
    'GRB_CALCULATORS',
    'FireballExpansionCalculator',
    'AfterglowSynchrotronCalculator',
    'ChirpMassCalculator',
    'RingdownQNMCalculator',
    'BinaryPulsarOrbitDecayCalculator',
    'PeriastronAdvanceCalculator',
    'KilonovaLightCurveCalculator',
    'InspiralFrequencyEvolutionCalculator',
    'MergerTimeCalculator',
    'GRBCalculator',
    'GRB_SYSTEMS',
    
    # Dark Matter Halos Module
    'DARK_MATTER_CALCULATORS',
    'NFWProfileCalculator',
    'RotationCurveCalculator',
    'SIDMCoreFormationCalculator',
    'FirstHaloCalculator',
    'StarFormationEfficiencyCalculator',
    'FeedbackEnergyInjectionCalculator',
    'VirialEquilibriumCalculator',
    'DarkMatterHaloCalculator',
    'HALO_SYSTEMS',
    
    # Stellar Evolution Module
    'STELLAR_EVOLUTION_CALCULATORS',
    'MainSequenceLifetimeCalculator',
    'MassLuminosityCalculator',
    'ConvectiveTurnoverCalculator',
    'TypeIaSupernovaCalculator',
    'CoreCollapseSupernovaCalculator',
    'NucleosynthesisYieldCalculator',
    'PlanetaryNebulaCalculator',
    'StellarWindCalculator',
    'StellarEvolutionCalculator',
    'STELLAR_SYSTEMS',
    
    # MHD Dynamo Module
    'MHD_DYNAMO_CALCULATORS',
    'KazantsevDynamoCalculator',
    'AlfvenMachNumberCalculator',
    'FieldReversalCalculator',
    'MeanFieldDynamoCalculator',
    'ISMTurbulenceCascadeCalculator',
    'MagneticFluxFreezeInCalculator',
    'MHDDynamoCalculator',
    'MHD_SYSTEMS',
    
    # Black Hole Thermodynamics Module
    'BH_THERMODYNAMICS_CALCULATORS',
    'HawkingTemperatureCalculator',
    'BekensteinHawkingEntropyCalculator',
    'BlackHoleEvaporationCalculator',
    'LQCBounceCalculator',
    'PrimordialBlackHoleCalculator',
    'HolographicPrincipleCalculator',
    'BlackHoleThermodynamicsCalculator',
    'BH_SYSTEMS',
]
