#!/usr/bin/env python3
"""
CondensedPhysicsAggregator.py - Unified UQFF Calculator Import Surface
======================================================================

Master aggregation module that provides a unified import surface for all
CondensedPhysics calculator modules. This enables scalable file clustering
while maintaining a single-import API.

ARCHITECTURE:
    CondensedPhysics.py      → Foundation (1,227 base classes, 168,784 lines)
    CondensedPhysics2.py     → Extension 1 (600 classes, 45,991 lines: Orb Analysis 10/11+ and Grok thread extensions)
    CondensedPhysics3.py     → Extension 2 (219 classes, 13,944 lines: 15+ categories, Sessions 41-96)
    CondensedPhysics4.py     → Extension 3 (94 classes, 7,275 lines: Phase 4, Sessions 97-120, 2026-03-22)
    Last updated: Session 120 v4.90 (2026-03-22) — 15 root-level UQFF C++ module pairs added (grok_share_dc707f5d3); CP4 84→94; PAPER_447–455 integrated
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
# IMPORT FROM CONDENSEDPHYSICS3.PY (EXTENSION 2 — 118 classes, 15+ categories)
# Source: Grok threads ba4c0789 (Session 41), 7f9068 (Session 47), 381a8fe7 (Session 48), 7514fe (Session 50)
# ═══════════════════════════════════════════════════════════════════════════════

from CondensedPhysics3 import (
    # Solar System
    SolarWindBubbleVerificationCalculator,
    HeliopausalBoundaryStepFunctionCalculator,
    # Stars
    StellarClusterUg3DiskTurbulenceCalculator,
    StellarUg1DipoleDefectCalculator,
    # Exoplanets
    ExoplanetAtmosphericMassLossUbCalculator,
    PlanetaryCoreUg3PenetrationScalingCalculator,
    # White Dwarf
    WhiteDwarfUQFFGravitationalDecayCalculator,
    WhiteDwarfDegenerateElectronUiCalculator,
    # Supernova
    KilonovaTransientQWaveParameterCalculator,
    SupernovaProgenitorNegativeTimeZoneCalculator,
    # Neutron Star
    NeutronStarCRPIceCubeFluxVerificationCalculator,
    NeutronStarMergerUbOutflowF_UCalculator,
    # Black Hole
    BlackHoleJetFluidAsymmetryRatioCalculator,
    BlackHoleUg4GalacticFeedbackCalculator,
    # Super Massive BH
    GaiaSgrADistanceErrorAnalysisCalculator,
    QuasarBlazerLuminosityEreactVerificationCalculator,
    # Milky Way Galaxy
    GalacticCenterUg4KappaDecayCalibrationCalculator,
    MilkyWayGalacticSpinUb_iCouplingCalculator,
    # Galaxy
    GalaxyIMFNucleosynthesisIndexCalculator,
    GalaxyEquationOfStateUCFCalculator,
    # Quasar
    QuasarJetAsymmetryCosRatioCalculator,
    QuasarEddingtonExcessJetVelocityCalculator,
    # Galaxy Cluster
    GalaxyClusterPSZ2UmTurbulenceCalculator,
    GalaxyClusterPLCKDoubleRelicShearCalculator,
    # Cosmological
    TwentySixLevelPolynomialHierarchyFullCalculator,
    CosmologicalLineFluximeSFRIntegralCalculator,
    PDGNuclearPolynomialFitVerificationCalculator,
    # Deep Field
    DeepFieldShearDeltaTauConstraintCalculator,
    HighRedshiftJWSTQWaveDeepFieldCalculator,
    DeepFieldG359ShearNISPConstraintCalculator,
    # Miscellaneous
    QScopeFrequencyResonanceUQFFCalculator,
    ATLASLHCQuarkEnergyLowNLevelCalculator,
    VacuumEnergyComponentRatioCalculator,
    UQFFIPCChainStatusCalculator,
    # Session 47 — Solar System + Hybrid MUGE + GW + HE (grok_share_7f9068)
    SolarSystemFUValidatorCalculator,
    HybridMUGEBlendingCalculator,
    WormholeMUGE13thTermCalculator,
    J1610QuasarRelativisticSCmCalculator,
    StressEnergyAMunuCouplingCalculator,
    GW231123MassGapUQFFCalculator,
    HighEnergyDatasetValidationCalculator,
    # Session 48 — CoAnQi UQFF+3D+Plugin Integration (grok_share_381a8fe7)
    CoAnQiCelestialBodyFUCalculator,
    CoAnQiModularCompressedMUGECalculator,
    CoAnQiModularResonanceMUGECalculator,
    CoAnQi26LevelEnergyDensityCalculator,
    CoAnQiQuasarJetFluidCalculator,
    CoAnQiArchitectureCalculator,
    DiPseudoMonopoleDPMTheoryCalculator,
    # Session 50 — PAPER_196-215 (grok_share_7514fe)
    TriadicMasterEquationCalculator,
    FUBiiExtendedIntegralCalculator,
    FUBiiTaxonomyCompactObjectCalculator,
    FUBiiTaxonomyCosmologicalCalculator,
    UmUniversalMagnetismTaxonomyCalculator,
    UQFFGravitationalWaveChirpQNMCalculator,
    UQFFReionizationBBNCalculator,
    UQFFCMBStructureGrowthCalculator,
    UQFFDarkMatterNFWSIDMCalculator,
    RamanujanPolynomialsQ26Calculator,
    MagnetarVortexAvalancheCalculator,
    QuTiPQuantumEntanglementCalculator,
    UQFFVariableCalibrationCalculator,
    UQFFvsLambdaCDMComparisonCalculator,
    UQFFvsMONDComparisonCalculator,
    UQFF99SystemCompressionCalculator,
    UQFF48ScaleMolecularRotorCIACalculator,
    HResDUniverseMasterCalculator,
    MHDClustersJetsAccretionCalculator,
    CosmicRaysWHIMFermiCalculator,
    # Session 52 — grok_share_7514fe unique physics extraction
    UQFFCompressedFriedmannCalculator,
    UQFFMultiFactorEvolutionMergerCalculator,
    UQFFVelocityStarFormationCollisionCalculator,
    UQFFSupernovaFeedbackMassLossCalculator,
    HydrogenNuclearShellResonanceCalculator,
    UQFFUniverseDiameterEstimationCalculator,
    TriadicSSqFeedbackEnhancedCalculator,
    DPMHarmonicBuoyancySeriesCalculator,
    DipoleVortexPrimeEncodingCalculator,
    UQFFRelativisticHierarchyDecayIntegralCalculator,
    # Session 53 — grok_share_7514fe second-pass unique physics
    SgrAStarSpinDragUQFFCalculator,
    UQFFLensingModulationRingsCalculator,
    HydrogenAtomUQFFGravityCalculator,
    FUBiiFullDPMPolynomialIntegralCalculator,
    UQFFNeutrinoDecayRateCouplingCalculator,
    MagnetarSGR1745DynamicModulationCalculator,
    # Session 54 — grok_share_7514fe third-pass unique physics
    UQFFBuoyancyMasterIntegralCalculator,
    UQFFCGMSSqMetallicityCalculator,
    # Session 55 — grok_share_7514fe fourth-pass unique physics
    NGC3603StellarPressureModulationCalculator,
    M16EagleNebulaRadiationSFRCalculator,
    CrabPWNUQFFCalculator,
    UQFFSombreroDustIntegratedCalculator,
    # Session 56 — grok_share_7514fe fifth-pass unique physics
    BubbleNebulaExpansionEnhancementCalculator,
    HorseheadNebulaPradBlackbodyCalculator,
    NGC1275PerseusAGNFilamentCalculator,
    SaturnDualGravityRingTensionCalculator,
    # Session 57 — grok_share_7514fe sixth-pass (final): early-universe (v/c)^2·L_UV
    UQFFEarlyUniverseRelativisticUVCalculator,
    # Session 58 — grok_share_8d951e12.txt: 10 new classes (PAPER_226–235)
    MagnetarSGR0501MUGEFullCalculator,
    StarbirthTapestryLMCUQFFCalculator,
    Westerlund2MUGEStellarWindCalculator,
    PillarsOfCreationErosionMUGECalculator,
    GalaxyNGC2525SNMassLossCalculator,
    HUDFGalaxiesCosmicFieldCalculator,
    GalaxyNGC1792StarburstForgeCalculator,
    SGR1745BHProximityMagEnergyCalculator,
    SgrAStarAccretionPrecessionCalculator,
    AntennaeGalaxiesMergerInteractionCalculator,
    # Session 59 — grok_share_8d951e12.txt second-pass: 5 new classes (PAPER_236–241)
    UQFFLearningAdvancementCalculator,
    UQFFSource10CatalogueCalculator,
    UQFFVacuumRepulsionCalculator,
    UQFFTHzConduitShockCalculator,
    UQFFSpookyActionDPMCalculator,
    # Session 60 — grok_share_8d951e12.txt third-pass: 2 new classes (PAPER_242–243)
    RingsOfRelativityEinsteinLensingMUGECalculator,
    NGC3603FullMUGECavityPressureCalculator,
    # Session 62 — grok_share_8d951e12.txt fourth-pass: 6 new classes (PAPER_244–249)
    MUGEQuantumUncertaintyTermCalculator,
    MUGEFluidSelfGravityTermCalculator,
    MUGEDualModeOscillatoryGravityCalculator,
    MUGEMergerInteractionModulationCalculator,
    UQFFSource10BatchProfiledCalculator,
    UQFFCUDAGPUOptimizationPatternCalculator,
    # Sessions 63-74 — C++ UQFF 2.0 upgrades + ALMA Cycle 12 (PAPER_250–272)
    SN1006TypeIaSNRFUBiCalculator,
    EtaCarinaeHomuculusFUBiCalculator,
    ChandraArchiveMultiSystemFUBiCalculator,
    SgrACenterNegativeBuoyancyCalculator,
    KeplerSNR1604FUBiCalculator,
    PSRJ0030NeutronStarFUBiCalculator,
    CrabNebulaM1FUBiCalculator,
    CassiopeiaASNRFUBiCalculator,
    MultiMessengerUQFFValidator,
    HUDFTRZCPTPhaseCalculator,
    HUDFInteractionCascadeBuoyancyCalculator,
    HUDFGravitationalMeissnerCalculator,
    NGC1792StarburstBuoyancyCoherenceCalculator,
    NGC1792HubbleSlowModeOscillatorCalculator,
    NGC1792RamPressureDegeneracyCalculator,
    Source10DPMResonanceAmplificationCalculator,
    Source10GravitationalVacuumDragCalculator,
    Source10THzDoubleGateConduitCalculator,
    AndromedaBlueshiftApproachAmplifierCalculator,
    AndromedaHI21cmUQFFResonanceCalculator,
    AndromedaDMShellPartitionCalculator,
    AndromedaFriedmannHzExpansionCalculator,
    SombreroRecessionDampingKappaCalculator,
    SombreroRingResonatorDustRingCalculator,
    SombreroSMBHDominanceRatioCalculator,
    SaturnSolarTidalPerturbationCalculator,
    SaturnRingTidalGravityResonanceCalculator,
    SaturnAtmosphericWindKineticPressureCalculator,
    SaturnSolarTidalHubbleExpansionCalculator,
    M16DualMassCoActionProductCalculator,
    M16ErosionSaturationHalfTimeCalculator,
    M16NebularFriedmannRedshiftCalculator,
    ResonanceSCDPMTHzCascadeCalculator,
    ResonanceSCCosmicAgeStandingWaveCalculator,
    ResonanceSCCooperDPMFreqSynthesisCalculator,
    CrabSNRDPMDilutionCalculator,
    CrabFilamentSpectralTriadCalculator,
    CrabPulsarOscResonanceWindowCalculator,
    CR24DualChannelArchitectureCalculator,
    CR24VacuumDifferentialHarmonicCalculator,
    CR24CompressedCooperSuperSeedingCalculator,
    UniverseDiameterLambdaVacuumAccelerationCalculator,
    UniverseDiameterSuperluminalHubbleRatioCalculator,
    UniverseDiameterGRCurvatureDominanceCalculator,
    HydrogenAtomLorentzEMDominanceCalculator,
    HydrogenAtomLymanCosmosBridgeCalculator,
    HydrogenAtomProtonGRSpectralMinimumCalculator,
    HydrogenPToEUg4iResonanceBridgeCalculator,
    HydrogenPToETHzQuantumDegeneracyCalculator,
    HydrogenPToEAetherGravitationalDominanceCalculator,
    LagoonNebulaSFRMassRunawayCalculator,
    LagoonNebulaHerschelRadiationErosionCalculator,
    LagoonNebulaDualRadiationEMBarrierCalculator,
    SpiralArmTorqueGravitationalAmplifierCalculator,
    SNIaHubbleTensionImprintCalculator,
    SpiralDMVisiblePartitionRotationCalculator,
    BiPolarPNWindShockGravitationalDominanceCalculator,
    BiPolarPNUVRadiationPressureCalculator,
    EquatorialTorusMagneticConfinementCalculator,
    BipolarPNLobeResonanceDPMMacroAntennaCalculator,
    ResonanceVacDiffTHzCrossoverRadiusCalculator,
    CooperDPMf1THz_AscConfirmationCalculator,
    OrionTrapeziumWindRamPressureDominanceCalculator,
    OrionTrapeziumOBUVRadiationChampagneFlowCalculator,
    OrionCompactHIISFRBindingCrossoverCalculator,
    CR34DPMForceDensitySpectralAtlasCalculator,
    CR34CrossChannelDominanceCrossoverCalculator,
    CR34HiIRegionTHzGeometricDifferentialCalculator,
    CR34bVacuumAetherFrequencyModeCalculator,
    CR34bSaturnFirstPlanetaryDualChannelCalculator,
    CR34bRhoISMFluidDensityCouplingCalculator,
    # Sessions 91-96 — gok_share_31b5c807a4.txt UQFF assimilation + supplemental
    TriadicMasterFUg1R26StateRamanujanCalculator,
    QWave47NonGaussianDistributionCalculator,
    AlphaBECNuclearLENREnhancementCalculator,
    UmBilinearHeavisideNeutrinoVacuumCascadeCalculator,
    HResNuclear6EquationDipolekNucCalculator,
    MUGE26StateFrequencyBasisProofIdentitiesCalculator,
    FUBi12TermExplicitIntegrandCalculator,
    BSMUQFFMultiExperimentCouplingCalculator,
    UiComplexSuperconductiveVacuumDensityCalculator,
    kkREBTrdicRamanujanFUBiBuoyancyKernelCalculator,
    gCompressedAllForcesR26ComponentCalculator,
    QWave81PhaseSeparationValidationCalculator,
    NineSystemSepAstroParameterCatalogueCalculator,
    UmRotorStringTorqueIntegrationCalculator,
    EDMSO10BSMRefinedFuCalculator,
    UQFFSupplementCalibration3VarCalculator,
    MagnetarDPMTHzFrequencyFormCalculator,
    SGR17452900SCmLxFreqFormCalculator,
    SgrAStarGWPrecessionSquaredCalculator,
    TapestryStarbirthDPMTHzFreqCalculator,
    M87JetBZModelFUBiCalculator,
    CentaurusAFUBiJetVshapeCalculator,
    StephansQuintetShockRidgeFUBiCalculator,
    SPTClJ2215CoolCoreStarburstCalculator,
    ElGordoACTCLJ0102MergerFUBiCalculator,
    ASASSN14liTDEOutflowFUBiCalculator,
    RAquariiSymbioticBinaryFUBiCalculator,
    DecayRateVacuumRhoRatioDoubleExpCalculator,
    DUniverseSpatialCurvatureFifthFactorCalculator,
)

CP3_CALCULATORS = {
    'SolarWindBubbleVerificationCalculator': SolarWindBubbleVerificationCalculator,
    'HeliopausalBoundaryStepFunctionCalculator': HeliopausalBoundaryStepFunctionCalculator,
    'StellarClusterUg3DiskTurbulenceCalculator': StellarClusterUg3DiskTurbulenceCalculator,
    'StellarUg1DipoleDefectCalculator': StellarUg1DipoleDefectCalculator,
    'ExoplanetAtmosphericMassLossUbCalculator': ExoplanetAtmosphericMassLossUbCalculator,
    'PlanetaryCoreUg3PenetrationScalingCalculator': PlanetaryCoreUg3PenetrationScalingCalculator,
    'WhiteDwarfUQFFGravitationalDecayCalculator': WhiteDwarfUQFFGravitationalDecayCalculator,
    'WhiteDwarfDegenerateElectronUiCalculator': WhiteDwarfDegenerateElectronUiCalculator,
    'KilonovaTransientQWaveParameterCalculator': KilonovaTransientQWaveParameterCalculator,
    'SupernovaProgenitorNegativeTimeZoneCalculator': SupernovaProgenitorNegativeTimeZoneCalculator,
    'NeutronStarCRPIceCubeFluxVerificationCalculator': NeutronStarCRPIceCubeFluxVerificationCalculator,
    'NeutronStarMergerUbOutflowF_UCalculator': NeutronStarMergerUbOutflowF_UCalculator,
    'BlackHoleJetFluidAsymmetryRatioCalculator': BlackHoleJetFluidAsymmetryRatioCalculator,
    'BlackHoleUg4GalacticFeedbackCalculator': BlackHoleUg4GalacticFeedbackCalculator,
    'GaiaSgrADistanceErrorAnalysisCalculator': GaiaSgrADistanceErrorAnalysisCalculator,
    'QuasarBlazerLuminosityEreactVerificationCalculator': QuasarBlazerLuminosityEreactVerificationCalculator,
    'GalacticCenterUg4KappaDecayCalibrationCalculator': GalacticCenterUg4KappaDecayCalibrationCalculator,
    'MilkyWayGalacticSpinUb_iCouplingCalculator': MilkyWayGalacticSpinUb_iCouplingCalculator,
    'GalaxyIMFNucleosynthesisIndexCalculator': GalaxyIMFNucleosynthesisIndexCalculator,
    'GalaxyEquationOfStateUCFCalculator': GalaxyEquationOfStateUCFCalculator,
    'QuasarJetAsymmetryCosRatioCalculator': QuasarJetAsymmetryCosRatioCalculator,
    'QuasarEddingtonExcessJetVelocityCalculator': QuasarEddingtonExcessJetVelocityCalculator,
    'GalaxyClusterPSZ2UmTurbulenceCalculator': GalaxyClusterPSZ2UmTurbulenceCalculator,
    'GalaxyClusterPLCKDoubleRelicShearCalculator': GalaxyClusterPLCKDoubleRelicShearCalculator,
    'TwentySixLevelPolynomialHierarchyFullCalculator': TwentySixLevelPolynomialHierarchyFullCalculator,
    'CosmologicalLineFluximeSFRIntegralCalculator': CosmologicalLineFluximeSFRIntegralCalculator,
    'PDGNuclearPolynomialFitVerificationCalculator': PDGNuclearPolynomialFitVerificationCalculator,
    'DeepFieldShearDeltaTauConstraintCalculator': DeepFieldShearDeltaTauConstraintCalculator,
    'HighRedshiftJWSTQWaveDeepFieldCalculator': HighRedshiftJWSTQWaveDeepFieldCalculator,
    'DeepFieldG359ShearNISPConstraintCalculator': DeepFieldG359ShearNISPConstraintCalculator,
    'QScopeFrequencyResonanceUQFFCalculator': QScopeFrequencyResonanceUQFFCalculator,
    'ATLASLHCQuarkEnergyLowNLevelCalculator': ATLASLHCQuarkEnergyLowNLevelCalculator,
    'VacuumEnergyComponentRatioCalculator': VacuumEnergyComponentRatioCalculator,
    'UQFFIPCChainStatusCalculator': UQFFIPCChainStatusCalculator,
    # Session 47 — grok_share_7f9068
    'SolarSystemFUValidatorCalculator': SolarSystemFUValidatorCalculator,
    'HybridMUGEBlendingCalculator': HybridMUGEBlendingCalculator,
    'WormholeMUGE13thTermCalculator': WormholeMUGE13thTermCalculator,
    'J1610QuasarRelativisticSCmCalculator': J1610QuasarRelativisticSCmCalculator,
    'StressEnergyAMunuCouplingCalculator': StressEnergyAMunuCouplingCalculator,
    'GW231123MassGapUQFFCalculator': GW231123MassGapUQFFCalculator,
    'HighEnergyDatasetValidationCalculator': HighEnergyDatasetValidationCalculator,
    # Session 48 — grok_share_381a8fe7
    'CoAnQiCelestialBodyFUCalculator': CoAnQiCelestialBodyFUCalculator,
    'CoAnQiModularCompressedMUGECalculator': CoAnQiModularCompressedMUGECalculator,
    'CoAnQiModularResonanceMUGECalculator': CoAnQiModularResonanceMUGECalculator,
    'CoAnQi26LevelEnergyDensityCalculator': CoAnQi26LevelEnergyDensityCalculator,
    'CoAnQiQuasarJetFluidCalculator': CoAnQiQuasarJetFluidCalculator,
    'CoAnQiArchitectureCalculator': CoAnQiArchitectureCalculator,
    'DiPseudoMonopoleDPMTheoryCalculator': DiPseudoMonopoleDPMTheoryCalculator,
    # Session 50 — PAPER_196-215 (grok_share_7514fe)
    'TriadicMasterEquationCalculator': TriadicMasterEquationCalculator,
    'FUBiiExtendedIntegralCalculator': FUBiiExtendedIntegralCalculator,
    'FUBiiTaxonomyCompactObjectCalculator': FUBiiTaxonomyCompactObjectCalculator,
    'FUBiiTaxonomyCosmologicalCalculator': FUBiiTaxonomyCosmologicalCalculator,
    'UmUniversalMagnetismTaxonomyCalculator': UmUniversalMagnetismTaxonomyCalculator,
    'UQFFGravitationalWaveChirpQNMCalculator': UQFFGravitationalWaveChirpQNMCalculator,
    'UQFFReionizationBBNCalculator': UQFFReionizationBBNCalculator,
    'UQFFCMBStructureGrowthCalculator': UQFFCMBStructureGrowthCalculator,
    'UQFFDarkMatterNFWSIDMCalculator': UQFFDarkMatterNFWSIDMCalculator,
    'RamanujanPolynomialsQ26Calculator': RamanujanPolynomialsQ26Calculator,
    'MagnetarVortexAvalancheCalculator': MagnetarVortexAvalancheCalculator,
    'QuTiPQuantumEntanglementCalculator': QuTiPQuantumEntanglementCalculator,
    'UQFFVariableCalibrationCalculator': UQFFVariableCalibrationCalculator,
    'UQFFvsLambdaCDMComparisonCalculator': UQFFvsLambdaCDMComparisonCalculator,
    'UQFFvsMONDComparisonCalculator': UQFFvsMONDComparisonCalculator,
    'UQFF99SystemCompressionCalculator': UQFF99SystemCompressionCalculator,
    'UQFF48ScaleMolecularRotorCIACalculator': UQFF48ScaleMolecularRotorCIACalculator,
    'HResDUniverseMasterCalculator': HResDUniverseMasterCalculator,
    'MHDClustersJetsAccretionCalculator': MHDClustersJetsAccretionCalculator,
    'CosmicRaysWHIMFermiCalculator': CosmicRaysWHIMFermiCalculator,
    # Session 52 — grok_share_7514fe unique physics extraction
    'UQFFCompressedFriedmannCalculator': UQFFCompressedFriedmannCalculator,
    'UQFFMultiFactorEvolutionMergerCalculator': UQFFMultiFactorEvolutionMergerCalculator,
    'UQFFVelocityStarFormationCollisionCalculator': UQFFVelocityStarFormationCollisionCalculator,
    'UQFFSupernovaFeedbackMassLossCalculator': UQFFSupernovaFeedbackMassLossCalculator,
    'HydrogenNuclearShellResonanceCalculator': HydrogenNuclearShellResonanceCalculator,
    'UQFFUniverseDiameterEstimationCalculator': UQFFUniverseDiameterEstimationCalculator,
    'TriadicSSqFeedbackEnhancedCalculator': TriadicSSqFeedbackEnhancedCalculator,
    'DPMHarmonicBuoyancySeriesCalculator': DPMHarmonicBuoyancySeriesCalculator,
    'DipoleVortexPrimeEncodingCalculator': DipoleVortexPrimeEncodingCalculator,
    'UQFFRelativisticHierarchyDecayIntegralCalculator': UQFFRelativisticHierarchyDecayIntegralCalculator,
    # Session 53 — grok_share_7514fe second-pass unique physics
    'SgrAStarSpinDragUQFFCalculator': SgrAStarSpinDragUQFFCalculator,
    'UQFFLensingModulationRingsCalculator': UQFFLensingModulationRingsCalculator,
    'HydrogenAtomUQFFGravityCalculator': HydrogenAtomUQFFGravityCalculator,
    'FUBiiFullDPMPolynomialIntegralCalculator': FUBiiFullDPMPolynomialIntegralCalculator,
    'UQFFNeutrinoDecayRateCouplingCalculator': UQFFNeutrinoDecayRateCouplingCalculator,
    'MagnetarSGR1745DynamicModulationCalculator': MagnetarSGR1745DynamicModulationCalculator,
    # Session 54 — grok_share_7514fe third-pass unique physics
    'UQFFBuoyancyMasterIntegralCalculator': UQFFBuoyancyMasterIntegralCalculator,
    'UQFFCGMSSqMetallicityCalculator': UQFFCGMSSqMetallicityCalculator,
    # Session 55 — grok_share_7514fe fourth-pass unique physics
    'NGC3603StellarPressureModulationCalculator': NGC3603StellarPressureModulationCalculator,
    'M16EagleNebulaRadiationSFRCalculator': M16EagleNebulaRadiationSFRCalculator,
    'CrabPWNUQFFCalculator': CrabPWNUQFFCalculator,
    'UQFFSombreroDustIntegratedCalculator': UQFFSombreroDustIntegratedCalculator,
    # Session 56 — grok_share_7514fe fifth-pass unique physics
    'BubbleNebulaExpansionEnhancementCalculator': BubbleNebulaExpansionEnhancementCalculator,
    'HorseheadNebulaPradBlackbodyCalculator': HorseheadNebulaPradBlackbodyCalculator,
    'NGC1275PerseusAGNFilamentCalculator': NGC1275PerseusAGNFilamentCalculator,
    'SaturnDualGravityRingTensionCalculator': SaturnDualGravityRingTensionCalculator,
    # Session 57 — grok_share_7514fe sixth-pass (final): early-universe (v/c)^2·L_UV
    'UQFFEarlyUniverseRelativisticUVCalculator': UQFFEarlyUniverseRelativisticUVCalculator,
    # Session 58 — grok_share_8d951e12.txt: 10 new classes (PAPER_226–235)
    'MagnetarSGR0501MUGEFullCalculator': MagnetarSGR0501MUGEFullCalculator,
    'StarbirthTapestryLMCUQFFCalculator': StarbirthTapestryLMCUQFFCalculator,
    'Westerlund2MUGEStellarWindCalculator': Westerlund2MUGEStellarWindCalculator,
    'PillarsOfCreationErosionMUGECalculator': PillarsOfCreationErosionMUGECalculator,
    'GalaxyNGC2525SNMassLossCalculator': GalaxyNGC2525SNMassLossCalculator,
    'HUDFGalaxiesCosmicFieldCalculator': HUDFGalaxiesCosmicFieldCalculator,
    'GalaxyNGC1792StarburstForgeCalculator': GalaxyNGC1792StarburstForgeCalculator,
    'SGR1745BHProximityMagEnergyCalculator': SGR1745BHProximityMagEnergyCalculator,
    'SgrAStarAccretionPrecessionCalculator': SgrAStarAccretionPrecessionCalculator,
    'AntennaeGalaxiesMergerInteractionCalculator': AntennaeGalaxiesMergerInteractionCalculator,
    # Session 59 — grok_share_8d951e12.txt second-pass: 5 new classes (PAPER_236–241)
    'UQFFLearningAdvancementCalculator': UQFFLearningAdvancementCalculator,
    'UQFFSource10CatalogueCalculator': UQFFSource10CatalogueCalculator,
    'UQFFVacuumRepulsionCalculator': UQFFVacuumRepulsionCalculator,
    'UQFFTHzConduitShockCalculator': UQFFTHzConduitShockCalculator,
    'UQFFSpookyActionDPMCalculator': UQFFSpookyActionDPMCalculator,
    # Session 60 — grok_share_8d951e12.txt third-pass: 2 new classes (PAPER_242–243)
    'RingsOfRelativityEinsteinLensingMUGECalculator': RingsOfRelativityEinsteinLensingMUGECalculator,
    'NGC3603FullMUGECavityPressureCalculator': NGC3603FullMUGECavityPressureCalculator,
    # Session 62 — grok_share_8d951e12.txt fourth-pass: 6 new classes (PAPER_244–249)
    'MUGEQuantumUncertaintyTermCalculator': MUGEQuantumUncertaintyTermCalculator,
    'MUGEFluidSelfGravityTermCalculator': MUGEFluidSelfGravityTermCalculator,
    'MUGEDualModeOscillatoryGravityCalculator': MUGEDualModeOscillatoryGravityCalculator,
    'MUGEMergerInteractionModulationCalculator': MUGEMergerInteractionModulationCalculator,
    'UQFFSource10BatchProfiledCalculator': UQFFSource10BatchProfiledCalculator,
    'UQFFCUDAGPUOptimizationPatternCalculator': UQFFCUDAGPUOptimizationPatternCalculator,
    # Sessions 63-74
    'SN1006TypeIaSNRFUBiCalculator': SN1006TypeIaSNRFUBiCalculator,
    'EtaCarinaeHomuculusFUBiCalculator': EtaCarinaeHomuculusFUBiCalculator,
    'ChandraArchiveMultiSystemFUBiCalculator': ChandraArchiveMultiSystemFUBiCalculator,
    'SgrACenterNegativeBuoyancyCalculator': SgrACenterNegativeBuoyancyCalculator,
    'KeplerSNR1604FUBiCalculator': KeplerSNR1604FUBiCalculator,
    'PSRJ0030NeutronStarFUBiCalculator': PSRJ0030NeutronStarFUBiCalculator,
    'CrabNebulaM1FUBiCalculator': CrabNebulaM1FUBiCalculator,
    'CassiopeiaASNRFUBiCalculator': CassiopeiaASNRFUBiCalculator,
    'MultiMessengerUQFFValidator': MultiMessengerUQFFValidator,
    'HUDFTRZCPTPhaseCalculator': HUDFTRZCPTPhaseCalculator,
    'HUDFInteractionCascadeBuoyancyCalculator': HUDFInteractionCascadeBuoyancyCalculator,
    'HUDFGravitationalMeissnerCalculator': HUDFGravitationalMeissnerCalculator,
    'NGC1792StarburstBuoyancyCoherenceCalculator': NGC1792StarburstBuoyancyCoherenceCalculator,
    'NGC1792HubbleSlowModeOscillatorCalculator': NGC1792HubbleSlowModeOscillatorCalculator,
    'NGC1792RamPressureDegeneracyCalculator': NGC1792RamPressureDegeneracyCalculator,
    'Source10DPMResonanceAmplificationCalculator': Source10DPMResonanceAmplificationCalculator,
    'Source10GravitationalVacuumDragCalculator': Source10GravitationalVacuumDragCalculator,
    'Source10THzDoubleGateConduitCalculator': Source10THzDoubleGateConduitCalculator,
    'AndromedaBlueshiftApproachAmplifierCalculator': AndromedaBlueshiftApproachAmplifierCalculator,
    'AndromedaHI21cmUQFFResonanceCalculator': AndromedaHI21cmUQFFResonanceCalculator,
    'AndromedaDMShellPartitionCalculator': AndromedaDMShellPartitionCalculator,
    'AndromedaFriedmannHzExpansionCalculator': AndromedaFriedmannHzExpansionCalculator,
    'SombreroRecessionDampingKappaCalculator': SombreroRecessionDampingKappaCalculator,
    'SombreroRingResonatorDustRingCalculator': SombreroRingResonatorDustRingCalculator,
    'SombreroSMBHDominanceRatioCalculator': SombreroSMBHDominanceRatioCalculator,
    'SaturnSolarTidalPerturbationCalculator': SaturnSolarTidalPerturbationCalculator,
    'SaturnRingTidalGravityResonanceCalculator': SaturnRingTidalGravityResonanceCalculator,
    'SaturnAtmosphericWindKineticPressureCalculator': SaturnAtmosphericWindKineticPressureCalculator,
    'SaturnSolarTidalHubbleExpansionCalculator': SaturnSolarTidalHubbleExpansionCalculator,
    'M16DualMassCoActionProductCalculator': M16DualMassCoActionProductCalculator,
    'M16ErosionSaturationHalfTimeCalculator': M16ErosionSaturationHalfTimeCalculator,
    'M16NebularFriedmannRedshiftCalculator': M16NebularFriedmannRedshiftCalculator,
    'ResonanceSCDPMTHzCascadeCalculator': ResonanceSCDPMTHzCascadeCalculator,
    'ResonanceSCCosmicAgeStandingWaveCalculator': ResonanceSCCosmicAgeStandingWaveCalculator,
    'ResonanceSCCooperDPMFreqSynthesisCalculator': ResonanceSCCooperDPMFreqSynthesisCalculator,
    'CrabSNRDPMDilutionCalculator': CrabSNRDPMDilutionCalculator,
    'CrabFilamentSpectralTriadCalculator': CrabFilamentSpectralTriadCalculator,
    'CrabPulsarOscResonanceWindowCalculator': CrabPulsarOscResonanceWindowCalculator,
    'CR24DualChannelArchitectureCalculator': CR24DualChannelArchitectureCalculator,
    'CR24VacuumDifferentialHarmonicCalculator': CR24VacuumDifferentialHarmonicCalculator,
    'CR24CompressedCooperSuperSeedingCalculator': CR24CompressedCooperSuperSeedingCalculator,
    'UniverseDiameterLambdaVacuumAccelerationCalculator': UniverseDiameterLambdaVacuumAccelerationCalculator,
    'UniverseDiameterSuperluminalHubbleRatioCalculator': UniverseDiameterSuperluminalHubbleRatioCalculator,
    'UniverseDiameterGRCurvatureDominanceCalculator': UniverseDiameterGRCurvatureDominanceCalculator,
    'HydrogenAtomLorentzEMDominanceCalculator': HydrogenAtomLorentzEMDominanceCalculator,
    'HydrogenAtomLymanCosmosBridgeCalculator': HydrogenAtomLymanCosmosBridgeCalculator,
    'HydrogenAtomProtonGRSpectralMinimumCalculator': HydrogenAtomProtonGRSpectralMinimumCalculator,
    'HydrogenPToEUg4iResonanceBridgeCalculator': HydrogenPToEUg4iResonanceBridgeCalculator,
    'HydrogenPToETHzQuantumDegeneracyCalculator': HydrogenPToETHzQuantumDegeneracyCalculator,
    'HydrogenPToEAetherGravitationalDominanceCalculator': HydrogenPToEAetherGravitationalDominanceCalculator,
    'LagoonNebulaSFRMassRunawayCalculator': LagoonNebulaSFRMassRunawayCalculator,
    'LagoonNebulaHerschelRadiationErosionCalculator': LagoonNebulaHerschelRadiationErosionCalculator,
    'LagoonNebulaDualRadiationEMBarrierCalculator': LagoonNebulaDualRadiationEMBarrierCalculator,
    'SpiralArmTorqueGravitationalAmplifierCalculator': SpiralArmTorqueGravitationalAmplifierCalculator,
    'SNIaHubbleTensionImprintCalculator': SNIaHubbleTensionImprintCalculator,
    'SpiralDMVisiblePartitionRotationCalculator': SpiralDMVisiblePartitionRotationCalculator,
    'BiPolarPNWindShockGravitationalDominanceCalculator': BiPolarPNWindShockGravitationalDominanceCalculator,
    'BiPolarPNUVRadiationPressureCalculator': BiPolarPNUVRadiationPressureCalculator,
    'EquatorialTorusMagneticConfinementCalculator': EquatorialTorusMagneticConfinementCalculator,
    'BipolarPNLobeResonanceDPMMacroAntennaCalculator': BipolarPNLobeResonanceDPMMacroAntennaCalculator,
    'ResonanceVacDiffTHzCrossoverRadiusCalculator': ResonanceVacDiffTHzCrossoverRadiusCalculator,
    'CooperDPMf1THz_AscConfirmationCalculator': CooperDPMf1THz_AscConfirmationCalculator,
    'OrionTrapeziumWindRamPressureDominanceCalculator': OrionTrapeziumWindRamPressureDominanceCalculator,
    'OrionTrapeziumOBUVRadiationChampagneFlowCalculator': OrionTrapeziumOBUVRadiationChampagneFlowCalculator,
    'OrionCompactHIISFRBindingCrossoverCalculator': OrionCompactHIISFRBindingCrossoverCalculator,
    'CR34DPMForceDensitySpectralAtlasCalculator': CR34DPMForceDensitySpectralAtlasCalculator,
    'CR34CrossChannelDominanceCrossoverCalculator': CR34CrossChannelDominanceCrossoverCalculator,
    'CR34HiIRegionTHzGeometricDifferentialCalculator': CR34HiIRegionTHzGeometricDifferentialCalculator,
    'CR34bVacuumAetherFrequencyModeCalculator': CR34bVacuumAetherFrequencyModeCalculator,
    'CR34bSaturnFirstPlanetaryDualChannelCalculator': CR34bSaturnFirstPlanetaryDualChannelCalculator,
    'CR34bRhoISMFluidDensityCouplingCalculator': CR34bRhoISMFluidDensityCouplingCalculator,
    # Sessions 91-96
    'TriadicMasterFUg1R26StateRamanujanCalculator': TriadicMasterFUg1R26StateRamanujanCalculator,
    'QWave47NonGaussianDistributionCalculator': QWave47NonGaussianDistributionCalculator,
    'AlphaBECNuclearLENREnhancementCalculator': AlphaBECNuclearLENREnhancementCalculator,
    'UmBilinearHeavisideNeutrinoVacuumCascadeCalculator': UmBilinearHeavisideNeutrinoVacuumCascadeCalculator,
    'HResNuclear6EquationDipolekNucCalculator': HResNuclear6EquationDipolekNucCalculator,
    'MUGE26StateFrequencyBasisProofIdentitiesCalculator': MUGE26StateFrequencyBasisProofIdentitiesCalculator,
    'FUBi12TermExplicitIntegrandCalculator': FUBi12TermExplicitIntegrandCalculator,
    'BSMUQFFMultiExperimentCouplingCalculator': BSMUQFFMultiExperimentCouplingCalculator,
    'UiComplexSuperconductiveVacuumDensityCalculator': UiComplexSuperconductiveVacuumDensityCalculator,
    'kkREBTrdicRamanujanFUBiBuoyancyKernelCalculator': kkREBTrdicRamanujanFUBiBuoyancyKernelCalculator,
    'gCompressedAllForcesR26ComponentCalculator': gCompressedAllForcesR26ComponentCalculator,
    'QWave81PhaseSeparationValidationCalculator': QWave81PhaseSeparationValidationCalculator,
    'NineSystemSepAstroParameterCatalogueCalculator': NineSystemSepAstroParameterCatalogueCalculator,
    'UmRotorStringTorqueIntegrationCalculator': UmRotorStringTorqueIntegrationCalculator,
    'EDMSO10BSMRefinedFuCalculator': EDMSO10BSMRefinedFuCalculator,
    'UQFFSupplementCalibration3VarCalculator': UQFFSupplementCalibration3VarCalculator,
    'MagnetarDPMTHzFrequencyFormCalculator': MagnetarDPMTHzFrequencyFormCalculator,
    'SGR17452900SCmLxFreqFormCalculator': SGR17452900SCmLxFreqFormCalculator,
    'SgrAStarGWPrecessionSquaredCalculator': SgrAStarGWPrecessionSquaredCalculator,
    'TapestryStarbirthDPMTHzFreqCalculator': TapestryStarbirthDPMTHzFreqCalculator,
    'M87JetBZModelFUBiCalculator': M87JetBZModelFUBiCalculator,
    'CentaurusAFUBiJetVshapeCalculator': CentaurusAFUBiJetVshapeCalculator,
    'StephansQuintetShockRidgeFUBiCalculator': StephansQuintetShockRidgeFUBiCalculator,
    'SPTClJ2215CoolCoreStarburstCalculator': SPTClJ2215CoolCoreStarburstCalculator,
    'ElGordoACTCLJ0102MergerFUBiCalculator': ElGordoACTCLJ0102MergerFUBiCalculator,
    'ASASSN14liTDEOutflowFUBiCalculator': ASASSN14liTDEOutflowFUBiCalculator,
    'RAquariiSymbioticBinaryFUBiCalculator': RAquariiSymbioticBinaryFUBiCalculator,
    'DecayRateVacuumRhoRatioDoubleExpCalculator': DecayRateVacuumRhoRatioDoubleExpCalculator,
    'DUniverseSpatialCurvatureFifthFactorCalculator': DUniverseSpatialCurvatureFifthFactorCalculator,
}

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
# IMPORT FROM GROKTHREADUQFFEXTENSIONS.PY (Thread 9c3666463ac14753b4f3bea869caaf01)
# 8 UQFF extension calculators + master: 13-term g_res, asymmetrical capacitor,
# variable light speed, fractal time, vacuum probability, 26-layer energies,
# compressed gravity (8 terms), 17 buoyancy proof variants
# ═══════════════════════════════════════════════════════════════════════════════
from GrokThreadUQFFExtensions import (
    GROK_THREAD_UQFF_CALCULATORS,
    UQFFConstants,
    SystemParams,
    ResonanceGravityCalculator,
    AsymmetricalCapacitorCalculator,
    UniversalMagnetismCalculator,
    AetherMetricTensor,
    UnifiedFieldCalculator,
    VariableLightSpeedCalculator,
    FractalTimeCalculator,
    VacuumFluctuationProbability,
    QuantumLevelEnergiesCalculator,
    CompressedGravityCalculator,
    BuoyancyForceProofCalculator,
    GrokThreadUQFFMasterCalculator,
)

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORT FROM CONDENSEDPHYSICS4.PY (EXTENSION 3 — 94 classes, Sessions 97-120)
# Source: gok_share_31b5c807a4 + grok_share_11254865 + grok_share_cfdcad2f5 + grok_share_755feea7
#         + grok_share_c020496d9e + grok_share_5fa36e4e035 + grok_share_dc707f5d3
# ═══════════════════════════════════════════════════════════════════════════════
from CondensedPhysics4 import (
    # Sessions 97-98 — CP4 creation + PSZ2 gap fill (PAPER_355–367)
    PLCKClusterG287MergerRelicTriadicCalculator,
    ASKAPUltraLongPeriodTransientFUBiCalculator,
    TOI1227bYoungNeptuneExoplanetFUBiCalculator,
    AT2024tvdWanderingMBHTDECalculator,
    G359FilamentGalacticCenterFUBiCalculator,
    J1610HighZQuasarJetFUBiCalculator,
    BubbleNebulaPositiveExpansionFUBiCalculator,
    H2OH2RotorPhillipsCSCrossSectionCalculator,
    NOMADMonophotonNeutrinoVacuumCouplingCalculator,
    ALICEMultiplicityCentralityRhoVacRatioCalculator,
    MagnetarMmagOutburstTimescaleCalculator,
    SgrAStarJWST2025FlareOmegaActDerivationCalculator,
    PSZ2G181MergerRelicTriadicFUBiCalculator,
    # Sessions 100-101 — grok_share_11254865 UQFF 2.0 Integration (PAPER_368–375)
    Ug4VacuumEnergyLambdaCDMGalacticBHCouplingCalculator,
    NavierStokesStableFluidUQFFQuasarJetCalculator,
    MultiBodySolarPcorePlanetaryScalingCalculator,
    StarMagic09SeptUQFFMultiBodyNSCalculator,
    MUGESuperconductive12TermResonanceCalculator,
    CompressedUQFFBcritSuperconductivityCalculator,
    MorrisThorneWormholeNullGeodesicsCalculator,
    J1610RelativisticQuasarJetUQFFNSCalculator,
    UQFFWormholeMeissnerRelativisticGammaCalculator,
    # Sessions 102-104 — grok_share_11254865 deep analysis (PAPER_376–386)
    StarMagic11254865MUGESessionHubCalculator,
    UQFFResonanceFormalProofSetCalculator,
    WormholeMUGETermImplSafetyCalculator,
    StarMagic11254865Session102HubCalculator,
    CohesiveUQFFIntegrationCalculator,
    DualModelMUGEComparisonCalculator,
    UQFFSolvableEquationSetCalculator,
    StarMagic11254865Session103HubCalculator,
    SGR1745CompressedMUGESpectralTermDecompositionCalculator,
    UQFF12TermSpectralLadderSGR1745Calculator,
    Ug4iTransientAgeDecayLawCalculator,
    SagAStarFullResonanceTermDecompositionCalculator,
    Canonical7SystemUQFFParameterRegistryCalculator,
    LaTeXDualBlockUQFFMasterEquationCalculator,
    # Sessions 106-107 — grok_share_cfdcad2f5 analysis (PAPER_387–399)
    vSCmRelativisticParameterUpdateCalculator,
    YangMillsMassGapVacuumDensityEvolutionCalculator,
    GalacticOmegaSVelocityDispersionCalibrationCalculator,
    SMBHMassSigmaDispersionRelationUQFFAnchorCalculator,
    HybridMUGEMeissnerBlendingModelCalculator,
    AetherMetricTensorPerturbationCalculator,
    SCmReactorEfficiencyDecayCalculator,
    FUThreeTermStarMagicMasterCalculator,
    WormholeUQFFResonanceAccelerationCalculator,
    HiggsEmergentLevel18UQFFStratumCalculator,
    Session107CfdcAd2f5HubCalculator,
    # Session 108 — grok_share_cfdcad2f5 C++ construction (PAPER_400–408)
    Ug2HeliosphereBubbleChargeCoupledEreactCalculator,
    Ug3MagneticStringsDiskPcoreCalculator,
    Ug4VacuumBHFeedbackCconcentrationCalculator,
    Ubi4TermSolarWindBuoyancyEpsilonSwCalculator,
    MusSCmAugmentedMagneticDipoleOmegaCCalculator,
    SCmDensityPlanetaryScalingLawCalculator,
    Ts00TwoComponentStressEnergyDecompositionCalculator,
    FU4BodySolarSystemNumericalVerificationCalculator,
    ResonanceMUGE14TermCompleteWormholeSumCalculator,
    Session108CfdcAd2f5OctConstructionFileHubCalculator,
    # Session 109 — grok_share_cfdcad2f5 refactoring exhaustion
    Session109CfdcAd2f5RefactoringSectionExhaustionHubCalculator,
    # Sessions 110-111 — grok_share_755feea7 Star Magic book + exhaustive re-analysis (PAPER_410–421)
    SCmHiddenElementUndetectableQsQuasarIgnitionCalculator,
    Ug1DPMDiPseudoMonopoleSolarCalibrationCalculator,
    HeliosphereHydrogenComplexSCmStellarAgeCalculator,
    Ug3CCWCWDifferentialRotationSCmPlanetaryCoreCalculator,
    QuasarJetNavierStokesUQFFFluidSolverBodyForceCalculator,
    EreactSCmReactivityAetherDensityReactorEfficiencyCalculator,
    TsUniverse5ComponentStressEnergyDecompositionCalculator,
    PiCyclesNegativeTimeCosineTemporalReversalCalculator,
    FUSunCompleteSCmSolarCycleFinalCalibrationCalculator,
    HamiltonianPlanetaryCoreHUg3HSCmHUAYangMillsMassGapCalculator,
    Session110Grok755feea7StarMagicBookPhysicsHubCalculator,
    FUCompleteLambdaI4thDissipationSumCalculator,
    UmHeavisideQuasiPeriodicSCmPhaseTransitionAmplifierCalculator,
    Session111Grok755feea7ExhaustiveReanalysisHubCalculator,
    UQFF29SystemCrossValidationMatrixCalculator,
    Session112GrokC020496d9ExhaustiveAuditHubCalculator,
    UmCompleteSSqVacuumThermalDampingCalculator,
    Session113GrokC020496d9ReAnalysisHubCalculator,
    # Session 114 — PAPER_424–429
    FUBiiUmUniversalCompanionCatalogCalculator,
    DPMFourComponentCorrelationCalculator,
    UAScmJWSTALMACERNValidationTableCalculator,
    TwentySixDResonanceLayerAmplitudeFrequencyCalculator,
    HResPeriodicTableUniversalNuclearCorrelationCalculator,
    ThreeNewNumberSystemsVacuumDipoleBuoyancyCalculator,
    Session114GrokC020496d9DeepPhysicsHubCalculator,
    # Session 115 — grok_share_5fa36e4e035.txt UQFF module library (PAPER_447–455)
    OrionNebulaHAlphaUQFFCalculator,
    MultiSystemUQFFCoreCalculator,
    YoungStarsOutflowsPressureCalculator,
    EagleNebulaWindRadiationCalculator,
    BigBangCosmicQGDMGWCalculator,
    CompressedUQFFEnvModularCalculator,
    MagnetarDualModeUQFFCalculator,
    MultiSystemCompressionCycle2Calculator,
    UQFFExpandedSystemRegistryCalculator,
    Session115GrokShare5fa36e4eHubCalculator,
)

CP4_CALCULATORS = {
    'PLCKClusterG287MergerRelicTriadicCalculator': PLCKClusterG287MergerRelicTriadicCalculator,
    'ASKAPUltraLongPeriodTransientFUBiCalculator': ASKAPUltraLongPeriodTransientFUBiCalculator,
    'TOI1227bYoungNeptuneExoplanetFUBiCalculator': TOI1227bYoungNeptuneExoplanetFUBiCalculator,
    'AT2024tvdWanderingMBHTDECalculator': AT2024tvdWanderingMBHTDECalculator,
    'G359FilamentGalacticCenterFUBiCalculator': G359FilamentGalacticCenterFUBiCalculator,
    'J1610HighZQuasarJetFUBiCalculator': J1610HighZQuasarJetFUBiCalculator,
    'BubbleNebulaPositiveExpansionFUBiCalculator': BubbleNebulaPositiveExpansionFUBiCalculator,
    'H2OH2RotorPhillipsCSCrossSectionCalculator': H2OH2RotorPhillipsCSCrossSectionCalculator,
    'NOMADMonophotonNeutrinoVacuumCouplingCalculator': NOMADMonophotonNeutrinoVacuumCouplingCalculator,
    'ALICEMultiplicityCentralityRhoVacRatioCalculator': ALICEMultiplicityCentralityRhoVacRatioCalculator,
    'MagnetarMmagOutburstTimescaleCalculator': MagnetarMmagOutburstTimescaleCalculator,
    'SgrAStarJWST2025FlareOmegaActDerivationCalculator': SgrAStarJWST2025FlareOmegaActDerivationCalculator,
    'PSZ2G181MergerRelicTriadicFUBiCalculator': PSZ2G181MergerRelicTriadicFUBiCalculator,
    'Ug4VacuumEnergyLambdaCDMGalacticBHCouplingCalculator': Ug4VacuumEnergyLambdaCDMGalacticBHCouplingCalculator,
    'NavierStokesStableFluidUQFFQuasarJetCalculator': NavierStokesStableFluidUQFFQuasarJetCalculator,
    'MultiBodySolarPcorePlanetaryScalingCalculator': MultiBodySolarPcorePlanetaryScalingCalculator,
    'StarMagic09SeptUQFFMultiBodyNSCalculator': StarMagic09SeptUQFFMultiBodyNSCalculator,
    'MUGESuperconductive12TermResonanceCalculator': MUGESuperconductive12TermResonanceCalculator,
    'CompressedUQFFBcritSuperconductivityCalculator': CompressedUQFFBcritSuperconductivityCalculator,
    'MorrisThorneWormholeNullGeodesicsCalculator': MorrisThorneWormholeNullGeodesicsCalculator,
    'J1610RelativisticQuasarJetUQFFNSCalculator': J1610RelativisticQuasarJetUQFFNSCalculator,
    'UQFFWormholeMeissnerRelativisticGammaCalculator': UQFFWormholeMeissnerRelativisticGammaCalculator,
    'StarMagic11254865MUGESessionHubCalculator': StarMagic11254865MUGESessionHubCalculator,
    'UQFFResonanceFormalProofSetCalculator': UQFFResonanceFormalProofSetCalculator,
    'WormholeMUGETermImplSafetyCalculator': WormholeMUGETermImplSafetyCalculator,
    'StarMagic11254865Session102HubCalculator': StarMagic11254865Session102HubCalculator,
    'CohesiveUQFFIntegrationCalculator': CohesiveUQFFIntegrationCalculator,
    'DualModelMUGEComparisonCalculator': DualModelMUGEComparisonCalculator,
    'UQFFSolvableEquationSetCalculator': UQFFSolvableEquationSetCalculator,
    'StarMagic11254865Session103HubCalculator': StarMagic11254865Session103HubCalculator,
    'SGR1745CompressedMUGESpectralTermDecompositionCalculator': SGR1745CompressedMUGESpectralTermDecompositionCalculator,
    'UQFF12TermSpectralLadderSGR1745Calculator': UQFF12TermSpectralLadderSGR1745Calculator,
    'Ug4iTransientAgeDecayLawCalculator': Ug4iTransientAgeDecayLawCalculator,
    'SagAStarFullResonanceTermDecompositionCalculator': SagAStarFullResonanceTermDecompositionCalculator,
    'Canonical7SystemUQFFParameterRegistryCalculator': Canonical7SystemUQFFParameterRegistryCalculator,
    'LaTeXDualBlockUQFFMasterEquationCalculator': LaTeXDualBlockUQFFMasterEquationCalculator,
    'vSCmRelativisticParameterUpdateCalculator': vSCmRelativisticParameterUpdateCalculator,
    'YangMillsMassGapVacuumDensityEvolutionCalculator': YangMillsMassGapVacuumDensityEvolutionCalculator,
    'GalacticOmegaSVelocityDispersionCalibrationCalculator': GalacticOmegaSVelocityDispersionCalibrationCalculator,
    'SMBHMassSigmaDispersionRelationUQFFAnchorCalculator': SMBHMassSigmaDispersionRelationUQFFAnchorCalculator,
    'HybridMUGEMeissnerBlendingModelCalculator': HybridMUGEMeissnerBlendingModelCalculator,
    'AetherMetricTensorPerturbationCalculator': AetherMetricTensorPerturbationCalculator,
    'SCmReactorEfficiencyDecayCalculator': SCmReactorEfficiencyDecayCalculator,
    'FUThreeTermStarMagicMasterCalculator': FUThreeTermStarMagicMasterCalculator,
    'WormholeUQFFResonanceAccelerationCalculator': WormholeUQFFResonanceAccelerationCalculator,
    'HiggsEmergentLevel18UQFFStratumCalculator': HiggsEmergentLevel18UQFFStratumCalculator,
    'Session107CfdcAd2f5HubCalculator': Session107CfdcAd2f5HubCalculator,
    'Ug2HeliosphereBubbleChargeCoupledEreactCalculator': Ug2HeliosphereBubbleChargeCoupledEreactCalculator,
    'Ug3MagneticStringsDiskPcoreCalculator': Ug3MagneticStringsDiskPcoreCalculator,
    'Ug4VacuumBHFeedbackCconcentrationCalculator': Ug4VacuumBHFeedbackCconcentrationCalculator,
    'Ubi4TermSolarWindBuoyancyEpsilonSwCalculator': Ubi4TermSolarWindBuoyancyEpsilonSwCalculator,
    'MusSCmAugmentedMagneticDipoleOmegaCCalculator': MusSCmAugmentedMagneticDipoleOmegaCCalculator,
    'SCmDensityPlanetaryScalingLawCalculator': SCmDensityPlanetaryScalingLawCalculator,
    'Ts00TwoComponentStressEnergyDecompositionCalculator': Ts00TwoComponentStressEnergyDecompositionCalculator,
    'FU4BodySolarSystemNumericalVerificationCalculator': FU4BodySolarSystemNumericalVerificationCalculator,
    'ResonanceMUGE14TermCompleteWormholeSumCalculator': ResonanceMUGE14TermCompleteWormholeSumCalculator,
    'Session108CfdcAd2f5OctConstructionFileHubCalculator': Session108CfdcAd2f5OctConstructionFileHubCalculator,
    'Session109CfdcAd2f5RefactoringSectionExhaustionHubCalculator': Session109CfdcAd2f5RefactoringSectionExhaustionHubCalculator,
    'SCmHiddenElementUndetectableQsQuasarIgnitionCalculator': SCmHiddenElementUndetectableQsQuasarIgnitionCalculator,
    'Ug1DPMDiPseudoMonopoleSolarCalibrationCalculator': Ug1DPMDiPseudoMonopoleSolarCalibrationCalculator,
    'HeliosphereHydrogenComplexSCmStellarAgeCalculator': HeliosphereHydrogenComplexSCmStellarAgeCalculator,
    'Ug3CCWCWDifferentialRotationSCmPlanetaryCoreCalculator': Ug3CCWCWDifferentialRotationSCmPlanetaryCoreCalculator,
    'QuasarJetNavierStokesUQFFFluidSolverBodyForceCalculator': QuasarJetNavierStokesUQFFFluidSolverBodyForceCalculator,
    'EreactSCmReactivityAetherDensityReactorEfficiencyCalculator': EreactSCmReactivityAetherDensityReactorEfficiencyCalculator,
    'TsUniverse5ComponentStressEnergyDecompositionCalculator': TsUniverse5ComponentStressEnergyDecompositionCalculator,
    'PiCyclesNegativeTimeCosineTemporalReversalCalculator': PiCyclesNegativeTimeCosineTemporalReversalCalculator,
    'FUSunCompleteSCmSolarCycleFinalCalibrationCalculator': FUSunCompleteSCmSolarCycleFinalCalibrationCalculator,
    'HamiltonianPlanetaryCoreHUg3HSCmHUAYangMillsMassGapCalculator': HamiltonianPlanetaryCoreHUg3HSCmHUAYangMillsMassGapCalculator,
    'Session110Grok755feea7StarMagicBookPhysicsHubCalculator': Session110Grok755feea7StarMagicBookPhysicsHubCalculator,
    'FUCompleteLambdaI4thDissipationSumCalculator': FUCompleteLambdaI4thDissipationSumCalculator,
    'UmHeavisideQuasiPeriodicSCmPhaseTransitionAmplifierCalculator': UmHeavisideQuasiPeriodicSCmPhaseTransitionAmplifierCalculator,
    'Session111Grok755feea7ExhaustiveReanalysisHubCalculator': Session111Grok755feea7ExhaustiveReanalysisHubCalculator,
    'UQFF29SystemCrossValidationMatrixCalculator': UQFF29SystemCrossValidationMatrixCalculator,
    'Session112GrokC020496d9ExhaustiveAuditHubCalculator': Session112GrokC020496d9ExhaustiveAuditHubCalculator,
    'UmCompleteSSqVacuumThermalDampingCalculator': UmCompleteSSqVacuumThermalDampingCalculator,
    'Session113GrokC020496d9ReAnalysisHubCalculator': Session113GrokC020496d9ReAnalysisHubCalculator,
    # Session 114 — PAPER_424–429
    'FUBiiUmUniversalCompanionCatalogCalculator': FUBiiUmUniversalCompanionCatalogCalculator,
    'DPMFourComponentCorrelationCalculator': DPMFourComponentCorrelationCalculator,
    'UAScmJWSTALMACERNValidationTableCalculator': UAScmJWSTALMACERNValidationTableCalculator,
    'TwentySixDResonanceLayerAmplitudeFrequencyCalculator': TwentySixDResonanceLayerAmplitudeFrequencyCalculator,
    'HResPeriodicTableUniversalNuclearCorrelationCalculator': HResPeriodicTableUniversalNuclearCorrelationCalculator,
    'ThreeNewNumberSystemsVacuumDipoleBuoyancyCalculator': ThreeNewNumberSystemsVacuumDipoleBuoyancyCalculator,
    'Session114GrokC020496d9DeepPhysicsHubCalculator': Session114GrokC020496d9DeepPhysicsHubCalculator,
    # Session 115 — PAPER_447–455 (grok_share_5fa36e4e035.txt)
    'OrionNebulaHAlphaUQFFCalculator': OrionNebulaHAlphaUQFFCalculator,
    'MultiSystemUQFFCoreCalculator': MultiSystemUQFFCoreCalculator,
    'YoungStarsOutflowsPressureCalculator': YoungStarsOutflowsPressureCalculator,
    'EagleNebulaWindRadiationCalculator': EagleNebulaWindRadiationCalculator,
    'BigBangCosmicQGDMGWCalculator': BigBangCosmicQGDMGWCalculator,
    'CompressedUQFFEnvModularCalculator': CompressedUQFFEnvModularCalculator,
    'MagnetarDualModeUQFFCalculator': MagnetarDualModeUQFFCalculator,
    'MultiSystemCompressionCycle2Calculator': MultiSystemCompressionCycle2Calculator,
    'UQFFExpandedSystemRegistryCalculator': UQFFExpandedSystemRegistryCalculator,
    'Session115GrokShare5fa36e4eHubCalculator': Session115GrokShare5fa36e4eHubCalculator,
}

# ═══════════════════════════════════════════════════════════════════════════════
# AGGREGATED MASTER REGISTRY
# ═══════════════════════════════════════════════════════════════════════════════

# Aggregate all calculators from all CP modules
ALL_CALCULATORS = {
    **CP2_CALCULATORS,
    # CP3 Extension 2 (219 classes, 15+ categories, Sessions 41-96 — 2026-03-20)
    **CP3_CALCULATORS,
    # CP4 Extension 3 (94 classes, Sessions 97-120 — 2026-03-22)
    **CP4_CALCULATORS,
    # Grok Thread UQFF Extensions (Thread 9c3666463ac14753b4f3bea869caaf01)
    **GROK_THREAD_UQFF_CALCULATORS,
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
AGGREGATOR_VERSION = "2.6.0"
TOTAL_MODULES = 22  # CP1 (1,227 classes, 168,784L), CP2 (600 classes, 45,991L), CP3 (219 classes, 13,944L, Sessions 41-96), CP4 (94 classes, 7,275L, Sessions 97-120), + 10 thread registries
# Updated: Session 120 v4.90 (2026-03-22) — 15 root-level UQFF C++ module pairs (grok_share_dc707f5d3); CP4 84→94 (PAPER_447–455)


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
        # CP3 Extension 2 (112 classes, 15+ categories, Sessions 41-60)
        'CP3_ALL': list(CP3_CALCULATORS.keys()),
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

    # CP3 exports (Extension 2 — 112 classes, Sessions 41-60)
    'CP3_CALCULATORS',
    'SolarWindBubbleVerificationCalculator',
    'HeliopausalBoundaryStepFunctionCalculator',
    'StellarClusterUg3DiskTurbulenceCalculator',
    'StellarUg1DipoleDefectCalculator',
    'ExoplanetAtmosphericMassLossUbCalculator',
    'PlanetaryCoreUg3PenetrationScalingCalculator',
    'WhiteDwarfUQFFGravitationalDecayCalculator',
    'WhiteDwarfDegenerateElectronUiCalculator',
    'KilonovaTransientQWaveParameterCalculator',
    'SupernovaProgenitorNegativeTimeZoneCalculator',
    'NeutronStarCRPIceCubeFluxVerificationCalculator',
    'NeutronStarMergerUbOutflowF_UCalculator',
    'BlackHoleJetFluidAsymmetryRatioCalculator',
    'BlackHoleUg4GalacticFeedbackCalculator',
    'GaiaSgrADistanceErrorAnalysisCalculator',
    'QuasarBlazerLuminosityEreactVerificationCalculator',
    'GalacticCenterUg4KappaDecayCalibrationCalculator',
    'MilkyWayGalacticSpinUb_iCouplingCalculator',
    'GalaxyIMFNucleosynthesisIndexCalculator',
    'GalaxyEquationOfStateUCFCalculator',
    'QuasarJetAsymmetryCosRatioCalculator',
    'QuasarEddingtonExcessJetVelocityCalculator',
    'GalaxyClusterPSZ2UmTurbulenceCalculator',
    'GalaxyClusterPLCKDoubleRelicShearCalculator',
    'TwentySixLevelPolynomialHierarchyFullCalculator',
    'CosmologicalLineFluximeSFRIntegralCalculator',
    'PDGNuclearPolynomialFitVerificationCalculator',
    'DeepFieldShearDeltaTauConstraintCalculator',
    'HighRedshiftJWSTQWaveDeepFieldCalculator',
    'DeepFieldG359ShearNISPConstraintCalculator',
    'QScopeFrequencyResonanceUQFFCalculator',
    'ATLASLHCQuarkEnergyLowNLevelCalculator',
    'VacuumEnergyComponentRatioCalculator',
    'UQFFIPCChainStatusCalculator',

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
