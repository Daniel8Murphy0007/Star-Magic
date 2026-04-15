#!/usr/bin/env python3
"""
CondensedPhysicsAggregator.py - Unified UQFF Calculator Import Surface
======================================================================

Master aggregation module that provides a unified import surface for all
CondensedPhysics calculator modules. This enables scalable file clustering
while maintaining a single-import API.

ARCHITECTURE:
    CondensedPhysics.py      → Foundation (1,264 base classes, 172,384 lines)
    CondensedPhysics2.py     → Extension 1 (680 classes, 50,893 lines: Orb Analysis 10/11+ + Grok thread extensions + Session 137 _84A767D3 + Session 138 SOURCE179 + Session 151 Millennium Prize)
    CondensedPhysics3.py     → Extension 2 (219 classes, 13,943 lines: 15+ categories, Sessions 41-96)
    CondensedPhysics4.py     → Extension 3 (540 classes, 40,597 lines, Sessions 97-225, v5.72)
    + 29+ standalone physics modules (Sessions 204-225)
    Last updated: Session 225 (2026-04-15) — tracking sync; 1,090 papers; 1,096 PDFs
    CondensedPhysicsAggregator.py → This file (unified API, v4.1.0)

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
# IMPORT FROM CONDENSEDPHYSICS.PY (FOUNDATION - 1,264 base classes)
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
    
    # Orb Analysis_10 registry
    ORB_ANALYSIS_10_PARAMS,
    ORB_ANALYSIS_10_CALCULATORS,
    
    # Orb Analysis_11 registry
    ORB_ANALYSIS_11_PARAMS,
    ORB_ANALYSIS_11_CALCULATORS,
    
    # Session 137 — Wolfram/_84A767D3 PhysicsTerm CP2 wrappers (PAPER_502–508)
    SOURCE_SESSION137_CP2,
    PIInfinityDecoderCalculator_84A767D3,
    WolframFieldUnityCalculator_84A767D3,
    SacredTimePhaseCalculator_84A767D3,
    HypergraphDimensionCalculator_84A767D3,
    BuoyantGravityHypergraphCalculator_84A767D3,
    WSTPBridgeValidationCalculator_84A767D3,

    # Session 138 — SOURCE179 PI Co-Resonance Field (PAPER_509–515)
    SOURCE_SOURCE179_CP2,
    GW150914PCRCalculator,
    PSRJ0437SacredOrbitCalculator,
    EtaCarinaBuoyantPCRCalculator,
    NGC1277HypergraphDimCalculator,
    TON618SacredPhaseCalculator,
    TXS0506PICoSumCalculator,

    # Aggregated CP2 registry (now includes Session 137 + 138 after merge)
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
# IMPORT FROM CONDENSEDPHYSICS4.PY (EXTENSION 3 — 580 classes, Sessions 97-220)
# Wildcard import: CP4 has grown from 103 to 580 classes across Sessions 97-220.
# Dynamic registry construction replaces manual per-class import.
# ═══════════════════════════════════════════════════════════════════════════════
from CondensedPhysics4 import *
import CondensedPhysics4 as _cp4_module
import inspect as _inspect

# Build CP4_CALCULATORS dynamically from all classes in the module
CP4_CALCULATORS = {}
for _name, _obj in _inspect.getmembers(_cp4_module, _inspect.isclass):
    if _name.startswith('_'):
        continue
    if hasattr(_obj, 'compute') or 'Calc' in _name or 'Calculator' in _name:
        CP4_CALCULATORS[_name] = _obj

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORT FROM STANDALONE PHYSICS MODULES (Sessions 204-220)
# Each module is guarded with try/except to avoid import failures breaking
# the aggregator if a module has unresolved dependencies.
# ═══════════════════════════════════════════════════════════════════════════════

STANDALONE_CALCULATORS = {}

def _safe_import(module_name, class_names):
    """Import calculator classes from a standalone module, silently skip on failure."""
    try:
        mod = __import__(module_name)
        for name in class_names:
            obj = getattr(mod, name, None)
            if obj is not None:
                globals()[name] = obj
                STANDALONE_CALCULATORS[name] = obj
    except Exception:
        pass

# Session 204 — Kozima-LENR + Ramanujan + WSTP modules
_safe_import('bcs_superconductivity_uqff', [
    'BCSGapEquation', 'BCSCriticalTemperature', 'CooperPairPhononCoupling',
])
_safe_import('spectral_ladder_26state', [
    'SpectralLadder26State', 'RamanujanAcceleration', 'SpectralLadderPhononMapping',
])
_safe_import('ramanujan_26d_summation', [
    'Ramanujan26DSummation', 'VDSPolylog26',
])
_safe_import('ramanujan_26d_expanded', [
    'ExpandedRamanujan26DCalculator',
])
_safe_import('qgp_ramanujan_application', [
    'QGPVacuumDensityCalculator', 'YangMillsMassGapCalculator',
    'ALICECentralityMultiplicityCalculator', 'ColorDeconfinementPhaseCalculator',
])

# Session 211 — SCm phonon gap
_safe_import('bh_phonon_interaction', [
    'PhononErgosphereSuperradiance', 'PhononModifiedHawkingTemperature',
    'QPOAccretionDiskPhononCoupling', 'PhononModifiedBHEntropy',
    'BlazarErgosphereCoupling',
])
_safe_import('scm_phonon_linewidth', [
    'LinewidthEnetEvolution', 'LinewidthNeutronDrop',
    'LinewidthBuoyancyReversal', 'LinewidthRegimeClassifier',
])
_safe_import('scm_phonon_resonance', [
    'ResonanceAcceleration', 'LinewidthGammaSweep', 'VacuumDensityCoupling',
    'ResonanceFrequencyScan', 'PhononDampingEvolution',
    'MultiLayerPhononGravityCoupling', 'KozimaNeutronDropLinewidth',
])

# Session 212 — Linewidth gap + jet modulation
_safe_import('blazar_jet_phonon', [
    'BlazarErgosphereResonance', 'BlazarJetPowerGammaCurve',
    'BlazarMultiMessengerPhononCorrelation', 'CentaurusAJetPhononCoupling',
    'TXS0506JetPhononCoupling', 'EtLinewidthModulation',
])
_safe_import('quasar_jet_phonon', [
    'JetModulationFactor', 'BlandfordZnajekJetPower',
    'JetPowerGammaSweep', 'JetPhononWSTPExporter',
])
_safe_import('nuclear_um_jwst_synthesis', ['UniversalMagnetismCalculator'])
_safe_import('linewidth_jet_modulation', [
    'LinewidthJetModulationSweep', 'CollimationPowerMapping', 'ReferenceSystemMatcher',
])

# Session 213 — AGN/blazar jet power curves + SMBH mergers
_safe_import('agn_jet_power_curves', [
    'MonteCarloJetPowerSampler', 'MultiAGNGammaSweep', 'MultiAGNMonteCarloBatch',
])
_safe_import('blazar_jet_power_curves_extended', [
    'CentaurusAJetPowerCurves', 'TXS0506JetPowerCurves', 'DualBlazarJetComparison',
])
_safe_import('smbh_binary_mergers', [
    'SMBHBinaryMergerPhonon', 'MergerStrainDamping',
    'MergerPhaseLag', 'MergerLagrangianVariation',
])

# Session 214 — BCS superconductivity + spectral ladder + E(t) modulation
_safe_import('triadic_solutions_next', [
    'CompressedGravityTriadic', 'ResonantGravityTriadic',
    'BuoyancyGravityTriadic', 'TriadicSolverNext',
])
_safe_import('triadic_validations_next', [
    'QGPTriadicValidator', 'NinetyNineSystemTriadicValidator', 'ALICETriadicCrossCheck',
])

# Session 215 — 26D Ramanujan + 3D MUGE + NS phonon + v11 scaling
_safe_import('muge_magnetar_3d_sim', [
    'SCmCoreModel', 'MagneticVortexModel', 'PhononResonanceShells', 'MUGEMagnetar3DSim',
])
_safe_import('muge_cluster_3d_sim', ['Galaxy3D', 'MUGECluster3DSim'])

# Session 216 — QGP + 99-system + cluster 3D
_safe_import('99system_master_equation', ['NinetyNineSystemMasterCalc', 'NinetyNineSystemAggregateCalc'])

# Session 217 — F_U_Bi_i master buoyancy
_safe_import('fubi_master_calculator', ['FUBiMasterCalculator'])

# Session 218 — F_U_Bi inside/outside + 99sys gamma sweep + v13
_safe_import('fubi_inside_outside', [
    'FUBiInsideOutsideCalc', 'FUBiDistinctionCalc', 'SolarCalibration147Calc',
])

# Session 219 — AGN/NS merger + SCm-QGP + revised curves + v14
_safe_import('fubi_agn_ns_mergers', [
    'AGNFUBiMergerCalc', 'NSMergerFUBiCalc', 'SMBHBinaryMergerFUBiCalc',
    'AGNAccretionBuoyancyCalc', 'SpectralLadderMergerCalc',
])
_safe_import('scm_qgp_dynamics', [
    'QGPVacuumDensityCalc', 'YangMillsMassGapCalc',
    'ALICEMultiplicityCalc', 'DeconfinementPhaseDiagramCalc',
])
_safe_import('fubi_i_curves_revised', [
    'CentaurusARevisedCurvesCalc', 'GW190425RevisedCurvesCalc',
    'TXS0506RevisedCurvesCalc', 'WormholeGeodesicBatchCalc',
])

# Session 220 — 3C273/TON618 AGN + GW170817 + SMBH merger + SCm DM halos + TXS0506 revised
_safe_import('fubi_i_curves_agn_ns_qgp', [
    'ThreeCTwoSevenThreeAGNCurvesCalc', 'TON618AGNCurvesCalc',
    'GW170817MergerCurvesCalc', 'GW190425UpgradedCurvesCalc',
    'QGPALICECentralityCurvesCalc',
])
_safe_import('fubi_smbh_mergers', [
    'SMBHInspiralFUBiCalc', 'SMBHCoalescenceFUBiCalc', 'SMBHRingdownFUBiCalc',
])
_safe_import('scm_dm_halos', [
    'SCmDMHaloDensityCalc', 'RotationCurveFlatteningCalc', 'HaloStabilizationCalc',
])
_safe_import('fubi_i_txs0506_revised', [
    'TXS0506ExtremeFlareCalc', 'TXS0506IceCubeCalc',
    'TXS0506SustainedEmissionCalc', 'TXS0506ThreeGammaProfileCalc',
])

# ═══════════════════════════════════════════════════════════════════════════════
# AGGREGATED MASTER REGISTRY
# ═══════════════════════════════════════════════════════════════════════════════

# Aggregate all calculators from all CP modules
ALL_CALCULATORS = {
    **CP2_CALCULATORS,
    # CP3 Extension 2 (219 classes, 15+ categories, Sessions 41-96 — 2026-03-20)
    **CP3_CALCULATORS,
    # CP4 Extension 3 (580 classes, Sessions 97-220 — 2026-04-13, dynamically built)
    **CP4_CALCULATORS,
    # Standalone physics modules (Sessions 204-220)
    **STANDALONE_CALCULATORS,
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
AGGREGATOR_VERSION = "4.0.0"
TOTAL_MODULES = 35  # CP1 (1,227 classes) + CP2 (668 classes) + CP3 (219 classes) + CP4 (580 classes, v5.75 Session 220) + 10 thread registries + 29 standalone physics modules
# Updated: Session 220 v5.75 (2026-04-13) — Housekeeping: aggregator catch-up; wildcard CP4 import (540→580 classes); 29 standalone modules added (Sessions 204-220); dynamic CP4_CALCULATORS registry; 1018/1000 papers; 1033 PDFs


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
        # CP3 Extension 2 (219 classes, 15+ categories, Sessions 41-96)
        'CP3_ALL': list(CP3_CALCULATORS.keys()),
        # CP4 Extension 3 (580 classes, Sessions 97-220, dynamically built)
        'CP4_ALL': list(CP4_CALCULATORS.keys()),
        # Standalone physics modules (Sessions 204-220)
        'STANDALONE': list(STANDALONE_CALCULATORS.keys()),
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
    'CP4_CALCULATORS',
    'STANDALONE_CALCULATORS',
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

    # Session 137 CP2 wrappers (_84A767D3 series)
    'SOURCE_SESSION137_CP2',
    'PIInfinityDecoderCalculator_84A767D3',
    'WolframFieldUnityCalculator_84A767D3',
    'SacredTimePhaseCalculator_84A767D3',
    'HypergraphDimensionCalculator_84A767D3',
    'BuoyantGravityHypergraphCalculator_84A767D3',
    'WSTPBridgeValidationCalculator_84A767D3',

    # Session 138 SOURCE179 PI Co-Resonance Field calculators
    'SOURCE_SOURCE179_CP2',
    'GW150914PCRCalculator',
    'PSRJ0437SacredOrbitCalculator',
    'EtaCarinaBuoyantPCRCalculator',
    'NGC1277HypergraphDimCalculator',
    'TON618SacredPhaseCalculator',
    'TXS0506PICoSumCalculator',
    # Orb Analysis_10
    'ORB_ANALYSIS_10_PARAMS',
    'ORB_ANALYSIS_10_CALCULATORS',
    
    # Orb Analysis_11
    'ORB_ANALYSIS_11_PARAMS',
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

    # Session 209 — CP4 #462-484 (PAPER_878–900)
    'SCmGaussianActivationBFieldSuppressionCalc',
    'BuoyancyKleinGordonScalarFieldEOMCalc',
    'PositiveEtBuoyancyExpansionMasterCalc',
    'KozimaExpansionNeutronDropCouplingCalc',
    'ExpansionLagrangianEulerLagrangeCalc',
    'NegativeEtBuoyancyErosionMasterCalc',
    'NetEnergyEplusEminusEvolutionCalc',
    'GWDampingErosion66PercentCalc',
    'ErosionLagrangianEulerLagrangeCalc',
    'UQFFVsStringTheory10AspectComparisonCalc',
    'EtFullLagrangianUnifiedDerivationCalc',
    'EtVsLambdaCDMDarkEnergyContrastCalc',
    'SCmVacuumDensityEvolutionCalc',
    'SCmNetEnergyBuoyancyRegimeCalc',
    'SCmKozimaPhononResonanceCouplingCalc',
    'SCmPhononModulatedEnergyPhiCalc',
    'SCmEtLagrangianVariationCalc',
    'EtVsQuintessenceScalarFieldContrastCalc',
    'PhononModulationFactor125THzGaussianCalc',
    'PhononModulatedEnergyEnetPhononCalc',
    'PhononLagrangianPhiS26DerivationCalc',
    'BuoyancyReversalSignFlipResonanceCalc',
    'EtVsKEssenceScherrerModelContrastCalc',
]

# Dynamically extend __all__ with CP4 and standalone calculator names
__all__ += list(CP4_CALCULATORS.keys()) + list(STANDALONE_CALCULATORS.keys())
