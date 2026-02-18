"""
Apply SelfSimulatingExpandingMixin inheritance to all 111 Model classes.

This script ONLY changes class declarations from:
    class XModel:
to:
    class XModel(SelfSimulatingExpandingMixin):

The init_self_simulating_expanding() call is done automatically via __init_subclass__.

Run with: python apply_mixin_to_models.py
"""

import re

# All 111 Model class names (from grep search results)
MODEL_CLASSES = [
    'DPMModel', 'HydrogenEvolutionModel', 'ProtoNucleusShellModel', 'UniversalGravityModel',
    'UniversalMagnetismModel', 'MagneticStringModel', 'HeavisideComponentModel',
    'HeliosphereThicknessModel', 'UniversalInertiaModel', 'UniversalInertiaVacuumModel',
    'SuperconductiveMaterialVacuumModel', 'StandardModelUQFFModel', 'TimeReversalZoneModel',
    'MSigmaRelationModel', 'FinalParsecProblemModel', 'CGMMetalRetentionModel',
    'CosmicDynamicsIntegrationModel', 'FastRadioBurstModel', 'WhittakerDecompositionModel',
    'BigBangOriginModel', 'CosmicEggHypergraphModel', 'PlasmaShieldCaptureModel',
    'BlackHolePhasesModel', 'TerahertzHolesModel', 'InertialOperatorModel',
    'CaduceusQuantumWaveModel', 'GlobularClusterStructureModel', 'HiggsSCmIntegrationModel',
    'DEVacuumPowerModel', 'MaxwellComponentFormModel', 'ProtonSaturationLevelsModel',
    'ERBridgeStateTransitionModel', 'MultiScaleGravityModel', 'AetherBlueQualitiesModel',
    'AGNFeedbackModel', 'CGMBaryonFractionModel', 'CompleteUnifiedFieldModel',
    'VacuumEnergyDensitySummaryModel', 'PseudoMonopoleModel', 'USPRModel',
    'UniversalBuoyancyInteractionModel', 'HillSphereModel', 'OortCloudBoundaryModel',
    'CrystallineGalaxyModel', 'BlackHoleTriangulationModel', 'StellarEquilibriumModel',
    'DensityWaveModel', 'DensityWaveCrystallineCouplingModel', 'USPRStellarConnectionModel',
    'SolarDomainModel', 'BlackHolePseudoMonopoleResonanceModel', 'CorotationResonanceModel',
    'GinzburgLandauFieldModel', 'BogoliubovDeGennesModel', 'QWaveResonanceModel',
    'TemporalDynamicsModel', 'AmplitudeStabilityModel', 'SuperconductingCoherenceModel',
    'GravitationalTimeDilationModel', 'GravitationalRedshiftModel', 'TidalForceModel',
    'BHMFEvolutionModel', 'BondiAccretionModel', 'EddingtonRatioModel',
    'TidalDisruptionEventModel', 'SMBHSpinEvolutionModel', 'SMBHUg1Model', 'SMBHUg2Model',
    'SMBHUg3Model', 'SMBHUg4Model', 'SMBHBulgeGravityModel', 'SMBHOmegaSGalacticModel',
    'SMBHCosmicTimeModel', 'VirgoClusterMassModel', 'VirgoClusterDarkMatterModel',
    'VirgoClusterVirialModel', 'VirgoClusterICMModel', 'VirgoClusterGravPotentialModel',
    'VirgoClusterM87JetModel', 'VirgoClusterTidalStrippingModel', 'VirgoClusterXRayModel',
    'VirgoClusterVelocityDispersionModel', 'NGC3603StarClusterModel', 'BubbleNebulaModel',
    'AntennaeGalaxiesModel', 'HorseheadNebulaModel', 'NGC1275Model', 'HUDFModel',
    'NGC1792Model', 'SombreroGalaxyModel', 'SaturnModel', 'M16EagleNebulaModel',
    'CrabNebulaModel', 'NGC2264Model', 'UGC10214Model', 'NGC4676Model',
    'RedSpiderNebulaModel', 'NGC3372Model', 'AGCarinaeModel', 'M42Model',
    'TarantulaNebulaModel', 'NGC2841Model', 'MysticMountainModel', 'NegativeTimeModel',
    'AetherVacuumEnergyModel', 'CosmicEggModel', 'SgrAStarGravityModel', 'RetrocausalModel',
    'TRZModel', 'VoidOscillationModel', 'TimeVaryingVacuumModel',
]

def apply_mixin_inheritance():
    """Apply SelfSimulatingExpandingMixin inheritance to all Model classes (declarations only)."""
    
    filepath = 'CondensedPhysics.py'
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    changes = []
    
    # Step 1: Change class declarations to inherit from SelfSimulatingExpandingMixin
    for model_name in MODEL_CLASSES:
        # Pattern: class ModelName:
        old_pattern = f'class {model_name}:'
        new_pattern = f'class {model_name}(SelfSimulatingExpandingMixin):'
        
        if old_pattern in content and new_pattern not in content:
            content = content.replace(old_pattern, new_pattern)
            changes.append(f'Added inheritance: {model_name}')
    
    # Write back
    if changes:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f'✅ Applied {len(changes)} inheritance changes:')
        for change in changes[:20]:
            print(f'   • {change}')
        if len(changes) > 20:
            print(f'   ... and {len(changes)-20} more')
    else:
        print('❌ No changes applied (may already have inheritance)')
    
    return changes


def verify_syntax():
    """Verify the file has valid Python syntax."""
    import ast
    try:
        with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
            code = f.read()
        ast.parse(code)
        print('✅ Syntax verification PASSED')
        return True
    except SyntaxError as e:
        print(f'❌ Syntax error: {e}')
        return False


if __name__ == '__main__':
    print('=' * 70)
    print('Applying SelfSimulatingExpandingMixin to all 111 Model classes')
    print('=' * 70)
    
    changes = apply_mixin_inheritance()
    
    print()
    print('=' * 70)
    print('Verifying syntax...')
    print('=' * 70)
    
    verify_syntax()
