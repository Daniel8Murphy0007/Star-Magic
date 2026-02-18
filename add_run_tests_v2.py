#!/usr/bin/env python3
"""
Targeted script to add run_tests() methods to CondensedPhysics.py classes.
Uses specific pattern matching for reliable insertion.
"""

import re

# All 115 classes that need run_tests()
CLASSES_NEEDED = [
    'UQFFScale', 'SystemParams', 'UQFF', 'UQFFMasterFramework', 'DPMModel', 
    'HydrogenEvolutionModel', 'AtomicModelUQFF', 'UniversalGravityModel', 
    'MagneticStringModel', 'UniversalInertiaModel', 'UniversalInertiaVacuumModel',
    'SuperconductiveMaterialVacuumModel', 'StandardModelUQFFModel', 'BuoyancyMassRelationship',
    'UnifiedFieldEquation', 'TimeReversalZoneModel', 'MSigmaRelationModel',
    'FinalParsecProblemModel', 'CGMMetalRetentionModel', 'FastRadioBurstModel',
    'WhittakerDecompositionModel', 'BigBangOriginModel', 'CosmicEggHypergraphModel',
    'PlasmaShieldCaptureModel', 'BlackHolePhasesModel', 'TerahertzHolesModel',
    'InertialOperatorModel', 'CaduceusQuantumWaveModel', 'GlobularClusterStructureModel',
    'HiggsSCmIntegrationModel', 'DEVacuumPowerModel', 'MaxwellComponentFormModel',
    'ProtonSaturationLevelsModel', 'ERBridgeStateTransitionModel', 'MultiScaleGravityModel',
    'AetherBlueQualitiesModel', 'AGNFeedbackModel', 'CGMBaryonFractionModel',
    'CompleteUnifiedFieldModel', 'PseudoMonopoleModel', 'USPRModel',
    'UniversalBuoyancyInteractionModel', 'HillSphereModel', 'OortCloudBoundaryModel',
    'CrystallineGalaxyModel', 'BlackHoleTriangulationModel', 'StellarEquilibriumModel',
    'DensityWaveModel', 'DensityWaveCrystallineCouplingModel', 'USPRStellarConnectionModel',
    'SolarDomainModel', 'BlackHolePseudoMonopoleResonanceModel', 'CorotationResonanceModel',
    'GinzburgLandauFieldModel', 'BogoliubovDeGennesModel', 'QWaveResonanceModel',
    'TemporalDynamicsModel', 'AmplitudeStabilityModel', 'SuperconductingCoherenceModel',
    'GravitationalTimeDilationModel', 'GravitationalRedshiftModel', 'TidalForceModel',
    'BHMFEvolutionModel', 'BondiAccretionModel', 'EddingtonRatioModel',
    'TidalDisruptionEventModel', 'SMBHSpinEvolutionModel', 'SMBHUg1Model',
    'SMBHUg2Model', 'SMBHUg3Model', 'SMBHUg4Model', 'SMBHBulgeGravityModel',
    'SMBHOmegaSGalacticModel', 'SMBHCosmicTimeModel', 'VirgoClusterMassModel',
    'VirgoClusterDarkMatterModel', 'VirgoClusterVirialModel', 'VirgoClusterICMModel',
    'VirgoClusterGravPotentialModel', 'VirgoClusterM87JetModel', 'VirgoClusterTidalStrippingModel',
    'VirgoClusterXRayModel', 'VirgoClusterVelocityDispersionModel', 'NavierStokesVortexModel',
    'LangevinDynamicsModel', 'BrainWaveSubharmonicModel', 'QScopeDataModel',
    'QScopeCalibrationModel', 'ComputationalComplexityUQFFModel', 'Ug1DefectFactorModel',
    'DiskPlaneUnitVectorModel', 'GalacticBlackHoleMassModel', 'SurfaceMagneticFieldModel',
    'SurfaceTemperatureModel', 'AetherVacuumEnergyModel', 'NegativeTimeModel',
    'CorePenetrationModel', 'OuterFieldBubbleModel', 'ReactorEfficiencyModel',
    'ReciprocatingDecayModel', 'SCmPenetrationModel', 'KappaDecayModel',
    'SolarWindModel', 'NavierStokesUQFFProofModel', 'AlphaBECModel',
    'DustYieldExtinctionModel', 'SgrAStarGravityModel', 'ShearChiSquaredModel',
    'TimeAsymmetryModel', 'RProcessOutflowModel', 'VacuumDensityAlignmentModel',
    'RelativisticJetAsymmetryModel', 'YangMillsMassGapModel', 'RiemannHypothesisModel',
    'PvsNPComplexityModel'
]


def generate_run_tests_method(class_name: str, is_enum: bool = False, is_dataclass: bool = False) -> str:
    """Generate a run_tests() staticmethod for a given class."""
    
    if is_enum:
        return f'''
    @staticmethod
    def run_tests() -> dict:
        """Test {class_name} enum values."""
        tests = []
        try:
            values = list({class_name})
            tests.append({{'name': 'Enum values exist', 'passed': len(values) > 0}})
        except Exception as e:
            tests.append({{'name': 'Enum values exist', 'passed': False, 'error': str(e)}})
        try:
            first = list({class_name})[0]
            tests.append({{'name': 'Value access', 'passed': first.value is not None}})
        except Exception as e:
            tests.append({{'name': 'Value access', 'passed': False, 'error': str(e)}})
        all_passed = all(t['passed'] for t in tests)
        return {{'class': '{class_name}', 'tests': tests, 'passed': sum(1 for t in tests if t['passed']), 'total': len(tests), 'all_passed': all_passed}}
'''
    elif is_dataclass:
        return f'''
    @staticmethod
    def run_tests() -> dict:
        """Test {class_name} dataclass."""
        tests = []
        try:
            params = {class_name}(name="Test", M=1.989e30, r=6.96e8)
            tests.append({{'name': 'Instantiation', 'passed': params is not None}})
        except Exception as e:
            tests.append({{'name': 'Instantiation', 'passed': False, 'error': str(e)}})
        try:
            params = {class_name}(name="Test", M=1.989e30, r=6.96e8)
            tests.append({{'name': 'Field access', 'passed': params.M == 1.989e30}})
        except Exception as e:
            tests.append({{'name': 'Field access', 'passed': False, 'error': str(e)}})
        all_passed = all(t['passed'] for t in tests)
        return {{'class': '{class_name}', 'tests': tests, 'passed': sum(1 for t in tests if t['passed']), 'total': len(tests), 'all_passed': all_passed}}
'''
    else:
        return f'''
    @staticmethod
    def run_tests() -> dict:
        """Test {class_name} functionality."""
        tests = []
        try:
            model = {class_name}()
            tests.append({{'name': 'Instantiation', 'passed': model is not None}})
        except Exception as e:
            tests.append({{'name': 'Instantiation', 'passed': False, 'error': str(e)}})
        try:
            model = {class_name}()
            has_method = hasattr(model, 'compute') or hasattr(model, 'calculate') or hasattr(model, 'get_equation_text')
            tests.append({{'name': 'Has compute method', 'passed': has_method}})
        except Exception as e:
            tests.append({{'name': 'Has compute method', 'passed': False, 'error': str(e)}})
        all_passed = all(t['passed'] for t in tests)
        return {{'class': '{class_name}', 'tests': tests, 'passed': sum(1 for t in tests if t['passed']), 'total': len(tests), 'all_passed': all_passed}}
'''


def process_file(filepath: str):
    """Add run_tests() methods to all classes in the file."""
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    lines = content.split('\n')
    
    # Track classes that need modification
    classes_updated = 0
    classes_skipped = 0
    classes_already_have = 0
    
    # Find each class and its boundaries
    class_info = []
    
    for class_name in CLASSES_NEEDED:
        # Find class definition
        pattern = rf'^class {class_name}[\[(:]'
        for i, line in enumerate(lines):
            if re.match(pattern, line):
                # Determine class type
                is_enum = '(Enum)' in line
                is_dataclass = i > 0 and '@dataclass' in lines[i-1]
                
                # Find class end (start of next class or global instance)
                class_start = i
                class_end = len(lines)
                
                for j in range(i + 1, len(lines)):
                    # Check for class ending markers
                    if (re.match(r'^class\s+\w+', lines[j]) or
                        re.match(r'^[A-Z][A-Z0-9_]+\s*=\s*\w+\(', lines[j]) or  # Global instance
                        (lines[j].startswith('# ') and '═' * 10 in lines[j])):  # Section comment
                        class_end = j
                        break
                
                # Check if run_tests already exists in this class
                class_block = '\n'.join(lines[class_start:class_end])
                if 'def run_tests' in class_block:
                    classes_already_have += 1
                    print(f"  SKIP {class_name}: already has run_tests()")
                else:
                    class_info.append({
                        'name': class_name,
                        'start': class_start,
                        'end': class_end,
                        'is_enum': is_enum,
                        'is_dataclass': is_dataclass
                    })
                break
    
    # Sort by position descending (process from end to start to avoid offset issues)
    class_info.sort(key=lambda x: x['start'], reverse=True)
    
    # Process each class
    for info in class_info:
        class_name = info['name']
        class_start = info['start']
        class_end = info['end']
        
        # Generate the run_tests method
        run_tests_code = generate_run_tests_method(
            class_name, 
            is_enum=info['is_enum'],
            is_dataclass=info['is_dataclass']
        )
        
        # Find insertion point - look for the last method or attribute in the class
        # We'll insert before the class ends but after the last actual class content
        
        # Find the proper indentation from the class content
        insert_line = class_end - 1
        
        # Back up past blank lines
        while insert_line > class_start and lines[insert_line].strip() == '':
            insert_line -= 1
        
        # If we're at a closing triple-quote, back up past it
        if insert_line > class_start and lines[insert_line].strip() in ['"""', "'''"]:
            insert_line -= 1
            # Back up to find the start of this docstring
            while insert_line > class_start and '"""' not in lines[insert_line] and "'''" not in lines[insert_line]:
                insert_line -= 1
        
        # Now we're at the last substantial line of the class
        # Insert the run_tests method after this line
        insert_line += 1
        
        # Insert the code
        new_lines = lines[:insert_line] + [run_tests_code] + lines[insert_line:]
        lines = new_lines
        
        classes_updated += 1
        print(f"  ADDED run_tests() to {class_name} at line {insert_line}")
    
    # Write back
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    return {
        'updated': classes_updated,
        'skipped': classes_skipped,
        'already_have': classes_already_have,
    }


if __name__ == '__main__':
    filepath = 'CondensedPhysics.py'
    print(f"Processing {filepath}...")
    print()
    
    result = process_file(filepath)
    
    print()
    print(f"=" * 60)
    print(f"Summary:")
    print(f"  Classes updated:        {result['updated']}")
    print(f"  Already had run_tests:  {result['already_have']}")
    print(f"  Skipped (not found):    {result['skipped']}")
    print(f"=" * 60)
