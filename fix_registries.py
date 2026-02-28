"""Fix registry references after removing duplicate classes"""
import re

# Load CP2
with open('CondensedPhysics2.py', 'r', encoding='utf-8') as f:
    cp2 = f.read()

# Classes that were removed and need registry entries removed
removed_classes = [
    'ACEDCEModulatedEnergyCalculator',
    'BulbDrivenPlasmaEnergeticsCalculator',
    'CounterClockwiseDiagonalCycleCalculator',
    'CyclicalConvectionPatternCalculator',
    'ExtendedCyclePatternAnalyzerCalculator',
    'FieldGeneratorResonanceCouplingCalculator',
    'Orb10RefinedFUCalculator',
    'Orb11RefinedFUCalculator',
    'QuadrantTransitionTrackerCalculator',
    'SpookyActionNonLocalTransferCalculator',
    'ThermalGradientDrivenDynamicsCalculator',
    'ThirtyNineFrameSequenceCalculator',
    'ThirtySixFrameSequenceCalculator',
    'WaxCapCoolingDynamicsCalculator',
]

# Remove registry entries
for class_name in removed_classes:
    # Remove from ORB_ANALYSIS_*_CALCULATORS registries
    pattern1 = rf"^\s*'{class_name}':\s*{class_name}\(\),?\n"
    cp2 = re.sub(pattern1, '', cp2, flags=re.MULTILINE)
    
    # Remove from __all__ exports
    pattern2 = rf"^\s*'{class_name}',?\n"
    cp2 = re.sub(pattern2, '', cp2, flags=re.MULTILINE)
    
    print(f"Removed registry/export refs: {class_name}")

# Write updated file
with open('CondensedPhysics2.py', 'w', encoding='utf-8') as f:
    f.write(cp2)

print("\nRegistry references cleaned up")
