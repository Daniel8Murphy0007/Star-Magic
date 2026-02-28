"""Remove duplicates and rename collisions between CP1 and CP2"""
import re

# Load CP2
with open('CondensedPhysics2.py', 'r', encoding='utf-8') as f:
    cp2 = f.read()

# 14 IDENTICAL duplicates to REMOVE (already exist in CP1)
identical_remove = [
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

# 6 DIFFERENT implementations to RENAME
rename_map = {
    'HydrogenResonanceCalculator': 'HydrogenResonanceOrb36Calculator',
    'IntelligentPlasmoidBehaviorCalculator': 'IntelligentPlasmoidBehaviorOrb11Calculator',
    'MagneticBubbleConfinementCalculator': 'MagneticBubbleConfinementOrb10Calculator',
    'NegativeTimeCalculator': 'NegativeTimeUPCalculator',
    'NegativeTimeOperatorCalculator': 'NegativeTimeOperatorExpCalculator',
    'TotalEnergyBudgetCalculator': 'TotalEnergyBudgetOrb11Calculator',
}

# Remove identical classes
removed_count = 0
for class_name in identical_remove:
    # Remove class definition
    pattern = rf'^class\s+{class_name}\b.*?(?=^class |\Z)'
    match = re.search(pattern, cp2, re.MULTILINE | re.DOTALL)
    if match:
        cp2 = cp2[:match.start()] + cp2[match.end():]
        removed_count += 1
        print(f"Removed: {class_name}")

# Rename different implementations
renamed_count = 0
for old_name, new_name in rename_map.items():
    if old_name in cp2:
        cp2 = cp2.replace(old_name, new_name)
        renamed_count += 1
        print(f"Renamed: {old_name} -> {new_name}")

# Write updated file
with open('CondensedPhysics2.py', 'w', encoding='utf-8') as f:
    f.write(cp2)

print()
print(f"Total removed: {removed_count}")
print(f"Total renamed: {renamed_count}")
print("CondensedPhysics2.py updated")
