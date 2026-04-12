#!/usr/bin/env python3
"""Verify CP4 class counts and registration blocks."""
import re

t = open('CondensedPhysics4.py', encoding='utf-8').read()

# Session 210 actual class names
s210_actual = [
    'PhononModifiedChristoffelGeodesicCalc',
    'MasterStellarWindPhononEtCalc',
    'RosetteNebulaNGC2237UQFFCalc',
    'NebulaObservationComparisonUQFFCalc',
    'PhononErgosphereSuperradianceCalc',
    'PhononQPOAccretionDiskCalc',
    'StellarWindBuoyancyLagrangianCalc',
    'PhononJetLaunchingM87SgrACalc',
    'PhononModulatedHawkingTemperatureCalc',
]
s210b_actual = [
    'BHJetModulationFactorLinewidthCalc',
    'JetCollimationLinewidthGammaCalc',
    'PhononNSSpinDownMagneticDipoleCalc',
    'MagnetarSpinDownPhononTimescaleCalc',
    'TidalDeformabilityPhononCorrectionCalc',
    'GW170817PhononStrainDampingCalc',
    'GW190425MassGapPhononSuppressionCalc',
]
s210c_actual = [
    'ExponentialStrainPhononEvolutionCalc',
    'MatchedFilterSNRPhononDampingCalc',
    'SgrAFlareContrastPhononGammaCalc',
    'MonteCarloJetPowerSamplingCalc',
    'InspiralPhaseLagPhononIntegralCalc',
    'M87JetPowerCurveGammaMatchCalc',
]

all_classes = re.findall(r'^class\s+(\w+)', t, re.MULTILINE)
print(f"Total CP4 class definitions: {len(all_classes)}")

for label, lst in [("S210", s210_actual), ("S210b", s210b_actual), ("S210c", s210c_actual)]:
    found = [n for n in lst if n in all_classes]
    missing = [n for n in lst if n not in all_classes]
    print(f"\n{label}: {len(found)}/{len(lst)} classes present")
    for m in missing:
        print(f"  MISSING: {m}")

# Check registration blocks
print("\n--- Registration blocks ---")
blocks = re.findall(r'(_SESSION_\w+_CLASSES)\s*=\s*\[([^\]]*)\]', t)
for name, content in blocks:
    entries = re.findall(r'"(\w+)"', content)
    status = f"{len(entries)} entries" if entries else "EMPTY"
    print(f"  {name}: {status}")

# Summary
print(f"\n--- Summary ---")
print(f"Class definitions: {len(all_classes)}")
print(f"Expected (claimed): 506")
print(f"Gap: {506 - len(all_classes)}")
