"""Validate all 23 new Session 209 CP4 classes (#462-#484)."""
import sys
sys.modules['CondensedPhysics2'] = type(sys)('CondensedPhysics2')
sys.modules['CondensedPhysics3'] = type(sys)('CondensedPhysics3')
exec(open('CondensedPhysics4.py', encoding='utf-8').read())

classes = [
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

passed = 0
failed = 0
for name in classes:
    try:
        cls = eval(name)
        r = cls().compute({})
        assert isinstance(r, dict), "compute() must return dict"
        assert 'primary_equations' in r, "missing primary_equations"
        passed += 1
        print(f'  PASS #{462 + classes.index(name)} {name}')
    except Exception as e:
        failed += 1
        print(f'  FAIL #{462 + classes.index(name)} {name}: {e}')

print(f'\n=== RESULTS: {passed}/23 PASS, {failed}/23 FAIL ===')
