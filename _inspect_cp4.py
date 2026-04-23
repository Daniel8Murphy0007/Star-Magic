import ast

with open('CondensedPhysics4.py', encoding='utf-8-sig', errors='ignore') as f:
    content = f.read()
lines = content.splitlines()

targets = [
    '_CP4Calculator',
    'ThreeDIPONonLinearProgressionCalculator',
    'FUBi26thGaussianTruncatedPolynomialBoundCalculator',
    'AldersOlbersParadoxDPMShellFluxCalculator',
    'UQFFCompEigenvalueQuantumGravityLinkageCalculator',
    'UQFFAllFormsEvolutionCatalogueCalculator',
    'UQFFGWAmplitudeLambdaCDMEmergenceCalculator',
    'UQFFMagneticGatewayCosmicFluxCalculator',
    'UQFF26thOrderFactorialBoundsCalculator',
    'UQFFUg26DPolynomialDefectExpansionCalculator',
    'UQFFUbDensityGradient26thDerivativeCalculator',
    'UQFFCompTensorFull26D13DCrossCalculator',
    'UQFFZeroMassAetherVacuumGradientReformulationCalculator',
    'UQFFExoticPocketedShellQuantumFrequencyCalculator',
    'UQFFNGC6278DwarfGalaxyVoidPocketShellCalculator',
    'UQFFMultiSystemJetHypergraphComparisonCalculator',
    'UQFFUniversalInertialOperatorCalculator',
    'NGC1316MergerEvolutionCalculator',
    'AetherResistanceFullUQFFCalculator',
    'GalacticDarkMatterNFWCouplingCalculator',
]

tree = ast.parse(content)
class_nodes = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.ClassDef)}

for cls_name in targets:
    node = class_nodes.get(cls_name)
    if not node:
        print(f'=== {cls_name}: NOT FOUND ===')
        continue
    print(f'=== {cls_name} (line {node.lineno}) ===')
    for item in node.body:
        if isinstance(item, ast.FunctionDef) and item.name == 'compute':
            end = min(item.end_lineno, item.lineno + 20)
            print(f'  compute() sig+body (lines {item.lineno}-{min(item.end_lineno, item.lineno+20)}):')
            for l in lines[item.lineno-1:end]:
                print(f'    {l}')
    print()
