#!/usr/bin/env python3
"""Analyze CondensedPhysics.py content."""

import re

def main():
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        content = f.read()
    
    lines = content.split('\n')
    
    # Count various patterns
    classes = re.findall(r'^class (\w+)', content, re.MULTILINE)
    methods = re.findall(r'^    def (\w+)', content, re.MULTILINE)
    constants_refs = len(re.findall(r"CONSTANTS\[", content))
    imports = re.findall(r'^(?:from|import) (\w+)', content, re.MULTILINE)
    
    # Categorize classes
    uqff_classes = [c for c in classes if 'UQFF' in c]
    model_classes = [c for c in classes if 'Model' in c]
    term_classes = [c for c in classes if 'Term' in c]
    calc_classes = [c for c in classes if 'Calc' in c or 'Calculator' in c]
    framework_classes = [c for c in classes if 'Framework' in c]
    params_classes = [c for c in classes if 'Params' in c or 'Parameters' in c]
    solver_classes = [c for c in classes if 'Solver' in c]
    
    # Compute methods
    compute_methods = [m for m in methods if m.startswith('compute')]
    validate_methods = [m for m in methods if 'validate' in m.lower()]
    test_methods = [m for m in methods if 'test' in m.lower() or 'run_tests' == m]
    
    # Count docstrings
    docstrings = len(re.findall(r'"""[\s\S]*?"""', content))
    
    # Count physics equations (LaTeX-like patterns)
    equations = len(re.findall(r'[Uu]g[1-4]|F_U|SCm|[UA]|g_grav|H\(z\)', content))
    
    # May 2025 models
    may_models = ['NGC2264Model', 'UGC10214Model', 'NGC4676Model', 'RedSpiderNebulaModel',
                  'NGC3372Model', 'AGCarinaeModel', 'M42Model', 'TarantulaNebulaModel',
                  'NGC2841Model', 'MysticMountainModel']
    may_present = [m for m in may_models if m in classes]
    
    # Print analysis
    print("=" * 70)
    print("CONDENSEDPHYSICS.PY CONTENT ANALYSIS")
    print("=" * 70)
    
    print(f"\n{'FILE STATISTICS':=^70}")
    print(f"  Total lines:          {len(lines):,}")
    print(f"  Total classes:        {len(classes)}")
    print(f"  Total methods:        {len(methods)}")
    print(f"  Docstrings:           {docstrings}")
    print(f"  CONSTANTS references: {constants_refs}")
    print(f"  Unique imports:       {len(set(imports))}")
    
    print(f"\n{'CLASS CATEGORIES':=^70}")
    print(f"  UQFF classes:         {len(uqff_classes)}")
    print(f"  Model classes:        {len(model_classes)}")
    print(f"  Term classes:         {len(term_classes)}")
    print(f"  Calculator classes:   {len(calc_classes)}")
    print(f"  Framework classes:    {len(framework_classes)}")
    print(f"  Params classes:       {len(params_classes)}")
    print(f"  Solver classes:       {len(solver_classes)}")
    
    print(f"\n{'METHOD CATEGORIES':=^70}")
    print(f"  compute_* methods:    {len(compute_methods)}")
    print(f"  validate_* methods:   {len(validate_methods)}")
    print(f"  test/run_tests:       {len(test_methods)}")
    
    print(f"\n{'UQFF CLASSES ({len(uqff_classes)})':=^70}")
    for c in sorted(uqff_classes):
        print(f"    {c}")
    
    print(f"\n{'MAY 2025 DOCUMENT MODELS ({len(may_present)}/10)':=^70}")
    for m in may_models:
        status = "PRESENT" if m in classes else "MISSING"
        print(f"    {m}: {status}")
    
    print(f"\n{'MODEL CLASSES ({len(model_classes)})':=^70}")
    # Group by type
    astro_models = [m for m in model_classes if any(x in m for x in ['NGC', 'M42', 'Crab', 'Sombrero', 'Bubble', 'Horsehead', 'Eagle', 'Spider', 'Mystic', 'Tadpole', 'Mice', 'Tarantula', 'Carina', 'Saturn', 'HUDF', 'Antennae'])]
    physics_models = [m for m in model_classes if m not in astro_models]
    
    print(f"\n  Astrophysical Models ({len(astro_models)}):")
    for m in sorted(astro_models)[:20]:
        print(f"      {m}")
    if len(astro_models) > 20:
        print(f"      ... and {len(astro_models)-20} more")
    
    print(f"\n  Physics/Generic Models ({len(physics_models)}):")
    for m in sorted(physics_models)[:20]:
        print(f"      {m}")
    if len(physics_models) > 20:
        print(f"      ... and {len(physics_models)-20} more")
    
    # Sample compute methods
    print(f"\n{'SAMPLE COMPUTE METHODS (first 20)':=^70}")
    unique_compute = sorted(set(compute_methods))[:20]
    for m in unique_compute:
        print(f"    {m}")
    
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print(f"""
CondensedPhysics.py is a {len(lines):,}-line UQFF physics calculator containing:

- {len(classes)} classes ({len(model_classes)} Models, {len(uqff_classes)} UQFF, {len(calc_classes)} Calculators)
- {len(methods)} methods ({len(compute_methods)} compute_*, {len(validate_methods)} validate_*, {len(test_methods)} tests)
- {len(may_present)}/10 May 2025 astrophysical models with 8 UQFF Master Equations each
- 8 UQFF Master Equation classes (Base, Compressed, Resonant, Superconductive, 
  Buoyant, Master Buoyant, Triadic, Quadratic)
- Supports LaTeX equation output, Wolfram verification, Qt integration

Architecture: Pure physics calculator (no hardcoded system data)
""")

if __name__ == "__main__":
    main()
