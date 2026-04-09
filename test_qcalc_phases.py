"""Quick test: QCalc.py Phase 5-8 integration validation."""
import sys
sys.path.insert(0, '.')

print("=" * 60)
print("QCalc Phase 5-8 Integration Test")
print("=" * 60)

# 1. Import
try:
    import QCalc
    print("[PASS] QCalc imported")
except Exception as e:
    print(f"[FAIL] QCalc import: {e}")
    sys.exit(1)

# 2. Check phase flags
for name in ['PHASE5_AVAILABLE', 'PHASE6_AVAILABLE', 'PHASE7_AVAILABLE', 'PHASE8_AVAILABLE']:
    val = getattr(QCalc, name, 'MISSING')
    status = "PASS" if val is True else "WARN"
    print(f"[{status}] {name} = {val}")

# 3. Instantiate solver
try:
    solver = QCalc.UnifiedFieldSolver()
    print("[PASS] UnifiedFieldSolver instantiated")
except Exception as e:
    print(f"[FAIL] Solver init: {e}")
    sys.exit(1)

# 4. Check methods exist
for method in ['_compute_phase5_extraction_physics', '_compute_phase6_galaxy_physics',
               '_compute_phase7_cosmological_physics', '_compute_phase8_kozima_ramanujan']:
    exists = hasattr(solver, method) and callable(getattr(solver, method))
    status = "PASS" if exists else "FAIL"
    print(f"[{status}] {method} exists")

# 5. Test Phase 5 directly
print("\n--- Testing Phase 5 directly ---")
try:
    params = QCalc.ComputeParams()
    params.M = 1.989e30
    params.r = 6.96e8
    params.z = 0.57
    params.t = 1.0
    p5_results = solver._compute_phase5_extraction_physics(params)
    print(f"[PASS] Phase 5: {len(p5_results)} equations")
    for eq in p5_results[:5]:
        print(f"    {eq.name}: {eq.result}")
except Exception as e:
    print(f"[FAIL] Phase 5: {e}")
    import traceback
    traceback.print_exc()

# 6. Test Phase 8 directly
print("\n--- Testing Phase 8 directly ---")
try:
    p8_results = solver._compute_phase8_kozima_ramanujan(params)
    print(f"[PASS] Phase 8: {len(p8_results)} equations")
    for eq in p8_results:
        print(f"    {eq.name}: {eq.result}")
except Exception as e:
    print(f"[FAIL] Phase 8: {e}")
    import traceback
    traceback.print_exc()

# 7. Test full solve
print("\n--- Testing full solve() ---")
try:
    result = solver.solve(params)
    equations = result.get('equations', [])
    long_form = result.get('long_form_equations', [])
    solutions = result.get('solutions', {})
    warnings = {k: v for k, v in solutions.items() if k.startswith('_')}
    non_warnings = {k: v for k, v in solutions.items() if not k.startswith('_')}
    
    print(f"  result keys: {list(result.keys())}")
    print(f"  'equations' count: {len(equations)}")
    print(f"  'long_form_equations' count: {len(long_form)}")
    print(f"  'solutions' count: {len(solutions)} ({len(non_warnings)} real + {len(warnings)} warnings)")
    
    # Show equation types
    if equations:
        print(f"\n  First 5 equations:")
        for eq in equations[:5]:
            t = type(eq).__name__
            if hasattr(eq, 'name'):
                print(f"    [{t}] {eq.name}: {eq.result}")
            else:
                print(f"    [{t}] {eq}")
    
    if long_form:
        print(f"\n  First 5 long_form_equations:")
        for eq in long_form[:5]:
            t = type(eq).__name__
            if hasattr(eq, 'name'):
                print(f"    [{t}] {eq.name}: {eq.result}")
            else:
                print(f"    [{t}] {eq}")
    
    # Show warnings
    if warnings:
        print(f"\n  Warnings:")
        for k, v in warnings.items():
            print(f"    {k}: {v}")
    
    # Show first 10 solution keys
    print(f"\n  First 20 solution keys:")
    for i, (k, v) in enumerate(non_warnings.items()):
        if i >= 20:
            print(f"    ... and {len(non_warnings) - 20} more")
            break
        print(f"    {k}: {v}")
    
    # Categorize ALL long_form_equations by source
    categories = {
        'Phase5': [], 'Phase6': [], 'Phase7': [], 'Phase8': [],
        'UQFF_Master': [], 'Enhanced_Gravity': [], 'Ug_Basic': [],
        'BATCH20_21': [], 'Phase3': [], 'Phase4': [], 'Phase1': [],
        'Wolfram': [], 'Other': []
    }
    
    for eq in long_form:
        name = eq.get('name', '') if isinstance(eq, dict) else getattr(eq, 'name', '')
        categorized = False
        
        if 'Phase5' in name or 'S52' in name or 'S57' in name or 'S60' in name or 'nebular' in name.lower() or 'source52' in name or 'source57' in name or 'source60' in name:
            categories['Phase5'].append(name); categorized = True
        elif 'M51' in name or 'NGC1316' in name or 'SMBH_Binary' in name:
            categories['Phase6'].append(name); categorized = True
        elif 'Andromeda' in name or 'Aether' in name or 'NGC346' in name or 'Source82' in name or 'Source86' in name or 'Source87' in name or 'Source88' in name or 'Source89' in name:
            categories['Phase7'].append(name); categorized = True
        elif 'Kozima' in name or 'Ramanujan' in name or 'MockTheta' in name or 'S26' in name or 'Fneutron' in name:
            categories['Phase8'].append(name); categorized = True
        elif 'UQFF' in name:
            categories['UQFF_Master'].append(name); categorized = True
        elif 'enhanced' in name.lower():
            categories['Enhanced_Gravity'].append(name); categorized = True
        elif name.startswith('Ug') or name == 'Um' or name == 'Ub':
            categories['Ug_Basic'].append(name); categorized = True
        elif 'MUGE' in name or 'BSM' in name or 'Universal' in name:
            categories['BATCH20_21'].append(name); categorized = True
        elif 'Phase3' in name or 'phase3' in name:
            categories['Phase3'].append(name); categorized = True
        elif 'Phase4' in name or 'Aether' in name.lower():
            categories['Phase4'].append(name); categorized = True
        elif 'E_' in name or 'Spectrum' in name or 'Level' in name or 'Vacuum' in name:
            categories['Phase1'].append(name); categorized = True
        elif 'Wolfram' in name or 'wolfram' in name:
            categories['Wolfram'].append(name); categorized = True
        
        if not categorized:
            categories['Other'].append(name)
    
    print(f"\n  === EQUATION CATEGORIZATION ({len(long_form)} total) ===")
    for cat, names in categories.items():
        if names:
            print(f"    {cat}: {len(names)} equations")
            for n in names[:5]:
                print(f"      - {n}")
            if len(names) > 5:
                print(f"      ... +{len(names)-5} more")
    
    # Check: solutions that have NO corresponding equation
    eq_names = set()
    for eq in long_form:
        n = eq.get('name', '') if isinstance(eq, dict) else getattr(eq, 'name', '')
        eq_names.add(n)
    
    orphan_solutions = {k: v for k, v in non_warnings.items() if k not in eq_names}
    if orphan_solutions:
        print(f"\n  === ORPHAN SOLUTIONS (in solutions but NOT in equations): {len(orphan_solutions)} ===")
        for k, v in list(orphan_solutions.items())[:20]:
            print(f"    {k}: {v}")
        if len(orphan_solutions) > 20:
            print(f"    ... +{len(orphan_solutions)-20} more")
    
    # Check: equations that have NO corresponding solution
    orphan_eqs = [n for n in eq_names if n and n not in solutions]
    if orphan_eqs:
        print(f"\n  === ORPHAN EQUATIONS (in equations but NOT in solutions): {len(orphan_eqs)} ===")
        for n in orphan_eqs[:20]:
            print(f"    {n}")

except Exception as e:
    print(f"[FAIL] solve(): {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 60)
print("Test complete")
