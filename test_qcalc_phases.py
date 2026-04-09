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
    solutions = result.get('solutions', {})
    warnings = {k: v for k, v in solutions.items() if k.startswith('_')}
    print(f"[PASS] solve(): {len(equations)} equations, {len(solutions)} solutions")
    if warnings:
        for k, v in warnings.items():
            print(f"  {k}: {v}")
except Exception as e:
    print(f"[FAIL] solve(): {e}")

print("\n" + "=" * 60)
print("Test complete")
