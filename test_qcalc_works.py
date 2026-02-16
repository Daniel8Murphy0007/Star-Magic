#!/usr/bin/env python3
"""Quick test to verify QCalc.py actually works"""

from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS
from IPData import create_manual_input

print("=" * 70)
print("QCalc.py FUNCTIONALITY TEST")
print("=" * 70)

# Test solar mass parameters
params = ComputeParams(
    query_name='Solar Test',
    M=1.989e30,  # Solar mass
    r=6.96e8     # Solar radius
)

solver = UnifiedFieldSolver()

# Test solver with solar parameters - it should compute ALL applicable equations
try:
    result = solver.solve(params)
    num_equations = len(result.get('long_form_equations', []))
    solutions = result.get('solutions', {})
    
    print(f"✓ UnifiedFieldSolver: WORKS ({num_equations} equations computed)")
    print(f"  → Solution keys: {list(solutions.keys())}")
    
    # Check specific solutions
    if 'F_U' in solutions:
        print(f"  → UQFF Base (F_U): {solutions['F_U']:.3e} N")
    if 'g_compressed' in solutions:
        print(f"  → UQFF Compressed (g): {solutions['g_compressed']:.3e} m/s²")
    if 'F_U_Bi' in solutions:
        print(f"  → UQFF Buoyant (F_U_Bi): {solutions['F_U_Bi']:.3e} N")
    if 'F_U_Bi_i' in solutions:
        print(f"  → UQFF Master Buoyant (F_U_Bi_i): {solutions['F_U_Bi_i']:.3e} N")
    if 'g_triadic' in solutions:
        print(f"  → UQFF Triadic (26-layer): {solutions['g_triadic']:.3e} m/s²")
    
    print(f"  → Total solutions: {len(solutions)}")
    
except Exception as e:
    print(f"✗ UnifiedFieldSolver: BROKEN - {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("COMPANION FILES TEST")
print("=" * 70)

# Test IPData
try:
    from IPData import InputParameters, recall_input
    print("✓ IPData.py: WORKS")
except Exception as e:
    print(f"✗ IPData.py: BROKEN - {e}")

# Test OPData
try:
    from OPData import OutputDataStore, QUERY_RESULTS
    print("✓ OPData.py: WORKS")
except Exception as e:
    print(f"✗ OPData.py: BROKEN - {e}")

# Test QCalc_validation
try:
    from QCalc_validation import ValidationDataFetcher
    print("✓ QCalc_validation.py: WORKS")
except Exception as e:
    print(f"✗ QCalc_validation.py: BROKEN - {e}")

# Test QCalc_test
try:
    import QCalc_test
    print("✓ QCalc_test.py: WORKS (pytest test file)")
except Exception as e:
    print(f"✗ QCalc_test.py: BROKEN - {e}")

# Test APIFetch
try:
    import APIFetch
    print("✓ APIFetch.py: WORKS")
except Exception as e:
    print(f"✗ APIFetch.py: BROKEN - {e}")

print("\n" + "=" * 70)
print("PHASE FILES STATUS (ARCHITECTURAL VIOLATIONS)")
print("=" * 70)

try:
    import Phase5_Consolidated
    print("⚠ Phase5_Consolidated.py: EXISTS (VIOLATION - should be deleted)")
except ImportError:
    print("✓ Phase5_Consolidated.py: NOT FOUND (correct)")

try:
    import Phase6_Consolidated
    print("⚠ Phase6_Consolidated.py: EXISTS (VIOLATION - should be deleted)")
except ImportError:
    print("✓ Phase6_Consolidated.py: NOT FOUND (correct)")

try:
    import Phase6_Enhanced
    print("⚠ Phase6_Enhanced.py: EXISTS (VIOLATION - should be deleted)")
except ImportError:
    print("✓ Phase6_Enhanced.py: NOT FOUND (correct)")

try:
    import Phase7_Consolidated
    print("⚠ Phase7_Consolidated.py: EXISTS (VIOLATION - should be deleted)")
except ImportError:
    print("✓ Phase7_Consolidated.py: NOT FOUND (correct)")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("QCalc.py core physics: Check results above")
print("Companion files: Check results above")
print("Phase files: Should NOT exist (architectural violations)")
print("=" * 70)
