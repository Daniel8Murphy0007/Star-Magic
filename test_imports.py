#!/usr/bin/env python3
"""Quick import test to diagnose CondensedPhysics import issues"""

import time
import sys

print("=" * 70)
print("IMPORT DIAGNOSTIC TEST")
print("=" * 70)

# Test 1: Basic imports
print("\n[1] Testing basic imports...")
start = time.time()
try:
    import numpy as np
    import json
    print(f"✓ Basic imports successful ({time.time()-start:.2f}s)")
except Exception as e:
    print(f"✗ Basic imports failed: {e}")
    sys.exit(1)

# Test 2: Grok modules
print("\n[2] Testing Grok modules...")
start = time.time()
try:
    from grok_100_equations_module import PhysicsConstants
    print(f"✓ grok_100_equations_module loaded ({time.time()-start:.2f}s)")
except Exception as e:
    print(f"⚠ grok_100_equations_module not available: {e}")

start = time.time()
try:
    from grok_100_equations_module_part2 import MUGECalculator
    print(f"✓ grok_100_equations_module_part2 loaded ({time.time()-start:.2f}s)")
except Exception as e:
    print(f"⚠ grok_100_equations_module_part2 not available: {e}")

# Test 3: CondensedPhysics
print("\n[3] Testing CondensedPhysics import...")
start = time.time()
try:
    from CondensedPhysics import UnifiedFieldSolver
    elapsed = time.time() - start
    print(f"✓ CondensedPhysics.UnifiedFieldSolver loaded ({elapsed:.2f}s)")
    
    if elapsed > 10:
        print(f"⚠ WARNING: Import took {elapsed:.2f}s (very slow!)")
except Exception as e:
    print(f"✗ CondensedPhysics import failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 4: Instantiate solver
print("\n[4] Testing UnifiedFieldSolver instantiation...")
start = time.time()
try:
    solver = UnifiedFieldSolver()
    print(f"✓ Solver instantiated ({time.time()-start:.2f}s)")
except Exception as e:
    print(f"✗ Solver instantiation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 70)
print("✓ All diagnostics passed")
print("=" * 70)
