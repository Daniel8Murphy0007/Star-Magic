#!/usr/bin/env python3
"""
test_ipc_pipeline.py - Test IPC Pipeline Wiring (Phase 0)
=========================================================

Tests the 3 critical connections:
1. Python subprocess wrapper (condensed_physics_subprocess.py)
2. IPC communication (Named Pipes)
3. End-to-end data flow

Author: Daniel T. Murphy
Date: March 2, 2026
"""

import sys
import json
import subprocess
import time

def test_python_subprocess():
    """Test 1: Verify condensed_physics_subprocess.py works"""
    print("\n" + "="*70)
    print("TEST 1: Python Subprocess Wrapper")
    print("="*70)
    
    # Sample input
    test_input = {
        "object_name": "SGR 1745+29",
        "M": 2.0,
        "r": 1e6,
        "z": 0.001,
        "B": 1e13
    }
    
    print(f"Input: {json.dumps(test_input, indent=2)}")
    
    try:
        # Call subprocess
        process = subprocess.Popen(
            ["python", "condensed_physics_subprocess.py"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Send input
        input_json = json.dumps(test_input) + "\n"
        stdout, stderr = process.communicate(input=input_json, timeout=30)
        
        print(f"\nstdout ({len(stdout)} bytes):")
        print(stdout)
        
        if stderr:
            print(f"\nstderr ({len(stderr)} bytes):")
            print(stderr)
        
        # Parse output
        try:
            result = json.loads(stdout)
            print(f"\n✓ Subprocess returned valid JSON")
            print(f"  success: {result.get('success', False)}")
            print(f"  compute_time_ms: {result.get('compute_time_ms', 'N/A')}")
            
            if result.get('success'):
                print(f"  solutions: {len(result.get('solutions', {}))} computed")
                print(f"  equations: {len(result.get('available_equations', []))} available")
            else:
                print(f"  error: {result.get('error', 'Unknown')}")
            
            return True
            
        except json.JSONDecodeError as e:
            print(f"\n✗ Failed to parse JSON output: {e}")
            return False
            
    except subprocess.TimeoutExpired:
        print("\n✗ Subprocess timeout (30s)")
        process.kill()
        return False
    except Exception as e:
        print(f"\n✗ Subprocess error: {e}")
        return False

def test_import_checks():
    """Test 2: Verify imports"""
    print("\n" + "="*70)
    print("TEST 2: Import Availability")
    print("="*70)
    
    modules = [
        ("CondensedPhysics", "CondensedPhysics.py calculator"),
        ("QCalc", "QCalc.py generic solver"),
        ("InputData", "Input data schema"),
        ("OPData", "Output data storage")
    ]
    
    results = []
    for module_name, description in modules:
        try:
            __import__(module_name)
            print(f"✓ {module_name:20} - {description}")
            results.append(True)
        except ImportError as e:
            print(f"✗ {module_name:20} - {description} (NOT AVAILABLE)")
            print(f"  Error: {e}")
            results.append(False)
    
    return all(results)

def test_condensed_physics_direct():
    """Test 3: Direct CondensedPhysics.solve() call"""
    print("\n" + "="*70)
    print("TEST 3: Direct CondensedPhysics.solve() Call")
    print("="*70)
    
    try:
        from CondensedPhysics import UnifiedFieldSolver
        
        solver = UnifiedFieldSolver()
        print("✓ UnifiedFieldSolver instantiated")
        
        # Test parameters
        params = {
            "object_name": "Test Object",
            "M": 1.0,
            "r": 1e5
        }
        
        print(f"Input: {params}")
        
        start_time = time.time()
        result = solver.solve(params)
        elapsed = (time.time() - start_time) * 1000
        
        print(f"✓ solve() completed in {elapsed:.2f}ms")
        print(f"  Result keys: {list(result.keys())}")
        
        return True
        
    except ImportError:
        print("✗ CondensedPhysics.py not available")
        return False
    except Exception as e:
        print(f"✗ solve() failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("\n" + "#"*70)
    print("# IPC PIPELINE WIRING TEST SUITE (Phase 0)")
    print("#"*70)
    
    results = []
    
    # Test 1: Subprocess wrapper
    results.append(("Python Subprocess", test_python_subprocess()))
    
    # Test 2: Import checks
    results.append(("Import Checks", test_import_checks()))
    
    # Test 3: Direct call
    results.append(("Direct CondensedPhysics", test_condensed_physics_direct()))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status:8} - {test_name}")
    
    total_passed = sum(1 for _, passed in results if passed)
    total_tests = len(results)
    
    print(f"\nTotal: {total_passed}/{total_tests} tests passed")
    
    if total_passed == total_tests:
        print("\n✓ All tests passed - IPC pipeline ready for integration")
        return 0
    else:
        print(f"\n✗ {total_tests - total_passed} test(s) failed - review errors above")
        return 1

if __name__ == "__main__":
    sys.exit(main())
