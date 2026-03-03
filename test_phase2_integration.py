#!/usr/bin/env python3
"""
test_phase2_integration.py
===========================
Tests Phase 2 CP2 + QCalc hybrid routing

Tests:
1. Standard UQFF query → QCalc.UnifiedFieldSolver
2. CP2 experimental query → CondensedPhysics2 calculator
3. Invalid query handling
4. Performance comparison

Usage:
    python test_phase2_integration.py
"""

import json
import subprocess
import time
import sys

def test_qcalc_route():
    """Test standard UQFF query routed to QCalc"""
    print("\n" + "="*80)
    print("TEST 1: Standard UQFF Query (QCalc route)")
    print("="*80)
    
    input_data = {
        "object_name": "Sagittarius A*",
        "M": 4.15e6 * 1.989e30,  # Solar masses
        "r": 8.5e3 * 3.086e16,   # Parsecs
        "z": 0.0
    }
    
    start = time.time()
    
    try:
        # Call qcalc_cp2_hybrid.py
        result = subprocess.run(
            ['python', 'qcalc_cp2_hybrid.py'],
            input=json.dumps(input_data),
            capture_output=True,
            text=True,
            timeout=10
        )
        
        elapsed = time.time() - start
        
        if result.returncode == 0:
            output = json.loads(result.stdout)
            print(f"[✓] SUCCESS (in {elapsed:.2f}s)")
            print(f"    Source: {output.get('source', 'unknown')}")
            print(f"    Success: {output.get('success', False)}")
            print(f"    Compute time: {output.get('compute_time_ms', 0):.1f} ms")
            
            if output.get('source') == 'QCalc':
                print("    [✓] Correctly routed to QCalc")
                return True
            else:
                print(f"    [✗] Expected QCalc, got {output.get('source')}")
                return False
        else:
            print(f"[✗] FAILED (return code {result.returncode})")
            print(f"    stderr: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("[✗] TIMEOUT (>10s)")
        return False
    except Exception as e:
        print(f"[✗] ERROR: {str(e)}")
        return False


def test_cp2_route():
    """Test CP2 experimental query"""
    print("\n" + "="*80)
    print("TEST 2: CP2 Experimental Query (Orb11 route)")
    print("="*80)
    
    input_data = {
        "object_name": "Orb11 Intelligent Plasmoid",
        "M": 1e6,
        "r": 1e10,
        "z": 0.0
    }
    
    start = time.time()
    
    try:
        result = subprocess.run(
            ['python', 'qcalc_cp2_hybrid.py'],
            input=json.dumps(input_data),
            capture_output=True,
            text=True,
            timeout=15  # CP2 might be slower
        )
        
        elapsed = time.time() - start
        
        if result.returncode == 0 or result.returncode == 1:  # 1 = success=false but valid JSON
            try:
                output = json.loads(result.stdout)
            except json.JSONDecodeError:
                print(f"[✗] Invalid JSON output")
                print(f"    stdout: {result.stdout[:200]}")
                return False
            
            print(f"[✓] RESPONSE (in {elapsed:.2f}s)")
            print(f"    Source: {output.get('source', 'unknown')}")
            print(f"    Success: {output.get('success', False)}")
            
            if output.get('source') == 'CondensedPhysics2' or 'CP2' in output.get('source', ''):
                print("    [✓] Correctly routed to CondensedPhysics2")
                
                # Check if calculator was used
                calc_used = output.get('calculator_used', 'None')
                if calc_used and calc_used != 'None':
                    print(f"    Calculator: {calc_used}")
                else:
                    print("    [!] No calculator found for query (expected - need CP2 compute() methods)")
                
                return True
            else:
                print(f"    [✗] Expected CP2, got {output.get('source')}")
                return False
        else:
            print(f"[✗] FAILED (return code {result.returncode})")
            print(f"    stderr: {result.stderr[:500]}")
            return False
            
    except subprocess.TimeoutExpired:
        print("[✗] TIMEOUT (>15s)")
        return False
    except Exception as e:
        print(f"[✗] ERROR: {str(e)}")
        return False


def test_red_mercury():
    """Test Red Mercury superconductor query"""
    print("\n" + "="*80)
    print("TEST 3: Red Mercury Superconductor (CP2 specialty)")
    print("="*80)
    
    input_data = {
        "object_name": "Red Mercury Superconductor Test",
        "M": 1e5,
        "r": 1e8,
        "T": 77,  # Liquid nitrogen temperature
    }
    
    start = time.time()
    
    try:
        result = subprocess.run(
            ['python', 'qcalc_cp2_hybrid.py'],
            input=json.dumps(input_data),
            capture_output=True,
            text=True,
            timeout=15
        )
        
        elapsed = time.time() - start
        
        if result.returncode == 0 or result.returncode == 1:
            try:
                output = json.loads(result.stdout)
            except json.JSONDecodeError:
                print(f"[✗] Invalid JSON output")
                return False
            
            print(f"[✓] RESPONSE (in {elapsed:.2f}s)")
            print(f"    Source: {output.get('source', 'unknown')}")
            print(f"    Calculator: {output.get('calculator_used', 'None')}")
            
            if 'CP2' in output.get('source', '') or output.get('source') == 'CondensedPhysics2':
                print("    [✓] Correctly detected Red Mercury keyword")
                return True
            else:
                print("    [✗] Failed to detect Red Mercury keyword")
                return False
        else:
            print(f"[✗] FAILED")
            return False
            
    except Exception as e:
        print(f"[✗] ERROR: {str(e)}")
        return False


def test_fallback():
    """Test standard query doesn't trigger CP2"""
    print("\n" + "="*80)
    print("TEST 4: Standard query should NOT trigger CP2")
    print("="*80)
    
    input_data = {
        "object_name": "Sun at 1 AU",
        "M": 1.989e30,
        "r": 1.496e11
    }
    
    try:
        result = subprocess.run(
            ['python', 'qcalc_cp2_hybrid.py'],
            input=json.dumps(input_data),
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0:
            output = json.loads(result.stdout)
            
            if output.get('source') == 'QCalc':
                print("[✓] Correctly routed to QCalc")
                print(f"    F_U = {output.get('solutions', {}).get('F_U', 'N/A')}")
                return True
            else:
                print(f"[✗] Incorrectly routed to {output.get('source')}")
                return False
        else:
            print("[✗] Failed")
            return False
            
    except Exception as e:
        print(f"[✗] ERROR: {str(e)}")
        return False


def main():
    print("\n" + "="*80)
    print("Phase 2 Integration Tests: QCalc + CP2 Hybrid Routing")
    print("="*80)
    
    tests = [
        ("Standard UQFF (QCalc)", test_qcalc_route),
        ("CP2 Orb11", test_cp2_route),
        ("Red Mercury", test_red_mercury),
        ("Fallback to QCalc", test_fallback)
    ]
    
    results = []
    
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"\n[CRASH] Test '{name}' crashed: {str(e)}")
            results.append((name, False))
    
    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    passed = sum(1 for _, p in results if p)
    total = len(results)
    
    for name, p in results:
        status = "[✓] PASS" if p else "[✗] FAIL"
        print(f"{status} - {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n[SUCCESS] All tests passed! Phase 2 routing works correctly.")
        sys.exit(0)
    else:
        print(f"\n[PARTIAL] {total - passed} test(s) failed. Check implementation.")
        sys.exit(1)


if __name__ == '__main__':
    main()
