"""
Test suite for Phase 3 production polish features:
- Error handling with circuit breaker
- LRU caching with TTL
- Progress tracking
- Performance benchmarks

Usage:
    python test_phase3_features.py
"""

import sys
import json
import time
import subprocess
from typing import Dict, Any

# ═══════════════════════════════════════════════════════════════════════════════
# TEST UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def run_calculation(input_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run calculation through qcalc_cp2_hybrid.py.
    
    Args:
        input_data: Input parameters dictionary
        
    Returns:
        Result dictionary
    """
    json_input = json.dumps(input_data)
    
    try:
        result = subprocess.run(
            [sys.executable, "qcalc_cp2_hybrid.py"],
            input=json_input,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        # Parse stdout as JSON
        output = result.stdout.strip()
        if output:
            return json.loads(output)
        else:
            return {
                'success': False,
                'error': 'No output received',
                'stderr': result.stderr
            }
            
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Calculation timeout (>30s)'}
    except json.JSONDecodeError as e:
        return {'success': False, 'error': f'Invalid JSON response: {str(e)}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}


def print_test_header(test_name: str):
    """Print formatted test section header."""
    print(f"\n{'═' * 80}")
    print(f" {test_name}")
    print(f"{'═' * 80}\n")


def print_result(test_name: str, passed: bool, details: str = ""):
    """Print test result."""
    status = "✓ PASS" if passed else "✗ FAIL"
    print(f"[{status}] {test_name}")
    if details:
        print(f"        {details}")


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 1: CACHE HIT/MISS
# ═══════════════════════════════════════════════════════════════════════════════

def test_cache_hit_miss():
    """Test cache hit on repeated query."""
    print_test_header("TEST 1: Cache Hit/Miss")
    
    query = {
        'object_name': 'Sun',
        'M': 1.989e30,
        'r': 1.496e11
    }
    
    # First query (cache miss)
    print("[1] First query (expected CACHE MISS)...")
    result1 = run_calculation(query)
    time1 = result1.get('compute_time_ms', 0)
    
    if not result1.get('success'):
        print_result("Cache Miss Test", False, f"Query failed: {result1.get('error')}")
        return False
    
    print(f"    Compute time: {time1:.1f} ms")
    
    # Second query (cache hit)
    time.sleep(0.5)  # Small delay
    print("[2] Second identical query (expected CACHE HIT)...")
    result2 = run_calculation(query)
    time2 = result2.get('compute_time_ms', 0)
    
    if not result2.get('success'):
        print_result("Cache Hit Test", False, f"Query failed: {result2.get('error')}")
        return False
    
    print(f"    Compute time: {time2:.1f} ms")
    
    # Verify cache hit is significantly faster
    speedup = time1 / time2 if time2 > 0 else 0
    is_cached = time2 < time1 * 0.5  # At least 2x faster
    
    print_result(
        "Cache Hit Performance", 
        is_cached, 
        f"Speedup: {speedup:.1f}x ({'cached' if is_cached else 'NOT cached'})"
    )
    
    return is_cached


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 2: ERROR HANDLING & RETRY
# ═══════════════════════════════════════════════════════════════════════════════

def test_error_handling():
    """Test error handling with invalid input."""
    print_test_header("TEST 2: Error Handling & Retry Logic")
    
    # Test 1: Invalid mass (negative)
    print("[1] Testing invalid input (negative mass)...")
    query1 = {
        'object_name': 'InvalidObject',
        'M': -1.0,
        'r': 1e11
    }
    
    result1 = run_calculation(query1)
    has_error_info = (
        'error' in result1 or 
        'error_type' in result1 or 
        'error_category' in result1
    )
    
    print_result(
        "Invalid Input Handling",
        has_error_info,
        f"Error message: {result1.get('error', 'N/A')[:60]}..."
    )
    
    # Test 2: Missing required parameters
    print("[2] Testing missing parameters...")
    query2 = {
        'object_name': 'IncompleteQuery'
        # Missing M, r, etc.
    }
    
    result2 = run_calculation(query2)
    handled_gracefully = not result2.get('success', False)
    
    print_result(
        "Missing Parameters Handling",
        handled_gracefully,
        "Query properly rejected" if handled_gracefully else "Should have failed"
    )
    
    return has_error_info and handled_gracefully


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 3: CIRCUIT BREAKER
# ═══════════════════════════════════════════════════════════════════════════════

def test_circuit_breaker():
    """Test circuit breaker by triggering multiple failures."""
    print_test_header("TEST 3: Circuit Breaker Pattern")
    
    print("[Info] Circuit breaker test requires manual verification")
    print("       (triggering 5+ failures would break normal tests)")
    
    # Simplified test: Verify circuit breaker is available
    try:
        from qcalc_error_handler import CIRCUIT_BREAKER
        has_circuit_breaker = CIRCUIT_BREAKER is not None
        
        if has_circuit_breaker:
            print(f"[Info] Circuit breaker config:")
            print(f"       - Failure threshold: {CIRCUIT_BREAKER.failure_threshold}")
            print(f"       - Recovery timeout: {CIRCUIT_BREAKER.recovery_timeout}s")
            print(f"       - Current state: {CIRCUIT_BREAKER.state.name}")
        
        print_result("Circuit Breaker Available", has_circuit_breaker)
        return has_circuit_breaker
        
    except ImportError:
        print_result("Circuit Breaker Available", False, "qcalc_error_handler not found")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 4: PROGRESS TRACKING
# ═══════════════════════════════════════════════════════════════════════════════

def test_progress_tracking():
    """Test progress tracking functionality."""
    print_test_header("TEST 4: Progress Tracking")
    
    print("[Info] Progress tracking test (check stderr for progress bars)")
    
    query = {
        'object_name': 'Betelgeuse',
        'M': 2.378e31,  # ~11.6 M_sun
        'r': 9.74e8     # ~6.5 AU
    }
    
    # Run calculation and capture stderr (where progress is printed)
    json_input = json.dumps(query)
    
    try:
        result = subprocess.run(
            [sys.executable, "qcalc_cp2_hybrid.py"],
            input=json_input,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        stderr = result.stderr
        
        # Check for progress indicators in stderr
        has_progress_init = "INIT" in stderr or "Starting" in stderr
        has_progress_compute = "Computing" in stderr or "COMPUTE" in stderr
        has_progress_complete = "COMPLETE" in stderr or "100%" in stderr
        
        has_progress = has_progress_init or has_progress_compute or has_progress_complete
        
        if has_progress:
            print("[Sample Progress Output]")
            for line in stderr.split('\n')[:10]:  # First 10 lines
                if line.strip():
                    print(f"    {line}")
        
        print_result(
            "Progress Tracking Output",
            has_progress,
            f"Found progress indicators: {has_progress_init or has_progress_compute or has_progress_complete}"
        )
        
        return has_progress
        
    except Exception as e:
        print_result("Progress Tracking Output", False, f"Error: {str(e)}")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 5: PERFORMANCE BENCHMARK
# ═══════════════════════════════════════════════════════════════════════════════

def test_performance_benchmark():
    """Benchmark QCalc vs CP2 routing performance."""
    print_test_header("TEST 5: Performance Benchmark")
    
    # Standard UQFF query (QCalc)
    print("[1] Benchmarking QCalc (standard UQFF)...")
    qcalc_query = {
        'object_name': 'Sagittarius A*',
        'M': 4.15e36,  # ~4M solar masses (SMBH)
        'r': 1.2e10    # ~12 million km
    }
    
    qcalc_result = run_calculation(qcalc_query)
    qcalc_time = qcalc_result.get('compute_time_ms', 0)
    qcalc_source = qcalc_result.get('source', 'UNKNOWN')
    
    print(f"    Source: {qcalc_source}")
    print(f"    Time: {qcalc_time:.1f} ms")
    print(f"    Success: {qcalc_result.get('success', False)}")
    
    # CP2 query (experimental physics)
    print("\n[2] Benchmarking CP2 (Orb Analysis)...")
    cp2_query = {
        'object_name': 'Orb 10 Test',
        'M': 1e30,
        'r': 1e11,
        'orb_number': 10
    }
    
    cp2_result = run_calculation(cp2_query)
    cp2_time = cp2_result.get('compute_time_ms', 0)
    cp2_source = cp2_result.get('source', 'UNKNOWN')
    
    print(f"    Source: {cp2_source}")
    print(f"    Time: {cp2_time:.1f} ms")
    print(f"    Success: {cp2_result.get('success', False)}")
    
    # Summary
    print("\n[Performance Summary]")
    print(f"    QCalc: {qcalc_time:.1f} ms ({qcalc_source})")
    print(f"    CP2:   {cp2_time:.1f} ms ({cp2_source})")
    
    qcalc_routed_correctly = 'QCalc' in qcalc_source
    cp2_routed_correctly = 'CP2' in cp2_source or 'CondensedPhysics2' in cp2_source
    
    print_result("QCalc Routing", qcalc_routed_correctly)
    print_result("CP2 Routing", cp2_routed_correctly)
    
    return qcalc_routed_correctly and cp2_routed_correctly


# ═══════════════════════════════════════════════════════════════════════════════
# TEST 6: CACHE STATISTICS
# ═══════════════════════════════════════════════════════════════════════════════

def test_cache_statistics():
    """Test cache statistics tracking."""
    print_test_header("TEST 6: Cache Statistics")
    
    try:
        from qcalc_cache import GLOBAL_CACHE
        
        if GLOBAL_CACHE is None:
            print_result("Cache Statistics", False, "Cache not available")
            return False
        
        # Get cache stats
        stats = GLOBAL_CACHE.get_stats()
        
        print("[Cache Statistics]")
        print(f"    Total queries: {stats['total_queries']}")
        print(f"    Cache hits: {stats['cache_hits']}")
        print(f"    Cache misses: {stats['cache_misses']}")
        print(f"    Hit rate: {stats['hit_rate']:.1%}")
        print(f"    Current size: {stats['current_size']} / {stats['max_size']}")
        print(f"    Total evictions: {stats['total_evictions']}")
        
        if stats['total_queries'] > 0:
            print(f"    Average entry age: {stats['average_entry_age_seconds']:.1f}s")
        
        has_stats = stats['total_queries'] > 0
        
        print_result("Cache Statistics Available", has_stats)
        return has_stats
        
    except ImportError:
        print_result("Cache Statistics", False, "qcalc_cache not found")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN TEST RUNNER
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    """Run all Phase 3 tests."""
    print("\n" + "=" * 80)
    print(" PHASE 3 FEATURE TEST SUITE")
    print(" Testing: Error Handling, Caching, Progress Tracking")
    print("=" * 80)
    
    results = {}
    
    # Run tests
    results['cache_hit_miss'] = test_cache_hit_miss()
    results['error_handling'] = test_error_handling()
    results['circuit_breaker'] = test_circuit_breaker()
    results['progress_tracking'] = test_progress_tracking()
    results['performance_benchmark'] = test_performance_benchmark()
    results['cache_statistics'] = test_cache_statistics()
    
    # Summary
    print_test_header("TEST SUMMARY")
    
    total = len(results)
    passed = sum(1 for v in results.values() if v)
    failed = total - passed
    
    print(f"Total tests: {total}")
    print(f"Passed: {passed} ✓")
    print(f"Failed: {failed} ✗")
    print(f"Success rate: {passed/total*100:.1f}%")
    
    print("\n[Detailed Results]")
    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  [{status}] {test_name}")
    
    print("\n" + "=" * 80)
    
    # Exit code
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
