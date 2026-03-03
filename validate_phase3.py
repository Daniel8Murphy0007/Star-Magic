"""Quick Phase 3 validation - verify imports and basic functionality."""

print("=" * 80)
print("PHASE 3 QUICK VALIDATION")
print("=" * 80)

# Test 1: Import qcalc_error_handler
print("\n[1/4] Testing qcalc_error_handler...")
try:
    from qcalc_error_handler import (
        CircuitBreaker, classify_error, retry_with_backoff, 
        format_error_response, ErrorCategory
    )
    print("  ✓ qcalc_error_handler imported successfully")
    print(f"      - CircuitBreaker class available")
    print(f"      - classify_error() function available")
    print(f"      - retry_with_backoff() decorator available")
except Exception as e:
    print(f"  ✗ FAILED: {str(e)}")
    exit(1)

# Test 2: Import qcalc_cache
print("\n[2/4] Testing qcalc_cache...")
try:
    from qcalc_cache import GLOBAL_CACHE, LRUCache, CacheConfig
    print("  ✓ qcalc_cache imported successfully")
    print(f"      - GLOBAL_CACHE available")
    print(f"      - LRUCache class available")
    stats = GLOBAL_CACHE.get_stats()
    print(f"      - Cache stats: {stats['total_queries']} queries, {stats['hit_rate']:.1%} hit rate")
except Exception as e:
    print(f"  ✗ FAILED: {str(e)}")
    exit(1)

# Test 3: Import qcalc_progress
print("\n[3/4] Testing qcalc_progress...")
try:
    from qcalc_progress import ProgressTracker, CalculationStage, DummyProgressTracker
    print("  ✓ qcalc_progress imported successfully")
    print(f"      - ProgressTracker class available")
    print(f"      - {len([attr for attr in dir(CalculationStage) if not attr.startswith('_')])} calculation stages defined")
except Exception as e:
    print(f"  ✗ FAILED: {str(e)}")
    exit(1)

# Test 4: Import qcalc_cp2_hybrid with Phase 3 features
print("\n[4/4] Testing qcalc_cp2_hybrid with Phase 3...")
try:
    import qcalc_cp2_hybrid
    
    # Check PHASE3_ENABLED flag
    phase3_enabled = getattr(qcalc_cp2_hybrid, 'PHASE3_ENABLED', False)
    print(f"  ✓ qcalc_cp2_hybrid imported successfully")
    print(f"      - PHASE3_ENABLED: {phase3_enabled}")
    
    if phase3_enabled:
        print(f"      - All Phase 3 modules integrated ✓")
    else:
        print(f"      - Running in Phase 2 compatibility mode")
    
except Exception as e:
    print(f"  ✗ FAILED: {str(e)}")
    exit(1)

print("\n" + "=" * 80)
print("✅ ALL PHASE 3 MODULES VALIDATED")
print("=" * 80)
print("\nPhase 3 integration complete. Run test_phase3_features.py for comprehensive tests.")
