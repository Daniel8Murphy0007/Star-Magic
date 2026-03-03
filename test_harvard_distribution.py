"""
Test suite for Harvard BH Mass Distribution Integration
GROK_THREAD_B6D9BC22_ANALYSIS.md - Priority 1 Implementation
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from grok_url_calculators import EPSBHMassFunctionCalculator, CONST

def test_load_harvard_distribution():
    """Test loading Harvard energy_distributions.json"""
    print("\n" + "="*80)
    print("TEST 1: Load Harvard Distribution")
    print("="*80)
    
    calc = EPSBHMassFunctionCalculator()
    dist = calc.load_harvard_distribution()
    
    if dist is None:
        print("❌ FAILED: Could not load Harvard distribution")
        return False
    
    # Check bins
    expected_bins = [2.75, 3.5, 4.5, 5.5, 6.5, 7.5, 8.0]
    if dist['bins'] != expected_bins:
        print(f"❌ FAILED: Expected bins {expected_bins}, got {dist['bins']}")
        return False
    print(f"✅ PASSED: Bins = {dist['bins']}")
    
    # Check peaks
    expected_peaks = [3.9, 5.6, 6.5, 7.0]
    actual_peaks = sorted(list(dist['peaks'].keys()))
    if actual_peaks != expected_peaks:
        print(f"❌ FAILED: Expected peaks {expected_peaks}, got {actual_peaks}")
        return False
    print(f"✅ PASSED: Peaks = {actual_peaks}")
    
    # Check integral
    if abs(dist['integral'] - 0.064) > 1e-6:
        print(f"❌ FAILED: Expected integral 0.064, got {dist['integral']}")
        return False
    print(f"✅ PASSED: Integral = {dist['integral']}")
    
    # Print peak interpretations
    print("\n--- UQFF Peak Interpretations ---")
    for log_M, interp in dist['peaks'].items():
        M_sun = 10**log_M
        print(f"  log(M/M_☉) = {log_M:.1f} (M = {M_sun:.2e} M_☉)")
        print(f"    {interp}")
    
    print("\n✅ TEST 1 PASSED: Harvard distribution loaded successfully")
    return True

def test_peak_interpretation():
    """Test UQFF peak interpretation matching"""
    print("\n" + "="*80)
    print("TEST 2: Peak Interpretation Matching")
    print("="*80)
    
    calc = EPSBHMassFunctionCalculator()
    calc.load_harvard_distribution()
    
    # Test exact matches
    test_cases = [
        (3.9, "Neutron drop"),
        (5.6, "THz coherence"),
        (6.5, "Vacuum feedback"),
        (7.0, "Supermassive BH")
    ]
    
    for log_M, expected_substr in test_cases:
        interp = calc.get_peak_interpretation(log_M)
        if expected_substr not in interp:
            print(f"❌ FAILED: log(M) = {log_M}, expected '{expected_substr}' in '{interp}'")
            return False
        print(f"✅ PASSED: log(M) = {log_M} → {interp[:50]}...")
    
    # Test near-miss (should match within tolerance)
    interp = calc.get_peak_interpretation(3.8, tolerance=0.3)
    if "Neutron drop" not in interp:
        print(f"❌ FAILED: log(M) = 3.8 should match 3.9 peak, got '{interp}'")
        return False
    print(f"✅ PASSED: log(M) = 3.8 matches 3.9 peak (within tolerance)")
    
    # Test far miss (should not match)
    interp = calc.get_peak_interpretation(4.2, tolerance=0.3)
    if "No peak match" not in interp:
        print(f"❌ FAILED: log(M) = 4.2 should not match any peak, got '{interp}'")
        return False
    print(f"✅ PASSED: log(M) = 4.2 correctly returns 'No peak match'")
    
    print("\n✅ TEST 2 PASSED: Peak interpretation matching works correctly")
    return True

def test_eps_integration():
    """Test EPS BH mass function with Harvard integration"""
    print("\n" + "="*80)
    print("TEST 3: EPS BH Mass Function with Harvard Integration")
    print("="*80)
    
    calc = EPSBHMassFunctionCalculator()
    
    # Test dataset for Sgr A* (M ~ 4.3e6 M_sun, log M ~ 6.63)
    dataset = {
        'M_BH': 4.3e6 * CONST['M_sun'],
        'z': 0.0,
        'sigma': 1.0
    }
    
    result = calc.compute(dataset)
    
    # Check basic EPS calculation
    if result['N_cumulative'] <= 0:
        print(f"❌ FAILED: N_cumulative should be positive, got {result['N_cumulative']}")
        return False
    print(f"✅ PASSED: N_cumulative = {result['N_cumulative']:.6e}")
    
    # Check log_M_sun calculation
    expected_log_M = 6.633  # log10(4.3e6)
    if abs(result['log_M_sun'] - expected_log_M) > 0.01:
        print(f"❌ FAILED: Expected log(M) ~ {expected_log_M}, got {result['log_M_sun']}")
        return False
    print(f"✅ PASSED: log(M/M_☉) = {result['log_M_sun']:.3f}")
    
    # Check UQFF interpretation (should match 6.5 or 7.0 peak)
    if result['uqff_interpretation'] == 'Harvard data not loaded':
        print(f"❌ FAILED: Harvard data should be auto-loaded")
        return False
    if "Vacuum feedback" not in result['uqff_interpretation'] and "SMBH" not in result['uqff_interpretation']:
        print(f"❌ FAILED: Expected 6.5/7.0 peak interpretation, got '{result['uqff_interpretation']}'")
        return False
    print(f"✅ PASSED: UQFF interpretation = {result['uqff_interpretation'][:60]}...")
    
    # Check Harvard peaks included
    if 'harvard_peaks' not in result:
        print(f"❌ FAILED: Harvard peaks should be included in result")
        return False
    print(f"✅ PASSED: Harvard peaks included ({len(result['harvard_peaks'])} peaks)")
    
    # Check Harvard integral
    if abs(result['harvard_integral'] - 0.064) > 1e-6:
        print(f"❌ FAILED: Expected integral 0.064, got {result['harvard_integral']}")
        return False
    print(f"✅ PASSED: Harvard integral = {result['harvard_integral']}")
    
    print("\n✅ TEST 3 PASSED: EPS integration with Harvard data successful")
    return True

def test_multi_modal_peaks():
    """Test all 4 Harvard peaks with EPS calculations"""
    print("\n" + "="*80)
    print("TEST 4: Multi-Modal Peak Validation")
    print("="*80)
    
    calc = EPSBHMassFunctionCalculator()
    
    # Test each peak
    peaks = [
        (3.9, 7943.0, "Neutron drop"),       # log(M) = 3.9 → M = 10^3.9 M_sun
        (5.6, 398107.0, "THz coherence"),    # log(M) = 5.6 → M = 10^5.6 M_sun
        (6.5, 3162277.66, "Vacuum feedback"), # log(M) = 6.5 → M = 10^6.5 M_sun
        (7.0, 10000000.0, "Supermassive BH") # log(M) = 7.0 → M = 10^7.0 M_sun
    ]
    
    print("\n--- Peak Validation ---")
    for log_M, M_sun_expected, name in peaks:
        dataset = {
            'M_BH': M_sun_expected * CONST['M_sun'],
            'z': 0.0,
            'sigma': 1.0
        }
        
        result = calc.compute(dataset)
        
        print(f"\nPeak: {name} (log M = {log_M})")
        print(f"  M_BH = {M_sun_expected:.2e} M_☉")
        print(f"  N(>M,z=0) = {result['N_cumulative']:.6e} Mpc^-3")
        print(f"  UQFF: {result['uqff_interpretation'][:70]}...")
        
        # Verify interpretation matches
        if name not in result['uqff_interpretation']:
            print(f"  ❌ FAILED: Expected '{name}' in interpretation")
            return False
        print(f"  ✅ PASSED")
    
    print("\n✅ TEST 4 PASSED: All 4 multi-modal peaks validated")
    return True

def run_all_tests():
    """Run all Harvard distribution tests"""
    print("\n" + "="*80)
    print("HARVARD ENERGY DISTRIBUTION TEST SUITE")
    print("Thread b6d9bc22 Priority 1 Implementation")
    print("="*80)
    
    tests = [
        test_load_harvard_distribution,
        test_peak_interpretation,
        test_eps_integration,
        test_multi_modal_peaks
    ]
    
    passed = 0
    failed = 0
    
    for test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"\n❌ TEST FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*80)
    print("FINAL RESULTS")
    print("="*80)
    print(f"✅ PASSED: {passed}/{len(tests)} tests")
    print(f"❌ FAILED: {failed}/{len(tests)} tests")
    
    if failed == 0:
        print("\n🎉 ALL TESTS PASSED! Harvard distribution integration successful.")
        print("\nNext steps:")
        print("  1. Git commit: 'feat: Add Harvard energy_distributions.json (thread b6d9bc22 Priority 1)'")
        print("  2. Update GROK_THREAD_B6D9BC22_ANALYSIS.md status: ✅ Priority 1 COMPLETE")
        print("  3. Proceed to Priority 2: C++ Advanced Calculator extraction")
    else:
        print("\n⚠️  SOME TESTS FAILED. Please review errors above.")
    
    return failed == 0

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
