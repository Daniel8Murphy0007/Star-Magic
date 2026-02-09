#!/usr/bin/env python3
"""
Test Document 11: IceCube pp/pγ SED < 0.1 PeV Verification
"""

import sys
import numpy as np
sys.path.insert(0, '.')
from CondensedPhysics import CRP_TERM_MODEL

def test_icecube_sed():
    """Run all tests for Document 11: IceCube pp/pγ SED."""
    print("=" * 70)
    print("ICECUBE pp/pγ SED < 0.1 PeV - VERIFICATION TESTS")
    print("=" * 70)
    
    # Run built-in IceCube tests
    results = CRP_TERM_MODEL.run_IceCube_tests()
    
    for result in results['results']:
        print(result)
    
    print("=" * 70)
    print(f"RESULTS: {results['tests_passed']}/{results['tests_total']} TESTS PASSED")
    if results['all_passed']:
        print("✓ ALL TESTS PASSED - Document 11 implementation verified")
    else:
        print(f"✗ {results['tests_total'] - results['tests_passed']} TEST(S) FAILED")
    print("=" * 70)
    
    return results['all_passed']

if __name__ == "__main__":
    success = test_icecube_sed()
    sys.exit(0 if success else 1)
