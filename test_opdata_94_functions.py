"""
Test OPData.py compatibility with 94-function outputs
"""

from OPData import OutputDataStore
import json

def test_opdata_94_functions():
    """
    Verify OPData.py can handle 94-equation outputs from QCalc.py
    """
    store = OutputDataStore()
    
    # Simulate 94-function output (3.5× more than previous 27)
    test_data = {
        'query': 'Sagittarius A*',
        'timestamp': '2026-02-13T17:00:00',
        'input_params': {
            'M': 4.154e6,
            'r': 26_000,
            'B': 1e10
        },
        'long_form_equations': [
            {'name': f'SOURCE{14 + i//10}_function_{i % 10}', 'result': i * 1e-8}
            for i in range(94)
        ],
        'solutions': {
            f'eq_{i}': i * 1e-8
            for i in range(94)
        },
        'available_equations': [f'eq_{i}' for i in range(94)],
        'simulation_set': {
            'time_steps': 1000,
            'dt': 0.01
        }
    }
    
    # Test storage
    query_id = store.store(test_data)
    print(f"✓ Stored query_id: {query_id}")
    
    # Test retrieval
    recalled = store.recall(query_id)
    assert recalled is not None, "Failed to recall data"
    print(f"✓ Recalled data successfully")
    
    # Verify equation count
    eq_count = len(recalled['long_form_equations'])
    assert eq_count == 94, f"Expected 94 equations, got {eq_count}"
    print(f"✓ Verified {eq_count} equations")
    
    # Verify solutions count
    sol_count = len(recalled['solutions'])
    assert sol_count == 94, f"Expected 94 solutions, got {sol_count}"
    print(f"✓ Verified {sol_count} solutions")
    
    # Test JSON serialization size
    json_str = json.dumps(recalled, indent=2)
    json_size = len(json_str)
    print(f"✓ JSON size: {json_size:,} bytes ({json_size / 1024:.2f} KB)")
    
    # Check if reasonable size (should be <1 MB for typical output)
    assert json_size < 1_000_000, f"JSON too large: {json_size:,} bytes"
    print(f"✓ JSON size within limits (<1 MB)")
    
    # Test search functionality
    search_results = store.search('Sagittarius')
    assert len(search_results) > 0, "Search failed"
    print(f"✓ Search found {len(search_results)} results")
    
    print("\n" + "=" * 70)
    print("✅ OPData.py VERIFIED: Can handle 94-function outputs!")
    print("=" * 70)
    print(f"\nCapacity Summary:")
    print(f"  - Equations stored: {eq_count}")
    print(f"  - Solutions stored: {sol_count}")
    print(f"  - JSON size: {json_size / 1024:.2f} KB")
    print(f"  - Storage overhead: {json_size / (94 * 8):.1f} bytes/result")
    print(f"  - Search functional: ✓")
    print(f"  - Recall functional: ✓")
    
    return True

if __name__ == "__main__":
    try:
        test_opdata_94_functions()
        print("\n✅ ALL TESTS PASSED - OPData.py is PRODUCTION READY!")
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
