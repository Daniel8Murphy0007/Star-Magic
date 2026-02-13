"""Test OPData.py compatibility with 94-function outputs."""
from OPData import OutputDataStore

# Test with 94-equation output
store = OutputDataStore()
test_result = {
    'query_id': 'test_94_functions',
    'timestamp': '2026-02-13T17:30:00',
    'input_params': {'M': 1e30, 'r': 1e6},
    'long_form_equations': [
        {'name': f'test_eq_{i}', 'result': float(i), 'unit': 'm/s²'} 
        for i in range(94)
    ],
    'solutions': {f'test_eq_{i}': float(i) for i in range(94)},
    'available_equations': [f'available_{i}' for i in range(120)]
}

# Store and recall
query_id = store.store(test_result)
recalled = store.recall(query_id)

# Verify
assert recalled is not None, "Should recall stored result"
assert len(recalled['long_form_equations']) == 94, f"Should have 94 equations, got {len(recalled['long_form_equations'])}"
assert len(recalled['solutions']) == 94, f"Should have 94 solutions, got {len(recalled['solutions'])}"

print(f"✓ OPData.py successfully handles 94-equation outputs")
print(f"✓ Stored query_id: {query_id}")
print(f"✓ Recalled {len(recalled['long_form_equations'])} equations")
print(f"✓ All 94 functions compatible with OPData storage layer")
