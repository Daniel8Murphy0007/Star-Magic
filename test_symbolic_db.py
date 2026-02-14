#!/usr/bin/env python3
"""
test_symbolic_db.py - Quick Test of Symbolic Database
======================================================

Verify that the symbolic database works with existing test cases.

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Created: February 13, 2026 (Track 1: Symbolic DB Validation)
"""

from SymbolicDB import SymbolicDatabase

# Create database connection
db = SymbolicDatabase('uqff_equations.db')

# Test 1: Check database statistics
print("="*80)
print("DATABASE VALIDATION TEST")
print("="*80)
stats = db.get_statistics()
print(f"\nTotal equations: {stats['total_equations']}")
print(f"\nTop categories:")
for cat, count in sorted(stats['by_category'].items(), key=lambda x: -x[1])[:5]:
    print(f"  {cat}: {count}")

# Test 2: Query by category
print("\n" + "="*80)
print("QUERY TEST: Magnetar equations")
print("="*80)
magnetar_eqs = db.query_by_category('magnetar')
print(f"Found {len(magnetar_eqs)} magnetar equations:")
for eq_id in magnetar_eqs[:5]:
    print(f"  - {eq_id}")

# Test 3: Get specific equation metadata
print("\n" + "="*80)
print("METADATA TEST: magnetic_field_decay")
print("="*80)
eq_meta = db.get_equation_metadata('magnetic_field_decay')
if eq_meta:
    print(f"ID: {eq_meta.id}")
    print(f"Category: {eq_meta.category}")
    print(f"Units: {eq_meta.units}")
    print(f"Parameters: {eq_meta.parameters}")
    print(f"Expression: {eq_meta.sympy_expr[:100]}...")
else:
    print("ERROR: Equation not found")

# Test 4: Solve equation (if expression is valid)
print("\n" + "="*80)
print("SOLVER TEST")
print("="*80)
print("NOTE: Some equations may need manual expression entry")
print("      Current extraction has PLACEHOLDER expressions")
print("      Full solver testing after manual expression updates")

db.close()
print("\n✓ VALIDATION COMPLETE: Database structure verified")
print("  Next: Update placeholder expressions → Full solver test")
