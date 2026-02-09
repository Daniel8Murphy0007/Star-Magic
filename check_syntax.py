#!/usr/bin/env python3
"""Find syntax issues in CondensedPhysics.py"""

with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Check for docstrings around line 2659
print("Lines around 2659:")
for i in range(2645, 2665):
    line = lines[i]
    has_triple = '"""' in line or "'''" in line
    print(f'{i+1}: {repr(line[:70])}  triple={has_triple}')

# Also check if there's something wrong with class definition
print("\nClass definition check:")
for i in range(2656, 2665):
    print(f'{i+1}: {lines[i].rstrip()}')
