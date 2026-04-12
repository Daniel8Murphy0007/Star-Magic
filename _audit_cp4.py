#!/usr/bin/env python3
"""Audit CP4 class definitions vs registration blocks."""
import re

t = open('CondensedPhysics4.py', encoding='utf-8').read()

# Find all class definitions
classes = re.findall(r'^class\s+(\w+)', t, re.MULTILINE)
print(f'Total class definitions: {len(classes)}')

# Find all _SESSION registration blocks
sessions = re.findall(r'(_SESSION_\w+_CLASSES)\s*=\s*\[(.*?)\]', t, re.DOTALL)
total_reg = 0
reg_names = set()
for name, block in sessions:
    entries = re.findall(r'"(\w+)"', block)
    total_reg += len(entries)
    reg_names.update(entries)
    print(f'  {name}: {len(entries)} entries')

print(f'Total registered: {total_reg}')

# Classes defined but not registered
def_names = set(classes)
missing_reg = def_names - reg_names
if missing_reg:
    print(f'\nDEFINED but NOT in _SESSION blocks: {len(missing_reg)}')
    for m in sorted(missing_reg)[:20]:
        print(f'  {m}')

# Registered but not defined
missing_def = reg_names - def_names
if missing_def:
    print(f'\nREGISTERED but NOT DEFINED (phantom classes): {len(missing_def)}')
    for m in sorted(missing_def):
        print(f'  {m}')

# Check last 30 classes (should include #501-#506)
print(f'\nLast 10 class definitions:')
for c in classes[-10:]:
    print(f'  {c}')
