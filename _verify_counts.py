#!/usr/bin/env python3
"""Final verification of counts."""
import re

t = open('CondensedPhysics4.py', encoding='utf-8').read()
blocks = re.findall(r'_SESSION_\w+_CLASSES\s*=\s*\[([^\]]*)\]', t)
total_reg = 0
for b in blocks:
    entries = re.findall(r"'(\w+)'", b)
    total_reg += len(entries)
classes = re.findall(r'^class\s+(\w+)', t, re.MULTILINE)
print(f'Class definitions: {len(classes)}')
print(f'Total registered: {total_reg}')
print(f'Diff: {len(classes) - total_reg}')
