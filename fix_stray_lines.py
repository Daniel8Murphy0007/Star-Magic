#!/usr/bin/env python3
"""Fix remaining uncommented lines in MAIN_1_final.cpp"""

with open('MAIN_1_final.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

fixed_lines = []
for line in lines:
    # Comment lines that start with "This element", "This integration", etc.
    if line.strip().startswith('This element') or \
       line.strip().startswith('This integration') or \
       line.strip().startswith('This note') or \
       line.strip().startswith('This marker'):
        fixed_lines.append('// ' + line)
    else:
        fixed_lines.append(line)

with open('MAIN_1_final.cpp', 'w', encoding='utf-8') as f:
    f.writelines(fixed_lines)

print("✓ Commented all 'This element/integration/note/marker' lines")
