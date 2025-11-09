#!/usr/bin/env python3
"""Remove excessive indentation from struct SystemParams"""

with open('MAIN_1_final.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

fixed_lines = []
in_struct = False

for line in lines:
    # Check if this is the struct SystemParams line or its comment
    if 'struct SystemParams' in line or ('// Struct for system params' in line and '                                 ' in line):
        # Remove excessive indentation
        fixed_lines.append(line.lstrip())
        if 'struct SystemParams' in line:
            in_struct = True
    else:
        fixed_lines.append(line)

with open('MAIN_1_final.cpp', 'w', encoding='utf-8') as f:
    f.writelines(fixed_lines)

print("✓ Fixed struct indentation")
