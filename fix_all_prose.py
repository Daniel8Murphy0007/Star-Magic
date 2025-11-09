#!/usr/bin/env python3
"""Comprehensive fix for all remaining uncommented prose in MAIN_1_final.cpp"""

with open('MAIN_1_final.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

fixed_lines = []
in_block_comment = False

for i, line in enumerate(lines):
    stripped = line.strip()
    
    # Check if we're in a block comment
    if '/*' in line:
        in_block_comment = True
    if '*/' in line:
        in_block_comment = False
        fixed_lines.append(line)
        continue
    
    # Skip if already a comment
    if stripped.startswith('//') or stripped.startswith('/*') or stripped.startswith('*'):
        fixed_lines.append(line)
        continue
    
    # Skip preprocessor directives and using statements
    if stripped.startswith('#') or stripped.startswith('using'):
        fixed_lines.append(line)
        continue
    
    # Check for uncommented prose (long sentences without semicolons)
    is_prose = (
        len(stripped) > 80 and
        not in_block_comment and
        not stripped.endswith(';') and
        not stripped.endswith('{') and
        not stripped.endswith('}') and
        not stripped.endswith(',') and
        not '(' in stripped[:20] and  # Not a function call
        not 'int ' in stripped and
        not 'double ' in stripped and
        not 'const ' in stripped and
        not 'std::' in stripped and
        (
            ' the ' in stripped.lower() or
            ' is ' in stripped.lower() or
            ' are ' in stripped.lower() or
            ' this ' in stripped.lower() or
            'each of' in stripped.lower() or
            'full assumed' in stripped.lower() or
            'buoyancy' in stripped.lower() or
            'experimental' in stripped.lower() or
            'challenges' in stripped.lower() or
            'explanation:' in stripped.lower()
        )
    )
    
    if is_prose:
        fixed_lines.append('// ' + line)
    else:
        fixed_lines.append(line)

with open('MAIN_1_final.cpp', 'w', encoding='utf-8') as f:
    f.writelines(fixed_lines)

print("✓ Commented all remaining prose lines")
