#!/usr/bin/env python3
"""Fix all remaining uncommented prose in MAIN_1_fixed.cpp"""

with open('MAIN_1_fixed.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

fixed_lines = []
in_block_comment = False

for i, line in enumerate(lines):
    stripped = line.strip()
    
    # Track block comment state
    if '/*' in line and '*/' not in line:
        in_block_comment = True
        fixed_lines.append(line)
        continue
    if '*/' in line and '/*' not in line:
        in_block_comment = False
        fixed_lines.append(line)
        continue
    
    # Skip if already commented or in block comment
    if in_block_comment or stripped.startswith('//') or stripped.startswith('*') or stripped.startswith('#') or stripped.startswith('using'):
        fixed_lines.append(line)
        continue
    
    # Check if line starts with prose that needs commenting
    needs_comment = (
        stripped.startswith('This ') or
        stripped.startswith('- ') or
        stripped.startswith('For i =') or
        (len(stripped) > 50 and 
         not stripped.endswith(';') and 
         not stripped.endswith('{') and
         not stripped.endswith('}') and
         not stripped.endswith(',') and
         not '=' in stripped[:20] and
         not '(' in stripped[:15] and
         ('comment in the C++' in stripped or
          'This equation' in stripped or
          'This element' in stripped or
          'This marker' in stripped or
          'This note' in stripped or
          'This integration' in stripped))
    )
    
    if needs_comment:
        fixed_lines.append('// ' + line)
    else:
        fixed_lines.append(line)

with open('MAIN_1_fixed.cpp', 'w', encoding='utf-8') as f:
    f.writelines(fixed_lines)

print("✓ Fixed all uncommented prose in MAIN_1_fixed.cpp")
