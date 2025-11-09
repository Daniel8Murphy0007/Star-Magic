#!/usr/bin/env python3
"""Comment all bullet point lists outside of block comments"""

with open('MAIN_1_final.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

fixed_lines = []
in_block_comment = False

for i, line in enumerate(lines):
    stripped = line.strip()
    
    # Track block comment state
    if '/*' in line:
        in_block_comment = True
    if '*/' in line:
        in_block_comment = False
        fixed_lines.append(line)
        continue
    
    # If already commented or in block comment, keep as is
    if in_block_comment or stripped.startswith('//') or stripped.startswith('*'):
        fixed_lines.append(line)
        continue
    
    # Skip preprocessor directives
    if stripped.startswith('#'):
        fixed_lines.append(line)
        continue
    
    # Check if this is a bullet point or markdown-style list item
    is_bullet = (
        stripped.startswith('- **') or  # Markdown bullet with bold
        stripped.startswith('- ') or    # Plain markdown bullet
        (stripped.startswith('**') and '**:' in stripped)  # Bold heading with colon
    )
    
    # Also catch prose lines that look like definitions
    starts_with_definition = (
        stripped.startswith('Full ') or
        stripped.startswith('Each ') or
        stripped.startswith('Where ')
    )
    
    if is_bullet or starts_with_definition:
        fixed_lines.append('// ' + line)
    else:
        fixed_lines.append(line)

with open('MAIN_1_final.cpp', 'w', encoding='utf-8') as f:
    f.writelines(fixed_lines)

print("✓ Commented all bullet points and definition lines")
