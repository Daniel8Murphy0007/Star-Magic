#!/usr/bin/env python3
"""Final targeted fix for remaining uncommented lines"""

import re

with open('MAIN_1_final.cpp', 'r', encoding='utf-8') as f:
    content = f.read()

# Find lines that start with prose keywords outside of comments
# These lines need to be commented out
lines = content.split('\n')
fixed_lines = []

in_block_comment = False
in_string = False

for line in lines:
    stripped = line.strip()
    
    # Track block comment state
    if '/*' in line:
        in_block_comment = True
    if '*/' in line:
        in_block_comment = False
        fixed_lines.append(line)
        continue
    
    # Skip if in comment
    if in_block_comment or stripped.startswith('//') or stripped.startswith('*'):
        fixed_lines.append(line)
        continue
    
    # Skip preprocessor and C++ statements
    if (stripped.startswith('#') or 
        stripped.startswith('using') or
        stripped.startswith('namespace') or
        stripped.startswith('int ') or
        stripped.startswith('double ') or
        stripped.startswith('const ') or
        stripped.startswith('std::') or
        stripped.startswith('return') or
        stripped.endswith(';') or
        stripped.endswith('{') or
        stripped.endswith('}')):
        fixed_lines.append(line)
        continue
    
    # Check if line starts with prose keywords
    starts_with_prose = (
        re.match(r'^(Key|The|This|Each|Full|Assumed|Equation|Learning:|Advancements:|Discoveries:)', stripped)
    )
    
    # Check if it's a standalone sentence (no code syntax)
    looks_like_prose = (
        len(stripped) > 20 and
        not '=' in stripped and
        not '(' in stripped[:10] and
        not '{' in stripped and
        not '}' in stripped and
        (starts_with_prose or
         ' the ' in stripped.lower() or
         ' is ' in stripped.lower() or
         ' are ' in stripped.lower())
    )
    
    if looks_like_prose and not stripped.endswith(','):
        fixed_lines.append('// ' + line)
    else:
        fixed_lines.append(line)

with open('MAIN_1_final.cpp', 'w', encoding='utf-8') as f:
    f.write('\n'.join(fixed_lines))

print("✓ Final prose fix applied")
