#!/usr/bin/env python3
"""Find unclosed docstrings in CondensedPhysics.py"""

with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    content = f.read()
    lines = content.split('\n')

# Count triple quote balance
in_docstring = False
docstring_start = None

for i, line in enumerate(lines, 1):
    count = line.count('"""')
    if count % 2 == 1:  # Odd number toggles state
        if in_docstring:
            in_docstring = False
            print(f'Docstring closed at line {i}')
        else:
            in_docstring = True
            docstring_start = i
            if i < 100 or (2600 < i < 2700):
                print(f'Docstring opened at line {i}')
    
    # Check for issues
    if i > 2650 and i < 2660 and in_docstring:
        print(f'  Line {i} (in docstring from {docstring_start}): {line[:50]}')

if in_docstring:
    print(f'\n*** UNCLOSED DOCSTRING starting at line {docstring_start} ***')
else:
    print(f'\nAll docstrings appear balanced')

# Also check for structure issues
print('\nChecking for structure issues before line 2657...')
# Find last complete class or function definition before line 2657
for i in range(2650, 0, -1):
    line = lines[i-1]
    if line.startswith('class ') or (line.startswith('def ') and not line.startswith('    def ')):
        print(f'Last top-level def/class before 2657: line {i}: {line[:60]}')
        break
