#!/usr/bin/env python3
"""Fix malformed math in Rossi file"""

import re

fpath = 'whitepapers/SCm_Rossi_ECat_Variants_Unified.md'

with open(fpath, 'r', encoding='utf-8') as f:
    content = f.read()

# Fix: S_{26^{(3) -> S_26^3
content = re.sub(r'S_\{26\^\{\(3\)', 'S_26^3', content)
content = re.sub(r'S_\{26\^', 'S_26^', content)

# Fix: E_{SCm-phonon = E_SCm-phonon
content = re.sub(r'E_\{SCm-', 'E_SCm-', content)
content = re.sub(r'E_\{phonon', 'E_phonon', content)

# Fix: times10^{ -> times10^
content = content.replace(r'\times10^{', r'\times10^')

# Remove unmatched closing braces after exponents
content = re.sub(r'(\^[\d]+)\}\s', r'\1 ', content)

# Fix dollar sign issues
lines = content.split('\n')
fixed_lines = []
for line in lines:
    # Count unmatched braces in math mode
    in_math = False
    brace_count = 0
    for i, char in enumerate(line):
        if char == '$':
            in_math = not in_math
        elif in_math:
            if char == '{':
                brace_count += 1
            elif char == '}':
                brace_count -= 1
    
    # If we end with unmatched braces, remove them
    if brace_count > 0:
        line = line.rstrip() + '}' * brace_count
    
    fixed_lines.append(line)

content = '\n'.join(fixed_lines)

with open(fpath, 'w', encoding='utf-8') as f:
    f.write(content)

print('[FIXED] Corrected math expressions')
