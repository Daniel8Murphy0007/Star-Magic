#!/usr/bin/env python3
"""Final fix for SCm_Rossi_ECat_Variants_Unified.md"""

import re

fpath = 'whitepapers/SCm_Rossi_ECat_Variants_Unified.md'

with open(fpath, 'r', encoding='utf-8') as f:
    content = f.read()

# Remove any remaining equation environments
content = re.sub(r'\\begin\{equation\}.*?\\end\{equation\}', '', content, flags=re.DOTALL)
content = re.sub(r'\\begin\{align\}.*?\\end\{align\}', '', content, flags=re.DOTALL)

# Remove problematic LaTeX commands
content = re.sub(r'\\@[a-zA-Z]+', '', content)
content = re.sub(r'\\[a-zA-Z]*box', '', content)
content = re.sub(r'\\textasciitilde', '~', content)
content = re.sub(r'\\checkmark', '', content)

# Fix broken \phi_T equation
content = re.sub(r'\\Phi_T\s*=\s*[^\\]*,', '(Gaussian resonance adjustment),', content)

# Remove any stray backslashes followed by special characters
content = re.sub(r'\\[^a-zA-Z\s{]', '', content)

with open(fpath, 'w', encoding='utf-8') as f:
    f.write(content)

print('[CLEANED] Removed all problematic LaTeX')
