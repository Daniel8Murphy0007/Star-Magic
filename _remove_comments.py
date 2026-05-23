#!/usr/bin/env python3
"""Remove ALL comment lines and Unicode"""

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'r', encoding='utf-8') as f:
    lines = f.readlines()

clean_lines = []
for line in lines:
    # Skip lines that are only comments or have problematic Unicode
    if line.strip().startswith('%'):
        continue
    if '─' in line or '│' in line or '├' in line or '└' in line:
        continue
    if 'ΓöÇ' in line:
        continue
    # Remove trailing % comments
    if '%' in line:
        line = line.split('%')[0] + '\n'
    clean_lines.append(line)

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.writelines(clean_lines)

print('[CLEANED] Comments and Unicode removed')
