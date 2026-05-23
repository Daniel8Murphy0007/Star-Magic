#!/usr/bin/env python3
"""Simple fix for Rossi file - just rebuild it cleanly"""

fpath = 'whitepapers/SCm_Rossi_ECat_Variants_Unified.md'

with open(fpath, 'r', encoding='utf-8') as f:
    lines = f.readlines()

cleaned_lines = []
for line in lines:
    # Remove LaTeX-specific patterns
    line = line.replace(r'\varepsilon', 'epsilon')
    line = line.replace(r'\epsilon', 'epsilon')
    line = line.replace(r'\kappa', 'kappa')
    line = line.replace(r'\text{', '')
    line = line.replace('}', '')
    line = line.replace(r'degree C', '°C')
    line = line.replace(r'\\', '')  # Remove double backslashes
    
    cleaned_lines.append(line)

with open(fpath, 'w', encoding='utf-8') as f:
    f.writelines(cleaned_lines)

print('[FIXED] Cleaned Rossi file')
