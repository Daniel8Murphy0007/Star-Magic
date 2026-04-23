#!/usr/bin/env python3
"""Deep audit of Wolfram Extensions."""
import re
with open('QCalc_Wolfram_Extensions.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

print('=== NEWTONIAN BASE ISSUES ===')
for i,l in enumerate(lines):
    if re.search(r'G \* M.*r\s*\*\*\s*2|G\*M.*r\*\*2|\(G \* M\).*r\s*\*\*\s*2', l):
        func_name = None
        for j in range(i, -1, -1):
            if lines[j].startswith('def '):
                func_name = lines[j].strip()[:80]
                break
        print(f'  Line {i+1} in {func_name}:')
        print(f'    {l.strip()}')

print('\n=== FUNCTIONS WITH Ug BUT NO E_react ===')
count = 0
for i,l in enumerate(lines):
    if l.startswith('def ') and 'calculate_' in l:
        block = lines[i:min(i+80, len(lines))]
        has_ereact = any('E_react' in bl or 'e_react' in bl.lower() for bl in block)
        has_ug = any(re.search(r'\bug1\b|\bUg1\b|\bug2\b|\bUg2\b', bl, re.I) for bl in block)
        if has_ug and not has_ereact:
            count += 1
            print(f'  Line {i+1}: {l.strip()[:80]}')

print(f'Total functions with Ug but no E_react: {count}')

print('\n=== FUNCTIONS TOTAL ===')
funcs = [l.strip() for l in lines if l.startswith('def ')]
print(f'Total functions: {len(funcs)}')
