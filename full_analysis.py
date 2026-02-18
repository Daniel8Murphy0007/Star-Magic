#!/usr/bin/env python
"""Full unabridged analysis of CondensedPhysics.py"""
import re

with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    content = f.read()
    lines = content.split('\n')

print('='*80)
print('FULL CONDENSEDPHYSICS.PY ANALYSIS')
print('='*80)

# Basic stats
print('\n1. FILE STATISTICS')
print('-'*80)
print(f'Total Lines: {len(lines)}')
print(f'Total Characters: {len(content)}')
print(f'File Size: {len(content.encode())} bytes')

# Find all classes
class_pattern = r'^class\s+(\w+)'
classes = []
for i, line in enumerate(lines, 1):
    m = re.match(class_pattern, line)
    if m:
        classes.append((i, m.group(1)))

print(f'\n2. CLASS INVENTORY ({len(classes)} total)')
print('-'*80)

# Categorize
models = [(ln, nm) for ln, nm in classes if nm.endswith('Model')]
frameworks = [(ln, nm) for ln, nm in classes if 'Framework' in nm or 'Calculator' in nm]
fetchers = [(ln, nm) for ln, nm in classes if 'Fetch' in nm]
equations = [(ln, nm) for ln, nm in classes if 'Equation' in nm and not nm.endswith('Model')]
other = [(ln, nm) for ln, nm in classes if not any([nm.endswith('Model'), 'Framework' in nm, 'Calculator' in nm, 'Fetch' in nm, ('Equation' in nm and not nm.endswith('Model'))])]

print(f'  Models: {len(models)}')
print(f'  Frameworks: {len(frameworks)}')
print(f'  Fetchers: {len(fetchers)}')
print(f'  Equations: {len(equations)}')
print(f'  Other/Utility: {len(other)}')

print(f'\n3. ALL {len(models)} MODELS (Line | Name | Size)')
print('-'*80)
for i, (ln, nm) in enumerate(models):
    next_ln = models[i+1][0] if i+1 < len(models) else len(lines)
    size = next_ln - ln
    print(f'  {ln:6} | {nm:45} | {size:5} lines')

print(f'\n4. FRAMEWORKS & CALCULATORS ({len(frameworks)})')
print('-'*80)
for ln, nm in frameworks:
    print(f'  {ln:6} | {nm}')

print(f'\n5. FETCHERS ({len(fetchers)})')
print('-'*80)
for ln, nm in fetchers:
    print(f'  {ln:6} | {nm}')

print(f'\n6. OTHER CLASSES ({len(other)})')
print('-'*80)
for ln, nm in other:
    print(f'  {ln:6} | {nm}')

# Duplicates check
print(f'\n7. DUPLICATE CHECK')
print('-'*80)
class_names = [nm for _, nm in classes]
seen = set()
dupes = []
for nm in class_names:
    if nm in seen:
        dupes.append(nm)
    seen.add(nm)
if dupes:
    print(f'  DUPLICATES FOUND: {dupes}')
else:
    print('  NO DUPLICATE CLASSES FOUND')

# Model method check
print(f'\n8. MODEL UPGRADE STATUS')
print('-'*80)
upgraded = []
needs_upgrade = []
may_2025 = []

may_2025_names = ['NGC2264Model','UGC10214Model','NGC4676Model','RedSpiderNebulaModel',
                  'NGC3372Model','AGCarinaeModel','M42Model','TarantulaNebulaModel',
                  'NGC2841Model','MysticMountainModel']

for i, (ln, nm) in enumerate(models):
    next_ln = models[i+1][0] if i+1 < len(models) else len(lines)
    model_code = '\n'.join(lines[ln-1:next_ln-1])
    has_clean = 'compute_clean_equation' in model_code
    has_res = 'compute_resonance' in model_code
    has_val = 'validate_model' in model_code
    has_test = 'run_tests' in model_code
    size = next_ln - ln
    
    if has_clean:
        upgraded.append((ln, nm, size))
    else:
        needs_upgrade.append((ln, nm, size, has_res, has_val, has_test))
    
    if nm in may_2025_names:
        may_2025.append((ln, nm, size, has_clean, has_res, has_val, has_test))

print(f'UPGRADED (have compute_clean_equation): {len(upgraded)}')
for ln, nm, sz in upgraded:
    print(f'  {ln:6} | {nm:42} | {sz:5} lines')

print(f'\nNEEDS UPGRADE: {len(needs_upgrade)}')
for ln, nm, sz, r, v, t in needs_upgrade:
    ry = 'Y' if r else 'N'
    vy = 'Y' if v else 'N'
    ty = 'Y' if t else 'N'
    print(f'  {ln:6} | {nm:42} | {sz:5} lines | Res:{ry} Val:{vy} Test:{ty}')

print(f'\n9. MAY 2025 MODELS STATUS')
print('-'*80)
for ln, nm, sz, has_clean, has_res, has_val, has_test in may_2025:
    status = 'UPGRADED' if has_clean else 'NEEDS compute_clean_equation'
    ry = 'Y' if has_res else 'N'
    vy = 'Y' if has_val else 'N'
    ty = 'Y' if has_test else 'N'
    print(f'  {ln:6} | {nm:30} | {sz:5} lines | {status} | Res:{ry} Val:{vy} Test:{ty}')

print(f'\n10. CONSTANTS DICTIONARY')
print('-'*80)
const_start = content.find('CONSTANTS = {')
if const_start >= 0:
    start_line = content[:const_start].count('\n') + 1
    # Count entries by finding key patterns
    const_section = content[const_start:const_start+500000]  # first 500k chars from CONSTANTS
    entries = len(re.findall(r"'[A-Z][A-Za-z0-9_]+':\s*{", const_section))
    print(f'  Start Line: {start_line}')
    print(f'  Approx Entries: {entries}')
else:
    print('  CONSTANTS dict not found')

# Top-level functions
print(f'\n11. TOP-LEVEL FUNCTIONS')
print('-'*80)
func_pattern = r'^def\s+(\w+)'
functions = []
for i, line in enumerate(lines, 1):
    m = re.match(func_pattern, line)
    if m:
        functions.append((i, m.group(1)))

print(f'Total top-level functions: {len(functions)}')
for ln, nm in functions:
    print(f'  {ln:6} | {nm}')

# Function duplicates
func_names = [nm for _, nm in functions]
seen = set()
func_dupes = []
for nm in func_names:
    if nm in seen:
        func_dupes.append(nm)
    seen.add(nm)
if func_dupes:
    print(f'  DUPLICATE FUNCTIONS: {func_dupes}')
else:
    print('  NO DUPLICATE FUNCTIONS')

print('\n' + '='*80)
print('FINAL SUMMARY')
print('='*80)
print(f'Total Lines:            {len(lines)}')
print(f'Total Characters:       {len(content)}')
print(f'Total Classes:          {len(classes)}')
print(f'Total Models:           {len(models)}')
print(f'Upgraded Models:        {len(upgraded)}')
print(f'Models Needing Upgrade: {len(needs_upgrade)}')
print(f'May 2025 Models:        {len(may_2025)}')
print(f'Top-Level Functions:    {len(functions)}')
print(f'Duplicate Classes:      {len(dupes)}')
print(f'Duplicate Functions:    {len(func_dupes)}')
print('='*80)
