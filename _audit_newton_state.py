"""Audit current Newton-as-seed violation state across all target files."""
import re

SEED_PAT = re.compile(
    r'(?:g_base|ug1|Ug1|Ug_base|Ug2|ug2|g_grav|base_gravity|g_proj)'
    r'\s*=\s*\(?.*?(?:self\.G|_G\b|G\b)\s*\*.*?/\s*\(?.*?\*\*\s*2'
)
RETURN_PAT = re.compile(
    r'return\s+\(?(?:self\.G|_G\b|G\b)\s*\*.*?/.*?\*\*\s*2'
)

FILES = [
    'CondensedPhysics.py',
    'CondensedPhysics2.py',
    'CondensedPhysics3.py',
    'CondensedPhysics4.py',
    'QCalc.py',
    'uqff_validation_test.py',
    'Phase6_Consolidated.py',
    'Phase7_Consolidated.py',
    'QCalc_Wolfram_Extensions.py',
    'add_uqff_to_8_models.py',
]

total = 0
for fname in FILES:
    try:
        with open(fname, encoding='utf-8-sig', errors='ignore') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f'{fname}: NOT FOUND')
        continue

    violations = []
    for i, line in enumerate(lines, 1):
        stripped = line.strip()
        # skip comments and string-only lines
        if stripped.startswith('#'):
            continue
        if SEED_PAT.search(stripped) or RETURN_PAT.search(stripped):
            violations.append((i, stripped[:100]))

    status = 'CLEAN' if not violations else f'{len(violations)} VIOLATIONS'
    print(f'{fname}: {status}')
    for lineno, text in violations[:8]:
        print(f'  L{lineno}: {text}')
    if len(violations) > 8:
        print(f'  ... +{len(violations)-8} more')
    total += len(violations)

print(f'\nTOTAL VIOLATIONS REMAINING: {total}')
