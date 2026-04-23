import re

files_to_check = [
    ('CondensedPhysics.py',  r'(?:self\.G|_G\b)\s*\*\s*(?:self\.M|[\w\.]+)\s*/\s*\(?[\w\.]+\s*\*\*\s*2'),
    ('CondensedPhysics2.py', r'(?:self\.G|_G\b)\s*\*\s*(?:self\.M|[\w\.]+)\s*/\s*\(?[\w\.]+\s*\*\*\s*2'),
    ('CondensedPhysics3.py', r'(?:self\.G|_G\b)\s*\*\s*(?:self\.M|[\w\.]+)\s*/\s*\(?[\w\.]+\s*\*\*\s*2'),
    ('QCalc.py',             r'_G\b\s*\*\s*[\w\.]+\s*/\s*\(?[\w\.]+\s*\*\*\s*2'),
    ('Phase7_Consolidated.py', r'\bG\b\s*\*\s*\w+\s*/\s*\(?\w+\s*\*\*\s*2'),
    ('QCalc_Wolfram_Extensions.py', r'\bG\b\s*\*\s*\w+\s*/\s*\(?\w+\s*\*\*\s*2'),
    ('add_uqff_to_8_models.py', r'\bG\b\s*\*\s*\w+\s*/\s*\(?\w+\s*\*\*\s*2'),
]

for fname, pat in files_to_check:
    try:
        with open(fname, encoding='utf-8-sig', errors='ignore') as f:
            lines = f.readlines()
        hits = []
        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            if '=' not in line:
                continue
            if not re.search(pat, line):
                continue
            if 'g_projection' in line or 'g_proj' in line:
                continue
            hits.append((i, line.rstrip()[:95]))
        print(f"{fname}: {len(hits)} potential violation(s)")
        for ln, text in hits[:15]:
            print(f"  L{ln}: {text}")
    except FileNotFoundError:
        print(f"{fname}: NOT FOUND")
