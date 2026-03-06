#!/usr/bin/env python3
"""
Analyze all duplicate compute_master_buoyant_equation definitions in CondensedPhysics.py.
Groups by unique body hash, then compares key equation components.
"""
import re
import hashlib
from collections import defaultdict

with open('CondensedPhysics.py', encoding='utf-8', errors='replace') as f:
    lines = f.read().splitlines()

results = []
for i, line in enumerate(lines):
    if re.match(r'\s+def compute_master_buoyant_equation', line):
        cls = 'UNKNOWN'
        for j in range(i - 1, -1, -1):
            m = re.match(r'^class (\w+)', lines[j])
            if m:
                cls = m.group(1)
                break
        indent = len(line) - len(line.lstrip())
        body = [line]
        for k in range(i + 1, min(i + 60, len(lines))):
            l2 = lines[k]
            if l2.strip() == '':
                body.append(l2)
                continue
            li2 = len(l2) - len(l2.lstrip())
            if li2 <= indent and re.match(r'\s*(def |class )', l2):
                break
            body.append(l2)
        norm = '\n'.join(l.rstrip() for l in body)
        h = hashlib.md5(norm.encode()).hexdigest()[:8]
        results.append({'line': i + 1, 'cls': cls, 'hash': h, 'body': norm})

# ── Group by hash ──────────────────────────────────────────────────────────────
by_hash = defaultdict(list)
for r in results:
    by_hash[r['hash']].append(r)

print(f'Total definitions   : {len(results)}')
print(f'Unique body variants: {len(by_hash)}')
print()

# ── Key equation fingerprints from each variant ─────────────────────────────
EQUATION_PATTERNS = [
    # formula form
    ('buoyancy_formula', re.compile(r'F_U_Bi_i\s*=.*')),
    ('F_Bi_formula',     re.compile(r'F_Bi\s*=.*')),
    ('galactic_factor',  re.compile(r'galactic_factor\s*=.*')),
    ('Omega_g_source',   re.compile(r'Omega_g\s*=.*|Omega_g\s*=.*CONSTANTS.*')),
    ('M_bh_source',      re.compile(r'M_bh\s*=.*')),
    ('d_g_source',       re.compile(r'd_g\s*=.*')),
    ('Ug_sum',           re.compile(r'Ug_sum\s*=.*|Ug\s*=.*')),
    ('f_TRZ',            re.compile(r'f_TRZ\s*=.*')),
    ('return_keys',      re.compile(r"return\s*\{")),
    ('calls_UQFF_base',  re.compile(r'compute_UQFF_base|compute_F_U_Bi_simple|compute_buoyant')),
]

sep = '=' * 72

for idx, (h, entries) in enumerate(sorted(by_hash.items(), key=lambda x: -len(x[1])), 1):
    first = entries[0]
    body_lines = first['body'].splitlines()
    class_list = [e['cls'] for e in entries]

    print(sep)
    print(f'VARIANT {idx:02d}  hash={h}  count={len(entries)}')
    print(f'  First line : {first["line"]}')
    print(f'  Classes    : {class_list[:6]}{"..." if len(class_list) > 6 else ""}')
    print()

    # Extract key equation lines
    print('  >>>> Key equation fingerprints <<<<')
    body_text = first['body']
    for label, pat in EQUATION_PATTERNS:
        m = pat.search(body_text)
        if m:
            print(f'    [{label}]  {m.group(0).strip()}')
    print()

    # Print full body
    print('  >>>> Full body <<<<')
    for bl in body_lines:
        print('    ' + bl)
    print()

# ── Intra-class duplicates ────────────────────────────────────────────────────
print(sep)
print('TRUE INTRA-CLASS DUPLICATES (same class, multiple definitions):')
by_class = defaultdict(list)
for r in results:
    by_class[r['cls']].append(r)
any_found = False
for cls, entries in sorted(by_class.items()):
    if len(entries) > 1:
        any_found = True
        hashes = [e['hash'] for e in entries]
        unique_h = set(hashes)
        same = 'ALL IDENTICAL' if len(unique_h) == 1 else f'DIFFERENT bodies ({len(unique_h)} variants)'
        print(f'  {cls}: {len(entries)} definitions — {same}')
        for e in entries:
            print(f'      line {e["line"]}  hash={e["hash"]}')
if not any_found:
    print('  None (all inter-class)')
print()

# ── Summary table ─────────────────────────────────────────────────────────────
print(sep)
print('SUMMARY TABLE: variant vs key formula')
print(f'  {"VAR":>3}  {"hash":>8}  {"count":>5}  {"F_U_Bi_i formula (truncated)"}')
print(f'  {"-"*3}  {"-"*8}  {"-"*5}  {"-"*45}')
for idx, (h, entries) in enumerate(sorted(by_hash.items(), key=lambda x: -len(x[1])), 1):
    body_text = entries[0]['body']
    fubi = re.search(r'F_U_Bi_i\s*=([^\n#]{0,60})', body_text)
    formula = fubi.group(1).strip() if fubi else '(no F_U_Bi_i assignment found)'
    print(f'  {idx:>3}  {h}  {len(entries):>5}  {formula}')
