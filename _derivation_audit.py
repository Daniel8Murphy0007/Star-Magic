"""
Honest audit: separate TRUE derivations from anchor lookups in uqff_pure_calculator.py.

A TRUE derivation = function whose return value is computed from PRIMITIVE constants
(BETA_I, RHO_UA, RHO_SCM, S_26, PHI_RESONANCE, S26_DPM, F_THZ, K_B, G_NEWTON, EV_J,
 D_BSFG, D_CRIT, TRZ, etc.) via arithmetic operations only.

An ANCHOR lookup = function whose return value comes from a literal table
(_SM_LITERAL_ANCHOR_SAT, CODATA constants, or hard-coded experimental targets).
"""
import re
from collections import defaultdict

with open('uqff_pure_calculator.py', encoding='utf-8') as f:
    src = f.read()

lines = src.splitlines()

# 1. Count all defined functions
all_defs = re.findall(r'^def\s+([A-Za-z_]\w*)\s*\(', src, re.MULTILINE)
print(f'Total `def` functions:          {len(all_defs)}')
print(f'  unique names:                 {len(set(all_defs))}')

# 2. Functions that return a literal anchor
anchor_fns = re.findall(r'^def\s+([A-Za-z_]\w*)\s*\([\s\S]*?return\s+_sm_literal_anchor\(', src, re.MULTILINE)
print(f'  -> route to _sm_literal_anchor (ANCHOR):  {len(set(anchor_fns))}')

# 3. Functions that route to _master_constant_primitive (predictive chain)
prim_fns = re.findall(r'^def\s+([A-Za-z_]\w*)\s*\([\s\S]*?return\s+_master_constant_primitive\(', src, re.MULTILINE)
print(f'  -> route to _master_constant_primitive:    {len(set(prim_fns))}')

# 4. _canonical_sat family (Layer 45)
canon = re.findall(r'^def\s+(_[a-zA-Z_]+_canonical_sat)\s*\(', src, re.MULTILINE)
print(f'  _*_canonical_sat (Layer 45):  {len(set(canon))}')

# 5. _primitive_sat family
prim_sat = re.findall(r'^def\s+(_[a-zA-Z_]+_primitive_sat)\s*\(', src, re.MULTILINE)
print(f'  _*_primitive_sat (predictive):{len(set(prim_sat))}')

# 6. Anchor table entries
m = re.search(r'_SM_LITERAL_ANCHOR_SAT\s*[:\w\[\],\s]*=\s*\{([\s\S]*?)\n\}', src)
if m:
    keys = re.findall(r"['\"]([A-Za-z0-9_]+)['\"]\s*:", m.group(1))
    print(f'  _SM_LITERAL_ANCHOR_SAT entries:           {len(keys)}')

# 7. MILLENNIUM_TARGETS entries
m = re.search(r'MILLENNIUM_TARGETS\s*:[^=]*=\s*\{([\s\S]*?)\n\}', src)
if m:
    keys = re.findall(r"['\"]([A-Za-z0-9_]+)['\"]\s*:", m.group(1))
    print(f'  MILLENNIUM_TARGETS entries:                {len(keys)}')

# 8. _MILLENNIUM_DERIVE entries
m = re.search(r'_MILLENNIUM_DERIVE\s*:[^=]*=\s*\{([\s\S]*?)\n\}', src)
if m:
    keys = re.findall(r"['\"]([A-Za-z0-9_]+)['\"]\s*:", m.group(1))
    print(f'  _MILLENNIUM_DERIVE entries:                {len(keys)}')

# 9. _LEDGER_PRIMITIVE entries
m = re.search(r'_LEDGER_PRIMITIVE\s*:[^=]*=\s*\{([\s\S]*?)\n\}', src)
if m:
    keys = re.findall(r"['\"]([A-Za-z0-9_]+)['\"]\s*:", m.group(1))
    print(f'  _LEDGER_PRIMITIVE entries (predictive):    {len(keys)}')

# 10. catalog references — "630 constants" / known catalog counts
print('\nCatalog phrases in source:')
for pat in [r'630\s+(?:scientific|physical|constants)', r'cataloged?\s+(\d+)', r'derived?\s+(\d+)', r'CATALOG_(\d+)']:
    for m in re.finditer(pat, src, re.IGNORECASE):
        line_no = src[:m.start()].count('\n') + 1
        print(f'  L{line_no}: {lines[line_no-1].strip()[:120]}')

# 11. Count `return <numeric literal>` in functions (hard-coded answers)
hardcoded = 0
in_def = False
cur_name = None
for ln in lines:
    md = re.match(r'^def\s+([A-Za-z_]\w*)', ln)
    if md:
        cur_name = md.group(1); in_def = True; continue
    if in_def and re.match(r'^def\s', ln):
        cur_name = None
    if cur_name and re.match(r'^\s+return\s+[-\d]', ln):
        hardcoded += 1
print(f'\nFunctions returning a bare numeric literal: {hardcoded}')

# 12. Dispatcher `_derive_constant` branch count
m = re.search(r'def _derive_constant\([\s\S]*?\n(?=def\s|\Z)', src)
if m:
    body = m.group(0)
    branches = len(re.findall(r'\bif\s+n\s*==|\belif\s+n\s*==|\bif\s+n\s+in\s*\(', body))
    print(f'_derive_constant branches: {branches}')
    print(f'_derive_constant length: {body.count(chr(10))} lines')
