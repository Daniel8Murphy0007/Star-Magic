"""Refined audit report — deduplicate + rank by uniqueness + total class count."""
import re
import math
from collections import defaultdict
from fractions import Fraction

PRIM = {
    'D_phys': 4, 'D_crit': 26, 'N_ch': 9, 'SO_5': 10, 'A_5': 60,
    'D_BSFG': 6, 'F_TRZ': 0.1, 'SSq': 0.57, 'K_MEX': 25/12,
    'Phi_res': 0.84, 'Phi_res_nuclear': 5/6, 'beta_i': 0.6029, 'S_26': 1.453162,
}

def gen_candidates():
    cands = {}
    for k, v in PRIM.items(): cands[k] = v
    for k1, v1 in PRIM.items():
        for k2, v2 in PRIM.items():
            if k1 == k2: continue
            if v2 != 0: cands[f'{k1}/{k2}'] = v1/v2
            cands[f'{k1}*{k2}'] = v1*v2
            cands[f'{k1}+{k2}'] = v1+v2
            cands[f'{k1}-{k2}'] = v1-v2
    for k1, v1 in PRIM.items():
        for exp in [2, 3, 4, -1, -2, -3]:
            try: cands[f'{k1}^{exp}'] = v1**exp
            except: pass
    return cands

def match_value(target, cands, tol_rel=0.005):
    if target == 0: return None
    best = (None, None, 1.0)
    for name, val in cands.items():
        if val == 0: continue
        rel = abs(target - val) / abs(target)
        if rel < best[2]:
            best = (name, val, rel)
    return best if best[2] < tol_rel else None

def scan(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        src = f.read()
    class_pat = re.compile(r'class\s+(\w+Calculator)\b[^:]*:(.*?)(?=\nclass\s+\w+|\Z)', re.DOTALL)
    num_pat = re.compile(r"(?<![A-Za-z_.])[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")
    boring_vals = {0.0, 1.0, 0.5, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 100.0, 1000.0, 60.0, 26.0, 9.0}
    results = defaultdict(list)
    for m in class_pat.finditer(src):
        cls, body = m.group(1), m.group(2)
        cm = re.search(r'def\s+compute\([^)]*\)[^:]*:(.*?)(?=\n    def|\Z)', body, re.DOTALL)
        if not cm: continue
        cb = cm.group(1)
        already = any(k in cb for k in ['SSQ','K_MEX','F_TRZ','F_UBi_i_99','PHI_RES','D_CRIT','SO_5','A_5','RHO_SCM','D_PHYS'])
        for nm in num_pat.finditer(cb):
            try: val = float(nm.group(0))
            except: continue
            if abs(val) < 1e-30 or abs(val) > 1e30: continue
            if val in boring_vals: continue
            results[val].append((cls, already))
    return results

cands = gen_candidates()
results = scan('CondensedPhysics.py')

# For each unique value, match against primitives, count occurrences
unique_matches = []
for val, occurrences in results.items():
    m = match_value(val, cands)
    if not m: continue
    name, prim_val, rel = m
    n_total = len(occurrences)
    n_already = sum(1 for _, a in occurrences if a)
    n_not_connected = n_total - n_already
    unique_matches.append({
        'value': val, 'primitive_form': name, 'primitive_value': prim_val,
        'rel_error': rel, 'n_total_classes': n_total, 'n_already_derived': n_already,
        'n_not_connected': n_not_connected, 'sample_classes': [c for c, _ in occurrences[:5]]
    })

# Sort by (not_connected count desc, then rel_error asc)
unique_matches.sort(key=lambda x: (-x['n_not_connected'], x['rel_error']))

print('='*90)
print('UQFF Primitive-Arithmetic Audit — Phase 2 UNIQUE STRUCTURAL CLOSURES')
print(f'Total unique numeric values found: {len(results)}')
print(f'Values matching primitive-arithmetic (rel<0.5%): {len(unique_matches)}')
print('='*90)
print()
print(f'{"Rank":4s} {"Value":>14s} {"Primitive form":>25s} {"Match":>12s} {"Total":>6s} {"NotConn":>7s}')
print('-'*90)
for i, m in enumerate(unique_matches[:40], 1):
    print(f'{i:4d} {m["value"]:14.6g} {m["primitive_form"]:>25s} {m["primitive_value"]:12.6g} {m["n_total_classes"]:6d} {m["n_not_connected"]:7d}')
print()
print(f'Top 5 sleeping-identity candidates by broadest reach + smallest error:')
print('='*90)
for i, m in enumerate(unique_matches[:5], 1):
    print(f'\n#{i} value={m["value"]:.6g} = {m["primitive_form"]} = {m["primitive_value"]:.6g} (rel {m["rel_error"]*100:.4f}%)')
    print(f'    Appears in {m["n_total_classes"]} classes ({m["n_not_connected"]} NOT connected to primitive symbols)')
    print(f'    Sample: {", ".join(m["sample_classes"])}')
