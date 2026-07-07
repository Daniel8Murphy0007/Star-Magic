"""
UQFF Primitive-Arithmetic Audit Script — Phase 2
Scans CondensedPhysics.py for numeric constants embedded in Round 1-47 stub upgrades.
Attempts to match each against primitive-arithmetic derivations from the 9 truly-independent
canonical primitives + simple combinations.

Triage:
  Category A — already primitive-derived (SSQ/K_MEX/F_TRZ/etc. in body) — pass
  Category B — structurally derivable but not connected — UPGRADE CANDIDATE
  Category C — per-system empirical (mass/T/distance) — no primitive form expected
"""

import re
import math
from itertools import product
from fractions import Fraction

PRIM = {
    'D_phys': 4,
    'D_crit': 26,
    'N_ch': 9,
    'SO_5': 10,
    'A_5': 60,
    'D_BSFG': 6,
    'F_TRZ': 0.1,
    'SSq': 0.57,
    'K_MEX': 25/12,
    'Phi_res': 0.84,
    'Phi_res_nuclear': 5/6,
    'beta_i': 0.6029,
    'S_26': 1.453162,
}

def gen_candidates(tol_rel=0.005):
    """Generate primitive-arithmetic candidates and their values."""
    cands = {}
    keys = list(PRIM.keys())
    for k, v in PRIM.items():
        cands[f'{k}'] = v
    for k1, v1 in PRIM.items():
        for k2, v2 in PRIM.items():
            if k1 == k2: continue
            if v2 != 0:
                cands[f'{k1}/{k2}'] = v1/v2
            cands[f'{k1}*{k2}'] = v1*v2
            cands[f'{k1}+{k2}'] = v1+v2
            cands[f'{k1}-{k2}'] = v1-v2
    for k1, v1 in PRIM.items():
        for exp in [2, 3, 4, 5, -1, -2, -3]:
            try:
                cands[f'{k1}^{exp}'] = v1**exp
            except:
                pass
    for k1, v1 in PRIM.items():
        for k2, v2 in PRIM.items():
            for k3, v3 in PRIM.items():
                if v2 != 0 and v3 != 0:
                    cands[f'{k1}/({k2}-{k3})' if v2>v3 else f'{k1}/({k3}-{k2})'] = v1/abs(v2-v3) if v2!=v3 else 0
                    cands[f'({k1}-{k2})/{k3}'] = (v1-v2)/v3
    fact_D_crit = math.factorial(26)
    cands['26!'] = fact_D_crit
    cands['26!*K_MEX'] = fact_D_crit * PRIM['K_MEX']
    cands['rho_SCm*26!*K_MEX'] = 7.09e-37 * fact_D_crit * PRIM['K_MEX']
    return cands

def match_value(target, cands, tol_rel=0.005):
    """Find primitive-arithmetic matches for a target value."""
    if target == 0: return []
    matches = []
    for name, val in cands.items():
        if val == 0: continue
        rel = abs(target - val) / abs(target) if target != 0 else 0
        if rel < tol_rel:
            matches.append((name, val, rel))
    matches.sort(key=lambda x: x[2])
    return matches[:5]

def scan_source(path):
    """Extract numeric constants from CondensedPhysics.py compute() bodies."""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        src = f.read()
    class_pat = re.compile(r'class\s+(\w+Calculator)\b[^:]*:(.*?)(?=\nclass\s+\w+|\Z)', re.DOTALL)
    number_pat = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+/\d+")
    results = []
    for m in class_pat.finditer(src):
        cls, body = m.group(1), m.group(2)
        compute_m = re.search(r'def\s+compute\([^)]*\)[^:]*:(.*?)(?=\n    def|\Z)', body, re.DOTALL)
        if not compute_m: continue
        compute_body = compute_m.group(1)
        already_derived = any(k in compute_body for k in ['SSQ', 'K_MEX', 'F_TRZ', 'F_UBi_i_99', 'PHI_RES', 'D_CRIT', 'SO_5', 'A_5', 'RHO_SCM', 'D_PHYS'])
        for numm in number_pat.finditer(compute_body):
            token = numm.group(0)
            try:
                val = float(Fraction(token)) if '/' in token else float(token)
            except:
                continue
            if abs(val) < 1e-40 or abs(val) > 1e40: continue
            if val in [0.0, 1.0, 2.0, 0.5, 10.0, 100.0, 1000.0, 60.0, 4.0, 6.0]: continue
            results.append({'class': cls, 'token': token, 'value': val, 'already_derived': already_derived})
    return results

if __name__ == '__main__':
    print('=' * 78)
    print('UQFF Primitive-Arithmetic Audit — Phase 2 CP1 P2 Rounds 1-47')
    print('=' * 78)
    cands = gen_candidates()
    print(f'\nGenerated {len(cands)} primitive-arithmetic candidates for matching')
    results = scan_source('CondensedPhysics.py')
    print(f'Scanned {len(results)} numeric constants across CP1 Calculator classes\n')

    cat_A = []
    cat_B = []
    cat_C = []
    seen = set()
    for r in results:
        key = (r['class'], r['token'])
        if key in seen: continue
        seen.add(key)
        matches = match_value(r['value'], cands, tol_rel=0.005)
        if matches and r['already_derived']:
            cat_A.append((r, matches[0]))
        elif matches:
            cat_B.append((r, matches[0]))
        else:
            cat_C.append(r)

    print(f'Category A (already primitive-derived): {len(cat_A)}')
    print(f'Category B (structurally derivable, NOT connected): {len(cat_B)}  <-- UPGRADE TARGETS')
    print(f'Category C (per-system empirical, no match): {len(cat_C)}\n')

    print('=' * 78)
    print(f'Category B — TOP 30 UPGRADE CANDIDATES (structurally derivable, not connected)')
    print('=' * 78)
    cat_B.sort(key=lambda x: x[1][2])
    for i, (r, (name, val, rel)) in enumerate(cat_B[:30], 1):
        print(f'  #{i:2d}  {r["class"]:60s}')
        print(f'       value={r["value"]:.6g} matches {name} = {val:.6g}  (rel={rel*100:.3f}%)')
