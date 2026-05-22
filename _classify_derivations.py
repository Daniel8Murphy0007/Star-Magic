"""Classify every legit closure script into one of four categories."""
import csv
import re
from pathlib import Path

rows = [r for r in csv.DictReader(open('master_closures.csv', encoding='utf-8'))
        if r['status'] in ('OK', 'OK_JSON') and r['ID'].isdigit() and int(r['ID']) <= 693]

UQFF_PRIMS = ['F_TRZ', 'Phi_res', 'SSq', 'K_Mex', 'D_phys', 'D_BSFG', 'D_crit',
              'N_ch', 'SO5', 'A_5', 'beta_i', 'phi_res', 'k_mex', 'd_phys', 'd_bsfg']
PHYS_CONSTS = ['hbar', 'c_light', 'G_N', 'k_B', 'alpha_em', 'm_e', 'm_p', 'm_n',
               'm_planck', 'rho_vac', 'rho_scm', 'RHO_SCM', 'V_SCM', 'kappa',
               'alpha_s', 'm_W', 'm_Z', 'm_H', 'M_sun', 'L_sun', 'G_F']
BOOKKEEP_SIGNALS = ['NULL extraction', 'null extraction', 'definitional null',
                    'cardinality of', 'session checkpoint', 'state-counter',
                    'state-counters', 'identity check', 'definitional ratio']

counts = {'PHYSICS': 0, 'PROOF': 0, 'DATA': 0, 'BOOKKEEP': 0}
samples = {k: [] for k in counts}
unclassified = []

for r in rows:
    sid = r['ID']
    script = (r.get('script') or '').strip()
    p = Path(script) if script else None
    if not (p and p.exists()):
        cands = list(Path('.').glob(f'_session{sid}*.py'))
        if not cands:
            unclassified.append(sid)
            continue
        name = r.get('name', '')
        match = [c for c in cands if name and name.lower() in c.name.lower()]
        p = match[0] if match else cands[0]
    try:
        t = p.read_text(encoding='utf-8', errors='ignore')
    except Exception:
        unclassified.append(sid)
        continue

    prim_hits = sum(1 for k in UQFF_PRIMS if re.search(r'\b' + re.escape(k) + r'\b', t))
    # Many scripts use single-letter aliases: F=Fraction(1,10), P=Fraction(5,6),
    # K=Fraction(25,12), SO=Fraction(10), S=Fraction(57,100), etc.
    # Detect the rational signatures of the locked primitives.
    PRIM_RATIONALS = ['Fraction(1,10)', 'Fraction(5,6)', 'Fraction(25,12)',
                      'Fraction(57,100)', 'Fraction(6029,10000)', 'Fraction(10)',
                      'Fraction(60)', 'Fraction(26)', 'Fraction(4,1)', 'Fraction(6,1)']
    prim_hits += sum(1 for k in PRIM_RATIONALS if k in t)
    # Aliased-primitive bundle pattern (one-liner UQFF derivations use this)
    ALIAS_PATTERNS = [
        r'FT\s*,\s*Ph\s*,\s*SQ\s*,\s*K',         # FT,Ph,SQ,K,...=...
        r'F\s*,\s*P\s*,\s*K\s*,\s*SO',           # F,P,K,SO=...
        r'F\s*=\s*Fraction\(1\s*,\s*10\)',       # F=Fraction(1,10)
        r'Fraction\s+as\s+F\b',                   # from fractions import Fraction as F
    ]
    if any(re.search(p, t) for p in ALIAS_PATTERNS):
        prim_hits += 5  # strong UQFF signature
    const_hits = sum(1 for k in PHYS_CONSTS if re.search(r'\b' + re.escape(k) + r'\b', t))
    bookkeep_hit = any(s in t for s in BOOKKEEP_SIGNALS)
    arith_lines = sum(1 for L in t.splitlines()
                      if re.search(r'[=+\-*/].*[+\-*/]', L) and not L.strip().startswith('#'))

    # Reclassify: any primitive reference + any math op = PHYSICS
    has_math = arith_lines >= 1 or 'math.' in t or 'np.' in t or 'sympy' in t or '**' in t
    if bookkeep_hit and prim_hits == 0 and arith_lines < 10:
        cat = 'BOOKKEEP'
    elif prim_hits >= 1 and has_math:
        cat = 'PHYSICS'
    elif const_hits >= 1 and has_math:
        cat = 'PROOF'
    elif arith_lines >= 10:
        cat = 'DATA'
    else:
        cat = 'BOOKKEEP'
    counts[cat] += 1
    if len(samples[cat]) < 10:
        samples[cat].append(f'  S{sid:>3s}  prims={prim_hits} consts={const_hits} arith={arith_lines:>3d}  {p.name}')

total = sum(counts.values())
print(f'Audited {total} closures.\n')
for cat in ['PHYSICS', 'PROOF', 'DATA', 'BOOKKEEP']:
    label = {
        'PHYSICS':  'PHYSICS  (UQFF-primitive derivations)',
        'PROOF':    'PROOF    (foundational-constant derivations)',
        'DATA':     'DATA     (observational-data closures)',
        'BOOKKEEP': 'BOOKKEEP (null/checkpoint, no derivation)',
    }[cat]
    print(f'  {label:<48s} {counts[cat]:>4d} / {total}')
print(f'  unclassified (no file):                          {len(unclassified):>4d}')
print()
for cat in ['PHYSICS', 'PROOF', 'DATA', 'BOOKKEEP']:
    print(f'\n--- {cat} samples ---')
    for s in samples[cat]:
        print(s)
