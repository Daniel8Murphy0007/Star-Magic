"""Produce final Job B v5.78 update list - grouped, numerically sorted, readable."""
import csv, re

rows = list(csv.DictReader(open('_job_b_categorization.csv', encoding='utf-8')))

def keyf(r):
    m = re.match(r'PAPER_(\d+)([a-z]?)', r['filename'])
    return (int(m.group(1)) if m else 99999, m.group(2) if m else '')

rows.sort(key=keyf)

BUCKETS = {
    'C': ('SI constants derivations (alpha, h, c, G)',
          'Add three-anchor SI closure banner (Sess 237-241); STRUCTURAL tier'),
    'A': ('Cosmology / Lambda / dark energy / vacuum / CMB / BBN',
          'Add closing section: v5.78 27-decade vacuum ledger (PAPER_1170) + xi=13/3 (PAPER_1171/1172)'),
    'B': ('Unified field / Lagrangian / master-equation DERIVATIONS',
          'Add v5.78 closed-Lagrangian cross-ref (PAPER_1159-1167 G1-G8) + CP4 #254'),
    'D': ('Kaluza-Klein / extra-dim / sub-mm / Yukawa',
          'Forward-pointer to PAPER_1171/1172 (xi=13/3 R26+KK lock) + PAPER_1173'),
    'E': ('Falsifiability / P-suite / DESI/Euclid/CMB-S4',
          'Forward-pointer to PAPER_1174 P1-P14 suite + PAPER_1177-1180'),
    'F': ('Index / catalog / KB / registry',
          'Update tables; add P6-P14, CP4 #254-#264, ref PAPER_1167/1170/1174'),
    'G': ('Millennium Prize (Yang-Mills, Navier-Stokes, Riemann, P-vs-NP, Hodge)',
          'Verify v5.78 closed-Lagrangian/ledger cross-ref is present (Sess 225 baseline)'),
    'H': ('Specific astrophysical system / framework application', 'NO UPDATE'),
    'I': ('LENR / nuclear / Kozima / Holmlid / Pons-Fleischmann', 'NO UPDATE unless ledger/xi enters; verify only'),
    'J': ('Audits / doctrine / methodology', 'Case-by-case'),
}

print('=' * 78)
print('Job B Phase B1: v5.78 Update Categorization')
print(f'Total papers scanned: {len(rows)}')
print('=' * 78)
print()
print('BUCKET COUNTS')
print('-' * 78)
counts = {}
for r in rows:
    counts[r['bucket']] = counts.get(r['bucket'], 0) + 1
for b in 'CABDEFGHIJ':
    desc = BUCKETS[b][0]
    print(f'  {b}  {desc:<60}  {counts.get(b,0):>5}')
print('-' * 78)
upd = sum(counts.get(b, 0) for b in 'ABCDEFGJ')
no = counts.get('H', 0) + counts.get('I', 0)
print(f'  NEEDS UPDATE (A,B,C,D,E,F,G,J): {upd}')
print(f'  NO UPDATE    (H, I):            {no}')
print()
print('=' * 78)
print('FULL LIST OF PAPERS NEEDING v5.78 UPDATE (144 papers)')
print('=' * 78)
for b in 'CABDEFGJ':
    sub = [r for r in rows if r['bucket'] == b]
    if not sub:
        continue
    desc, action = BUCKETS[b]
    print()
    print(f'BUCKET {b} - {desc}   ({len(sub)} papers)')
    print(f'  ACTION: {action}')
    print('-' * 78)
    for r in sub:
        print(f"  {r['filename']}")
