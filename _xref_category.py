import csv
from collections import Counter
rows = list(csv.DictReader(open('master_closures.csv', encoding='utf-8')))
prov = {r['script']: r['classification'] for r in csv.DictReader(open('provenance_audit_v2.csv', encoding='utf-8'))}

xtab = Counter()
for r in rows:
    cat = r.get('category') or '(empty)'
    cls = prov.get(r['script'], 'EXTRA_or_UNKNOWN')
    xtab[(cls, cat)] += 1

print('Cross-tab: provenance x master_closures.category')
print('-' * 80)
classes = sorted(set(c for c, _ in xtab))
cats = sorted(set(c for _, c in xtab))
print(f'{"class":<24s} ' + ' '.join(f'{c[:24]:>26s}' for c in cats))
for cls in classes:
    cells = ' '.join(f'{xtab[(cls, c)]:>26d}' for c in cats)
    print(f'{cls:<24s} {cells}')

print('\nATOM_SEARCH_FIT scripts and their categories:')
for r in rows:
    if prov.get(r['script']) == 'ATOM_SEARCH_FIT':
        print(f"  {r['script']:<50s} category={r['category']:<35s} err={r['error_pct']}")
