import json
with open('_test_results.json') as f:
    data = json.load(f)
pf = data['per_file']
for fname in ['CondensedPhysics2.py', 'CondensedPhysics4.py']:
    entry = pf.get(fname, {})
    failures = [r for r in entry.get('results', []) if not r.get('passed')]
    print(f'=== {fname}: {len(failures)} failures ===')
    for r in failures:
        cls = r.get('class', '?')
        err = r.get('error', '?')
        print(f'  {cls}: {err}')
    print()
