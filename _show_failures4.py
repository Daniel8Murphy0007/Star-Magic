import json
with open('_test_results.json') as f:
    data = json.load(f)
pf = data['per_file']
for fname in ['CP2', 'CP4']:
    entry = pf.get(fname, {})
    print(f'=== {fname}: {entry.get("fail", 0)} failures ===')
    for r in entry.get('failures', []):
        cls = r.get('class', '?')
        err = r.get('error', '?')
        print(f'  {cls}: {err}')
    print()
