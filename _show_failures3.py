import json
with open('_test_results.json') as f:
    data = json.load(f)
pf = data['per_file']
for fname in ['CP2', 'CP4']:
    entry = pf.get(fname, {})
    print(f'=== {fname}: {entry.get("fail", 0)} failures ===')
    print('Entry keys:', list(entry.keys()))
    results = entry.get('results', [])
    if results:
        print('First result keys:', list(results[0].keys()))
        failures = [r for r in results if not r.get('passed', True)]
        for r in failures:
            cls = r.get('class', r.get('name', '?'))
            err = r.get('error', r.get('err', '?'))
            print(f'  {cls}: {err}')
    else:
        # Try other keys
        for k, v in entry.items():
            if isinstance(v, list) and v:
                print(f'  {k}: {v[:2]}')
    print()
