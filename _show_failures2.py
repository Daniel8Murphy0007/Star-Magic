import json
with open('_test_results.json') as f:
    data = json.load(f)
pf = data['per_file']
print('Files in per_file:', list(pf.keys()))
gt = data['grand_total']
print('Grand total:', gt)
for fname, entry in pf.items():
    fails = entry.get('fail', 0)
    if fails > 0:
        print(f'\n=== {fname}: {fails} failures ===')
        for r in entry.get('results', []):
            if not r.get('passed'):
                cls = r.get('class', '?')
                err = r.get('error', '?')
                print(f'  {cls}: {err}')
