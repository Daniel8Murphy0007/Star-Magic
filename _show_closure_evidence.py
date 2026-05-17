import csv
rows = [r for r in csv.DictReader(open('master_closures.csv', encoding='utf-8')) if r['status'] == 'OK']
def err_key(r):
    try: return abs(float(r['error_pct']))
    except: return 9e9
rows.sort(key=err_key)
print(f'TOTAL OK CLOSURES: {len(rows)}')
print(f'UNIQUE PRODUCING SCRIPTS: {len(set(r["script"] for r in rows))}')
print()
print('=== 20 BEST (lowest |residual %|) ===')
for r in rows[:20]:
    print(f"{r['ID']:>4}  {r['label'][:38]:38}  pred={r['predicted'][:16]:16}  obs={r['observed'][:16]:16}  err={r['error_pct'][:10]:>10}%  {r['script']}")
print()
print('=== 15 LARGEST RESIDUALS (still OK / accepted) ===')
for r in rows[-15:]:
    print(f"{r['ID']:>4}  {r['label'][:38]:38}  pred={r['predicted'][:16]:16}  obs={r['observed'][:16]:16}  err={r['error_pct'][:10]:>10}%  {r['script']}")
print()
# residual bands
bands = {'EXACT (0%)': 0, '<0.01%': 0, '0.01-0.1%': 0, '0.1-1%': 0, '1-5%': 0, '>=5%': 0}
for r in rows:
    e = err_key(r)
    if e == 0: bands['EXACT (0%)'] += 1
    elif e < 0.01: bands['<0.01%'] += 1
    elif e < 0.1: bands['0.01-0.1%'] += 1
    elif e < 1: bands['0.1-1%'] += 1
    elif e < 5: bands['1-5%'] += 1
    else: bands['>=5%'] += 1
print('=== RESIDUAL DISTRIBUTION ===')
for k, v in bands.items():
    print(f'  {k:15}  {v:4}  ({100*v/len(rows):5.1f}%)')
