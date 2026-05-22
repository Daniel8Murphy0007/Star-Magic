import json
from collections import Counter
d = json.load(open('unified_closure_audit.json'))
print('--- All entries by status ---')
print(Counter(e['status'] for e in d))
print()
print('--- 3 CALIBRATED entries ---')
for e in d:
    if e['status']=='CALIBRATED':
        print(f"  id={e['id']}  claim={e.get('claim','')[:80]}")
        print(f"    chain={e.get('chain','')[:120]}")
print()
print('--- 8 AXIOM entries ---')
for e in d:
    if e['status']=='AXIOM':
        print(f"  id={e['id']}  claim={e.get('claim','')[:80]}")
print()
print('--- 5 POSTULATED entries ---')
for e in d:
    if e['status']=='POSTULATED':
        print(f"  id={e['id']}  claim={e.get('claim','')[:80]}")
