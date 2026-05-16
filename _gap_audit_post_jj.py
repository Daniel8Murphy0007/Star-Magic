import csv
rows = list(csv.DictReader(open('master_closures.csv',encoding='utf-8')))
fails_in_range = [r for r in rows if r['status']!='OK' and 257<=int(r['ID'])<=342]
ok_in_range    = [r for r in rows if r['status']=='OK'  and 257<=int(r['ID'])<=342]
fail_scripts = sorted(set(r['script'] for r in fails_in_range))
ok_scripts   = sorted(set(r['script'] for r in ok_in_range))
print(f"OK in S257-S342:         {len(ok_in_range)} closures from {len(ok_scripts)} scripts")
print(f"PARSE_FAIL in S257-S342: {len(fails_in_range)} rows from {len(fail_scripts)} scripts")
# scripts that have ZERO OK rows
fail_only = [s for s in fail_scripts if s not in ok_scripts]
print(f"Scripts with no OK closure at all: {len(fail_only)}")
for s in fail_only[:25]:
    print(' -', s)
print("\nGlobal totals:")
print(f"  Total rows: {len(rows)}")
print(f"  OK total:   {sum(1 for r in rows if r['status']=='OK')}")
print(f"  PARSE_FAIL: {sum(1 for r in rows if r['status']!='OK')}")
