import csv
rows = list(csv.DictReader(open('master_closures.csv',encoding='utf-8')))
fails = [r for r in rows if r['status']!='OK' and 257<=int(r['ID'])<=342]
print(f"Sample PARSE_FAIL outputs ({len(fails)} total):")
for r in fails[:8]:
    print(f"--- S{r['ID']}  {r['script']}")
    print(r['raw_output'][-500:].replace('\\n','\n'))
    print()
