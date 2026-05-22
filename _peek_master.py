import csv
rdr = csv.DictReader(open('master_closures.csv', encoding='utf-8'))
print(rdr.fieldnames)
rows = list(rdr)
print(f'rows={len(rows)}')
print('first row:', rows[0])
cats = set(r.get('category', '') for r in rows)
print('categories:', cats)
def isexact(r):
    try:
        return abs(float(r.get('error_pct') or 'nan')) < 1e-9
    except Exception:
        return False
print('rows with error_pct == 0:', sum(1 for r in rows if isexact(r)))
print('rows with category EXACT:', sum(1 for r in rows if 'EXACT' in (r.get('category') or '')))
