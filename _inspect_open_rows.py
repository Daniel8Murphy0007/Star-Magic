import csv
rows = list(csv.DictReader(open('master_closures.csv', encoding='utf-8')))
for i, r in enumerate(rows, start=2):
    if i in (608, 614, 621, 627):
        print(f"L{i}:")
        for k, v in r.items():
            print(f"   {k!r:15s} = {v!r}")
        print()
