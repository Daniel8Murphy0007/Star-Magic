import csv
from collections import defaultdict

# Read CSV
data = list(csv.DictReader(open('COMPLETE_PHYSICS_CLASS_INVENTORY.csv', 'r', encoding='utf-8')))

# Group by file
by_file = defaultdict(list)
for r in data:
    by_file[r['SourceFile']].append(f"{r['ClassName']} ({r['PhysicsType']})")

# Write output
with open('QUICK_CLASS_LOOKUP.txt', 'w', encoding='utf-8') as f:
    f.write('QUICK PHYSICS CLASS LOOKUP BY SOURCE FILE\n')
    f.write('='*80 + '\n\n')
    
    for file, classes in sorted(by_file.items(), key=lambda x: int(x[0].replace('source','').replace('.cpp',''))):
        f.write(f'{file}:\n')
        for cls in classes:
            f.write(f'  {cls}\n')
        f.write('\n')

print("Quick lookup file created: QUICK_CLASS_LOOKUP.txt")
