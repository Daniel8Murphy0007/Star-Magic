"""Tag the 92 S694-S785 rows in master_closures.csv as CALIBRATION_FROM_PARAMETERS.

- Keep all 570 rows.
- Add new 'category' column.
- For rows whose ID is S694-S785 (or numeric >693): set category and status to CALIBRATION_FROM_PARAMETERS.
- For all other 478 rows: category = DERIVATION_FIRST_PRINCIPLES.
- Backup original to master_closures.csv.bak.
"""
import csv
import shutil
from pathlib import Path
from collections import Counter

ROOT = Path('.')
src = ROOT / 'master_closures.csv'
bak = ROOT / 'master_closures.csv.bak'

if not bak.exists():
    shutil.copy(src, bak)
    print(f"Backup created: {bak.name}")
else:
    print(f"Backup already exists: {bak.name} (not overwritten)")

with open(src, encoding='utf-8', newline='') as f:
    reader = csv.DictReader(f)
    fieldnames = list(reader.fieldnames)
    rows = list(reader)

if 'category' not in fieldnames:
    fieldnames.append('category')


def is_calibration(row_id: str) -> bool:
    rid = row_id.strip()
    if rid.startswith('S'):
        try:
            n = int(rid[1:])
            return 694 <= n <= 785
        except ValueError:
            return False
    if rid.isdigit():
        return int(rid) > 693
    return False


tagged_cal = 0
tagged_deriv = 0
for r in rows:
    if is_calibration(r['ID']):
        r['category'] = 'CALIBRATION_FROM_PARAMETERS'
        r['status'] = 'CALIBRATION_FROM_PARAMETERS'
        tagged_cal += 1
    else:
        r['category'] = 'DERIVATION_FIRST_PRINCIPLES'
        tagged_deriv += 1

with open(src, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print()
print(f"Total rows kept: {len(rows)}")
print(f"  DERIVATION_FIRST_PRINCIPLES:    {tagged_deriv}")
print(f"  CALIBRATION_FROM_PARAMETERS:    {tagged_cal}")
print()
print(f"New columns: {fieldnames}")
print()
print("Category distribution by status:")
print(Counter((r['status'], r['category']) for r in rows))
