"""Stamp batch 4 commit (0bfa00ee) on CSV rows for PAPER_022-029 + b-variants."""
from pathlib import Path
import csv

CSV = Path('_job_b_categorization.csv')
COMMIT = '0bfa00ee'
SESSION = 'Session 262'

TARGET_IDS = {
    'PAPER_022', 'PAPER_023', 'PAPER_024', 'PAPER_025', 'PAPER_025b',
    'PAPER_026', 'PAPER_026b', 'PAPER_027', 'PAPER_028', 'PAPER_029',
}
TARGET_FN_HINT = {
    'PAPER_026': 'Sterile_Neutrino_Mass_Generation',  # disambiguate dup
}

rows = list(csv.DictReader(CSV.open(encoding='utf-8')))
stamped = 0
for r in rows:
    pid = r.get('paper_id', '').strip()
    if pid not in TARGET_IDS:
        continue
    if pid == 'PAPER_026':
        if 'Mass_Generation' not in r.get('filename', ''):
            continue
    if r.get('status', '').startswith('DONE'):
        continue
    r['status'] = 'DONE v5.78'
    r['session'] = SESSION
    r['commit'] = COMMIT
    if r.get('bucket', '') and not r['bucket'].endswith('*'):
        r['bucket'] = r['bucket'] + '*'
    stamped += 1

with CSV.open('w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print(f'Stamped {stamped} rows with {COMMIT}')
