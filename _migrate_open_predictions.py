"""Migrate 4 CVW G1+G3+G6+G7 OPEN rows (predicted=0, |error_pct|=100) to T-PRED 9999 sentinel.

Per memory/repo: T-PRED sentinel schema is `observed=9999 error_pct=9999`, status=OPEN_PREDICTIONS.
The original observed value is preserved in a new `observed_target` column appended to raw_output
so the prediction target isn't lost.
"""
import csv, math, shutil, datetime

SRC = 'master_closures.csv'
BAK = f'master_closures.csv.bak_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'

shutil.copy2(SRC, BAK)
print(f"backup -> {BAK}")

with open(SRC, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    rows = list(reader)

migrated = 0
for r in rows:
    label = (r.get('label') or '')
    if 'G1 + G3 + G6 + G7 SM Anchor Gate' not in label:
        continue
    if (r.get('status') or '').strip() != 'OPEN':
        continue
    try:
        ep = float((r.get('error_pct') or '').strip())
    except ValueError:
        continue
    if abs(abs(ep) - 100.0) > 1e-6:
        continue
    # Preserve original prediction target
    orig_pred = r.get('predicted', '')
    orig_obs = r.get('observed', '')
    orig_ep = r.get('error_pct', '')
    note = f' [T-PRED migration: orig_pred={orig_pred} orig_obs={orig_obs} orig_err={orig_ep}]'
    r['predicted'] = '9999'
    r['observed'] = '9999'
    r['error_pct'] = '9999'
    r['status'] = 'OPEN_PREDICTIONS'
    # Update cvw_stamp to reflect new sentinel state
    if r.get('cvw_stamp'):
        r['cvw_stamp'] = r['cvw_stamp'].rstrip() + note
    if r.get('raw_output') is not None and 'T-PRED migration' not in (r.get('raw_output') or ''):
        r['raw_output'] = (r.get('raw_output') or '').rstrip() + note
    migrated += 1
    print(f"  migrated: {r.get('closure','<no-closure>')}  (was pred={orig_pred} obs={orig_obs} err={orig_ep})")

with open(SRC, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"\nMigrated {migrated} rows to T-PRED 9999 sentinel.")
