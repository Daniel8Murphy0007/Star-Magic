"""Level-1 desk audit of all 570 rows in master_closures.csv.

For each row, check:
- script_exists: does the referenced .py file exist on disk?
- has_predicted: is predicted field populated?
- has_observed: is observed field populated?
- error_pct_parsed: can error_pct be parsed as float?
- residual_tier: EXACT / <0.01% / 0.01-0.1% / 0.1-1% / >=1%

Emit:
- confirmed_derivations.csv  (full per-row audit table)
- console summary
"""
import csv
from pathlib import Path
from collections import Counter

ROOT = Path('.')
src = ROOT / 'master_closures.csv'
out = ROOT / 'confirmed_derivations.csv'

with open(src, encoding='utf-8', newline='') as f:
    reader = csv.DictReader(f)
    rows = list(reader)


def tier_of(err_str):
    try:
        e = float(err_str)
    except (TypeError, ValueError):
        return 'UNPARSED'
    if e == 0.0:
        return 'EXACT'
    if e < 0.01:
        return '<0.01%'
    if e < 0.1:
        return '0.01-0.1%'
    if e < 1.0:
        return '0.1-1%'
    return '>=1%'


out_rows = []
counters = {
    'script_exists': 0, 'script_missing': 0,
    'has_predicted': 0, 'has_observed': 0,
    'parseable_error': 0, 'tier': Counter(),
    'category': Counter(), 'status': Counter(),
    'all_checks_passed': 0,
}

for r in rows:
    rid = r['ID'].strip()
    script_ref = (r.get('script') or '').strip()
    script_path = ROOT / script_ref if script_ref else None
    script_exists = bool(script_path and script_path.exists())
    has_pred = bool((r.get('predicted') or '').strip())
    has_obs = bool((r.get('observed') or '').strip())
    err = (r.get('error_pct') or '').strip()
    try:
        float(err)
        parseable = True
    except (TypeError, ValueError):
        parseable = False
    tier = tier_of(err)
    cat = r.get('category', '')
    status = r.get('status', '')

    all_pass = script_exists and has_pred and has_obs and parseable

    counters['script_exists'] += int(script_exists)
    counters['script_missing'] += int(not script_exists)
    counters['has_predicted'] += int(has_pred)
    counters['has_observed'] += int(has_obs)
    counters['parseable_error'] += int(parseable)
    counters['tier'][tier] += 1
    counters['category'][cat] += 1
    counters['status'][status] += 1
    counters['all_checks_passed'] += int(all_pass)

    out_rows.append({
        'ID': rid,
        'name': r.get('name', ''),
        'category': cat,
        'status': status,
        'script_exists': 'YES' if script_exists else 'NO',
        'has_predicted': 'YES' if has_pred else 'NO',
        'has_observed': 'YES' if has_obs else 'NO',
        'predicted': (r.get('predicted') or '')[:30],
        'observed': (r.get('observed') or '')[:30],
        'error_pct': err,
        'tier': tier,
        'confirmed_level1': 'YES' if all_pass else 'NO',
    })

with open(out, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)

print(f"=== LEVEL-1 DESK AUDIT: {len(rows)} rows ===")
print()
print(f"Category breakdown:           {dict(counters['category'])}")
print(f"Status breakdown:             {dict(counters['status'])}")
print()
print(f"Script-on-disk:               {counters['script_exists']}/{len(rows)}")
print(f"Script MISSING:               {counters['script_missing']}/{len(rows)}")
print(f"Has predicted value:          {counters['has_predicted']}/{len(rows)}")
print(f"Has observed value:           {counters['has_observed']}/{len(rows)}")
print(f"Parseable error_pct:          {counters['parseable_error']}/{len(rows)}")
print()
print(f"Residual tier distribution:")
ORDER = ['EXACT', '<0.01%', '0.01-0.1%', '0.1-1%', '>=1%', 'UNPARSED']
for t in ORDER:
    cnt = counters['tier'].get(t, 0)
    if cnt:
        print(f"  {t:<12} {cnt}")
print()
print(f"ALL LEVEL-1 CHECKS PASSED:    {counters['all_checks_passed']}/{len(rows)}")
print()
print(f"Full per-row table written to: {out.name}  ({len(out_rows)} rows)")
