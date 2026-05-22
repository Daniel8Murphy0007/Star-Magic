"""Emit the canonical gap list: rows the calibration scan needs to address."""
import csv, math
rows = list(csv.DictReader(open('master_closures.csv', encoding='utf-8')))
gaps = []
for r in rows:
    try:
        e = float((r.get('error_pct') or '').strip())
    except (ValueError, TypeError):
        continue
    if math.isnan(e):
        continue
    if r.get('status', '').strip() == 'OPEN_PREDICTIONS':
        continue
    a = abs(e)
    if a < 1.0:  # within 1% = locked
        continue
    if abs(a - 9999) < 1e-6:
        continue
    gaps.append((a, r))

gaps.sort(key=lambda x: -x[0])
print(f"=== Out-of-sync UQFF systems (|error_pct| >= 1%): {len(gaps)} ===\n")
print(f"{'rank':>4} {'|err%|':>12} {'ID':>5} {'status':10} {'closure / label':50} {'script':25}")
print("-" * 120)
for i, (a, r) in enumerate(gaps[:30], 1):
    closure = (r.get('closure') or r.get('label') or '')[:50]
    print(f"{i:4d} {a:12.4e} {r.get('ID','') or '':>5} {r['status']:10} {closure:50} {(r.get('script') or '')[:25]:25}")

# Save full list to CSV for subagent
with open('_calibration_gap_targets.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['rank', 'abs_error_pct', 'ID', 'status', 'closure', 'predicted', 'observed', 'script', 'label'])
    for i, (a, r) in enumerate(gaps, 1):
        w.writerow([i, f"{a:.6e}", r.get('ID', ''), r['status'],
                    r.get('closure', ''), r.get('predicted', ''), r.get('observed', ''),
                    r.get('script', ''), r.get('label', '')])
print(f"\nFull list -> _calibration_gap_targets.csv ({len(gaps)} rows)")
