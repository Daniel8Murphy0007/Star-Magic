"""One-shot: list rows with |error_pct|==100 (T-PRED sentinel candidates) and NaN error_pct rows."""
import csv, math
rows = list(csv.DictReader(open('master_closures.csv', encoding='utf-8')))
hundred, nans = [], []
for i, r in enumerate(rows, start=2):
    ep = (r.get('error_pct') or '').strip()
    if not ep:
        continue
    try:
        v = float(ep)
        if math.isnan(v):
            nans.append((i, r))
        elif abs(abs(v) - 100.0) < 1e-6:
            hundred.append((i, r))
    except ValueError:
        nans.append((i, r))

print(f"=== |error_pct|=100 rows: {len(hundred)} ===")
for l, r in hundred:
    print(f"L{l} status={r['status']:10s} closure={r['closure'][:60]:60s} pred={r['predicted']} obs={r['observed']}")

print(f"\n=== NaN error_pct rows: {len(nans)} ===")
for l, r in nans:
    print(f"L{l} status={r['status']:10s} closure={r['closure'][:60]:60s} pred={r['predicted']!r} obs={r['observed']!r} ep={r['error_pct']!r}")
