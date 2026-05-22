"""Find rows with non-numeric predicted/observed (the '18 NaN' from profiler)."""
import csv, math
rows = list(csv.DictReader(open('master_closures.csv', encoding='utf-8')))
bad = []
for i, r in enumerate(rows, start=2):
    p_raw = (r.get('predicted') or '').strip()
    o_raw = (r.get('observed') or '').strip()
    try:
        p = float(p_raw) if p_raw else float('nan')
        o = float(o_raw) if o_raw else float('nan')
        if math.isnan(p) or math.isnan(o):
            bad.append((i, r, 'blank_or_nan'))
    except ValueError:
        bad.append((i, r, 'parse_fail'))

print(f"rows with non-numeric predicted/observed: {len(bad)}")
for i, r, why in bad:
    print(f"L{i} [{why}] status={r['status']:12s} pred={(r.get('predicted') or '')!r:25s} obs={(r.get('observed') or '')!r:25s} closure={(r['closure'] or '')[:50]}")
