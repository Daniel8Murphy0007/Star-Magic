"""One-shot closure-report gap audit (Session 304+)."""
import uqff_pure_calculator as u

r = u._constant_closure_report()
items = []
for name, info in r.items():
    if name == '_summary':
        continue
    err = info.get('err_pct', None)
    status = info.get('status', '?')
    items.append((name, err, status, info.get('si_anchor'), info.get('uqff_sat')))

items_sorted = sorted(
    [x for x in items if isinstance(x[1], (int, float))],
    key=lambda t: (0 if t[1] != t[1] else -abs(t[1])),
)

print(f"{'name':<14} {'err%':>14} {'status':<14} {'si_anchor':>16} {'uqff_sat':>16}")
print('-' * 80)
for n, e, s, t, v in items_sorted:
    if isinstance(e, float) and e != e:
        estr = 'NaN'
    else:
        estr = f'{e:.3e}'
    tstr = f'{t:.6e}' if isinstance(t, (int, float)) else str(t)
    vstr = f'{v:.6e}' if isinstance(v, (int, float)) else str(v)
    print(f'{n:<14} {estr:>14} {s:<14} {tstr:>16} {vstr:>16}')

print()
print('--- _summary ---')
for k, v in r.get('_summary', {}).items():
    print(f'  {k}: {v}')

print()
print('--- BROKEN / ANCHOR / HARDCODED rows ---')
for n, info in r.items():
    if n == '_summary':
        continue
    if info.get('status') in ('broken', 'anchor', 'hardcoded'):
        print(f"  {n:<14} status={info['status']:<10}  err={info['err_pct']}  notes={info['notes']}")
