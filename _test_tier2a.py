"""Smoke-test particle-mass pack derived probe."""
import uqff_pure_calculator as u
r = u._l96_particle_mass_pack_probe()
print(r['summary'])
print()
for k in r['anchors']:
    d = r['derived'][k]
    a = r['anchors'][k]
    e = r['err_pct'][k]
    if 'kg' in k:
        print(f'  {k:<14} derived={d:>12.4e}  anchor={a:>12.4e}  err={e:.3f}%')
    elif d < 0.01:
        print(f'  {k:<14} derived={d:>12.6f}  anchor={a:>12.6f}  err={e:.3f}%')
    else:
        print(f'  {k:<14} derived={d:>12.4f}  anchor={a:>12.4f}  err={e:.3f}%')
print()
print('all_within_1pct =', r['all_within_1pct'])
print('all_within_5pct =', r['all_within_5pct'])
