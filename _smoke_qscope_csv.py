"""Verify qscope_thz_earth_core_20231003.csv P_ch1_W / P_ch2_W columns
match uqff_pure_calculator._oscilloscope_power_W(V_peak) for every row.
Also verify f_carrier_THz lies inside the catalog 1.2-1.3 THz q-scope band."""
import csv, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []
with open('qscope_thz_earth_core_20231003.csv', 'r', encoding='utf-8') as f:
    rows = [r for r in csv.DictReader(line for line in f if not line.startswith('#'))]

tests.append(('CSV row count = 30', len(rows) == 30))

for r in rows:
    sid = int(r['signal_id'])
    v1 = float(r['V_ch1_peak_V']); v2 = float(r['V_ch2_peak_V'])
    p1_csv = float(r['P_ch1_W']);  p2_csv = float(r['P_ch2_W'])
    p1_calc = m._oscilloscope_power_W(V_peak_V=v1)
    p2_calc = m._oscilloscope_power_W(V_peak_V=v2)
    tests.append((f'Signal {sid:02d} Ch1 P matches calc',
                  abs(p1_csv - p1_calc) < 1e-9))
    tests.append((f'Signal {sid:02d} Ch2 P matches calc',
                  abs(p2_csv - p2_calc) < 1e-9))
    f_hz = float(r['f_carrier_THz']) * 1.0e12
    tests.append((f'Signal {sid:02d} carrier in 1.2-1.3 THz band',
                  m._qscope_measured_in_band(f_hz=f_hz)))

# Aggregate: max Ch1 = 0.0128 W (signal 1, V=0.8); min Ch1 = 0.005 W (signals 15,18,19,20)
ch1_powers = [float(r['P_ch1_W']) for r in rows]
tests.append(('Max Ch1 power = 0.0128 W (signal 1, V=0.8)',
              abs(max(ch1_powers) - 0.0128) < 1e-9))
tests.append(('Min Ch1 power = 0.005 W (V=0.5)',
              abs(min(ch1_powers) - 0.005) < 1e-9))

# Timeline continuity: t_seconds_relative monotone non-decreasing across both batches
ts = [int(r['t_seconds_relative']) for r in rows]
tests.append(('t_seconds_relative monotone',
              all(ts[i] <= ts[i+1] for i in range(len(ts)-1))))
tests.append(('Total span = 385 s (16:39:35 -> 16:46:00)',
              ts[-1] - ts[0] == 385))

# Batch split: rows 1-10 are batch 1; rows 11-20 are batch 2; rows 21-30 are batch 3
b1 = [int(r['batch']) for r in rows[:10]]
b2 = [int(r['batch']) for r in rows[10:20]]
b3 = [int(r['batch']) for r in rows[20:]]
tests.append(('Rows 1-10 are batch 1', all(b == 1 for b in b1)))
tests.append(('Rows 11-20 are batch 2', all(b == 2 for b in b2)))
tests.append(('Rows 21-30 are batch 3', all(b == 3 for b in b3)))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
if passed != total:
    for name, ok in tests:
        if not ok:
            print(f"  FAIL {name}")
