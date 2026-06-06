"""Verify q-scope cycle / envelope / stability algorithms against the
50-row qscope_thz_earth_core_20231003.csv (5 batches, master clock 0-645s)."""
import csv, math, importlib
m = importlib.import_module('uqff_pure_calculator')

# Load CSV (strip # comment lines)
with open('qscope_thz_earth_core_20231003.csv', 'r', encoding='utf-8') as f:
    rows = list(csv.DictReader(line for line in f if not line.startswith('#')))

t       = [int(r['t_seconds_relative']) for r in rows]
batch   = [int(r['batch'])              for r in rows]
flow    = [r['flow_state']              for r in rows]
v_ch1   = [float(r['V_ch1_peak_V'])     for r in rows]
v_ch2   = [float(r['V_ch2_peak_V'])     for r in rows]
ch2_sh  = [r['ch2_shape']               for r in rows]

tests = []

# ---- _qscope_event_times_s ----
frc = m._qscope_event_times_s(t, flow, 'full_reversal_cycle')
tests.append(('full_reversal_cycle events at t=[68,190,320,450,580]',
              frc == [68.0, 190.0, 320.0, 450.0, 580.0]))

rcc = m._qscope_event_times_s(t, flow, 'reversal_cycle_complete')
tests.append(('reversal_cycle_complete events at t=[124,242,372,502,632]',
              rcc == [124.0, 242.0, 372.0, 502.0, 632.0]))

stab_evts = m._qscope_event_times_s(t, flow, 'stabilized')
# 'stabilized' appears once per batch in batches 2..5 (rows 20,30,40,50)
tests.append(('stabilized events at t=[255,385,515,645]',
              stab_evts == [255.0, 385.0, 515.0, 645.0]))

# Length mismatch raises
try:
    m._qscope_event_times_s([1, 2], ['a', 'b', 'c'], 'a'); ok = False
except ValueError:
    ok = True
tests.append(('event_times_s rejects length mismatch', ok))

# ---- _qscope_mean_period_s ----
# diffs of frc: [122, 130, 130, 130] -> mean = 512/4 = 128.0
tests.append(('mean period full_reversal_cycle = 128.0 s',
              m._qscope_mean_period_s(frc) == 128.0))
# diffs of rcc: [118, 130, 130, 130] -> mean = 508/4 = 127.0
tests.append(('mean period reversal_cycle_complete = 127.0 s',
              m._qscope_mean_period_s(rcc) == 127.0))
# Empty / single -> 0.0
tests.append(('mean period of [] = 0.0',     m._qscope_mean_period_s([])     == 0.0))
tests.append(('mean period of [5.0] = 0.0',  m._qscope_mean_period_s([5.0])  == 0.0))

# ---- _qscope_batch_envelope ----
env = m._qscope_batch_envelope(v_ch1, batch)
tests.append(('envelope has 5 batch entries', set(env.keys()) == {1, 2, 3, 4, 5}))
# Batch 1 Ch1: [0.80,0.75,0.70,0.65,0.60,0.60,0.70,0.65,0.60,0.60] max=0.80 min=0.60 mean=0.665
tests.append(('batch 1 Ch1 max = 0.80',  env[1]['max']      == 0.80))
tests.append(('batch 1 Ch1 min = 0.60',  env[1]['min']      == 0.60))
tests.append(('batch 1 Ch1 pp_range = 0.20',
              abs(env[1]['pp_range'] - 0.20) < 1e-12))
tests.append(('batch 1 Ch1 mean = 0.665',
              abs(env[1]['mean'] - 0.665) < 1e-12))
tests.append(('batch 1 Ch1 n = 10', env[1]['n'] == 10.0))
# Batches 2..5 Ch1 identical: [0.60,0.65,0.60,0.55,0.50,0.60,0.55,0.50,0.50,0.50]
# max=0.65 min=0.50 mean=0.555
for b in (2, 3, 4, 5):
    tests.append((f'batch {b} Ch1 max = 0.65', env[b]['max']      == 0.65))
    tests.append((f'batch {b} Ch1 min = 0.50', env[b]['min']      == 0.50))
    tests.append((f'batch {b} Ch1 mean = 0.555',
                  abs(env[b]['mean'] - 0.555) < 1e-12))
    tests.append((f'batch {b} Ch1 pp_range = 0.15',
                  abs(env[b]['pp_range'] - 0.15) < 1e-12))

# ---- _qscope_stability_index ----
stab = m._qscope_stability_index(v_ch1, batch, reference_batch=1)
# batch 1 std = sqrt(0.04525/10) = sqrt(0.004525) ~ 0.0672681...
# batch 2 std = sqrt(0.02725/10) = sqrt(0.002725) ~ 0.0521920...
# ratio batch2/batch1 = 0.0521920 / 0.0672681 ~ 0.77588
ref_std = math.sqrt(0.04525 / 10.0)
b2_std  = math.sqrt(0.02725 / 10.0)
expected_ratio = b2_std / ref_std
tests.append(('stability[1] == 1.0 (self-ratio)',
              abs(stab[1] - 1.0) < 1e-12))
for b in (2, 3, 4, 5):
    tests.append((f'stability[{b}] == sqrt(0.02725/0.04525) ~ 0.77589',
                  abs(stab[b] - expected_ratio) < 1e-12))
    tests.append((f'stability[{b}] < 1.0 (flow stabilizing)',  stab[b] < 1.0))

# Zero reference -> all zeros
stab_zero = m._qscope_stability_index([1, 1, 1, 1], [1, 1, 2, 2], reference_batch=1)
tests.append(('zero reference std -> all-zero map',
              all(v == 0.0 for v in stab_zero.values())))

# ---- _qscope_cycle_summary ----
summ = m._qscope_cycle_summary(t, flow)
tests.append(('summary full_reversal_cycle_count = 5',
              summ['full_reversal_cycle_count'] == 5))
tests.append(('summary full_reversal_cycle_mean_period_s = 128.0',
              summ['full_reversal_cycle_mean_period_s'] == 128.0))
tests.append(('summary reversal_cycle_complete_count = 5',
              summ['reversal_cycle_complete_count'] == 5))
tests.append(('summary reversal_cycle_complete_mean_period_s = 127.0',
              summ['reversal_cycle_complete_mean_period_s'] == 127.0))
tests.append(('summary full_reversal_cycle_times_s matches frc',
              summ['full_reversal_cycle_times_s'] == frc))

# ---- Cross-check Ch2 inverted_sinusoid via same primitive ----
# ch2_shape == 'inverted_sinusoid' rows: 6,10,15,16,19,25,26,29,35,36,39,45,46,49
# t = 68,124,190,203,242,320,333,372,450,463,502,580,593,632 (14 events)
# diffs (13 entries) sum to 564; mean period = 564/13 = 43.3846...
inv = m._qscope_event_times_s(t, ch2_sh, 'inverted_sinusoid')
tests.append(('Ch2 inverted_sinusoid 14 events',  len(inv) == 14))
tests.append(('Ch2 inverted_sinusoid first/last = 68 / 632',
              inv[0] == 68.0 and inv[-1] == 632.0))
tests.append(('Ch2 inverted_sinusoid mean period = 564/13 s',
              abs(m._qscope_mean_period_s(inv) - 564.0/13.0) < 1e-12))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
if passed != total:
    for name, ok in tests:
        if not ok:
            print(f"  FAIL {name}")
