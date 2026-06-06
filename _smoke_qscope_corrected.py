"""Verify q-scope corrected-interpretation primitives:
   V_peak = V_pp / 2
   V_eff  = V_pp / (2 sqrt(2))      [sinusoid RMS]
   P_rms  = V_eff^2 / Z              [Z = 50 ohm]
   QSCOPE_AMP_RANGE_A = 3.102 A; QSCOPE_DIFFERENTIAL_AMP_A = 6.205 A
Spec literal worked example: V_pp = 1.00 V -> V_eff = 0.35 (2dec), P_rms = 0.00245 W."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')
tests = []

# Constants
tests.append(('QSCOPE_AMP_RANGE_A = 3.102',
              m.QSCOPE_AMP_RANGE_A == 3.102))
tests.append(('QSCOPE_DIFFERENTIAL_AMP_A = 6.205',
              m.QSCOPE_DIFFERENTIAL_AMP_A == 6.205))
# Spec literals (independent): 3.102 and 6.205 are quoted constants; 2*3.102=6.204
# in IEEE-754, off spec by 1e-3 — match literal values, do not derive one from other.
tests.append(('dA approx 2 * range within spec precision (1e-3)',
              abs(m.QSCOPE_DIFFERENTIAL_AMP_A - 2*m.QSCOPE_AMP_RANGE_A) < 1.1e-3))

# V_peak from V_pp
tests.append(('V_pp=1.00 -> V_peak=0.50',
              m._v_peak_from_pp(1.00) == 0.50))
tests.append(('V_pp=0.95 -> V_peak=0.475',
              m._v_peak_from_pp(0.95) == 0.475))
tests.append(('V_pp=0.70 -> V_peak=0.35',
              m._v_peak_from_pp(0.70) == 0.35))

# V_eff (RMS) for sinusoid
exact_eff_1V = 1.0 / (2.0 * math.sqrt(2.0))   # ~0.35355
tests.append(('V_pp=1.00 -> V_eff = 1/(2 sqrt(2)) ~ 0.35355',
              abs(m._v_eff_sinusoid(1.00) - exact_eff_1V) < 1e-12))
tests.append(('V_pp=1.00 V_eff rounds to 0.35 V (spec 2-dec display)',
              round(m._v_eff_sinusoid(1.00), 2) == 0.35))
tests.append(('V_pp=0.95 V_eff rounds to 0.34 V (spec match)',
              round(m._v_eff_sinusoid(0.95), 2) == 0.34))
tests.append(('V_pp=0.90 V_eff rounds to 0.32 V (spec match)',
              round(m._v_eff_sinusoid(0.90), 2) == 0.32))
tests.append(('V_pp=0.70 V_eff rounds to 0.25 V (spec match)',
              round(m._v_eff_sinusoid(0.70), 2) == 0.25))
tests.append(('V_pp=0.65 V_eff rounds to 0.23 V (spec match)',
              round(m._v_eff_sinusoid(0.65), 2) == 0.23))
tests.append(('V_pp=0.60 V_eff rounds to 0.21 V (spec match)',
              round(m._v_eff_sinusoid(0.60), 2) == 0.21))
tests.append(('V_pp=0.50 V_eff rounds to 0.18 V (spec match)',
              round(m._v_eff_sinusoid(0.50), 2) == 0.18))

# Composition: V_eff = V_peak / sqrt(2) (consistency)
for v_pp in [1.00, 0.85, 0.70, 0.60]:
    vp = m._v_peak_from_pp(v_pp)
    ve = m._v_eff_sinusoid(v_pp)
    tests.append((f'V_pp={v_pp}: V_eff == V_peak/sqrt(2)',
                  abs(ve - vp/math.sqrt(2.0)) < 1e-12))

# RMS power
tests.append(('V_eff=0.35 V, Z=50 -> P = 0.00245 W (spec exact)',
              abs(m._oscilloscope_power_rms_W(0.35) - 0.00245) < 1e-12))
tests.append(('V_eff=0.25 V, Z=50 -> P = 0.00125 W',
              abs(m._oscilloscope_power_rms_W(0.25) - 0.00125) < 1e-12))
# Direct from V_pp via composition: P = V_pp^2 / (8 Z)
for v_pp in [1.00, 0.85, 0.70]:
    ve = m._v_eff_sinusoid(v_pp)
    p_rms = m._oscilloscope_power_rms_W(ve)
    p_direct = (v_pp * v_pp) / (8.0 * m.Z_OSCILLOSCOPE_OHM)
    tests.append((f'V_pp={v_pp}: P_rms == V_pp^2 / (8 Z)',
                  abs(p_rms - p_direct) < 1e-15))

# Cross-check vs existing peak helper: P_peak / P_rms == 2 for sinusoid
for v_pp in [1.00, 0.80, 0.65]:
    vp = m._v_peak_from_pp(v_pp)
    ve = m._v_eff_sinusoid(v_pp)
    p_peak = m._oscilloscope_power_W(vp)
    p_rms  = m._oscilloscope_power_rms_W(ve)
    tests.append((f'V_pp={v_pp}: P_peak / P_rms == 2.0',
                  abs(p_peak / p_rms - 2.0) < 1e-12))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
if passed != total:
    for name, ok in tests:
        if not ok:
            print(f"  FAIL {name}")
