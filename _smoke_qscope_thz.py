"""Smoke test for q-scope THz oscilloscope artifact (Earth-core 1.246 THz, 10-image bundle)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Named constants
tests.append(('F_QSCOPE_MEASURED_HZ = 1.246e12 Hz',
              abs(m.F_QSCOPE_MEASURED_HZ - 1.246e12) / 1.246e12 < 1e-12))
tests.append(('Z_OSCILLOSCOPE_OHM = 50.0 Ohm',
              abs(m.Z_OSCILLOSCOPE_OHM - 50.0) < 1e-12))

# Pre-existing OMEGA_SCM primitive remains 1.25e12 (rounded center)
tests.append(('OMEGA_SCM primitive = 1.25e12 Hz (rounded center)',
              abs(m.OMEGA_SCM - 1.25e12) / 1.25e12 < 1e-12))

# Measured value lies inside 1.2-1.3 THz band
tests.append(('1.246 THz is in [1.2, 1.3] THz band',
              m._qscope_measured_in_band() is True))
# Below band -> False
tests.append(('1.0 THz NOT in [1.2, 1.3] THz band',
              m._qscope_measured_in_band(f_hz=1.0e12) is False))
# Above band -> False
tests.append(('1.5 THz NOT in [1.2, 1.3] THz band',
              m._qscope_measured_in_band(f_hz=1.5e12) is False))
# Boundary inclusive
tests.append(('1.20 THz IS in band (lower boundary)',
              m._qscope_measured_in_band(f_hz=1.20e12) is True))
tests.append(('1.30 THz IS in band (upper boundary)',
              m._qscope_measured_in_band(f_hz=1.30e12) is True))

# Pre-existing _l24_qscope_in_band on primitive still True
tests.append(('Pre-existing _l24_qscope_in_band() still True',
              m._l24_qscope_in_band() is True))

# Oscilloscope power: P = V^2 / Z at spec defaults V=0.8 V, Z=50 -> 0.0128 W
P_spec = m._oscilloscope_power_W()
tests.append(('P(0.8 V, 50 Ohm) = 0.0128 W (spec match)',
              abs(P_spec - 0.0128) < 1e-12))

# Verify spec amplitude series produces sensible powers (10-image envelope check)
# Ch1 amplitudes (V) from spec: 0.8, 0.75, 0.7, 0.65, 0.6, 0.6, 0.7, 0.65, 0.6, 0.6
spec_ch1_amps_V = [0.8, 0.75, 0.7, 0.65, 0.6, 0.6, 0.7, 0.65, 0.6, 0.6]
powers = [m._oscilloscope_power_W(V_peak_V=v) for v in spec_ch1_amps_V]
# Max power = (0.8)^2 / 50 = 0.0128 W
tests.append(('Max Ch1 power = 0.0128 W (spec V=0.8 max)',
              abs(max(powers) - 0.0128) < 1e-12))
# Min power = (0.6)^2 / 50 = 0.0072 W
tests.append(('Min Ch1 power = 0.0072 W (spec V=0.6 min)',
              abs(min(powers) - 0.0072) < 1e-12))

# Frequency separation: measured anchor 1.246 vs primitive 1.25 THz = 4 GHz delta
delta_hz = m.OMEGA_SCM - m.F_QSCOPE_MEASURED_HZ
tests.append(('Measured-vs-primitive separation = 4 GHz',
              abs(delta_hz - 4.0e9) < 1e3))

# Quick consistency with bundle 9 (regression)
tests.append(('Bundle 9 RHO_VAC_AETHER still present',
              abs(m.RHO_VAC_AETHER - 1.0e-23) / 1.0e-23 < 1e-10))
tests.append(('Bundle 9 V_SCM_PROPAGATION still present',
              abs(m.V_SCM_PROPAGATION - 1.0e8) / 1.0e8 < 1e-10))
tests.append(('Bundle 8 _defect_factor_ug1 still present',
              hasattr(m, '_defect_factor_ug1')))

# Angular frequency relation (informational, not a stored derivation)
omega_measured = 2.0 * math.pi * m.F_QSCOPE_MEASURED_HZ
tests.append(('omega = 2 pi * 1.246e12 ~ 7.83e12 rad/s (informational)',
              abs(omega_measured - 7.83e12) / 7.83e12 < 1e-3))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
