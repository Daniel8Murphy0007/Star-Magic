"""Smoke test for quantum-variable bundle 5 (gamma, E_react, f_quasi, R_b)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Spec-named constants
tests.append(('GAMMA_DECAY_PER_DAY = 5e-5',    abs(m.GAMMA_DECAY_PER_DAY - 5.0e-5) < 1e-15))
tests.append(('KAPPA_REACT_PER_DAY = 5e-4',    abs(m.KAPPA_REACT_PER_DAY - 5.0e-4) < 1e-15))
tests.append(('E_REACT_AMPLITUDE = 1e46',      abs(m.E_REACT_AMPLITUDE - 1.0e46) / 1e46 < 1e-10))
tests.append(('R_B_FIELD_BUBBLE = 1.496e13',   abs(m.R_B_FIELD_BUBBLE - 1.496e13) < 1.0))

# E_react(t=0) = 1e46 (spec match eq 6)
e0 = m._e_react_full(t_days=0.0)
tests.append(('E_react(0) = 1e46', abs(e0 - 1.0e46) / 1.0e46 < 1e-10))

# E_react(t=1000) = 1e46 * exp(-0.5) = 6.0653e45
e1000 = m._e_react_full(t_days=1000.0)
expected = 1.0e46 * math.exp(-0.5)
tests.append(('E_react(1000d) = 1e46*exp(-0.5)', abs(e1000 - expected) / expected < 1e-10))

# E_react monotonic decay
tests.append(('E_react monotonic decay', e1000 < e0))

# U_m buildup at t=1000d, t_n=1000 (cos(1000 pi)=1) -> 1 - exp(-0.05) ~ 0.04877
buildup = m._u_m_buildup_factor(t_days=1000.0, t_0_days=0.0)
expected_buildup = 1.0 - math.exp(-5.0e-5 * 1000.0 * math.cos(math.pi * 1000.0))
tests.append(('U_m buildup at t=1000d (literal eq)', abs(buildup - expected_buildup) < 1e-15))
# cos(1000 pi) = 1 numerically -> 1 - exp(-0.05) = 0.0487705...
tests.append(('U_m buildup ~ 0.04877 (spec eq 2)', abs(buildup - 0.04877) < 1e-4))

# U_m buildup at t=0 -> 1 - exp(0) = 0
b0 = m._u_m_buildup_factor(t_days=0.0)
tests.append(('U_m buildup(t=0) = 0', b0 == 0.0))

# Negative-time test: at t_n = -1, cos(-pi) = -1, so 1 - exp(+gamma*t) is negative
b_neg = m._u_m_buildup_factor(t_days=1000.0, t_0_days=1001.0)  # t_n = -1
# t_n = -1 -> cos(-pi) = -1, so 1 - exp(+0.05) < 0
tests.append(('U_m buildup negative for t_n=-1', b_neg < 0.0))

# Cross-check: existing _l23_e_react_decay (dimensionless envelope) * 1e46 = full E_react
env = m._l23_e_react_decay(1000.0)
full = m._e_react_full(t_days=1000.0,
                        kappa_per_day=m._L23_KAPPA_PER_DAY)
tests.append(('_l23_e_react_decay * 1e46 = _e_react_full (matching kappa)',
              abs(full - env * 1.0e46) / max(full, 1e-300) < 1e-10))

# f_quasi already captured from bundle 3
tests.append(('F_QUASI_DEFAULT = 0.01 (bundle 3)', abs(m.F_QUASI_DEFAULT - 0.01) < 1e-15))

# R_b numerically equals R_J_MAGNETIC (semantic alias)
tests.append(('R_B_FIELD_BUBBLE == R_J_MAGNETIC',
              m.R_B_FIELD_BUBBLE == m.R_J_MAGNETIC))

# Existing Heaviside S(r-R_b) Ug2 form using R_b alias works
ug2_in  = m._u_g2_heliosphere_uqff(r=1.0e12, R_b=m.R_B_FIELD_BUBBLE)
ug2_out = m._u_g2_heliosphere_uqff(r=m.R_B_FIELD_BUBBLE, R_b=m.R_B_FIELD_BUBBLE)
tests.append(('Ug2(r<R_b) = 0', ug2_in == 0.0))
tests.append(('Ug2(r=R_b) > 0', ug2_out > 0.0))

# Regression: bundle 4 still loads
tests.append(('M_BH_SGR_A still present', abs(m.M_BH_SGR_A - 8.15e36) / 8.15e36 < 1e-10))
tests.append(('_mu_j_time_dependent still present', hasattr(m, '_mu_j_time_dependent')))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
