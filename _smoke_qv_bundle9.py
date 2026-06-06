"""Smoke test for quantum-variable bundle 9 (rho_UA, rho_Ui, v_SCm, rho_A, rho_SCm)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Pre-existing constants
tests.append(('RHO_UA  = 7.09e-36 J/m^3', abs(m.RHO_UA  - 7.09e-36) / 7.09e-36 < 1e-10))
tests.append(('RHO_SCM = 7.09e-37 J/m^3', abs(m.RHO_SCM - 7.09e-37) / 7.09e-37 < 1e-10))

# Existing G-lock relationship: RHO_UA = 10 * RHO_SCM
tests.append(('RHO_UA / RHO_SCM = 10 (G-lock per Map sec 2)',
              abs(m.RHO_UA / m.RHO_SCM - 10.0) < 1e-10))

# Bundle-9 new named constants
tests.append(('RHO_VAC_UI = 2.84e-36 J/m^3',
              abs(m.RHO_VAC_UI - 2.84e-36) / 2.84e-36 < 1e-10))
tests.append(('RHO_VAC_AETHER = 1.0e-23 J/m^3',
              abs(m.RHO_VAC_AETHER - 1.0e-23) / 1.0e-23 < 1e-10))
tests.append(('V_SCM_PROPAGATION = 1.0e8 m/s',
              abs(m.V_SCM_PROPAGATION - 1.0e8) / 1.0e8 < 1e-10))

# v_SCm propagation ~ c/3 check
c_over_3 = 299792458.0 / 3.0
tests.append(('V_SCM_PROPAGATION ~ c/3 (within 10%)',
              abs(m.V_SCM_PROPAGATION - c_over_3) / c_over_3 < 0.10))

# Numerical rho_Ui relationship: rho_Ui / rho_SCm = 4.0056e0 (level-13 amplification)
ratio_ui_scm = m.RHO_VAC_UI / m.RHO_SCM
tests.append(('RHO_VAC_UI / RHO_SCM ~ 4.006 (level-13 amplification)',
              abs(ratio_ui_scm - 4.0056) < 0.01))
# rho_Ui / rho_UA = 0.4006
ratio_ui_ua = m.RHO_VAC_UI / m.RHO_UA
tests.append(('RHO_VAC_UI / RHO_UA ~ 0.4006',
              abs(ratio_ui_ua - 0.4006) < 0.001))

# _e_react_from_scm_velocity literal eval at t=0 -> 709.0 (spec literal arithmetic)
e_react_lit = m._e_react_from_scm_velocity(t_days=0.0)
expected    = (7.09e-37 * (1.0e8) ** 2) / 1.0e-23
tests.append(('E_react literal eq 7 at t=0 = 709.0',
              abs(e_react_lit - expected) / expected < 1e-12))
tests.append(('E_react literal eq 7 at t=0 ~ 709',
              abs(e_react_lit - 709.0) < 1e-6))

# Decay at t=1000 days -> 709 * exp(-0.5) = 430.06
e_react_t1000 = m._e_react_from_scm_velocity(t_days=1000.0)
expected_decay = expected * math.exp(-0.5)
tests.append(('E_react literal eq 7 at t=1000d ~ 430.06',
              abs(e_react_t1000 - expected_decay) / expected_decay < 1e-12))

# Decay ratio check: e_react(1000)/e_react(0) = exp(-0.5) ~ 0.6065
ratio_decay = e_react_t1000 / e_react_lit
tests.append(('E_react decay ratio exp(-0.5) at t=1000d',
              abs(ratio_decay - math.exp(-0.5)) < 1e-12))

# Aether perturbation magnitude at default (eta=1e-22, T_smunu=1.123e7)
pert = m._aether_perturbation_eta_T()
tests.append(('eta * T_smunu = 1.123e-15 (spec eq 6/7)',
              abs(pert - 1.123e-15) / 1.123e-15 < 1e-12))

# Aether perturbation propagates into _aether_metric_a_munu (bundle 1)
metric = m._aether_metric_a_munu()
# Diagonal: [1 + pert, -1 + pert, -1 + pert, -1 + pert]
tests.append(('A_munu[0] = 1 + 1.123e-15', abs(metric[0] - (1.0 + 1.123e-15)) < 1e-20))
tests.append(('A_munu[1] = -1 + 1.123e-15', abs(metric[1] - (-1.0 + 1.123e-15)) < 1e-20))
tests.append(('A_munu[2] = -1 + 1.123e-15', abs(metric[2] - (-1.0 + 1.123e-15)) < 1e-20))
tests.append(('A_munu[3] = -1 + 1.123e-15', abs(metric[3] - (-1.0 + 1.123e-15)) < 1e-20))

# Bundle-5 _e_react_full (A_0=1e46) still produces ~1e46 at t=0 (parallel composer)
e_react_full_0 = m._e_react_full(t_days=0.0)
tests.append(('Bundle 5 _e_react_full(0) = 1e46 (parallel form)',
              abs(e_react_full_0 - 1.0e46) / 1.0e46 < 1e-10))

# Both composers decay at same kappa rate (consistency)
e_full_t1000 = m._e_react_full(t_days=1000.0)
ratio_full   = e_full_t1000 / e_react_full_0
tests.append(('Bundle 5 and Bundle 9 E_react share kappa decay',
              abs(ratio_full - ratio_decay) < 1e-12))

# Regression: bundle 8 still loads
tests.append(('Bundle 8 _defect_factor_ug1 still present',
              hasattr(m, '_defect_factor_ug1')))
tests.append(('Bundle 8 T_S_SUN_K still present',
              abs(m.T_S_SUN_K - 5778.0) < 1e-10))
tests.append(('Bundle 8 F_TRZ_DEFAULT still present',
              abs(m.F_TRZ_DEFAULT - 0.1) < 1e-15))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
