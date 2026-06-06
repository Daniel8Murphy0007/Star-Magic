"""Smoke test for quantum-variable bundle 3 (f_Heaviside, i, H_SCm, lambda_i, j)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []

# Named constants
tests.append(('F_HEAVISIDE_DEFAULT = 0.01', abs(m.F_HEAVISIDE_DEFAULT - 0.01) < 1e-15))
tests.append(('F_QUASI_DEFAULT = 0.01',     abs(m.F_QUASI_DEFAULT - 0.01) < 1e-15))
tests.append(('H_SCM_DEFAULT = 1.0',        abs(m.H_SCM_DEFAULT - 1.0) < 1e-15))
tests.append(('LAMBDA_I_DEFAULT = 1.0',     abs(m.LAMBDA_I_DEFAULT - 1.0) < 1e-15))

# Ug2 heliospheric: literal eq 6 evaluation with Sun defaults
# k_2*(rho_UA+rho_SCm)*M_s/r^2 * S * (1+delta_sw*v_sw) * H_SCm * E_react
# = 1.2 * 7.80e-36 * 1.989e30 / 2.238e26 * 1 * 5001 * 1 * 1e46 ~ 4.16e18 J/m^3
ug2 = m._u_g2_heliosphere_uqff()
expected_ug2 = (1.2 * (7.09e-36 + 7.09e-37) * 1.989e30 / (1.496e13 ** 2)
                * 1.0 * (1.0 + 0.01 * 5.0e5) * 1.0 * 1.0e46)
tests.append(('U_g2 heliospheric matches literal eq 6',
              abs(ug2 - expected_ug2) / expected_ug2 < 1e-10))

# Ug2 verify (1 + delta_sw * v_sw) factor: delta=0.01, v=5e5 -> factor = 1 + 5000 = 5001
ug2_factor1 = m._u_g2_heliosphere_uqff(delta_sw=0.0)
ug2_factor2 = m._u_g2_heliosphere_uqff(delta_sw=0.01, v_sw=5.0e5)
ratio = ug2_factor2 / ug2_factor1
tests.append(('(1 + delta_sw*v_sw) factor = 5001 (Sun spec)', abs(ratio - 5001.0) < 1e-6))

# Ug2 with H_SCm=1.1 should be 10% higher
ug2_high = m._u_g2_heliosphere_uqff(H_SCm=1.1)
tests.append(('H_SCm=1.1 multiplies Ug2 by 1.1', abs(ug2_high / ug2 - 1.1) < 1e-10))

# Ug2 Heaviside step: r < R_b -> S = 0 -> Ug2 = 0
ug2_inside = m._u_g2_heliosphere_uqff(r=1.0e12, R_b=1.0e13)
tests.append(('Heaviside S(r<R_b) = 0 -> U_g2 = 0', ug2_inside == 0.0)

)
# Ug3 string-summed magnetism: Sun defaults t=0, B_j=1e3, k_3=1.8, E=1e46 -> ~1.8e49
ug3 = m._u_g3_string_magnetism_uqff()
tests.append(('U_g3 string ~ 1.8e49', abs(ug3 - 1.8e49) / 1.8e49 < 0.01))

# Ug3 cos factor: t=0 -> cos(0) = 1, scaling linear in B_j_sum and k_3
ug3_2x = m._u_g3_string_magnetism_uqff(B_j_sum=2.0e3)
tests.append(('U_g3 linear in B_j_sum', abs(ug3_2x / ug3 - 2.0) < 1e-12))

# Regression: existing _l96_bearden_Um_trz (f_heaviside scaling 1 + 1e13 * f_heaviside)
import inspect
sig = inspect.signature(m._l96_bearden_Um_trz)
tests.append(('_l96_bearden_Um_trz has f_heaviside param',
              'f_heaviside' in sig.parameters))
# Magnetic Um with f_heaviside=0.01 -> scaling (1 + 1e11)
um_with = m._l96_bearden_Um_trz(t=1.0, t_n=0.0, mu_sum_over_r=2.26e10,
                                  gamma=0.00005, P_SCm=1.0, E_react=1.0e46,
                                  f_heaviside=0.01, f_quasi=0.01)
um_without = m._l96_bearden_Um_trz(t=1.0, t_n=0.0, mu_sum_over_r=2.26e10,
                                     gamma=0.00005, P_SCm=1.0, E_react=1.0e46,
                                     f_heaviside=0.0, f_quasi=0.01)
ratio_h = um_with / um_without
tests.append(('Heaviside scaling (1 + 1e11)', abs(ratio_h - (1.0 + 1.0e11)) / 1.0e11 < 1e-6))

# Regression: existing _l96_bearden_Ui_trz (lambda_i parameter)
sig2 = inspect.signature(m._l96_bearden_Ui_trz)
tests.append(('_l96_bearden_Ui_trz has lambda_i param',
              'lambda_i' in sig2.parameters))
# Inertia spec at Sun: lambda_i=1, rho_SCm=7.09e-37, rho_UA=7.09e-36,
# omega_s=2.5e-6, t_n=0 -> cos=1, f_TRZ=0.1 -> 1.1 -> result ~ 1.38e-47
ui = m._l96_bearden_Ui_trz(t=0.0, t_n=0.0, lambda_i=1.0,
                            omega_s=2.5e-6, f_TRZ=0.1)
expected_ui = 1.0 * 7.09e-37 * 7.09e-36 * 2.5e-6 * 1.0 * 1.1
tests.append(('U_i spec value ~ 1.38e-47', abs(ui - expected_ui) / expected_ui < 1e-10))

# Bundle 2 regressions
tests.append(('R_J_MAGNETIC = 1.496e13 still present', abs(m.R_J_MAGNETIC - 1.496e13) < 1.0))
tests.append(('F_FEEDBACK_DEFAULT = 0.1 still present', abs(m.F_FEEDBACK_DEFAULT - 0.1) < 1e-15))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
for name, ok in tests:
    print(f"  {'OK  ' if ok else 'FAIL'} {name}")
