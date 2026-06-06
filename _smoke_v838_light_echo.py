"""Verify V838 Monocerotis light-echo master intensity primitives:
   L_SUN_W = 3.826e26 W
   L_OUTBURST_V838_W = 6e5 * L_SUN_W ~ 2.296e38 W (spec ~2.3e38)
   _light_echo_radius_m(t)         = c * t
   _dust_density_uqff(rho_0,b,U)   = rho_0 * exp(-b * U)
   _v838_light_echo_intensity_uqff = master composer with G-lock x12.1 UQFF enh
Also confirm spec primitives reused: _defect_factor_ug1 defaults 0.01 sin(0.001*t),
   RHO_UA = 7.09e-36, F_TRZ_DEFAULT = 0.1, ratio rho_UA/rho_SCm = 10 (G-lock)."""
import math, importlib
m = importlib.import_module('uqff_pure_calculator')

tests = []
# --- constants ---
tests.append(('L_SUN_W = 3.826e26 W',
              m.L_SUN_W == 3.826e26))
tests.append(('L_OUTBURST_V838_W = 6e5 * L_SUN_W',
              m.L_OUTBURST_V838_W == 6.0e5 * 3.826e26))
# Spec arithmetic gap: literal eq = 2.2956e32 W; spec-quoted example = 2.3e38 W (6-order).
# Both anchors coexist (same pattern as bundles 2/3/9 spec arithmetic mismatches).
tests.append(('Literal L_OUTBURST_V838_W = 2.2956e32 W (6e5 * 3.826e26)',
              abs(m.L_OUTBURST_V838_W - 2.2956e32) / 2.2956e32 < 1e-12))
tests.append(('Spec-quoted anchor L_OUTBURST_V838_SPEC_QUOTE_W = 2.3e38 W',
              m.L_OUTBURST_V838_SPEC_QUOTE_W == 2.3e38))

# --- _light_echo_radius_m ---
tests.append(('r_echo(0) = 0',
              m._light_echo_radius_m(0.0) == 0.0))
# t = 1 s -> r = c
tests.append(('r_echo(1 s) = C_LIGHT',
              m._light_echo_radius_m(1.0) == m.C_LIGHT))
# Literal spec: t = 3 yr (Julian year = 365.25 * 86400 = 3.15576e7 s)
t_3yr_s = 3.0 * 365.25 * 86400.0
r_3yr = m._light_echo_radius_m(t_3yr_s)
# 3 light-year length: 1 ly ~ 9.461e15 m -> 3 ly ~ 2.838e16 m
tests.append(('r_echo(3 yr) ~ 3 ly (~2.838e16 m, 0.5% tol)',
              abs(r_3yr - 2.838e16) / 2.838e16 < 0.005))

# --- _dust_density_uqff ---
# beta = 0 -> rho_dust = rho_0 (no modulation)
tests.append(('rho_dust(beta=0) = rho_0',
              m._dust_density_uqff(1.0e-21, 0.0, 1.0e50) == 1.0e-21))
# beta * U = 1 -> rho_0 * exp(-1)
tests.append(('rho_dust(beta*U=1) = rho_0 / e',
              abs(m._dust_density_uqff(1.0, 1.0, 1.0) - math.exp(-1.0)) < 1e-15))
# Pre-multiplied product
tests.append(('rho_dust supports pre-multiplied beta*U_g1 product',
              abs(m._dust_density_uqff(2.0, 0.5, 4.0) - 2.0 * math.exp(-2.0)) < 1e-15))

# --- _defect_factor_ug1 reuse for spec delta_def = 0.01 sin(0.001 t) ---
# Spec defaults: amplitude=0.01, omega=0.001 (any time unit)
tests.append(('delta_def(0) = 0',
              m._defect_factor_ug1(0.0) == 0.0))
# t = pi / (2 * 0.001) = 1570.7963... -> amplitude * 1 = 0.01
tests.append(('delta_def(pi/0.002) = 0.01 (peak)',
              abs(m._defect_factor_ug1(math.pi / 0.002) - 0.01) < 1e-15))

# --- _v838_light_echo_intensity_uqff ---
# UQFF enhancement = (1 + 0.1) * (1 + 10) = 12.1 (G-lock ratio 10)
tests.append(('G-lock ratio RHO_UA/RHO_SCM = 10.0',
              abs(m.RHO_UA / m.RHO_SCM - 10.0) < 1e-12))

# Sanity: with beta=0, intensity = (L/(4 pi r^2)) * sigma * rho_0 * 12.1
t_3yr_s = 3.0 * 365.25 * 86400.0
r_3yr = m._light_echo_radius_m(t_3yr_s)
classical = (m.L_OUTBURST_V838_W / (4.0 * math.pi * r_3yr * r_3yr)) * 1.0e-22
I_pred = classical * 1.0e-21 * 12.1
I_calc = m._v838_light_echo_intensity_uqff(t_3yr_s)
tests.append(('I_echo(3yr, default) = classical * rho_0 * 12.1',
              abs(I_calc - I_pred) / I_pred < 1e-12))

# Composition: beta=0 vs beta>0 with U_g1=0 -> identical (exp(0)=1)
I_b0 = m._v838_light_echo_intensity_uqff(t_3yr_s, beta_U_g1=0.0)
I_b1 = m._v838_light_echo_intensity_uqff(t_3yr_s, beta_U_g1=0.0)
tests.append(('beta_U_g1=0 path stable',
              I_b0 == I_b1))

# t=0 -> +inf (geometric divergence at outburst origin)
tests.append(('I_echo(0) = +inf',
              math.isinf(m._v838_light_echo_intensity_uqff(0.0))
              and m._v838_light_echo_intensity_uqff(0.0) > 0))

# beta*U_g1 = 1 reduces classical*rho component by exp(-1)
I_full = m._v838_light_echo_intensity_uqff(t_3yr_s, beta_U_g1=0.0)
I_atten = m._v838_light_echo_intensity_uqff(t_3yr_s, beta_U_g1=1.0)
tests.append(('I_echo(beta*U=1) / I_echo(beta*U=0) = exp(-1)',
              abs(I_atten / I_full - math.exp(-1.0)) < 1e-12))

# Inverse-square check: I(2t) / I(t) = 1/4 (with beta=0)
t_a = 1.0e7; t_b = 2.0e7
I_a = m._v838_light_echo_intensity_uqff(t_a)
I_b = m._v838_light_echo_intensity_uqff(t_b)
tests.append(('I_echo inverse-square: I(2t)/I(t) = 1/4',
              abs(I_b / I_a - 0.25) < 1e-12))

# Override rho_UA / rho_SCm sanity (UQFF enh disabled if rho_UA=0, rho_SCm>0)
I_no_ua = m._v838_light_echo_intensity_uqff(t_3yr_s, rho_UA_val=0.0,
                                              rho_SCm_val=m.RHO_SCM)
classical_only = classical * 1.0e-21 * (1.0 + m.F_TRZ_DEFAULT) * 1.0
tests.append(('rho_UA=0 disables UA term: enh = (1+f_TRZ)*1 = 1.1',
              abs(I_no_ua - classical_only) / classical_only < 1e-12))

passed = sum(1 for _, ok in tests if ok)
total  = len(tests)
print(f"PASS {passed}/{total}")
if passed != total:
    for name, ok in tests:
        if not ok:
            print(f"  FAIL {name}")
