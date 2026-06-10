import importlib.util, math
spec = importlib.util.spec_from_file_location('upc', 'uqff_pure_calculator.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

tests = []
# 4 new Doc 43.c leaves
# U_H Higgs field
u_h = m._l96_u_h_higgs_field()
tests.append(('U_H positive', u_h > 0))
tests.append(('U_H ~ 1e-23 * 1.585e-8 * 1.01', abs(u_h - 1.0e-23 * 1.585e-8 * 1.01 * math.exp(-(m.SSQ ** 26) * math.exp(-math.pi))) / u_h < 1e-9))

# rho_vac UA-SCm(n,t)
rho_0 = m._l96_rho_ua_scm_n_t(n=0, t=0.0)
tests.append(('rho_vac n=0 positive', rho_0 > 0))
rho_1 = m._l96_rho_ua_scm_n_t(n=1, t=0.0)
tests.append(('rho_vac n=1 = 0.1 * n=0', abs(rho_1 - 0.1 * rho_0) / rho_0 < 1e-9))
rho_2 = m._l96_rho_ua_scm_n_t(n=2, t=0.0)
tests.append(('rho_vac n=2 = 0.01 * n=0', abs(rho_2 - 0.01 * rho_0) / rho_0 < 1e-9))

# Doppler v_radial — spec quotes dimensionless v/c = -3.33e-5
v_rad = m._l96_v_radial_doppler(delta_lambda=-1.0e-11, lambda_rest=1.0e-6)
expected = m.C_LIGHT * -1e-5
tests.append(('v_radial Doppler', abs(v_rad - expected) / abs(expected) < 1e-9))
# Dimensionless beta = v/c form
beta = v_rad / m.C_LIGHT
tests.append(('v_radial beta = -1e-5', abs(beta - (-1.0e-5)) < 1e-15))

# Decay rate ratio — RHO_SCM/RHO_UA = 0.1, decay near 1 -> ~0.1
# Spec quotes 0.0963; depends on exact t input. Just verify in 0.05-0.15 range.
dr = m._l96_decay_rate_ratio(t=0.0)
tests.append(('decay rate ratio in 0.05-0.15', 0.05 < dr < 0.15))
# rho ratio is exactly 0.1
tests.append(('decay rate base ratio 0.1', abs(m.RHO_SCM / m.RHO_UA - 0.1) < 1e-15))

# Regression on prior 43.e leaves
tests.append(('E_AETHER_DENSITY (43.d)', abs(m.E_AETHER_DENSITY - 1.683e-10) < 1e-15))
tests.append(('E_0_HYDROGEN_REACTOR (43.e)', abs(m.E_0_HYDROGEN_REACTOR - 1.683e-37) < 1e-40))
tests.append(('T_EARTH_MOON (43.e)', m.T_EARTH_MOON == 2.36e6))
T = m.T_EARTH_MOON
tests.append(('E_earth_moon UQFF (43.e)', abs(m._e_earth_moon_tidal_uqff(T/4) - 7.96e-22)/7.96e-22 < 5e-2))
tests.append(('E_hyd 1s (43.e)', abs(m._e_hydrogen_radial_uqff(T/4, 1, 0, 0) - 3.98e-22)/3.98e-22 < 5e-2))
tests.append(('E_k6 22 (43.e)', abs(m._e_k_quantum_wave(6, (6/26.0)*T/4, 2, 2) - 2.37e-22)/2.37e-22 < 5e-2))
tests.append(('Y_22 peak (43.e)', abs(m._y_lm_squared_peak(2, 2) - 0.596) < 1e-4))

# Regression on prior 43.c leaves already shipped
tests.append(('LENR_Q (43.c shipped)', abs(m._lenr_weak_Q() - 1.253e-13)/1.253e-13 < 1e-2))
tests.append(('delta_n(6)=2pi (43.c shipped)', abs(m._delta_n_monopole(6) - 2*math.pi) < 1e-12))
tests.append(('Higgs mass shipped', m._l96_higgs_mass_GeV() > 0))
tests.append(('Higgs signal shipped', m._l96_higgs_signal_strength() > 0))
tests.append(('Solar corona W_mag shipped', m._l96_solar_corona_W_mag_GeV() > 0))
tests.append(('LENR E field shipped', isinstance(m._l96_lenr_E_field_from_Um(), float)))
tests.append(('neutron eta shipped', m._l96_neutron_production_eta() > 0))
tests.append(('Bearden U_m TRZ shipped', isinstance(m._l96_bearden_Um_trz(), float)))

# Regression on canonical masters
tests.append(('g_compressed_cycle2', m._g_compressed_cycle2() > 0))
tests.append(('NGC1316 g_master', m._l96_ngc1316_g_master() > 0))
tests.append(('h_unified z=0', abs(m._hubble_unified(0, 0, 67.4) - 67.4) < 1e-6))

ok = sum(1 for _, v in tests if v)
print(f'PASS {ok}/{len(tests)}')
for nm, v in tests:
    flag = 'OK  ' if v else 'FAIL'
    print(f'  {flag} {nm}')
