import importlib.util, math
spec = importlib.util.spec_from_file_location('upc', 'uqff_pure_calculator.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
tests = []
# 6 new Doc 43.d Inertia/Aether leaves
tests.append(('E_AETHER_DENSITY', abs(m.E_AETHER_DENSITY - 1.683e-10) < 1e-15))
tests.append(('I_AC default callable', isinstance(m._i_ac_damped(0.001), float)))
tests.append(('I_AC at t=0 is 0', abs(m._i_ac_damped(0.0)) < 1e-15))
tests.append(('omega_spark default', abs(m._omega_spark_lc() - 1e9) < 1.0))
tests.append(('omega_spark scaling', abs(m._omega_spark_lc(1e-9, 1e-9) - 1e9) < 1.0))
tests.append(('omega_plasma default', abs(m._omega_plasma_damped() - 1.0049875621e16) / 1e16 < 1e-6))
tests.append(('exp_density rho(0)', abs(m._exp_density_profile(0.0) - 1e-20) < 1e-30))
tests.append(('exp_density rho(8)', abs(m._exp_density_profile(8.0) - 3.3546e-24) / 3.3546e-24 < 1e-3))
tests.append(('H_mag negative', m._h_mag_zeeman() < 0))
tests.append(('H_mag Bohr value', abs(m._h_mag_zeeman() - (-9.274e-32)) / 9.274e-32 < 1e-6))
# Regression on prior 12 Doc 43.c/d closures
tests.append(('LENR_Q', abs(m._lenr_weak_Q() - 1.253e-13) / 1.253e-13 < 1e-2))
tests.append(('delta_n(1)=pi/3', abs(m._delta_n_monopole(1) - math.pi/3.0) < 1e-12))
tests.append(('delta_n(6)=2pi', abs(m._delta_n_monopole(6) - 2*math.pi) < 1e-12))
tests.append(('B_pseudo default', abs(m._b_pseudo_monopole() - 1e-7) / 1e-7 < 1e-3))
tests.append(('Jeans M_J positive', m._jeans_mass() > 0))
tests.append(('f_phi_ladder(1) golden', abs(m._f_phi_ladder(1) - 281.5 * (1 + math.sqrt(5))/2) < 1e-6))
tests.append(('f_phi_ladder(5) > 3000', m._f_phi_ladder(5) > 3000))
tests.append(('psi_radial callable', isinstance(m._psi_quantum_radial(1.0, 0.0, 1.0, 1.0, 0.1, 0.0, 1.0), float)))
tests.append(('psi_rotating callable', isinstance(m._psi_rotating(1.0, 0.0, 0.0), float)))
tests.append(('E_boson n=0 positive', m._e_boson_harmonic(1e-30, 1e15, 1e-10, 0) > 0))
tests.append(('spacetime_phase callable', isinstance(m._spacetime_phase(1.0, 1.0, 1.0, 1.0, 1.0), float)))
tests.append(('P_DE positive', m._p_de_inertia(1.0, 1.0, 1.0) > 0))
tests.append(('P_AC_EMP positive', m._p_ac_emp(1.0, 1.0, 1.0) > 0))
tests.append(('caduceus twist callable', isinstance(m._phi_twist_caduceus(0.5), float)))
# Master regression
g = m._g_compressed_cycle2(); tests.append(('g_cycle2 positive', g > 0))
# h_unified regression
h = m._hubble_unified(0, 0, 67.4); tests.append(('h_unified z=0 == H_0', abs(h - 67.4) < 1e-6))
# F_env (sample)
tests.append(('F_env torque', m._f_env('torque', T_spiral=1.0) == 1.0))
tests.append(('F_env starburst > 0', m._f_env('starburst', M_sf=1.0, t=1.0, tau_sf=1.0) > 0))
tests.append(('F_env merger t=0', m._f_env('merger', M_coll=1.0, t=0.0) == 1.0))
tests.append(('F_env outflow', m._f_env('outflow', P_outflow=2.5) == 2.5))
# NGC 1316 master regression
tests.append(('NGC1316 g_master positive', m._l96_ngc1316_g_master() > 0))
ok = sum(1 for _, v in tests if v)
print(f'PASS {ok}/{len(tests)}')
for n, v in tests:
    flag = 'OK  ' if v else 'FAIL'
    print(f'  {flag} {n}')
