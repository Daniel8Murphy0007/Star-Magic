import importlib.util, math
spec = importlib.util.spec_from_file_location('upc', 'uqff_pure_calculator.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

tests = []
# 4 new constants
tests.append(('EPSILON_SW = 0.001', m.EPSILON_SW == 0.001))
tests.append(('RHO_VAC_SW = 8e-21', m.RHO_VAC_SW == 8e-21))
tests.append(('ETA_AETHER_COUPLING = 1e-22', m.ETA_AETHER_COUPLING == 1e-22))
tests.append(('K_GRAV_LADDER = (1.5,1.2,1.8,1.0)', m.K_GRAV_LADDER == (1.5, 1.2, 1.8, 1.0)))
tests.append(('BETA_I = 0.6 (existing)', m.BETA_I == 0.6))

# Solar wind buoyancy modulation
# eps*rho = 0.001 * 8e-21 = 8e-24 (sub-epsilon vs 1.0 in IEEE-754 -> verify product in isolation)
eps_rho = m.EPSILON_SW * m.RHO_VAC_SW
tests.append(('eps_sw * rho_sw = 8e-24', abs(eps_rho - 8e-24) < 1e-35))
modulation = 1.0 + eps_rho
tests.append(('modulation factor ≈ 1 (sub-eps perturbation)', abs(modulation - 1.0) < 1e-15))

# U_bi solar wind sample with spec values
# Eq 2: U_b1 = -0.6 * 1.39e26 * 7.3e-16 * 3.20e16 * 1 * 1 * 1 ≈ -1.94e27
u_b1 = m._u_bi_solar_wind_modulated(U_gi=1.39e26, Omega_g=7.3e-16,
                                     M_bh=3.20e16, d_g=1.0,
                                     U_UA=1.0, t_n=0.0)
# = -0.6 * 1.39e26 * 7.3e-16 * 3.20e16 * 1 * 1 = -1.948e27
tests.append(('U_bi solar wind spec ≈ -1.94e27', abs(u_b1 - (-1.948416e27))/1.948416e27 < 1e-3))

# With epsilon_sw=0 -> same value (spec: no change)
u_b1_no_sw = m._u_bi_solar_wind_modulated(U_gi=1.39e26, Omega_g=7.3e-16,
                                          M_bh=3.20e16, d_g=1.0,
                                          U_UA=1.0, t_n=0.0,
                                          epsilon_sw=0.0)
tests.append(('U_bi epsilon_sw=0 unchanged', abs(u_b1_no_sw - u_b1)/abs(u_b1) < 1e-15))

# t_n=1 -> cos(pi) = -1, flips sign
u_b1_tn1 = m._u_bi_solar_wind_modulated(U_gi=1.39e26, Omega_g=7.3e-16,
                                         M_bh=3.20e16, d_g=1.0,
                                         U_UA=1.0, t_n=1.0)
tests.append(('U_bi t_n=1 flips sign', u_b1_tn1 > 0 and abs(u_b1_tn1 + u_b1) < 1e10))

# Aether metric A_munu
# Eq 4: eta * T_smunu = 1e-22 * 1.123e7 = 1.123e-15
a = m._aether_metric_a_munu(T_smunu=1.123e7)
tests.append(('A_00 = 1 + 1.123e-15', abs(a[0] - (1.0 + 1.123e-15)) < 1e-30))
tests.append(('A_11 = -1 + 1.123e-15', abs(a[1] - (-1.0 + 1.123e-15)) < 1e-30))
tests.append(('A_22 = A_11', a[2] == a[1]))
tests.append(('A_33 = A_11', a[3] == a[1]))
tests.append(('A_munu length 4', len(a) == 4))

# Defaults — verify the eta*T product in isolation (FP rounding when added to 1.0 is non-zero but lossy)
eta_T = m.ETA_AETHER_COUPLING * 1.123e7
tests.append(('eta * T_smunu = 1.123e-15 (isolated)', abs(eta_T - 1.123e-15) < 1e-30))
a_def = m._aether_metric_a_munu()
tests.append(('A_munu default perturbation present (FP-aware)', abs(a_def[0] - 1.0) > 0.0 and abs(a_def[0] - 1.0) < 1e-14))

# eta=0 -> pure Minkowski
a_min = m._aether_metric_a_munu(T_smunu=1.123e7, eta=0.0)
tests.append(('Minkowski A_00 = 1', a_min[0] == 1.0))
tests.append(('Minkowski A_11 = -1', a_min[1] == -1.0))

# Eq 10 sum k_i * U_gi spec
# (1.5 * 1.39e26) + (1.2 * 1.18e53) + (1.8 * 1.8e49) + (1.0 * 2.50e-20)
sum_k_ugi = 1.5 * 1.39e26 + 1.2 * 1.18e53 + 1.8 * 1.8e49 + 1.0 * 2.50e-20
# Spec quotes ≈ 1.42e53
tests.append(('sum k_i U_gi ≈ 1.42e53', abs(sum_k_ugi - 1.416e53)/1.42e53 < 1e-2))
# K_GRAV_LADDER usage
k = m.K_GRAV_LADDER
sum_via_ladder = k[0]*1.39e26 + k[1]*1.18e53 + k[2]*1.8e49 + k[3]*2.50e-20
tests.append(('K_GRAV_LADDER matches spec sum', abs(sum_via_ladder - sum_k_ugi)/sum_k_ugi < 1e-15))

# Regression on prior leaves
tests.append(('Doc 43 _l96_universal_permanence', m._l96_universal_permanence() > 0))
tests.append(('Doc 43 _l96_nonlocal_jump_prob ≈ 0.313', abs(m._l96_nonlocal_jump_prob() - 0.313) < 1e-3))
tests.append(('Doc 43 _l96_reactive_energy_density ≈ 9.86e14', abs(m._l96_reactive_energy_density() - 9.86e14)/9.86e14 < 1e-2))
tests.append(('Doc 43.c _l96_u_h_higgs_field', m._l96_u_h_higgs_field() > 0))
tests.append(('Doc 43.c _l96_rho_ua_scm_n_t', m._l96_rho_ua_scm_n_t() > 0))
tests.append(('Doc 43.c _l96_v_radial_doppler', abs(m._l96_v_radial_doppler() - m.C_LIGHT*(-1e-5))/abs(m.C_LIGHT*1e-5) < 1e-9))
tests.append(('Doc 43.c _l96_decay_rate_ratio', m._l96_decay_rate_ratio() > 0))
tests.append(('Doc 43.d E_AETHER_DENSITY', abs(m.E_AETHER_DENSITY - 1.683e-10) < 1e-15))
tests.append(('Doc 43.d _i_ac_damped', abs(m._i_ac_damped(0.0)) < 1e-15))
tests.append(('Doc 43.e E_0_HYDROGEN_REACTOR', abs(m.E_0_HYDROGEN_REACTOR - 1.683e-37) < 1e-40))
tests.append(('Doc 43.e _e_k_quantum_wave k=6', m._e_k_quantum_wave(6, (6/26.0)*m.T_EARTH_MOON/4, 2, 2) > 0))
tests.append(('Doc 43.b _l96_agn_feedback_Ug4 default', m._l96_agn_feedback_Ug4() > 0))
tests.append(('Doc 43.b _l96_final_parsec_Ug4 default', m._l96_final_parsec_Ug4() > 0))
tests.append(('Master g_compressed_cycle2', m._g_compressed_cycle2() > 0))
tests.append(('Master NGC1316 g_master', m._l96_ngc1316_g_master() > 0))

ok = sum(1 for _, v in tests if v)
print(f'PASS {ok}/{len(tests)}')
for nm, v in tests:
    flag = 'OK  ' if v else 'FAIL'
    print(f'  {flag} {nm}')
