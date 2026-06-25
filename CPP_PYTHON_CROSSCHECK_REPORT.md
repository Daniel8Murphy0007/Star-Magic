# C++ vs Python cross-check report

**C++ source:** `uqff_exact_closures.cpp`
**Tolerance:** relative=1e-06, absolute floor=1e-12
**Total C++ functions:** 630
**Captured numeric values:** 630

## Summary

| Category | Count | Meaning |
|---|---:|---|
| MATCH      |   277 | Python and C++ agree within tolerance |
| DRIFT      |     0 | Both return scalars but they differ |
| MISSING    |   312 | C++ function has no matching Python closure key |
| UNCALLABLE |    41 | Python closure exists but doesn't return a comparable scalar |


## MISSING — 312 C++ functions with no Python counterpart

These are C++-only functions; no PARADOX_TO_CLOSURE key matches.

```
  F_U_alpha_decay_per_day                             = 0.001
  F_U_genesis_n_components                            = 4
  R_d_duality_range_exp                               = 7
  T_SCm_activation_K                                  = 60
  adult_height_cm                                     = 170
  aharonov_bohm_phase                                 = 6.28319
  aharonov_casher_phase                               = 6.28319
  alpha_inverse                                       = 137.04
  aluminum_density                                    = 2700
  amino_acids_20                                      = 20
  atiyah_singer_index                                 = 22
  atm_pressure_kPa                                    = 101.32
  atm_scale_height_km                                 = 8.56
  au_over_r_earth                                     = 23483
  avogadro_n_a_aa                                     = 6.0228
  baryon_fraction                                     = 0.502417
  bertrand_p_1_4                                      = 0.25
  bertrand_uniform_p                                  = 0.25
  bh_4_laws_prefactor                                 = 3.125
  bh_seed_mass_msun                                   = 56160
  blood_glucose_mg_dl                                 = 100
  blood_ph                                            = 7.4
  body_temp_celsius                                   = 37
  bohm_prefactor                                      = 0.0625
  bohr_magneton_lead                                  = 9.27333
  bohr_radius_a0_lead                                 = 5.29033
  bqp_bound                                           = 4
  br_mu_to_e_gamma                                    = 1.26844e-13
  braid_gate_max                                      = 26
  carbon_steel_density                                = 7850
  cdf_w_mass_shift_MeV                                = 74.2427
  cfl_gap_MeV                                         = 85.0812
  ckm_row1_unitarity                                  = 1
  clifford_bundle_qualia                              = 8192
  cnub_temp_K                                         = 1.95393
  co2_atmospheric_ppm                                 = 420
  codons_64                                           = 64
  compton_lambda_lead                                 = 2.42633
  concrete_density_kg_m3                              = 2400
  concrete_fc_MPa                                     = 30
  cosmic_filament_dim                                 = 2
  cosmological_constant                               = 5.95695e-10
  coulomb_ke_lead                                     = 8.983
  coulomb_log_lnL                                     = 16.98
  cr_ankle_eV                                         = 3.61624e+18
  cross_section_A2                                    = 10.5
  dT_pulse_ms                                         = 25
  d_crit_feedback_loops                               = 23
  d_crit_universal                                    = 26
  dark_flow_km_s                                      = 600
  ... and 262 more
```

## UNCALLABLE — Python closure exists but isn't a scalar

41 entries. Python closures returning dicts/None/structured data don't have a single number to compare.

_First 20:_

```
  F_TRZ_identity                            C++=0.1           python=dict
  tsirelson_bound                           C++=2.82843       python=dict
  hayflick_limit                            C++=60            python=dict
  V_little_V_big_ratio                      C++=0.030303      python=dict
  proton_core_density                       C++=2.14644e-36   python=dict
  twenty_six_pre_mass_states                C++=26            python=dict
  t_violation_meson                         C++=0.06029       python=dict
  pulsar_glitch                             C++=3.26419e-07   python=dict
  crab_TeV_cutoff                           C++=79.261        python=dict
  inflaton_n_s                              C++=0.964681      python=dict
  U_UA_coupling_constant                    C++=0.0001        python=dict
  kerr_ringdown_offset_coeff                C++=4.33333       python=dict
  gw170817_phonon_prefactor                 C++=0.666667      python=dict
  transcendental_ln_10                      C++=2.30267       python=dict
  transcendental_ln_2                       C++=0.693167      python=dict
  transcendental_pi_squared                 C++=9.87083       python=dict
  pcr_quantum_triadic                       C++=3             python=dict
  d_bsfg_derived_from_d_crit                C++=6             python=dict
  k_mex_derived_from_phi_5_6                C++=2.08333       python=dict
  f221_f220_qnm_ratio                       C++=0.983426      python=dict
```

## Verdict

Of 277 comparable entries, **100.0% match** within relative tolerance 1e-06.

DRIFT entries are the high-priority audit targets. Each represents a case where the C++ value encoded does not match the current Python closure value. Likely causes: (a) the Python closure formula changed after the C++ port was authored; (b) the C++ value was a target rather than a derived value; (c) a typo. Verify each manually and either update the C++ value or flag the closure as deprecated.
