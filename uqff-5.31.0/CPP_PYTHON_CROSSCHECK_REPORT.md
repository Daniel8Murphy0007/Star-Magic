# C++ vs Python cross-check report

**C++ source:** `uqff_exact_closures.cpp`
**Tolerance:** relative=1e-06, absolute floor=1e-12
**Total C++ functions:** 630
**Captured numeric values:** 630

## Summary

| Category | Count | Meaning |
|---|---:|---|
| MATCH      |   303 | Python and C++ agree within tolerance |
| DRIFT      |    15 | Both return scalars but they differ |
| MISSING    |   312 | C++ function has no matching Python closure key |
| UNCALLABLE |     0 | Python closure exists but doesn't return a comparable scalar |

## DRIFT — 15 value mismatches

These need investigation. C++ and Python disagree for the same name.

| Function | C++ value | Python value | Relative error |
|---|---|---|---:|
| `hubble_bubble` | -0.30145 | -30 | 9.900e-01 |
| `F_TRZ_identity` | 0.1 | 6 | 9.833e-01 |
| `hubble_tension` | 5.6 | 67.4 | 9.169e-01 |
| `pulsar_glitch` | 3.26419e-07 | 1e-06 | 6.736e-01 |
| `goldbach_weak` | 1 | 3 | 6.667e-01 |
| `crab_TeV_cutoff` | 79.261 | 80 | 9.238e-03 |
| `transcendental_pi_over_4` | 0.779167 | 0.785398 | 7.934e-03 |
| `transcendental_zeta_2` | 1.64733 | 1.64493 | 1.459e-03 |
| `transcendental_e` | 2.72083 | 2.71828 | 9.386e-04 |
| `transcendental_e_squared` | 7.39583 | 7.38906 | 9.172e-04 |
| `inflaton_n_s` | 0.964681 | 0.9655 | 8.484e-04 |
| `transcendental_catalan_g` | 0.916667 | 0.915966 | 7.654e-04 |
| `transcendental_pi_squared` | 9.87083 | 9.8696 | 1.245e-04 |
| `transcendental_ln_10` | 2.30267 | 2.30259 | 3.543e-05 |
| `transcendental_ln_2` | 0.693167 | 0.693147 | 2.811e-05 |

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

## Verdict

Of 318 comparable entries, **95.3% match** within relative tolerance 1e-06.

DRIFT entries are the high-priority audit targets. Each represents a case where the C++ value encoded does not match the current Python closure value. Likely causes: (a) the Python closure formula changed after the C++ port was authored; (b) the C++ value was a target rather than a derived value; (c) a typo. Verify each manually and either update the C++ value or flag the closure as deprecated.
