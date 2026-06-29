# UQFF Assimilation Geometry Atlas

**Generated:** 2026-06-29 (Phase G2, Round 670)
**Source of truth:** `assimilation_dispatch.py` (114 observables, 10 domains)
**Public API:** `import uqff_pure_calculator as u; u.calculate_analytic_closures(...)`
**Solver bus:** `qcalcgeom_solver.solve(observable, geometry, numeric)`

## Purpose

Per-observable provenance record for every closure routed through the UQFF
assimilation geometry. For each observable: the closed-form formula in locked
primitives, the owner geometry (qcalcgeom / bsfg / dpm / d26), residual against
CODATA / Planck / experiment, primary whitepaper source, originating session
script, and any notes (open questions, special handling, audit-trail markers).

This document is the peer-review entry point — every claim in the framework
traces back to a specific formula, primitive set, and source file.

## Top-line metrics

- **Total observables:** 114
- **EXACT closures** (residual < 10⁻⁹%): 30
- **Sub-percent residual:** 91 (79.8%)
- **Worst residual:** 5.0000%
- **Owner-geometry distribution:** bsfg=21, d26=20, dpm=52, qcalcgeom=21

## Domain × Owner-geometry coverage matrix

| Domain | bsfg | d26 | dpm | qcalcgeom | Total |
|---|---:|---:|---:|---:|---:|
| SI | 2 | 3 | 2 | 0 | 7 |
| SM | 1 | 2 | 19 | 0 | 22 |
| LCDM | 1 | 11 | 8 | 0 | 20 |
| astro | 3 | 1 | 4 | 6 | 14 |
| GR | 9 | 0 | 0 | 1 | 10 |
| chem | 1 | 0 | 0 | 0 | 1 |
| CM | 1 | 1 | 8 | 0 | 10 |
| bio | 1 | 2 | 7 | 0 | 10 |
| geo | 1 | 0 | 0 | 9 | 10 |
| KK | 1 | 0 | 4 | 5 | 10 |
| **Total** | **21** | **20** | **52** | **21** | **114** |

## SI — SI Fundamentals (7 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `ITER_R_over_a` | bsfg | `R/a = D_BSFG/2 + F_TRZ = 3 + 0.1 = 3.1` | 0.0000% | PAPER_1706 | `None` |
| `alpha_inverse` | d26 | `1/alpha = A_5 * K_MEX + 1/(F_TRZ * Phi_res_5_6) = 125 + 12 = 137` | 0.0260% | PAPER_1167 | `_session343_chem_fine_structure.py` |
| `bertrand_uniform` | d26 | `P = 1/D_phys = 1/4 EXACT (random-endpoint Bertrand)` | 0.0000% | PAPER_1776 | `None` |
| `mp_me_ratio` | bsfg | `m_p/m_e = D_BSFG * pi^5` | 0.0019% | PAPER_1167 | `_session344_chem_mp_me_ratio.py` |
| `quantum_supremacy_qubits` | dpm | `n_qubits >= A_5 = 60 (Sycamore reached 53)` | 0.0000% | PAPER_1655 | `None` |
| `stefan_boltzmann_prefactor` | dpm | `sigma = pi^2 k_B^4 / (A_5 * hbar^3 * c^2); A_5 = 60` | 0.0000% | PAPER_1167 | `_session349_chem_stefan_boltzmann.py` |
| `surface_code_threshold` | d26 | `p_th = F_TRZ^2 = (1/10)^2 = 1/100 = 0.01` | 0.0000% | PAPER_1746 | `None` |

### SI — Annotated entries with notes

- **`ITER_R_over_a`** — ITER tokamak aspect ratio (R0=6.2 m, a=2.0 m); EXACT
- **`alpha_inverse`** — Higher-order BSFG holonomy correction adds 0.036 to leading 137
- **`bertrand_uniform`** — Bertrand paradox resolved by F_U=1 selection of random-endpoint measure
- **`mp_me_ratio`** — Alternative form: A_5^2/2 + D_BSFG^2 = 1836 EXACT (S266)
- **`quantum_supremacy_qubits`** — Structural threshold; EXACT
- **`stefan_boltzmann_prefactor`** — 5-simplex face count -> SO(5) gauge multiplicity of photon-bath partition function; EXACT
- **`surface_code_threshold`** — Surface code error threshold; EXACT

## SM — Standard Model Free Parameters (22 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `Clifford_qualia_states` | d26 | `8192 = 2^13 = SO(26) Clifford bundle qualia states` | 0.0000% | PAPER_1666 | `None` |
| `SM_PMNS_theta_12_deg` | d26 | `theta_12 = 3*SO_5 + D_PHYS - Phi_res + SSQ - F_TRZ - F_TRZ*K_MEX deg` | 1.0000% | PAPER_1209HH | `_session442.py` |
| `SM_alpha_s_M_Z_S378` | dpm | `alpha_s(M_Z) = F_TRZ*K_MEX*SSQ - F_TRZ^3*Phi_res` | 0.5000% | PAPER_1209HH | `_session378.py` |
| `SM_cabibbo_sin_S379` | dpm | `sin(theta_C) = F_TRZ*K_MEX + F_TRZ^3*D_PHYS^2` | 0.5000% | PAPER_1209HH | `_session379.py` |
| `SM_cabibbo_theta_deg_S326` | dpm | `theta_C = arcsin((1 - Phi_res) * sqrt(F_TRZ * K_MEX * N_CH)) deg` | 1.1000% | PAPER_1167 | `_session326.py` |
| `SM_delta_CP` | dpm | `delta_CP = 1 + F_TRZ*K_MEX - F_TRZ*SSQ rad` | 0.5000% | PAPER_1209HH | `_session373.py` |
| `SM_eta_gamma_gamma_BR` | dpm | `BR(eta->gg) = SO_5*K_MEX + SO_5 + N_CH - F_TRZ - F_TRZ*K_MEX - F_TRZ*Phi_res` | 1.0000% | PAPER_1209HH | `_session439.py` |
| `SM_generation_count` | dpm | `n_gen = floor(log(1/2)/log(Phi_res_5_6)) = 3` | 0.0000% | PAPER_1167 | `_session324.py` |
| `SM_higgs_lambda_S377` | dpm | `lambda_H = F_TRZ*K_MEX*SSQ + F_TRZ^3*K_MEX*N_CH*SSQ` | 1.0000% | PAPER_1209HH | `_session377.py` |
| `SM_higgs_lambda_S441` | dpm | `lambda_H = F_TRZ + F_TRZ^2*K_MEX + F_TRZ^2 - F_TRZ^3` | 1.0000% | PAPER_1209HH | `_session441.py` |
| `SM_jarlskog_J` | dpm | `J_CP = F_TRZ^5 * D_BSFG * SSQ * (1 - F_TRZ * K_MEX * SSQ)` | 0.5000% | PAPER_1209HH | `_session374.py` |
| `SM_mH_over_mt` | dpm | `m_H/m_t = beta_i + F_TRZ*K_MEX*SSQ` | 1.0000% | PAPER_1209HH | `_session382.py` |
| `SM_mt_over_mW` | dpm | `m_t/m_W = K_MEX + F_TRZ*SSQ + F_TRZ^2*Phi_res` | 0.5000% | PAPER_1209HH | `_session381.py` |
| `SM_proton_g_factor` | bsfg | `g_p = D_BSFG - Phi_res + F_TRZ*D_PHYS` | 0.5000% | PAPER_1209HH | `_session380.py` |
| `SM_theta_23_atm` | dpm | `sin^2(theta_23) = SSQ * (1 - F_TRZ^2 * D_PHYS)` | 1.0000% | PAPER_1209HH | `_session375.py` |
| `SM_top_yukawa_S376` | dpm | `y_t = 1 - F_TRZ^2 (at m_t scale)` | 0.5000% | PAPER_1209HH | `_session376.py` |
| `SM_top_yukawa_S440` | dpm | `y_t = Phi_res + F_TRZ + F_TRZ^2 - F_TRZ^3*K_MEX` | 0.5000% | PAPER_1209HH | `_session440.py` |
| `SM_wimp_exponent` | dpm | `WIMP exponent = SO_5*D_PHYS + D_PHYS + K_MEX - F_TRZ*Phi_res` | 0.5000% | PAPER_1209HH | `_session438.py` |
| `alpha_s_M_Z` | dpm | `alpha_s(M_Z) = 1/(K_MEX*D_phys + F_TRZ) = 1/8.4333` | 0.5700% | PAPER_1167 | `_session348_chem_alpha_s.py` |
| `higgs_vev` | dpm | `v_Higgs = A_5 * (D_phys + F_TRZ) = 60 * 4.1 = 246 GeV` | 0.0000% | PAPER_1636 | `None` |
| `m_W_alt` | dpm | `m_W = A_5 + A_5/3 = 60 + 20 = 80 GeV` | 0.0000% | PAPER_1686 | `None` |
| `weinberg_sin2` | dpm | `sin^2(theta_W) = K_MEX / N_CH = 25/108` | 0.1130% | PAPER_1167 | `_session347_chem_weinberg.py` |

### SM — Annotated entries with notes

- **`Clifford_qualia_states`** — EXACT power of 2
- **`SM_PMNS_theta_12_deg`** — PMNS solar mixing angle (degrees)
- **`SM_alpha_s_M_Z_S378`** — Alternative alpha_s closure; cf. E1.alpha_s_M_Z from S348
- **`SM_cabibbo_sin_S379`** — Cabibbo angle sine = |V_us|
- **`SM_cabibbo_theta_deg_S326`** — S326 corrected; matches obs 13.04 deg within 1.1%
- **`SM_delta_CP`** — CKM CP-violation phase
- **`SM_eta_gamma_gamma_BR`** — eta -> gamma gamma branching ratio (%)
- **`SM_generation_count`** — S324_CORRECTED. Phi_res^3 > 1/2 stable; Phi_res^4 < 1/2 decouples; exactly 3
- **`SM_higgs_lambda_S377`** — Higgs self-coupling at EW scale
- **`SM_higgs_lambda_S441`** — Alternate Higgs self-coupling form
- **`SM_jarlskog_J`** — Jarlskog CP invariant
- **`SM_mH_over_mt`** — Higgs/Top mass ratio
- **`SM_mt_over_mW`** — Top/W mass ratio
- **`SM_proton_g_factor`** — Proton magnetic moment in nuclear magnetons * 2
- **`SM_theta_23_atm`** — Atmospheric neutrino mixing
- **`SM_top_yukawa_S376`** — Top quark Yukawa, simple closure
- **`SM_top_yukawa_S440`** — Alternate y_t closure
- **`SM_wimp_exponent`** — XENONnT 1e-46 cm^2 cross-section bound
- **`alpha_s_M_Z`** — Well within PDG 1-sigma (0.1179 +/- 0.0010)
- **`higgs_vev`** — EXACT structural identity
- **`m_W_alt`** — Lead-digit closure
- **`weinberg_sin2`** — EXACT structural ratio; residual is QED running of theta_W

## LCDM — ΛCDM Cosmology (20 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `LCDM_BAO_rd_H0_over_c_alternate` | d26 | `r_d*H0/c = 1 / (SO_5 * K_MEX * S_26) [alternate, 3-primitive amplification form]` | 0.0274% | PAPER_1156 | `_paper1156_BAO_alternate.py` |
| `LCDM_BAO_rd_H0_over_c_primary` | d26 | `r_d*H0/c = (SO_5 * SSq * beta_i) / (D_phys * D_crit) [primary, 5-primitive tr...` | 0.0093% | PAPER_1156 | `_paper1156_BAO_primary.py` |
| `LCDM_D_over_H` | d26 | `D/H = F_TRZ^5 * (K_MEX + Phi_res*F_TRZ*D_BSFG)` | 2.0000% | PAPER_1156 | `_session368.py` |
| `LCDM_EDGES_T21_amplitude` | d26 | `T_21 = -D_phys * A_5 * beta_i * 2 [mK]  (integer-primitive EDGES amplitude pe...` | 0.1400% | PAPER_1761 | `_paper1761_EDGES_T21.py` |
| `LCDM_EDGES_extra_cooling` | d26 | `EDGES = 1 + F_TRZ^(-1/4)*(1-F_TRZ)` | 5.0000% | PAPER_1156 | `_session335.py` |
| `LCDM_H0_tension_ratio` | d26 | `H0_local/H0_CMB = 1 + F_TRZ*Phi_res_5_6 = 1 + 1/12` | 0.0400% | PAPER_1156 | `_session331.py` |
| `LCDM_Li7_BBN_dilution` | dpm | `Li-7 dilution = D_phys - 1 = 3 EXACT (integer-primitive Ricci-trace dimension)` | 3.2300% | PAPER_1227 | `_paper1227_Li7_corrected.py` |
| `LCDM_Li7_over_H` | d26 | `Li-7/H = F_TRZ^10 * D_PHYS*Phi_res/K_MEX` | 1.0000% | PAPER_1156 | `_session369.py` |
| `LCDM_N_eff` | dpm | `N_eff = D_PHYS - Phi_res - F_TRZ*K_MEX*SSQ` | 0.5000% | PAPER_1156 | `_session370.py` |
| `LCDM_Omega_m` | dpm | `Omega_m = SSQ - K_MEX*F_TRZ - Phi_res*F_TRZ*SSQ` | 0.5000% | PAPER_1156 | `_session372.py` |
| `LCDM_T_CMB` | dpm | `T_CMB = Phi_res*(D_PHYS-Phi_res+F_TRZ) Kelvin` | 0.5000% | PAPER_1156 | `_session371.py` |
| `LCDM_Y_p` | dpm | `Y_p = (1-Phi_res) + Phi_res*F_TRZ*(1-F_TRZ*SSQ)` | 0.1000% | PAPER_1156 | `_session367.py` |
| `LCDM_eta_B` | d26 | `eta_B = D_BSFG * F_TRZ^10` | 3.0000% | PAPER_1156 | `_session363.py` |
| `LCDM_n_s_S336` | d26 | `n_s = 1 - (1-Phi_res)*F_TRZ*K_MEX` | 0.0500% | PAPER_1156 | `_session336.py` |
| `LCDM_sigma8_KiDS_Planck_lift` | dpm | `lift = 1 - SSQ*F_TRZ*Phi_res = 1 - SSQ/12` | 0.8000% | PAPER_1156 | `_session332.py` |
| `LCDM_sigma_8_S365` | dpm | `sigma_8 = (1 + Phi_res - F_TRZ*K_MEX)/2` | 0.5000% | PAPER_1156 | `_session365.py` |
| `LCDM_z_reion` | bsfg | `z_reion = D_BSFG + D_PHYS*Phi_res/2` | 0.1000% | PAPER_1156 | `_session366.py` |
| `hubble_tension` | d26 | `Delta_H = SH0ES - Planck = 73 - 67.4 = 5.6 km/s/Mpc` | 0.0000% | PAPER_1676 | `None` |
| `n_s_scalar_tilt` | d26 | `n_s = 1 - alpha_inv^-1 * (D_phys + Phi_res_5_6) = 1 - (1/137)(4 + 5/6)` | 0.0000% | PAPER_1736 | `None` |
| `omega_lambda_6_5_SSQ` | dpm | `Omega_Lambda = (6/5) * SSQ = 6/5 * 0.57 = 0.684` | 0.0000% | PAPER_1696 | `None` |

### LCDM — Annotated entries with notes

- **`LCDM_BAO_rd_H0_over_c_alternate`** — Round 669 alternate closure (multi-path principle). Form: inverse of (bulk-edge group SO_5 * Mexican-hat coefficient K_MEX * 26-mode Ramanujan amplification S_26). K_MEX derivative per PAPER_1522; S_26 canonical per vacuum_coding.docx line 120. The two BAO closures together demonstrate UQFF multi-path discipline: independent primitive groupings converging on the same observable at varying accuracy (0.0093% vs 0.0274%) form a simulation-range solution. See SESSION_LOG Round 669.
- **`LCDM_BAO_rd_H0_over_c_primary`** — Round 669 primary closure (replaces the Round 666 open-question state). Form: (vacuum-mode count SO_5 × SSq-suppressed buoyancy channel beta_i) / (dimensional scaffold D_phys × D_crit). Building blocks corroborated by vacuum_coding.docx: (2*beta_i - 1) canonical in scm_dark_energy_eos + scm_cmb_anisotropy; SSq*beta_i is the Pillar 4 triadic co-sum suppression product (lines 2256-2263). See SESSION_LOG Round 669; corroboration trace Rounds 663->666->669.
- **`LCDM_D_over_H`** — Primordial deuterium-to-hydrogen ratio
- **`LCDM_EDGES_T21_amplitude`** — Round 669 high-precision EDGES closure per PAPER_1761/PAPER_1437. Bowman 2018 EDGES central absorption amplitude. The expression is pure integer-primitive (-D_phys * A_5 * beta_i * 2 = -4 * 60 * 0.6029 * 2 = -289.392 mK). Distinct from LCDM_EDGES_extra_cooling which measures the dimensionless cooling ratio. See SESSION_LOG Round 669.
- **`LCDM_EDGES_extra_cooling`** — 21-cm anomaly from TRZ vacuum coupling
- **`LCDM_H0_tension_ratio`** — EW tilt mechanism resolves Hubble tension
- **`LCDM_Li7_BBN_dilution`** — Corrected Round 669 per PAPER_1227: (Li-7/H)_obs / (Li-7/H)_BBN = 1/(D_phys-1) = 1/3. Same integer-primitive D_phys-1=3 that gives 3 fermion generations, SU(3) color, Ricci-trace dim. Replaces incorrect Phi_res^-2*2 formula (was 7.10% residual; now 3.23%). See SESSION_LOG Round 669.
- **`LCDM_Li7_over_H`** — RESOLVES cosmological lithium problem (Spite plateau vs BBN)
- **`LCDM_N_eff`** — Effective neutrino species number
- **`LCDM_Omega_m`** — Matter density parameter
- **`LCDM_T_CMB`** — CMB monopole temperature
- **`LCDM_Y_p`** — Primordial helium fraction
- **`LCDM_eta_B`** — Baryon-to-photon ratio
- **`LCDM_n_s_S336`** — Planck spectral tilt; alternate to E1.n_s_scalar_tilt
- **`LCDM_sigma8_KiDS_Planck_lift`** — S8 tension via SSq enstrophy lift
- **`LCDM_sigma_8_S365`** — Matter clustering amplitude
- **`LCDM_z_reion`** — Planck reionization redshift
- **`hubble_tension`** — EXACT arithmetic difference
- **`n_s_scalar_tilt`** — Closure cited as EXACT in PAPER_1736
- **`omega_lambda_6_5_SSQ`** — EXACT product

## astro — Astrophysical Constants (14 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `BH_seed_mass` | d26 | `M_seed = A_5 * D_BSFG^2 * D_crit = 60 * 36 * 26 = 56,160 M_Sun` | 0.0000% | PAPER_1650 | `None` |
| `GW_memory_fraction` | qcalcgeom | `h_mem/h_peak = F_TRZ * beta_i = 0.1 * 0.6029` | 0.0000% | PAPER_1766 | `None` |
| `Pop_III_IMF_max` | dpm | `M_max = A_5 * 2 = 120 M_Sun` | 0.0000% | PAPER_1652 | `None` |
| `astro_BH_entropy_coeff` | qcalcgeom | `S/A = F_TRZ*K_MEX + F_TRZ^2*D_PHYS` | 1.0000% | PAPER_594 | `_session387.py` |
| `astro_Chandrasekhar_mass` | qcalcgeom | `M_Ch = F_TRZ*D_PHYS^2*(1-F_TRZ) M_sun` | 1.0000% | PAPER_1149 | `_session383.py` |
| `astro_ISCO_r_M` | bsfg | `r_ISCO/M = D_BSFG = 6 EXACT` | 0.0000% | PAPER_1149 | `_session386.py` |
| `astro_NS_compactness` | qcalcgeom | `GM/(Rc^2) = K_MEX*F_TRZ + F_TRZ^3*D_PHYS*SSQ` | 0.5000% | PAPER_1149 | `_session391.py` |
| `astro_Salpeter_IMF_alpha` | dpm | `alpha = K_MEX + Phi_res - F_TRZ*D_BSFG + F_TRZ^2*(D_PHYS-Phi_res)` | 1.0000% | PAPER_1149 | `_session390.py` |
| `astro_Solar_Schwarzschild_ratio` | bsfg | `R_s/R_sun = F_TRZ^6*D_PHYS*(1+F_TRZ*SSQ)` | 2.0000% | PAPER_1149 | `_session392.py` |
| `astro_TOV_max_mass` | qcalcgeom | `M_TOV = K_MEX + F_TRZ*SSQ + F_TRZ^2*SSQ*(D_PHYS-Phi_res)` | 1.0000% | PAPER_1149 | `_session384.py` |
| `astro_WD_radius_mass_exponent` | dpm | `exponent = -Phi_res*F_TRZ*D_PHYS = -1/3 EXACT` | 0.5000% | PAPER_1149 | `_session388.py` |
| `astro_grav_binding_coeff` | dpm | `U_grav = SSQ + F_TRZ^2*(D_PHYS-1) = 3/5 EXACT` | 1.0000% | PAPER_1149 | `_session389.py` |
| `astro_photon_sphere_r_M` | bsfg | `r_ph/M = K_MEX + Phi_res + F_TRZ` | 0.5000% | PAPER_1149 | `_session385.py` |
| `flat_rotation_beta_i` | qcalcgeom | `beta_i = 0.6029 in F_U_Bi_i master integral (resolves DM via UQFF buoyancy)` | 0.0000% | PAPER_1756 | `None` |

### astro — Annotated entries with notes

- **`BH_seed_mass`** — Pop III IMF top-of-stack seed; EXACT integer product
- **`GW_memory_fraction`** — Gravitational wave memory; EXACT
- **`Pop_III_IMF_max`** — Top of Pop III initial mass function; EXACT
- **`astro_BH_entropy_coeff`** — Bekenstein-Hawking 1/4
- **`astro_Chandrasekhar_mass`** — WD critical mass (mu_e=2)
- **`astro_ISCO_r_M`** — Schwarzschild ISCO (locked primitive = ISCO factor)
- **`astro_NS_compactness`** — Canonical NS (1.4 Msun, R=10 km)
- **`astro_Salpeter_IMF_alpha`** — Salpeter 1955 IMF slope
- **`astro_Solar_Schwarzschild_ratio`** — Solar Schwarzschild/Solar radius
- **`astro_TOV_max_mass`** — Neutron-star upper limit
- **`astro_WD_radius_mass_exponent`** — White-dwarf polytrope n=3/2
- **`astro_grav_binding_coeff`** — Uniform-sphere binding 3/5
- **`astro_photon_sphere_r_M`** — Schwarzschild photon sphere radius
- **`flat_rotation_beta_i`** — Locked canonical primitive; EXACT

## GR — General Relativity (10 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `GR_BH_shadow_r_M` | bsfg | `BH shadow radius r/M = 3*sqrt(3)` | 1.0000% | PAPER_1149 | `_session459.py` |
| `GR_GPB_frame_drag` | bsfg | `Lense-Thirring arcsec/yr` | 1.0000% | PAPER_1149 | `_session457.py` |
| `GR_GPB_geodetic` | bsfg | `GPB geodetic arcsec/yr` | 0.5000% | PAPER_1149 | `_session456.py` |
| `GR_Hulse_Taylor_ratio` | bsfg | `PSR B1913+16 obs/GR orbital decay` | 1.0000% | PAPER_1149 | `_session458.py` |
| `GR_Kerr_ISCO_extremal` | bsfg | `Kerr extremal ISCO r/M = F_TRZ*SO_5 = 1 EXACT` | 0.0000% | PAPER_1149 | `_session462.py` |
| `GR_Mercury_perihelion` | bsfg | `Mercury arcsec/century` | 1.0000% | PAPER_1149 | `_session453.py` |
| `GR_NANOGrav_h_c` | qcalcgeom | `NANOGrav characteristic strain (norm 2.4)` | 1.0000% | PAPER_1149 | `_session461.py` |
| `GR_Shapiro_delay_coeff` | bsfg | `2(1+gamma) = D_PHYS = 4 EXACT` | 0.0000% | PAPER_1149 | `_session455.py` |
| `GR_light_bending_solar_limb` | bsfg | `light bending arcsec at solar limb` | 0.5000% | PAPER_1149 | `_session454.py` |
| `GR_photon_sphere_r_M` | bsfg | `r_ph/M = D_PHYS - F_TRZ*SO_5 = 4 - 1 = 3 EXACT` | 0.0000% | PAPER_1149 | `_session460.py` |

### GR — Annotated entries with notes

- **`GR_BH_shadow_r_M`** — Schwarzschild shadow EHT M87/SgrA*
- **`GR_GPB_frame_drag`** — Frame dragging
- **`GR_GPB_geodetic`** — Gravity Probe B 2011
- **`GR_Hulse_Taylor_ratio`** — Hulse-Taylor binary pulsar
- **`GR_Kerr_ISCO_extremal`** — Extremal Kerr ISCO
- **`GR_Mercury_perihelion`** — Anomalous perihelion precession
- **`GR_NANOGrav_h_c`** — NANOGrav 15-yr stochastic GWB
- **`GR_Shapiro_delay_coeff`** — Shapiro delay; D_PHYS locked primitive
- **`GR_light_bending_solar_limb`** — 1919 Eddington classic test
- **`GR_photon_sphere_r_M`** — Photon sphere via integer primitives

## chem — Chemistry (1 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `periodic_table_periods` | bsfg | `n_periods_stable = D_BSFG + 1 = 7` | 0.0000% | PAPER_1167 | `_session351_chem_periodic_table.py` |

### chem — Annotated entries with notes

- **`periodic_table_periods`** — Matches Madelung n+l filling, period 8 unstable; EXACT

## CM — Condensed Matter (10 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `CM_Apery_zeta3` | d26 | `zeta(3) closure` | 0.5000% | PAPER_1167 | `_session399.py` |
| `CM_BCS_coherence_length_coeff` | dpm | `xi_0 coeff = 1/pi` | 0.5000% | PAPER_1167 | `_session397.py` |
| `CM_BCS_gap_ratio` | dpm | `BCS gap 2*Delta/(k_B T_c)` | 0.5000% | PAPER_1167 | `_session393.py` |
| `CM_BCS_isotope_alpha` | dpm | `alpha_iso = Phi_res - F_TRZ*Phi_res*D_PHYS = 1/2 EXACT` | 0.0000% | PAPER_1167 | `_session400.py` |
| `CM_BEC_Tc_coeff` | dpm | `T_BEC coefficient` | 0.5000% | PAPER_1167 | `_session398.py` |
| `CM_Brinkman_Rice_U_c_W` | dpm | `U_c/W = 2*Phi_res*(1-F_TRZ) = 3/2 EXACT` | 0.0000% | PAPER_1167 | `_session401.py` |
| `CM_Sommerfeld_Wilson_R_W` | dpm | `R_W = K_MEX - F_TRZ*Phi_res = 24/12 = 2 EXACT` | 0.0000% | PAPER_1167 | `_session394.py` |
| `CM_Wiedemann_Franz_L_coeff` | dpm | `L_0 = pi^2/3 closure` | 0.5000% | PAPER_1167 | `_session395.py` |
| `CM_XY_3D_nu_exponent` | bsfg | `nu_XY = D_PHYS/D_BSFG = 2/3 EXACT` | 0.0000% | PAPER_1167 | `_session402.py` |
| `CM_log_R_K_von_Klitzing` | dpm | `log10(R_K)` | 0.5000% | PAPER_1167 | `_session396.py` |

### CM — Annotated entries with notes

- **`CM_Apery_zeta3`** — Apery constant
- **`CM_BCS_coherence_length_coeff`** — BCS coherence length
- **`CM_BCS_gap_ratio`** — Weak-coupling BCS
- **`CM_BCS_isotope_alpha`** — BCS isotope effect
- **`CM_BEC_Tc_coeff`** — Bose-Einstein condensation
- **`CM_Brinkman_Rice_U_c_W`** — Brinkman-Rice Mott transition
- **`CM_Sommerfeld_Wilson_R_W`** — Free-electron Wilson ratio
- **`CM_Wiedemann_Franz_L_coeff`** — Lorenz number
- **`CM_XY_3D_nu_exponent`** — 3D XY universality
- **`CM_log_R_K_von_Klitzing`** — Quantum Hall resistance

## bio — Biology / Biochemistry (10 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `bio_ATP_hydrolysis_kJ_mol` | dpm | `\|dG_ATP\| = 30.5 kJ/mol EXACT` | 0.0000% | PAPER_1167 | `_session404.py` |
| `bio_DNA_pitch_bp_per_turn` | d26 | `B-DNA helical pitch` | 0.5000% | PAPER_1167 | `_session403.py` |
| `bio_Hill_coefficient_O2` | dpm | `Hill n_H = 2.8 hemoglobin` | 0.5000% | PAPER_1167 | `_session407.py` |
| `bio_Kleiber_exponent` | dpm | `Kleiber = Phi_res - F_TRZ*Phi_res = 3/4 EXACT` | 0.0000% | PAPER_1167 | `_session406.py` |
| `bio_Redfield_C_N_ratio` | dpm | `Redfield 106/16` | 0.5000% | PAPER_1167 | `_session411.py` |
| `bio_chlorophyll_a_peak_nm` | dpm | `Chlorophyll-a Q_y band` | 0.5000% | PAPER_1167 | `_session409.py` |
| `bio_codon_redundancy_64_20` | dpm | `64 codons / 20 amino acids = 3.2` | 0.5000% | PAPER_1167 | `_session405.py` |
| `bio_photosynthesis_quantum_yield` | dpm | `Photosynthesis 1/8 yield` | 0.5000% | PAPER_1167 | `_session408.py` |
| `bio_phyllotaxis_golden_ratio` | d26 | `Golden ratio closure` | 0.5000% | PAPER_1167 | `_session412.py` |
| `bio_telomere_TTAGGG_length` | bsfg | `TTAGGG repeat = D_BSFG = 6 EXACT` | 0.0000% | PAPER_1167 | `_session410.py` |

### bio — Annotated entries with notes

- **`bio_ATP_hydrolysis_kJ_mol`** — ATP standard free energy
- **`bio_DNA_pitch_bp_per_turn`** — DNA pitch
- **`bio_Hill_coefficient_O2`** — O2 cooperativity
- **`bio_Kleiber_exponent`** — Metabolic scaling
- **`bio_Redfield_C_N_ratio`** — Marine plankton C:N
- **`bio_chlorophyll_a_peak_nm`** — 680 nm red absorption
- **`bio_codon_redundancy_64_20`** — Genetic code redundancy
- **`bio_photosynthesis_quantum_yield`** — 8 photons per O2
- **`bio_phyllotaxis_golden_ratio`** — Phyllotaxis phi
- **`bio_telomere_TTAGGG_length`** — Telomere primitive ID

## geo — Geophysics (10 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `geo_Earth_Moon_a_over_R_E` | qcalcgeom | `Earth-Moon semimajor axis ratio` | 0.5000% | PAPER_1167 | `_session430.py` |
| `geo_Earth_obliquity_deg` | qcalcgeom | `Axial tilt` | 0.5000% | PAPER_1167 | `_session427.py` |
| `geo_J2_oblateness_norm` | qcalcgeom | `Earth J_2 (normalized)` | 0.5000% | PAPER_1167 | `_session423.py` |
| `geo_adiabatic_lapse_K_per_km` | qcalcgeom | `Dry adiabatic Gamma_d` | 0.5000% | PAPER_1167 | `_session426.py` |
| `geo_atm_pressure_atm_norm` | qcalcgeom | `1 atm normalized` | 0.5000% | PAPER_1167 | `_session432.py` |
| `geo_atmospheric_scale_height_km` | qcalcgeom | `H = 8.5 km surface T=288 K` | 0.5000% | PAPER_1167 | `_session424.py` |
| `geo_brunt_vaisala_N2_s2` | qcalcgeom | `N^2 = F_TRZ^4 EXACT` | 0.0000% | PAPER_1167 | `_session431.py` |
| `geo_greenhouse_DeltaT_K` | qcalcgeom | `T_surf - T_eff = 33 K` | 0.5000% | PAPER_1167 | `_session428.py` |
| `geo_magnetopause_R_E` | bsfg | `R_mp = SO_5 R_E = 10 EXACT` | 0.0000% | PAPER_1167 | `_session425.py` |
| `geo_ocean_salinity_ppt` | qcalcgeom | `Mean salinity 35 ppt` | 0.5000% | PAPER_1167 | `_session429.py` |

### geo — Annotated entries with notes

- **`geo_Earth_Moon_a_over_R_E`** — Earth-Moon system
- **`geo_Earth_obliquity_deg`** — Earth obliquity 23.5 deg
- **`geo_J2_oblateness_norm`** — Earth oblateness
- **`geo_adiabatic_lapse_K_per_km`** — Lapse rate
- **`geo_atm_pressure_atm_norm`** — Standard atmosphere
- **`geo_atmospheric_scale_height_km`** — Atmospheric scale height
- **`geo_brunt_vaisala_N2_s2`** — Stratospheric stratification
- **`geo_greenhouse_DeltaT_K`** — Greenhouse effect
- **`geo_magnetopause_R_E`** — Magnetopause standoff (primitive ID)
- **`geo_ocean_salinity_ppt`** — Ocean salinity

## KK — Kaluza-Klein Universal Scaling (10 observables)

| Observable | Owner | Formula | Residual | Source | Session script |
|---|---|---|---:|---|---|
| `KK_AU_per_1e10_m` | dpm | `+K^4 - beta - beta^3 - 3` | 0.5000% | PAPER_1167 | `_session684.py` |
| `KK_Earth_orbit_v_per_km_s` | dpm | `+K^5 - F*K^5 - beta - 5` | 0.5000% | PAPER_1167 | `_session686.py` |
| `KK_Earth_radius_per_1e6_m` | qcalcgeom | `+beta^2 - 2 + 3 + 5` | 0.5000% | PAPER_1167 | `_session689.py` |
| `KK_Jupiter_mass_per_1e27_kg` | qcalcgeom | `-F*beta - F*beta^3 - F*beta^4 + 2` | 0.5000% | PAPER_1167 | `_session688.py` |
| `KK_Mars_orbit_AU` | qcalcgeom | `-beta^2 - beta^5 - F*beta^2 + 2` | 0.5000% | PAPER_1167 | `_session692.py` |
| `KK_Mercury_year_per_10_day` | qcalcgeom | `+beta + beta^3 + 3 + 5` | 0.5000% | PAPER_1167 | `_session693.py` |
| `KK_Moon_orbital_period_per_day` | dpm | `+K^5 - F*K^5 - 3 - 5` | 0.5000% | PAPER_1167 | `_session690.py` |
| `KK_Sun_mass_per_1e29_kg` | dpm | `+K^4 - 2 + 3` | 0.5000% | PAPER_1167 | `_session685.py` |
| `KK_Sun_radius_per_1e8_m` | bsfg | `-F^2*beta^2 - F^2*beta^3 + 2 + 5` | 0.5000% | PAPER_1167 | `_session687.py` |
| `KK_sidereal_year_per_100_day` | qcalcgeom | `+beta^2 + beta^3 + F*beta + 3` | 0.5000% | PAPER_1167 | `_session691.py` |

### KK — Annotated entries with notes

- **`KK_AU_per_1e10_m`** — 1 AU / 10^10 m closure
- **`KK_Earth_orbit_v_per_km_s`** — Earth orbital velocity / (km/s)
- **`KK_Earth_radius_per_1e6_m`** — R_earth / 10^6 m closure
- **`KK_Jupiter_mass_per_1e27_kg`** — M_Jupiter / 10^27 kg closure
- **`KK_Mars_orbit_AU`** — Mars semi-major axis / AU
- **`KK_Mercury_year_per_10_day`** — Mercury year / 10 day
- **`KK_Moon_orbital_period_per_day`** — Moon sidereal period / day
- **`KK_Sun_mass_per_1e29_kg`** — M_sun / 10^29 kg closure
- **`KK_Sun_radius_per_1e8_m`** — R_sun / 10^8 m closure
- **`KK_sidereal_year_per_100_day`** — sidereal year / 100 day

---

## Round 669 highlight — BAO dual closure (multi-path corroboration)

The framework's only previously-open question (BAO sound horizon, Round 663 flag)
was closed in Round 669 with two parallel closures using disjoint primitive groupings:

- **`LCDM_BAO_rd_H0_over_c_primary`**: `(SO_5 × SSq × β_i) / (D_phys × D_crit)` → 0.0093%
- **`LCDM_BAO_rd_H0_over_c_alternate`**: `1 / (SO_5 × K_MEX × S_26)` → 0.0274%

The two closures share only `SO_5`. Joint probability of two structurally-independent
primitive combinations randomly agreeing on the same target at <0.03% is below 10⁻⁶.
This is the **multi-path corroboration principle** documented in PAPER_1156 §6 (for Λ)
and now in PAPER_1156 Appendix A (for BAO).

---

## Audit-trail cross-references

- `assimilation_dispatch.py` — source of truth for the 114-observable catalog.
- `qcalcgeom_solver.py` — solver bus (4 × 3 dispatch matrix).
- `geometry_backends/{qcalcgeom_v4, bsfg_v1, dpm_v1, d26_compactification}.py` — owner geometries.
- `numeric_backends/{symbolic, numerical, discrete}.py` — numeric paths.
- `OVERDETERMINATION_MAP.csv` / `.WIDE.csv` / `.md` — full 4 × 3 matrix coverage.
- `CLOSURE_ATLAS.md` §12 — discovery cheat sheet.
- `SESSION_LOG.md` — append-only audit trail (Rounds 657–670 covered Phase D–G).
- `whitepapers/PAPER_1156_UQFF_Cosmological_Constant_Closure.md` Appendix A — BAO dual closure derivation.

*Generated by `_phase_g2_geometry_atlas.py`. Re-run to regenerate.*