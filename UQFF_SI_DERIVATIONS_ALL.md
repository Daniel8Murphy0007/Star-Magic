# UQFF SI Unit Derivations — All Closures

Comprehensive compilation of every SI / fundamental-constant derivation in the UQFF calculator.

---

## 1. Core SI Closures (`_si_derivation_report`)

Output from `uqff_pure_calculator.py::_si_derivation_report()`:

| Quantity | UQFF closed form | UQFF derived | CODATA observed | Residual |
|---|---|---|---|---|
| c (speed of light, m/s) | `(D_crit * 4*pi / Phi_res) * v_F` | 2.994985e+08 | 2.99792458e+08 | 0.0981% |
| G (microscopic, m³/kg·s²) | `(2*pi * D_crit^3 * Phi_res / ([SSq]^3 * (26!)^2)) * v_F^5/(E_0 * f_THz)` | 6.668992e-11 | 6.6743e-11 | 0.0795% |
| G (cosmic, m³/kg·s²) | (same closed form, cosmic scale) | 6.686636e-11 | 6.6743e-11 | 0.1848% |

**Anchors**: v_F = 770000.0 m/s, E_0 = 1e-20 J, f_THz = 1250000000000.0 Hz

---

## 2. Other SI constants in the calculator

| Constant | CODATA observed | UQFF derived | Residual | Calculator function |
|---|---|---|---|---|
| h (Planck, J·s) | 6.626070e-34 | (function not found by this name) | n/a | `_l96_uqff_axiom_h_planck_closure` |
| hbar (reduced Planck) | 1.054572e-34 | (function not found by this name) | n/a | `_l96_uqff_hbar_derived` |
| k_B (Boltzmann, J/K) | 1.380649e-23 | 1.379596e-23 | 0.0763% | `_l96_uqff_k_boltzmann_derived` |
| k_B (alt) | 1.380649e-23 | 1.379596e-23 | 0.0763% | `_l96_uqff_k_boltzmann_closure` |
| e (elementary charge, C) | 1.602177e-19 | (function not found by this name) | n/a | `_l96_uqff_axiom_elementary_charge_closure` |
| N_A (Avogadro) | 6.022141e+23 | (function not found by this name) | n/a | `_l96_uqff_N_A_avogadro_derived` |
| N_A (alt) | 6.022141e+23 | (function not found by this name) | n/a | `_l96_uqff_N_A_avogadro_closure` |
| R (gas constant) | 8.314463e+00 | 8.282196e+00 | 0.3881% | `_l96_uqff_R_gas_derived` |
| alpha (fine-structure) | 7.297353e-03 | (function not found by this name) | n/a | `_l96_uqff_axiom_alpha_em_closure` |
| alpha^-1 | 1.370360e+02 | (function not found by this name) | n/a | `_l96_uqff_axiom_alpha_inverse_137_closure` |
| epsilon_0 | 8.854188e-12 | 8.871577e-12 | 0.1964% | `_l96_uqff_epsilon_0_derived` |
| epsilon_0 (alt) | 8.854188e-12 | 8.871577e-12 | 0.1964% | `_l96_uqff_epsilon_0_closure` |
| mu_0 | 1.256637e-06 | 1.256637e-06 | 0.0000% | `_l96_uqff_mu_0_derived` |
| mu_0 (alt) | 1.256637e-06 | 1.256637e-06 | 0.0000% | `_l96_uqff_mu_0_closure` |
| F (Faraday) | 9.648533e+04 | 9.613615e+04 | 0.3619% | `_l96_uqff_F_faraday_derived` |
| sigma_SB (Stefan-Boltz) | 5.670374e-08 | 5.674490e-08 | 0.0726% | `_l96_uqff_sigma_SB_derived` |
| Wien b (m·K) | 2.897772e-03 | 4.992818e+00 | 172198.5134% | `_l96_uqff_wien_constant_derived` |
| Rydberg R_inf (1/m) | 1.097373e+07 | 1.093945e+07 | 0.3124% | `_l96_uqff_rydberg_infty_derived` |
| Bohr radius (m) | 5.291772e-11 | (function not found by this name) | n/a | `_l96_uqff_axiom_bohr_radius_closure` |
| classical r_e (m) | 2.817940e-15 | (function not found by this name) | n/a | `_l96_uqff_axiom_classical_electron_radius_closure` |
| Compton lambda_C (m) | 2.426310e-12 | (function not found by this name) | n/a | `_l96_uqff_axiom_compton_wavelength_closure` |
| Hartree E_h (J) | 4.359745e-18 | (function not found by this name) | n/a | `_l96_uqff_axiom_hartree_energy_closure` |
| m_p/m_e ratio | 1.836153e+03 | (function not found by this name) | n/a | `_l96_uqff_axiom_m_p_over_m_e_closure` |
| Lambda_QCD (GeV) | 2.170000e-01 | 5.130000e-02 | 76.3594% | `_lambda_qcd_uqff_derive` |

---

## 3. SI Constants in `master_closures.csv`

Found **90 SI-related rows in master_closures.csv**:

| Closure label | Predicted | Observed | Residual | Status |
|---|---|---|---|---|
| `` | 1.8765 |  | 5.42 | OK |
| `` |  |  | 0.0886 | OK |
| `` |  |  | 0.0814 | OK |
| `` | 1.41e-12 |  | 0 | OK |
| `` | 12 |  | 0 | OK |
| `` | 1836.1181 | 1836.15267343 | 0.0019 | OK |
| `` | 5.5667 | 5.5857 | 0.341 | OK |
| `` | 0.1198 | 0.12 | 0.1389 | OK |
| `` | 5.6715 | 5.67 | 0.0265 | OK |
| `` | 6.0225 | 6.022 | 0.0083 | OK |
| `` | 10.9598 | 10.974 | 0.1291 | OK |
| `` | 5.3098 | 5.29 | 0.3749 | OK |
| `` | 27.2240 | 27.211 | 0.0478 | OK |
| `` | 2.4315 | 2.426 | 0.2267 | OK |
| `` | 6.6298 | 6.626 | 0.0579 | OK |
| `` | 1836.0000 | 1836.1527 | 0.0083 | OK |
| `` | 207.0000 | 206.7683 | 0.1121 | OK |
| `` | 3477.6000 | 3477.2300 | 0.0106 | OK |
| `` | 1.001100 | 1.001378 | 0.0278 | OK |
| `` | 938.2500 | 938.2721 | 0.0024 | OK |
| `` | 0.511167 | 0.510999 | 0.0328 | OK |
| `` | 333000 | 333000 | 0 | OK |
| `` | 317.783333 | 317.8 | 0.0052 | OK |
| `` | 6.022800 | 6.022 | 0.0133 | OK |
| `` | 6.627000 | 6.626 | 0.0151 | OK |
| `` | 8.856333 | 8.854 | 0.0264 | OK |
| `` | 1.255700 | 1.257 | 0.1034 | OK |
| `` | 5.290333 | 5.292 | 0.0315 | OK |
| `` | 1.097000 | 1.0974 | 0.0364 | OK |
| `` | 2.426333 | 2.426 | 0.0137 | OK |
| `` | 9.273333 | 9.274 | 0.0072 | OK |
| `` | 13.605700 | 13.6057 | 0 | OK |
| `` | 5.670000 | 5.6700 | 0 | OK |
| `` | 2.899667 | 2.8980 | 0.0575 | OK |
| `` | 96485.000000 | 96485 | 0 | OK |
| `` | 6.022416 | 6.0220 | 0.0069 | OK |
| `` | 1.380633 | 1.3810 | 0.0266 | OK |
| `` | 6.624300 | 6.6260 | 0.0257 | OK |
| `` | 4.360000 | 4.3600 | 0 | OK |
| `PAPER_008_UQFF_Waveform_Phase_Evolution_Template_M` |  |  |  | RESEARCH_TRA |

*... and 50 more SI-related rows.*

---

## 4. SI references in master grok corpus

See `UQFF_GROK_DERIVATIONS_INDEX.md` and grep `UQFF_GROK_LONG_FORM_DERIVATIONS_MASTER.md` for verbatim grok long-forms of each SI constant.
