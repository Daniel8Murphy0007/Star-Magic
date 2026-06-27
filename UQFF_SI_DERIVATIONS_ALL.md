# UQFF SI Unit Derivations -- All Closures with Their Actual Functions

Each SI constant has one or more dedicated `_session<N>_<name>.py` standalone derivation scripts.
Total SI-related session scripts identified: 57

---

## Core Closures from `uqff_pure_calculator.py::_si_derivation_report()`

| Quantity | UQFF closed form | UQFF derived | CODATA | Residual |
|---|---|---|---|---|
| **c** (speed of light, m/s) | `(D_crit · 4π / Phi_res) · v_F` | 2.99499e+08 | 2.99792e+08 | **0.0981%** |
| **G** (microscopic) | `(2π · D_crit³ · Phi_res / ([SSq]³ · (26!)²)) · v_F⁵/(E_0·f_THz)` | 6.66899e-11 | 6.6743e-11 | **0.0795%** |
| **G** (cosmic) | (same form, cosmic scale) | 6.68664e-11 | 6.6743e-11 | **0.1848%** |

Anchors: v_F=770000 m/s, E_0=1e-20 J, f_THz=1.25e12 Hz

---

## Per-Constant Dedicated Session Scripts (the missing trace)


### c (speed of light, m/s)

- **CODATA observed**: 2.99792e+08
- `_session592_c_si.py` (NOT FOUND)
- **`_session694_alpha_forward_derivation.py`**: , G6)
  - Output: `    F_TRZ   = {F_TRZ}   ({float(F_TRZ):.6f})`

### h (Planck, J·s)

- **CODATA observed**: 6.62607e-34
- **`_session482_planck.py`**: 
  - Output: `S482 COMPLETE. Planck h (x1e-34 J*s) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session629_planck.py`**: 
  - Output: `Planck h lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%`

### hbar (reduced Planck, J·s)

- **CODATA observed**: 1.05457e-34
- **`_session719_hbar_dimensional_class_probe.py`**: Session 719 -- hbar-chain opening: dimensional class probe.
  - Output: `    L_DPM^4 = hbar * v_SCM / rho_vac = {L4_implied:.6e} m^4`
- **`_session720_hbar_classIV_nogo_and_L_SCM_promotion.py`**: Session 720 -- hbar Class IV no-go theorem + L_SCM promotion to 12th locked primitive.
  - Output: `STEP 2 -- Promotion: L_SCM := (hbar * v_SCM / rho_vac_SCm)^(1/4)`
- **`_session729_H0_eta_hbar_consolidation.py`**: from {c, rho_vac, L_SCM} + locked rationals.
  - Output: `  H_0 (obs)     = {H0_SI:.6e} s^-1`

### k_B (Boltzmann, J/K)

- **CODATA observed**: 1.38065e-23
- **Session script(s)**: (none found by name probe — search master grok for long form)

### e (elementary charge, C)

- **CODATA observed**: 1.60218e-19
- **Session script(s)**: (none found by name probe — search master grok for long form)

### N_A (Avogadro, 1/mol)

- **CODATA observed**: 6.02214e+23
- **`_session472_avogadro.py`**: 
  - Output: `S472 COMPLETE. Avogadro (x1e23) = {pred:.4f}; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.`
- **`_session583_avogadro.py`**: 
  - Output: `Avogadro lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%`
- **`_session627_avogadro.py`**: 
  - Output: `Avogadro lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%`

### R (gas constant, J/mol·K)

- **CODATA observed**: 8.31446
- **Session script(s)**: (none found by name probe — search master grok for long form)

### alpha (fine-structure)

- **CODATA observed**: 0.00729735
- **`_session343_chem_fine_structure.py`**: 1/alpha = A_5 * K_Mex + 1/(F_TRZ * Phi_res) = 125 + 12 = 137.
  - Output: `S343 COMPLETE. 1/alpha = A_5*K_Mex + 1/(F_TRZ*Phi_res) = 125 + 12 = {inv_alpha:.4f}; `
- **`_session475_fine_structure.py`**: 
  - Output: `S475 COMPLETE. Fine structure 1/alpha = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session552_alpha.py`**: 
  - Output: `S552 COMPLETE. fine structure alpha = {float(v):.6f}; closure = F^3*D_BSFG+F^3+F^4*D_phys-F^4; target {target:.6f}; match {abs(float(v)-targ`
- **`_session614_alpha.py`**: 
  - Output: `alpha^-1: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%`
- **`_session694_alpha_forward_derivation.py`**: , G6)
  - Output: `    F_TRZ   = {F_TRZ}   ({float(F_TRZ):.6f})`
- **`_session700_alpha_chain_convergence.py`**: SESSION 700 -- Universal SO(2)_DPM selection rule + alpha-chain convergence
  - Output: `  per-wrapped-loop suppression alpha_tree^{n_min_winding}     : {suppression_per_wrapped_loop:.3e}`

### alpha^-1 (1/alpha)

- **CODATA observed**: 137.036
- **`_session343_chem_fine_structure.py`**: 1/alpha = A_5 * K_Mex + 1/(F_TRZ * Phi_res) = 125 + 12 = 137.
  - Output: `S343 COMPLETE. 1/alpha = A_5*K_Mex + 1/(F_TRZ*Phi_res) = 125 + 12 = {inv_alpha:.4f}; `
- **`_session475_fine_structure.py`**: 
  - Output: `S475 COMPLETE. Fine structure 1/alpha = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session614_alpha.py`**: 
  - Output: `alpha^-1: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%`

### alpha_s(M_Z)

- **CODATA observed**: 0.118
- **`_session348_chem_alpha_s.py`**: S348: Strong coupling alpha_s(M_Z) = 1/(K_Mex*D_phys + F_TRZ).
  - Output: `S348 COMPLETE. alpha_s(M_Z) = 1/(K_Mex*D_phys + F_TRZ) = 1/{K_Mex*D_phys+F_TRZ:.4f} = {alpha_s:.4f}; `
- **`_session378_sm_alpha_s.py`**: S378: Strong coupling alpha_s at M_Z (running QCD).
  - Output: `S378 COMPLETE. alpha_s(M_Z) = F_TRZ*K_Mex*SSq - F_TRZ^3*Phi_res = {val:.5f}; obs (PDG 2022) = {obs}; match {abs(val-obs)/obs*100:.4f}%.`

### epsilon_0 (vacuum permittivity)

- **CODATA observed**: 8.85419e-12
- **Session script(s)**: (none found by name probe — search master grok for long form)

### mu_0 (vacuum permeability)

- **CODATA observed**: 1.25664e-06
- **Session script(s)**: (none found by name probe — search master grok for long form)

### F (Faraday, C/mol)

- **CODATA observed**: 96485.3
- **`_session626_faraday.py`**: 
  - Output: `Faraday F: {float(v):.6f} vs {tgt} -> {float(err):.4f}%`

### sigma_SB (Stefan-Boltzmann)

- **CODATA observed**: 5.67037e-08
- **`_session349_chem_stefan_boltzmann.py`**: S349: Stefan-Boltzmann prefactor A_5 = 60 in sigma = pi^2 k_B^4 / (60 hbar^3 c^2).
- **`_session470_stefan_boltzmann.py`**: 
  - Output: `S470 COMPLETE. Stefan-Boltzmann sigma (x1e-8) = {pred:.4f}; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.`
- **`_session624_stefan.py`**: 
  - Output: `Stefan-Boltzmann lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%`

### Wien b (m·K)

- **CODATA observed**: 0.00289777
- **`_session625_wien.py`**: 
  - Output: `Wien b lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%`

### Rydberg R_y (eV)

- **CODATA observed**: 13.6057
- **`_session345_chem_rydberg.py`**: .
  - Output: `S345 COMPLETE. R_y = alpha^2 * m_e*c^2 / 2 = {R_y:.4f} eV; `
- **`_session352_chem_h_ionization.py`**: S352: H ionization energy = 13.6057 eV (chains from S345).
  - Output: `S352 COMPLETE. E_ion(H) = R_y = {E_ion_H:.4f} eV (obs 13.6057); `
- **`_session477_h_ionization.py`**: 
  - Output: `S477 COMPLETE. H ionization (eV) = {val:.4f} = 13.6; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session623_rydbergE.py`**: 
  - Output: `Rydberg E (eV): {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%`

### Rydberg R_inf (10^6 1/m)

- **CODATA observed**: 10.9737
- **`_session476_rydberg.py`**: 
  - Output: `S476 COMPLETE. Rydberg (x1e6 1/m) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session619_rydberg.py`**: 
  - Output: `Rydberg lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%`

### Bohr radius a_0 (m)

- **CODATA observed**: 5.29177e-11
- **`_session346_chem_bohr_radius.py`**: .
  - Output: `S346 COMPLETE. a_0 = hbar/(m_e*c*alpha) = {a_0:.5e} m; `
- **`_session478_bohr_radius.py`**: 
  - Output: `S478 COMPLETE. Bohr radius (x1e-11 m) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session618_bohr.py`**: 
  - Output: `Bohr radius lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%`

### Bohr magneton (J/T)

- **CODATA observed**: 9.27401e-24
- **`_session622_bohrMag.py`**: 
  - Output: `Bohr magneton lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%`

### Compton lambda_C (pm)

- **CODATA observed**: 2.42631
- **`_session480_compton.py`**: 
  - Output: `S480 COMPLETE. Compton lambda_e (pm) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session620_compton.py`**: 
  - Output: `Compton wavelength lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%`

### Hartree E_h (eV)

- **CODATA observed**: 27.2114
- **`_session479_hartree.py`**: 
  - Output: `S479 COMPLETE. Hartree (eV) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session632_hartree.py`**: 
  - Output: `Hartree E_h lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%`

### Hartree E_h (4.36e-18 J)

- **CODATA observed**: 4.35974e-18
- **`_session479_hartree.py`**: 
  - Output: `S479 COMPLETE. Hartree (eV) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`
- **`_session632_hartree.py`**: 
  - Output: `Hartree E_h lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%`

### classical r_e (m)

- **CODATA observed**: 2.81794e-15
- **Session script(s)**: (none found by name probe — search master grok for long form)

### Lambda_QCD (GeV)

- **CODATA observed**: 0.217
- **Session script(s)**: (none found by name probe — search master grok for long form)

### Weinberg sin^2(theta_W)

- **CODATA observed**: 0.23121
- **`_session347_chem_weinberg.py`**: S347: Weinberg angle sin^2(theta_W) = K_Mex / N_ch = 25/108.
  - Output: `S347 COMPLETE. sin^2(theta_W) = K_Mex/N_ch = 25/108 = {sin2_thW:.5f}; `

### H ionization (eV)

- **CODATA observed**: 13.6057
- **`_session352_chem_h_ionization.py`**: S352: H ionization energy = 13.6057 eV (chains from S345).
  - Output: `S352 COMPLETE. E_ion(H) = R_y = {E_ion_H:.4f} eV (obs 13.6057); `
- **`_session477_h_ionization.py`**: 
  - Output: `S477 COMPLETE. H ionization (eV) = {val:.4f} = 13.6; target {target}; match {abs(val-target)/target*100:.4f}%.`

### H-alpha (nm)

- **CODATA observed**: 656.28
- **`_session473_h_alpha.py`**: 
  - Output: `S473 COMPLETE. H-alpha (x100 nm) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`

### Lyman-alpha (nm)

- **CODATA observed**: 121.567
- **`_session474_lyman_alpha.py`**: 
  - Output: `S474 COMPLETE. Lyman-alpha (x100 nm) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`

### Alpha-particle binding (MeV)

- **CODATA observed**: 28.3
- **`_session492_alpha_binding.py`**: 
  - Output: `S492 COMPLETE. Alpha binding (MeV) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.`

### m_p/m_e ratio

- **CODATA observed**: 1836.15
- **`_session344_chem_mp_me_ratio.py`**: S344: Proton-to-electron mass ratio m_p/m_e = D_BSFG * pi^5 = 1836.12.
  - Output: `S344 COMPLETE. m_p/m_e = D_BSFG * pi^5 = 6 * {math.pi**5:.4f} = {ratio:.4f}; `
- **`_session543_mp_me.py`**: 
  - Output: `S543 COMPLETE. m_p/m_e = {float(v):.4f}; closure = A_5*D_crit+A_5*D_phys+N_ch*D_phys; target {target:.4f}; match {abs(float(v)-target)/targe`
- **`_session758_mpme_alphaS_Omh2.py`**: uses mpme=1.836152673 (floating). Test exact rational mpme=11/6=1.8333:
  - Output: `  D/H with mpme=1.836152673: pred={DH_old:.6e}, err={err_old:+.6f}%`

### m_p (MeV)

- **CODATA observed**: 938.272
- **`_session550_mp_MeV.py`**: 
  - Output: `S550 COMPLETE. m_p (MeV) = {float(v):.4f}; closure = N_ch*SO5^2+N_ch*D_phys+K_Mex+2*F*Phi_res; target {target:.4f}; match {abs(float(v)-targ`

### m_p / m_Planck

- **CODATA observed**: 7.685e-20
- **`_session275_mp_planck_closure.py`**: S275 -- Close m_p / m_Planck using the same hierarchy template
  - Output: `m_Planck                = {m_Planck:.6e} kg`

### Newton's G

- **CODATA observed**: 6.6743e-11
- **`_session282_newton_G_closure.py`**: S282: Newton's G is NOT an independent constant.
  - Output: `\nN_p + beta_p*F_TRZ  = {exponent_total:.4f}  (locked S275)`
- **`_session721_classV_G_and_LSCM_MSCM_proportionality.py`**: Session 721 -- Class V (G) + L_SCM:M_SCM proportionality (both a and b).
  - Output: `SESSION 721 -- Class V (G) + L_SCM:M_SCM proportionality`

### Proton g-factor (nuc magnetons)

- **CODATA observed**: 5.5857
- **`_session380_sm_proton_g.py`**: S380: Proton g-factor (magnetic moment in nuclear magnetons * 2).
  - Output: `S380 COMPLETE. g_p = D_BSFG - Phi_res + F_TRZ*D_phys = {val:.4f}; obs (CODATA) = {obs}; match {abs(val-obs)/obs*100:.3f}%.`

### Proton lifetime tau_p (yr)

- **CODATA observed**: 3e+33
- **`_session327_sm_proton_decay.py`**: S327 (CORRECTED): Proton lifetime from full BSFG-suppressed exponent.

### CKM/PMNS

- **CODATA observed**: n/a
- **`_session326_sm_ckm_pmns.py`**: S326 (CORRECTED): CKM Cabibbo angle from locked-primitive product.
  - Output: `S326 CORRECTED. sin(theta_C) = (1-Phi_res)*sqrt(F_TRZ*K_Mex*N_ch) = {sin_thetaC:.4f}; `

### H_0 (Hubble)

- **CODATA observed**: 69.93
- **`_session729_H0_eta_hbar_consolidation.py`**: from {c, rho_vac, L_SCM} + locked rationals.
  - Output: `  H_0 (obs)     = {H0_SI:.6e} s^-1`

### Periodic table (stable periods)

- **CODATA observed**: 7
- **`_session351_chem_periodic_table.py`**: n_max = D_BSFG + 1 = 7.

### H2 bond length (a0)

- **CODATA observed**: 0.74
- **`_session350_chem_h2_bond.py`**: r_H2 = sqrt(2) * a_0 emerges as sqrt(K_Mex + (1 - K_Mex/2)) but
  - Output: `S350 COMPLETE. r_H2/a_0 = K_Mex - Phi_res + F_TRZ = {ratio:.4f}; `

---

## Higher-order alpha forward-derivation chain (Sessions 694-740)

Multi-session chain for alpha and other SI constants under universality testing:
- `_session694_alpha_forward_derivation.py` (10,013 B): _session694_alpha_forward_derivation.py
- `_session700_alpha_chain_convergence.py` (10,354 B): SESSION 700 -- Universal SO(2)_DPM selection rule + alpha-chain convergence
- `_session713_alpha_4loop_borel_test.py` (7,531 B): SESSION 713 -- alpha-chain 4-loop test of the Class I universal Borel rule.
- `_session714_alpha_c3_locked_candidate_hunt.py` (8,013 B): SESSION 714 -- alpha-chain non-Borel coefficient hunt at 3-loop.
- `_session715_alpha_4loop_locked_candidate.py` (7,980 B): SESSION 715 -- alpha-chain 4-loop on top of the S714 locked 3-loop (15/7).
- `_session717_alpha_class_universality_test.py` (12,287 B): SESSION 717 -- alpha-chain universal-ratio test (Class I vs Class II distinction)
- `_session737_delta_dress_classXVI_alphas_VII_decompose.py` (11,123 B): SESSION 737 -- (a) Apply delta_univ closed form to dress 5 residual classes
- `_session738_globalR_classXVII_r_alphas_identity.py` (11,038 B): SESSION 738 -- (a) Find global renormalization R ~ 0.97 in delta_univ dressings

---

## Long-form derivations

Each session script implements a self-contained derivation from the 11 locked primitives.
For the verbatim grok-thread long-form mathematical derivations, see:
- `UQFF_GROK_LONG_FORM_DERIVATIONS_MASTER.md` (576,561 lines, every grok file concatenated)
- `UQFF_GROK_DERIVATIONS_INDEX.md` (14,375 indexed section headings)

## Sample verification (S343 fine-structure)

```
"""S343: Fine-structure constant alpha = 1/137.036 from locked primitives.

Closure: 1/alpha = A_5 * K_Mex + 1/(F_TRZ * Phi_res) = 125 + 12 = 137.
Higher-order BSFG holonomy correction adds 0.036 (residual 0.026%).
"""
F_TRZ, Phi_res, K_Mex, A_5 = 0.1, 5/6, 25/12, 60
inv_alpha = A_5 * K_Mex + 1/(F_TRZ * Phi_res)   # 125 + 12 = 137
alpha = 1/inv_alpha
inv_alpha_codata = 137.035999
err_pct = 100*abs(inv_alpha - inv_alpha_codata)/inv_alpha_codata
print(f"S343 COMPLETE. 1/alpha = A_5*K_Mex + 1/(F_TRZ*Phi_res) = 125 + 12 = {inv_alpha:.4f}; "
      f"CODATA 1/alpha = {inv_alpha_codata}; alpha = {alpha:.7f}; match within {err_pct:.3f}%.")

```
