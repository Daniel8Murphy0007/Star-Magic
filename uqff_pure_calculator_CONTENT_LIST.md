# uqff_pure_calculator.py — Content List / Inventory

**File:** `uqff_pure_calculator.py`  
**Version:** `1.2.0+pure_calculator_strip_07Jun2026`  
**Approx. size:** >24,500 lines (monolithic self-contained module)  
**Style:** Pure functional / procedural. Zero top-level `class` definitions. All logic in `_`-prefixed helper functions + large `Dict[str, Callable]` registries + a few high-level dispatchers.  
**Dependencies:** Only `math` + `typing` (stdlib). "Pure/stripped" — no imports from sibling modules in the Star-Magic workspace.  
**Purpose:** Standalone UQFF (Universal Quantum Field Framework) / MUGE / Star-Magic calculator + derivation engine + exhaustive closure/validation ledger. Computes "primitive saturation" predictions for SM + cosmological + astrophysical + LENR observables from a small locked set of UQFF parameters (26D ledger). Provides per-object g-force / resonance / buoyancy masters and hundreds of session/paper-specific falsifiable claims (Sxxxx / Lxx closures) with predicted vs. observed + % diff.

---

## 1. Header, Imports & Core Constants (≈ lines 1-160)
- `__version__`
- Base SI-ish: `PLANCK_H`, `EV_J`, `K_B`, `N_AVOGADRO`
- UQFF "first principles" seeds: `V_FERMI_PROXY`, `E_ZERO_UQFF`, `F_THZ_UQFF`, `_PHI_RES_INIT`, `_D_CRIT_INIT=26`, `_SSQ_INIT=0.57`, `_G8_26_INIT`
- **Derived fundamentals** (not CODATA):
  - `C_LIGHT` (from D_CRIT * pi / phi * V_fermi)
  - `G_NEWTON` (complex expression involving factorial(26), SSQ, E0, FTHz, Vfermi)
- UQFF densities & scales: `RHO_SCM=7.09e-37`, `RHO_UA=10*RHO_SCM`, `BETA_I≈0.603`, `OMEGA_SCM=1.25e12`, `PHI_RESONANCE=0.84`, `SSQ=0.57`, `KAPPA`, `TRZ=0.1`, `D_CRIT=26`, `D_BSFG=6`, `S_26≈1.453`, `A_26 = sum(i**6 for i=1..26)`, `S26_3`, `S26_DPM=1.4531e26`, `G1_K=5/6`, `G2_BETA_BASE=3/5`, `G3_RICCI_COEF=0.5`, `G4_BSFG_COEF=0.15`, `G5_KK_SUPPRESS=1.624e-37`, `G8_26_BARRIER=factorial(26)`
- Weights: `W_C=W_R=W_B=0.33...`
- Many system-specific constants (magnetar B decay times, radii for SgrA*, Crab, Pillars, Westerlund2, NGC objects, Saturn, M16, Horsehead, Antennae, Sombrero, HUDF, etc.)
- Small pure helpers: `_e_react_uqff`, `_e_crack_uqff_*`, `_mu_s_uqff`, `_a_dpm_uqff`, `_F_DPM_uqff`

## 2. Provenance / Cluster Registry (≈120-145)
- Long `PROV_*` string constants documenting historical source clusters (DaVinci handwritten files, Bearden scalar tech, UFE orb experiments, arXiv batches, 11Sept/11Oct astro systems, Lagrangian G1-G8, 14Sept 71-eq catalog, 99system, UA vacuum manifold, etc.).
- `CLUSTER_REGISTRY` dict with many alias keys → provenance text (used for traceability/auditing).

## 3. Millennium Problem Targets & Derivers (≈147-224)
- `MILLENNIUM_TARGETS`: 8 Clay + open problems (Yang-Mills gap, Riemann zeta, BSD, Navier-Stokes, Hodge, Poincaré, P-vs-NP, Black-hole info/Page curve) with target values, units, reference kind, descriptions.
- Individual `_millennium_*_derive()` functions (some call other primitives or `_bh_finite_bound_report`).
- `_MILLENNIUM_DERIVE` registry + `_millennium(name)` dispatcher with alias normalization + diff computation against ledger.

## 4. L96 / Session-276 Paper Verification Layer (early) (≈226-276)
- `_l96_verify_session276_image_papers()` — runs probes over a hardcoded list of ~30 PAPER_98x-101x titles, collects closures/transcribed/exact counts.
- `_l96_dump_session276_lambdas()`
- Sample lookups (S3851, S3876).

## 5. Master Chain / Literal Anchors + Spinor (≈277-355)
- `_master_chain_base()`, `_BASE_CHAIN`
- `_SM_LITERAL_ANCHOR_SAT` (alpha, proton_mass, YM gap, neutron lifetime, H0, t0 mapped to /_BASE_CHAIN)
- `_sm_literal_anchor` + `_master_constant_formula` alias
- Spinor bundle locks (`SPINOR_ANCHORS`, `_spinor_canonical_engine_derive`, `_spinor_closure` with % diffs to B9 anchors).

## 6. Primitive Saturation Functions — The Heart of "Zero-Param" Derivations (≈357-954 + many more)
Huge family of `_xxx_primitive_sat() -> float`:
- **Particle / SM constants**: alpha, m_p, m_e, m_mu/tau/t, m_W/Z/H, v_higgs, g_f, alpha_s, sin2_theta_w, CKM/PMNS elements, a_e/a_mu, g_e, E_hartree, hyperfine_cs, etc.
- **Cosmological**: Omega_m/b/dm, n_s, A_s, eta, Y_p, z_re, tau_reion, w_z05, f_NL, sigma_8, T_CMB, T_neutrino, BAO_rd, hudf_z, H0/t0 primitives, etc.
- **Astro g-values**: sgr_a_g, sgr_1745_g (magnetar), tapestry_g, westerlund_g, pillars_g, crab_g, m16_g, horsehead/antennae/sombrero/ngc*/hudf_g, rings_of_relativity, etc.
- **LENR / condensed matter**: multiple rossi, parkhomov, pons-fleischmann, mizuno, mckubre, stringham, brillouin, plus 630eV cluster.
- **QED / atomic / SI**: r_infinity, sigma_SB, b_wien, a0, lambda_c, r_e, mu_B, R_gas, Faraday, von Klitzing, Josephson, vacuum permittivity/permeability.
- **Planck / fundamental**: h_planck, c_light, N_A, planck_mass/length/time/temp, hbar, k_b, etc.
- **Inflation / GW / slow-roll**: r_tensor_scalar, dn_s, f_NL_equil/orth, epsilon/eta_slow_roll, N_efolds, t_reh, H_inflation, etc.
- **P1-P5 predictions**, quark masses, other exotics.
- `_LEDGER_PRIMITIVE: Dict[str, Callable]` — central map of 70+ short names → the sat functions (extended with more aliases later).

## 7. High-Level Constant Derivation Dispatch (≈1124-1500+)
- `_master_constant_primitive(name)` — giant alias normalizer (fine_structure_alpha → alpha, m_p_mev, yang_mills_gap_gev, hubble, t0, w, eht_sgra_shadow, jades_gs_z14_mass, gw170817_inspiral, sgr_a_g, lenr_rossi_ev, many Planck/SI/astro/LENR keys) then calls the primitive or falls back.
- `_derive_constant(name, *args)` — even higher-level public-ish entry: handles "h","c","g","alpha","muge","g_compressed","master_lagrangian","99system_catalog","rho_scm_energy", "all_modes", "si_units", special cases for 630eV, 1.25THz, G1..G8, ua_4layer, u_mi, f_u_bi_i_*, etc. Returns raw value or dicts for complex items.
- Related: `_g1_mexican_hat`, `_g2.._g8_*`, `_vacuum_ledger_4term`, prediction helpers (P1-P5), `_muge_dual_validate`, `_generate_99system_catalog`, `_si_unit_derivations`, `_f_u_bi_i_mode(s)`, `_rho_*_energy_density` / mass equiv.

## 8. Core UQFF Dynamics & Master Equations (≈2505-2600+ and L95 sections)
- `_cos_pi_tn(t_n=0)` helper (time modulation, implied).
- `_master_lagrangian(R_26=1, ...)` — PAPER_1167 6-term L_FU (T1..T6): GR-like, DPM gauge, interaction (beta_i * Ug*Ub), Um self, UA kinetic, Mexican-hat SCm term. Returns breakdown + L_total. Requires 4-element Ug/Ub vectors.
- `_g_compressed(M, r, t_n)` — 8-term master gravity (Ug1..4 + psi + phi + quantum_integral + buoyancy_fluid).
- `_g_resonance(...)` — DPM + THz + aether resonance modulated form.
- `_ua_4layer_explicit` — 4-layer UA DPM (prime/double/triple/quad).
- `_u_mi` — Universal Inertial Operator (DaVinci heritage).
- `_f_u_bi_i_steps` — 4-step Archimedes → LEP → Q-wave (1.25 THz) → compressed buoyancy chain.
- L95 family (scattered 27xx-34xx):
  - `_l95_l_uqff_master_synthesis`, `_l95_l_uqff_master_sum`, `_l95_fubi_closure_identity`, many F_UBi / F_U / stationarity-at-F_U=1 helpers, `_l95_g_uqff_compressed_master`, galactic bar, AGN BZ jet, etc.
  - Vacuum energy, stationarity residual, L_total → F_U_at_stationarity logic.
- Other: `_L24_F_UBI_HZ` etc., 4-mode (compressed/resonant/buoyant/superconductive) support.

## 9. Per-Astrophysical-System Master g Calculators (≈3507 — ~4500+)
Detailed, time-/B-field-/spin-aware implementations (many with decay, dOmega/dt, GW quadrupole, H0 modulation, Lambda, B_crit suppression, Ug2/Ug3 injection, cosm terms):
- `_magnetar_g_master_uqff` (and v2) — full B(t) exp decay, spin Omega + dOmega, GW term.
- `_sgr_a_g_master_uqff`
- `_tapestry_g_master_uqff`, `_westerlund2_g_master_uqff`
- `_pillars_g_master_uqff`, `_rings_g_master_uqff` (Einstein rings), `_ngc2525_g...`, `_ngc3603`, `_bubble_nebula`, `_antennae`, `_horsehead`, `_ngc1275`, `_hudf`, `_ngc1792`, `_sombrero`, `_saturn`, `_m16_eagle`, `_crab_g_master_uqff`
- `_l96_ngc1316_g_master` (Fornax merger elliptical, special)
- Plus earlier primitive g versions for many of the same objects.
- Some L95 galactic/AGN helpers.

## 10. L94 / L95 / L96 / Higher Session Layers & Massive Closure Registries (thousands of lines)
- L94: Higgs off-shell, polyakov, Calabi-Yau, inflation/CMB Lagrangians, Hubble modulated, swampland, QEC phonon, NCG, CGC, GW strain, BH shadow, GUP, QGP Alice, SCM/UA duality, cosmic egg partition, cos(pi tn) net zero, A26, etc.
- L95: A26 sextic, M_AMU_dpm, Lambda/H0 asymmetry, G6/G7 closed, many master L/FU synthesis + FUBi identities + stationarity + cross-validate vs SM/GR claims (Yukawa vs F_U=1 umbilicus, integrated forces, etc.).
- **L96 / S276+**: 
  - Paper probe/verification harnesses.
  - Taylor-Green (NS smoothness proxy via enstrophy boundedness with UA_scalar ledger saturation + phonon damping).
  - Dozens of `_l96_uqff_*_derived()` / `_report()` / `_closure()` for:
    - pH (-37 derivation), input power (27W), Hubble time, observable universe radius (c/H0 + expansion), CMB LSS distance.
    - Cosmological constant finetuning (rho_vac = rho_SCm * S26 * Phi resolves 120 orders).
    - Cusp/core + missing satellites (SSQ * Phi * K_MEX buoyancy smooths NFW → Burkert).
    - Quantum measurement (spinor bundle projection dim 2^(13), Born emergent, F_U=1 observer).
    - Arrow of time (entropy from low at BB via SCm-UA contact → 26^2 * S26_DPM).
    - Firewall/AMPS (ER=EPR + finite 26! bound + Page recovery).
    - Dark matter identity (7 keV sterile neutrino from D_phys chain + 3.5 keV line).
    - QG UV completion (26D compactification, 26! finite cutoff replaces renormalization).
    - Proton radius puzzle (reconciles muonic/electronic via UQFF form).
    - de Sitter swampland (static ledger + F_U=1 + Mexican-hat min stable, no quintessence).
    - Maxwell demon (Landauer erasure), and many more (from grep: S2815, S3215 etc.).
- **Closure Registries** (huge `Dict[str, tuple]`):
  - `_PA_S276_CLOSURE_REGISTRY`, `_PA_S272_CLOSURE_REGISTRY`, `_PA_S273_CLOSURE_REGISTRY`, ...
  - Format per entry: (paper_tag, label, optional_lambda_form, predicted, observed, unit)
  - Used by verification / dump / probe functions that report "closures", "transcribed_count", "exact_count", % diffs, honesty notes about SSQ lock conflicts etc.
- Later layers: L46 math helpers?, _L55_JWST_HIGHZ tension ledger (8-row z=7-13 massive galaxy abundance, photo vs spec split, inverse-variance, inter-kind tension), `_l55_*` evaluation/summary/anchor_validation/inventory functions. More Sxxxx (S2600+ DE_GAMMA, EOS, FUBII, INFL_LAG, ...).

## 11. Other / Supporting / End Sections (scattered + tail)
- BH finite bound / Page curve reports, vacuum ledger paper-1170 explicit, many cross-master comparisons (g_compressed vs resonance vs FUBi).
- 99-system catalog generator + triadic residual logic (links to 99system_master_equation.py heritage).
- Extensive string metadata inside dicts: "uqff_mechanism", "uqff_claim", "verdict", "primary_source", "paper_basis", "F_U_normalization":1.0, "F_U_eq_1_normalization":True, etc. for every major claim.
- Proof outlines, delta-S-delta-phi=0 at F_U=1, 26-factorial finite bounds, etc.
- (Tail) More L55 JWST high-z ledger processing, additional S/PAPER metadata lists, final summary dicts. No obvious `if __name__ == "__main__":` demo runner visible in samples (may be present or intended for external callers).

---

## Key Architectural Patterns
- **Zero free parameters post-locks**: All predictions flow from fixed UQFF seeds (RHO_SCM, 26, 0.57, 0.84, G-coeffs derived from 26D Lagrangian, 1.25 THz phonon, S26_DPM, etc.).
- **F_U = 1 stationarity / normalization**: Central "belly-button" / umbilicus condition where buoyancy + resonance + gravity terms cross; mass, observers, stability, etc. emerge here.
- **Triadic / Compressed / Resonance / Buoyancy modes**: Wc+Wr+Wb weighting, 8-term g_compressed, explicit 4-step F_UBi chain, 4 explicit modes.
- **Downward 26D projection + DPM grind + SCm/UA contact** (echoes 26D_DOWNWARD_PROJECTION.md).
- **Verifiability first**: Every important number has a closure entry with form (transcribed lambda or None), predicted, observed, unit + %diff machinery. "Honesty rule" comments about locked params vs paper values.
- **Layered sessions (L94/L95/L96 + S27x + PAPER_xxxx)**: Incremental addition of closures without breaking the pure calculator core.
- **Self-documenting**: Massive inline provenance + "primary_source" + "uqff_derivation_form" strings.

---

## Usage Entry Points (call these)
- `_derive_constant("alpha")`, `("hubble")`, `("sgr_a_g")`, `("lenr_rossi_ev")`, `("master_lagrangian")`, `("g_compressed")`, `("all_modes")`, `("muge")`, ...
- `_master_constant_primitive("yang_mills_gap_gev")`
- Specific: `_magnetar_g_master_uqff(...)`, `_crab_g_master_uqff(...)`, `_l96_ngc1316_g_master(...)`
- `_l95_l_uqff_master_synthesis(...)`, `_master_lagrangian(...)`
- Verification: `_l96_verify_session276_image_papers()`, `_l96_dump_session276_lambdas()`, `_l96_uqff_*_closure()` family, `_l55_jwst_highz_inventory()`
- Millennium: `_millennium("riemann")` or `"yang_mills"`
- 99-system: `_generate_99system_catalog()`

---

**Notes for maintainers / auditors:**
- This file is the "pure" distillation of years of session work (DaVinci → 2025 batches → 2026 L/S papers). It is intentionally large and repetitive in structure so a single `python uqff_pure_calculator.py` (or import) can stand alone for audits, cross-checks, or embedding.
- Many "closures" intentionally have `form=None` when they are direct observable anchors or conflict with locked SSQ etc.; the point is transparency.
- For production pipeline / batch astro work see earlier sprint artifacts (production_pipeline.py etc.); this is the pure physics kernel.

*Generated from direct source analysis on 2026-06-12. File continues with additional S/PAPER registries and L-layer evaluations beyond the sampled tail.*