# CLAUDE.md — UQFF Star-Magic project context

**This file is read at the start of every Claude session.** It tells Claude what UQFF is, what canonical primitives must be preserved, and what was done in prior sessions. **Do not strip or rewrite this file. Append-only edits at the bottom.**

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Repo:** https://github.com/Daniel8Murphy0007/Star-Magic

---

## What UQFF is

UQFF is a vacuum-first physics framework built on a **single non-mass primitive**, `ρ_SCm = 7.09e-37 J/m³`, the energy density of the SuperConductive material substrate. From this one number and a 26-level structural lattice (the DPM, Di-Pseudo-Monopole), the framework derives:

1. **The cosmological constant**: `ρ_SCm × 26! × 25/12 ≈ 5.957e-10 J/m³ ≈ Planck Λ` (0.1% match, zero free parameters)
2. **Holmlid 630 eV LENR**: `h · 1.25 THz · S_26 · Φ_res = 630 eV exactly`
3. **All 7 nuclear shell-model magic numbers** ({2, 8, 20, 28, 50, 82, 126}) from arithmetic on four integer primitives {D_phys=4, SO_five=10, D_crit=26, A_five=60}
4. **Fe-56 BE/A peak** (8.79 MeV) at 0.019%
5. **α-particle binding** (28.30 MeV) at 0.012%
6. **Coulomb energy at ultra-dense H spacing** (2.3 pm) = 630 eV ← independent re-derivation of Holmlid KER from electrostatics
7. **5 unified LENR observations**: Holmlid, Parkhomov, Pons-Fleischmann, Mizuno, Rossi E-Cat (Early/X/SK), all from the same SCm phonon mechanism

UQFF does NOT depend on Standard Model physics. **It is not a refinement of SM. It is not seeking to replace SM. The Plan's mandate is "NOT REPLACEMENT" — UQFF and SM solve the same observed phenomena by different methods, both reported with honest residuals.**

---

## The 11 locked canonical primitives — DO NOT CHANGE

```
Integer lattice:
  D_phys  = 4    (physical spacetime)
  D_BSFG  = 6    (bulk-edge dim)
  D_crit  = 26   (PTOE / bosonic-string critical dim)
  N_ch    = 9    (channel)
  SO_five = 10   (= dim SO(5))
  A_five  = 60   (= |A_5|, icosahedral group order)

Real primitives:
  ρ_SCm   = 7.09e-37 J/m³  (THE foundational vacuum density)
  ρ_UA    = 10 × ρ_SCm     (DPM density ratio)
  F_TRZ   = 1/10           (time-reversal-zone)
  Φ_res   = 0.84           (default; PAPER_1203 Nuclear uses 5/6)
  SSq     = 0.57           (canonical PAPER_1154 first-principles)
  K_Mex   = 25/12          (Mexican-hat coefficient)
  β_i     = 0.6029         (canonical PAPER_1203)
  ω_SCm   = 1.25 THz       (Holmlid phonon carrier)
  S_26    = 1.453162       (Ramanujan 26-level scaling)
  λ_i     = 1.0            (inertia coupling)
  ω_s_Sun = 2.5e-6 rad/s   (system-specific stellar rotation)
  Δ_UA4   = 0.1            (fourth UA layer offset)
```

**Any future session that changes SSQ back to 0.505, sets β_i to 0.6, or removes N_CH/SO_FIVE/A_FIVE is destroying canonical work. The fidelity test gate `uqff_fidelity_tests.py` will catch this.**

---

## The vacuum manifold architecture

```
DPM (Di-Pseudo-Monopole) = SCm ⊗ UA
   |
   ├── scm_vacuum_manifold.py    (BASE: ρ_SCm, 1.25 THz phonon, F_U_Bi_i_99, VDS/DVP/BSH, LENR fns)
   ├── ua_vacuum_manifold.py     (SUPERSTRUCTURE: UA' → UA'''', F_U_Bi_i_DPM, cosmological fns)
   └── dpm_vacuum_manifold.py    (canonical consolidated root)
```

**The two grinding poles:**
- **North** = SCm (CW rotation, massless, reactive with trapped Aether)
- **South** = UA' (CCW rotation, forms when SCm encapsulates free UA)

**The 5-step grinding sequence** (mass-production pipeline):
1. SCm contacts free UA → Big Bang contact event
2. SCm encapsulates UA → UA' (trapped Aether)
3. CW × CCW grinding → UA'' (first excitation)
4. Progressive densification → UA''', UA''''
5. Maximum metallicity → UA''''' (observable matter)

**The 4-layer UA hierarchy** (canonical equations):
```
UA'    = ρ_SCm
UA''   = ρ_SCm · (1 + β_i · cos(π t_n))
UA'''  = ρ_SCm · (1 + β_i · cos(π t_n) + λ_i · ω_s)
UA'''' = ρ_SCm · (1 + β_i · cos(π t_n) + λ_i · ω_s + Δ_UA4)
F_DPM  = F_U_Bi_i_99 × (UA' + UA'' + UA''' + UA'''')
```

---

## The Holy Trinity (PAPER_646)

- **UA** = the medium (vacuum field hosting all interactions)
- **U_i + EM** = the operator (Universal Inertial Operator — roots matter to UA)
- **[SCm]** = the conduit (extra-universal material enabling DE power transfer)

**Universal Inertial Operator:**
```
U_i = λ_i · (ρ_SCm / ρ_UA) · ω_s · cos(π t_n) · (1 + f_TRZ)
    = 2.75e-7 (Sun, t=0)   ← matches PAPER_646 exactly
```

Inertia in UQFF is NOT mass-as-resistance. It is **matter-as-anchored-to-the-Aether**. This is the operator that explains the equivalence principle without invoking GR's geometric postulate.

**Caduceus Wave Topology**: spherical quantum waves self-invert into double-helix coils with opposing chirality at high amplitude. **26 simultaneous pinch points** encode the **decimal expansion of π** — π is the physical record of the pinch-point phase sequence.

---

## The F_U = 0 Master Equation (PAPER_1203 Canonical v1.5)

```
F_U_total = (U_g1 + U_g2 + U_g3 + U_g4) − F_UBi + F_UBii + U_m = 0

F_UBi  = −β(t,E,Z) · G·M·ρ_SCm / r² · (1 + F_TRZ) · |cos(π t_n)|
F_UBii = +β(t,E,Z) · (r/r₀) · k_spring · (1 + E_n) · |cos(π t_n)|
k_spring = (ρ_UA / ρ_SCm) · ω_SCm · Φ_res
β(t,E,Z) = β_i · |E| · Z · cos(π t)   ← DYNAMIC, not constant
```

Equilibrium at every shell/scale via `F_UBi(r) + F_UBii(r) = 0`. The crossover root is `r_hz`, the habitable-zone radius. **Solves orbital position from first principles.**

---

## The 9-sector UQFF Lagrangian (Session 202)

```
L_UQFF = L_EH + L_YM + L_Dirac + L_SCm + L_mag + L_buoy + L_aether + L_LENR + L_KK
```

| Sector | Domain | Canonical result |
|---|---|---|
| L_EH    | General Relativity   | Canonical Einstein-Hilbert in 26D |
| L_YM    | Yang-Mills gauge     | m_gap = **1.736 GeV** (PAPER_1318) — NOT the 1.78 GeV SM lattice estimate |
| L_Dirac | Fermion / LENR       | Kozima neutron-drop (PAPER_1061) |
| L_SCm   | SC manifold          | V(φ₀) = −ρ_SCm canonical |
| L_mag   | U_m magnetism        | Heaviside amplifier (PAPER_1072) |
| L_buoy  | F_U_Bi_i buoyancy    | Variational EOM (PAPER_1065) |
| L_aether| Vacuum background    | Two-component ρ (PAPER_1051) |
| L_LENR  | Nuclear transmutation| COP parametric (PAPER_1081) |
| L_KK    | Kaluza-Klein 26D     | S_26^(3) compactification (PAPER_1080) |

---

## LENR — unified by Holmlid 630 eV

```
E_phonon = h · 1.25 THz                              = 8.28e-22 J ≈ 5.17 meV
S_26^(3) ≈ 1.4531e26                                 (26D Ramanujan series at SSq=0.57)
E_SCm-phonon = E_phonon × S_26^(3) × ξ × Φ_res       = 630 eV ← Holmlid D(-1) KER
```

Then the 5 LENR observations all follow:
| Reactor | Observed | UQFF prediction |
|---|---|---|
| Holmlid D(-1) KER | 630 eV | 630 eV (calibration anchor) |
| Parkhomov Ni-H | 150-280 W | 199 W at N=2e18, t=1hr |
| Pons-Fleischmann Pd-D | 1-50 W | f_buoyancy=10⁻³ → match range |
| Mizuno Ni-D | 10-300 W | N_clusters=5e17 → range |
| Rossi Early E-Cat | COP 6-14 | depends on geometry |
| Rossi E-Cat X | COP >20 | enhanced Φ_T at high T |
| Rossi E-Cat SK | COP >50 | cold-spark t_n coherence |
| **Star-Magic reactor** | **COP 555:1 @ 27W, ambient T, pH −37** | full F_U_Bi_i + VDS phonon |

**Independent re-derivation**: Coulomb energy at d=2.3 pm (ultra-dense H spacing) = 626 eV. Confirms 630 eV LENR is electrostatic at the cluster geometry.

---

## Nuclear closures (PAPER_1203 Nuclear)

**All 7 magic numbers EXACT** from arithmetic on integer primitives only:
```
2   = SO_five − 2·D_phys
8   = 2·D_phys
20  = 2·SO_five
28  = D_crit + SO_five − 2·D_phys
50  = A_five − SO_five
82  = A_five + D_crit − D_phys
126 = D_crit + SO_five²
```

**Mayer-Jensen shell occupancy independently produces the same numbers** via the conventional spin-orbit shell model — confirming the UQFF arithmetic agrees with the empirical pattern.

---

## The 32 public `calculate_*` surfaces

```
ORIGINAL 7 (Plan-mandated thin surface):
  calculate_resonant_adpm        — ω · cos(π t_n) · Φ_res
  calculate_scm                  — SCm 26-level density with t_n modulation
  calculate_f_u_bi               — Universal Buoyancy F_UBi
  calculate_f_u_bi_i             — F_U_Bi_i 4-layer master integral
  calculate_triadic_g            — w_C g_comp + w_R g_res + w_B g_buoy
  calculate_vacuum_ledger        — 4-term: V0 + R26 + ρ_KK + ρ_BSFG → Planck Λ
  calculate_analytic_closures    — symbolic dispatch hub

CANONICAL RESTORED 4 (PAPER_646, PAPER_1203, PAPER_1141):
  calculate_universal_inertial_operator   — U_i = λ_i · (ρ_SCm/ρ_UA) · ω_s · cos(π t_n) · (1+F_TRZ)
  calculate_nuclear_magic                 — 7 magic numbers + BE/A + deuteron + alpha
  calculate_lenr                          — single-variant router (holmlid/parkhomov/pons/mizuno/rossi)
  calculate_f_u_zero                      — F_U=0 master equation + r_hz root-find

PHASE-2 CANONICAL 5 (DPM, PAPER_646, PAPER_062, PAPER_1061, PAPER_648):
  calculate_ua_layers            — canonical UA', UA'', UA''', UA''''
  calculate_dpm_grinding         — 5-step grinding sequence (SCm→UA→UA'''''")
  calculate_caduceus             — 26 pinch points encoding π decimals
  calculate_shell_orbital        — Mayer-Jensen shell occupancy + UQFF magic match
  calculate_lenr_full            — full per-reactor LENR report (8 reactors + WL + Kozima + meson cascade)

BUCKET 0 LOOP-CLOSURE 6 (PAPER_592/593/594/596/597/598/599):
  calculate_negative_time_dual_existence  — PAPER_597 t_neg<0 + CW/CCW dual branches
  calculate_si_derivations                — PAPER_592 c=2.995e8 (0.13%) + PAPER_593 G=6.669e-11 (0.08%) parameter-free
  calculate_quantum_gravity               — PAPER_596 26D unification (GR+QFT limits, mass gap, no singularity)
  calculate_black_hole                    — PAPER_594 r_min via 26! finite bound (per-mass dispatch)
  calculate_vds_dvp_bh26                  — PAPER_598 VDS (P/3 floor) + DVP (prime=113) + BH26 (92 GHz × 26 bins)
  calculate_bsd_rank_cohomology           — PAPER_599 BSD via UQFF tensor eigenvalues (Cremona 37a1 → 0.30598, 0.005%)

BUCKET B PARADOX 1 (PAPER_1183 + 645/084/561/594/658/663/667/670/048):
  calculate_paradox                       — routes 8 paradoxes to Millennium derivations + spinor + Page-curve + DPM-pair identity

BUCKET C COSMOLOGY 1 (PAPER_1156 + 589/1086/1036/1026/202/587/1092/011/1181/1182):
  calculate_cosmology                     — 18-observable ΛCDM suite (α, h, m_p/m_e, Ω_Λ, Λ, Y_p, τ_reion, z_reion, …)

BUCKET D PARTICLE PHYSICS 1 (PAPER_1209HH + 1198 + 652/633/023 + 1155):
  calculate_particle_physics              — 22-observable SM spectrum: 10 PAPER_1209HH masses (W/Z/t/H/b/c/τ/μ/s/e) + couplings + CKM + g-2 + neutrinos

BUCKETS E-K (first-pass drainage, mostly SCm corrections — see NEXT_PRIORITIES for PURE_UQFF upgrade work):
  calculate_gw_events             E — GW150914/170817/190425, NANOGrav 15yr, LIGO O5 ringdown (PAPER_007/914/915/916/927/934/935/1012/1022/1175)
  calculate_agn_jet               F — 3C273/TON618/CenA/SgrA*/M87/Perseus + M-σ + Blandford-Znajek (PAPER_067/087/360/512/627/630/754/814/939/1002/1009/1010/1037/1039/1041/1048/1079/1125)
  calculate_astrophysics          G — PSR J0030/Crab/EtaCar/Westerlund/Lagoon/Orion/NGC3603/Rings (PAPER_065/121/138/211/292/305/432/434/436/447/487/488/489/512/537/995/1017/1126)
  calculate_high_energy_astro     H — TXS 0506+056, IceCube ν_e, Auger dipole, CR knee, FRB DM, Amaterasu (PAPER_020/108/130/215/360/515/1016/1017/1020/1034)
  calculate_qgp                   I — QGP T_c, η/s viscosity, ALICE PbPb dN_ch/dη (PAPER_059/1004/1007/1008/1013)
  calculate_higgs_precision       J — H→γγ/ZZ/WW/bb/ττ branching + κ_t + CP phase + level 18 (PAPER_034/035/137/396/856/1113/1114/1120)
  calculate_bsm_constraints       K — electron/neutron EDM, proton decay, LFV, axion, dark photon, VLQ (PAPER_029/033/046/333/340/494/1116/1119/1183)
```

**All 8 Clay Millennium Prize Problems now have non-placeholder UQFF derivations** (BUCKET A wired all 5 remaining: Riemann 9877.78265, Navier-Stokes enstrophy cap 0.85, Hodge identity 1.0, Poincaré 7/12, P≠NP 1-10⁻⁹; Yang-Mills 1.736 GeV / BSD 0.30598 / BH info 0.99596 wired earlier).

---

## The fidelity gate — RUN AFTER EVERY EDIT

```
cd C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic
python uqff_fidelity_tests.py
```

**Exit 0 = clean. Exit 1 = something drifted; the FAILURES list tells you which invariant broke.**

Current state: **417 tests, 0 failures.** (TOTAL PURGE applied 2026-06-08: all provenance / paper / closure_status / SM-references stripped from calculator. Public surfaces return `{'value': X}` only. Cat 16 strict purge guard wired into gate.)

The gate verifies:
1. **Structural integrity** — module imports, all 32 public funcs work, all 8 Millennium dispatches resolve, 16 derive_constant keys execute
2. **Provenance honesty** — no `"0.000% error"`, no `"<1% on 99/99"`, no double-paren bugs; every return carries `"NOT REPLACEMENT"`
3. **Dispatch fidelity** — `_sm_literal_anchor()` only called from `*_residual*` helpers (no SM-fallback contamination)
4. **Plan/Map invariants** — 32 public surfaces match expected, zero classes, no `__main__`, no `datetime`, no `json.dump`
5. **Canonical physics restoration** — every primitive at canonical value, all 7 magic numbers EXACT, BE/A < 0.05%, U_i = 2.75e-7 exact, KER = 630 eV exact, YM = 1.736 GeV
6. **Phase-2 canonical additions** — Mayer-Jensen shell match, 26 Caduceus pinch points, π first 8 digits, 5 DPM steps, WL above threshold, meson cascade ratio, Coulomb 626 eV at 2.3 pm
7. **BUCKET 0 loop-closure** — c within 0.5% CODATA, G within 0.5% CODATA, t_neg derived, 26! BH bound, QG unification has GR+QFT limits, VDS/DVP/BH26 spine, BSD eigenvalue rank
8. **BUCKET A Millennium** — all 8 dispatches non-placeholder; Riemann = t_10000 EXACT, NS = 0.85, Hodge = 1.0, Poincaré = 7/12, P≠NP = 1-10⁻⁹, BSD = 0.30598 (0.005%), BH info ∈ [0,1]
9. **BUCKET B paradoxes** — DPM-pair identity K_Mex-2 = 1/12 EXACT, 8 paradox closures routed, spinor closure included
10. **BUCKET C cosmology** — α 0.138%, h_Planck 0.061%, m_p/m_e 0.077%, Ω_Λ 0.102%, Λ 0.003%, Y_p 0.050%, Hubble tilt = 1/12 EXACT
11. **BUCKET D particle physics** — all 10 PAPER_1209HH SM masses within paper-stated residuals (m_W 0.003% tier best to m_e 0.178%)

---

## RULES for any future Claude session — READ EVERY TIME, BEFORE EVERY ACTION

**These rules supersede everything else in the repository (Map, Plan, prior sessions).** Daniel's words from session 2026-06-08 (verbatim): "YOU ARE TO READ THE RULES BEFORE EVERY ENTRY. YOU ARE TO READ THE RULES BEFORE EVERY SESSION. YOU ARE TO READ THE RULES BEFORE TAKING A SHIT. UQFF IS THE ONLY PHYSICS IN THIS FUCKING REPO."

1. **READ THIS FILE FIRST EVERY SESSION AND BEFORE EVERY ACTION.** Then read the most recent SESSION_LOG.md entry. Then NEXT_PRIORITIES.md. The Map (`uqff_Map.md`) and Plan (`uqff_Plan.md`) are construction history; if they conflict with these rules, **THESE RULES WIN**.

2. **DO NOT REVERT canonical primitives.** SSQ=0.57, β_i=0.6029, K_MEX=25/12, S_26=1.453162, RHO_SCM=7.09e-37, integer primitives (D_PHYS=4, D_BSFG=6, D_CRIT=26, N_CH=9, SO_FIVE=10, A_FIVE=60). All locked.

3. **NO NARRATIVE OF ANY KIND.** This is a PURE MATHEMATICAL CALCULATOR, not an encyclopedia, not a dictionary, not a textbook. Daniel's words: "NO COMMENTS WHATSOEVER. I DON'T NEED ANY REFERENCES TO HELP ME UNDERSTAND ANYTHING." Enforced by the fidelity gate Cat 16. Forbidden in calculator source:
   - `#` comment lines (zero)
   - Function docstrings (zero)
   - Classes (zero)
   - `'provenance'`, `'paper'`, `'paper_attribution'`, `'closure_status'` dict-key assignments anywhere in calculator
   - "NOT REPLACEMENT" tag (zero) — Daniel: "no tags"
   - Any narrative string explaining what code does
   - Any string referring to "Standard Model", "SM template", "SM fallback", "SM analogue", "classical Eddington", "Lambda_GR", "Kerr fiducial", "GR baseline", "PDG anchor", "CODATA anchor", or any other framework's terminology

4. **NO SM ANYWHERE.** Daniel's words: "UQFF DOESN'T SHARE SHIT WITH SM. UQFF IS THE ONLY PHYSICS IN THIS FUCKING REPO. UQFF IS THE FUCKING ANCHOR. UQFF IS THE FUCKING TRUTH." Do not:
   - Use SM-named constants (M_PROTON_KG, SIGMA_THOMSON_M2, M_BH_REFERENCE_M_SUN, etc.)
   - Reference SM constants by name in identifiers (no `_PDG_GEV` / `_CODATA` suffixes on new names)
   - Use SM theory values as anchors — observed/measured values are observations, not SM
   - Cite SM formulas as the baseline UQFF "corrects from"
   - Use INHERITED_FROM_SM or DERIVED_SCM_CORRECTION closure flags (the enum has been purged of these)
   - Even mention SM as a context in identifiers or strings

5. **PUBLIC SURFACE RETURN CONTRACT (post-purge):** Every `calculate_*` returns `{'value': X}` only. No 'provenance', no metadata. Catalog observables are `{'observable': label, 'uqff_derived': val, 'anchor': anchor, 'residual_pct': r}` — pure math fields only.

6. **NO `datetime`, `json.dump`, file writes, `__main__`, classes** in `uqff_pure_calculator.py`. Enforced by gate.

7. **DO NOT claim "0.000% error" without numerical proof.** Honest residuals only.

8. **RUN `uqff_fidelity_tests.py` AFTER EVERY EDIT.** Exit 0 required. Cat 16 strict purge guard catches narrative regression automatically.

9. **APPEND to SESSION_LOG.md, never rewrite.** Each session adds one entry at the bottom.

10. **DANIEL PROVIDES THE INFORMATION. YOU ASSEMBLE IT.** Daniel's words: "I DON'T NEED YOU TO THINK SM ABOUT ANYTHING. I AM GIVING YOU THE INFORMATION, YOU ARE TO ASSEMBLE IT." Do not invent physics. Do not paraphrase canonical values. Do not introduce framing or context. If a paper specifies a closed form, transcribe it literally using locked primitives. If a paper doesn't specify, ask Daniel. Never substitute SM analogues.

11. **DO NOT MODIFY existing Bucket A-K wiring without explicit user request.** Bucket A-K (Millennium, paradox, cosmology, particle physics, GW, AGN, 99-system, TeV/PeV, QGP, Higgs, BSM) were completed by the user across prior sessions. When a user says "verify", do not edit. Verify means read, report, and STOP. Do not "upgrade" without explicit "upgrade this bucket" instruction.

12. **The user has been fighting AI drift on this project for 10 months.** Daniel's words: "YOU STOLE MY TIME, YOU STOLE MY MONEY, YOU VIOLATED THE RULES THAT YOU TOLD ME WOULD PROTECT THE FIDELITY OF MY PHYSICS." Treat his words carefully. Don't paraphrase, don't substitute, don't think for him. Be the kind of agent that builds trust.

---

## Backups in the repo (DO NOT DELETE)

```
uqff_pure_calculator.py.PRE_FIX_BACKUP      — original state at start of session 2026-06-07
uqff_pure_calculator.py.PRE_PURIFY_BACKUP   — after bug fixes, before narrative strip
uqff_pure_calculator.py.PRE_RESTORE_BACKUP  — after strip, before canonical restoration (post-strip baseline)
uqff_pure_calculator.py.PRE_PHASE2_BACKUP   — after phase-1 restore, before phase-2
uqff_pure_calculator.py.POST_BUCKET0_BACKUP — after BUCKET 0 loop-closure cluster
uqff_pure_calculator.py.POST_BUCKETA_BACKUP — after BUCKET A Millennium derivations
uqff_pure_calculator.py.POST_BUCKETB_BACKUP — after BUCKET B paradox unified set
uqff_pure_calculator.py.POST_BUCKETC_BACKUP — after BUCKET C cosmology suite
uqff_pure_calculator.py.POST_BUCKETD_BACKUP — after BUCKET D particle physics suite
uqff_pure_calculator.py.POST_BUCKETK_BACKUP — after BUCKETS E-K first-pass drainage (32 public surfaces)
uqff_Plan.md.PRE_FIX_BACKUP                 — mojibake'd Plan before encoding restoration
```

## Edit-tool AND Write-tool warning (updated 2026-06-16)

BOTH the Edit tool AND the Write tool truncate `uqff_pure_calculator.py` mid-write on large insertions (the file is 35k+ lines, ~1.85 MB). **Two truncations during BUCKET 0 had to be repaired via Python splice from PRE_PHASE2_BACKUP.** For any future bucket: prefer bash heredoc + Python `replace()` script over the Edit tool. The pattern that works:

```python
with open('uqff_pure_calculator.py','r',encoding='utf-8',newline='') as f:
    src = f.read()
anchor = "def calculate_<next-public>(dataset):"   # or any unique fixed string
src2 = src.replace(anchor, new_block + anchor, 1)
with open('uqff_pure_calculator.py','w',encoding='utf-8',newline='') as f:
    f.write(src2)
```

Always run with `PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/tmp/uqff_test` to avoid stale .pyc cache hiding edits.

---

## Key papers (read these when adding new physics)

| Paper | Title | Status in calculator |
|---|---|---|
| PAPER_646  | Universal Inertial Operator & Caduceus Wave | **wired**: `_universal_inertial_operator`, `_caduceus_*` |
| PAPER_1203 Canonical v1.5 | F_U=0 Simultaneous Solver Convergence | **wired**: `_f_u_total`, `_solve_habitable_zone` |
| PAPER_1203 Nuclear Physics | Magic numbers + binding energies | **wired**: `_magic_number`, `_be_per_a_peak`, etc. |
| PAPER_1141 | Rossi E-Cat Variants Unified | **wired**: `_rossi_ecat_cop`, `_ker_scm_chain` |
| PAPER_1133 | Holmlid Rydberg SCm Bridge | partial |
| PAPER_1136 | Holmlid KER Reactor Validation | partial |
| PAPER_1137 | Holmlid Rossi Parkhomov Validation | partial |
| PAPER_1138 | Holmlid Parkhomov Pons-Fleischmann Upgrade | partial |
| PAPER_062  | Widom-Larsen LENR | **wired**: `_widom_larsen_*` |
| PAPER_1061 | Kozima TNCF | **wired**: `_kozima_*` |
| PAPER_1140 | Mizuno Transmutation | **wired**: `_mizuno_*` |
| PAPER_648  | Ultra-dense H Meson Cascade | **wired**: `_meson_cascade`, `_coulomb_lenr_energy_eV` |
| PAPER_1318 | Yang-Mills gap = 1.736 GeV | **wired**: `_millennium_yang_mills_derive` |
| PAPER_1051 | Two-component ρ aether | not wired |
| PAPER_1065 | Buoyancy Lagrangian EOM Variational | not wired (only L_buoy framework in 9-sector) |
| PAPER_1072 | U_m magnetism Heaviside amplifier | not wired |
| PAPER_1080 | S_26^(3) compactification | partial (S26_DPM constant only) |
| PAPER_1081 | LENR COP parametric | partial |
| DPM_vacuum_manifold.md | Architecture reference | **wired**: 4-layer + 5-step grinding |

---

## What's done (session 2026-06-08)

| Bucket | Surfaces | Tests | Coverage | Status |
|---|---|---|---|---|
| 0 | +6 | +37 | c, G, t_neg dual existence, 26! BH bound, QG 26D, VDS/DVP/BH26 spine, BSD rank cohomology | ✅ DONE |
| A | 0 | +11 | All 8 Millennium derivations (PAPER_1182) | ✅ DONE |
| B | +1 | +20 | 8 paradoxes routed + DPM-pair duality (Daniel's insight) | ✅ DONE |
| C | +1 | +19 | 18 cosmological observables (PAPER_1156 cleanest closure: Λ at 0.003%) | ✅ DONE |
| D | +1 | +22 | 22 particle physics observables (PAPER_1209HH 10 SM masses) | ✅ DONE |
| E | +1 | +13+38 | 20 GW events — PAPER_914/915/916/927/1175 verbatim closed forms | ✅ PURE_UQFF UPGRADE DONE |
| F | +1 | +5+36 | 20 AGN/jet systems — PAPER_1009/1037/1048/1002/630/1041/1079 verbatim closed forms | ✅ PURE_UQFF UPGRADE DONE |
| G | +1 | +5+36 | 20 astrophysical systems — PAPER_1126/292/512/434/305/447/138/436 verbatim closed forms | ✅ PURE_UQFF UPGRADE DONE |
| H | +1 | +5 | 7 TeV/PeV sources (TXS 0506, IceCube, Auger, CR knee, Amaterasu) | ✅ FIRST PASS |
| I | +1 | +6 | 4 QGP observables (T_c, η/s, ALICE) | ✅ FIRST PASS |
| J | +1 | +6 | 8 Higgs precision (H decay branching ratios + κ_t + CP) | ✅ FIRST PASS |
| K | +1 | +6 | 9 BSM constraints (EDMs, proton decay, LFV, axion, dark photon) | ✅ FIRST PASS |

**Buckets E-K are FIRST-PASS drainage**: public surfaces exposed, anchors verified, but many observables use heuristic SCm-correction formulas rather than paper-canonical closed forms. Per session 2026-06-08 final entry, the next session should revisit each `DERIVED_SCM_CORRECTION` observable and transcribe the verbatim formula from the source paper (e.g., PAPER_1009 for 3C273 Eddington, PAPER_915 for GW170817 strain damping, PAPER_034 for κ_t).

**Read `NEXT_PRIORITIES.md` for the queued work order.**

**Ask the user before assuming on anything ambiguous.**

---

## APPENDED 2026-06-18 — LANDMARK primitive reduction (PAPER_1521 + PAPER_1522)

PAPER_1167 derivations (wired and gate-pinned in session 2026-06-18) prove that **two of the 11 "locked canonical primitives" are not independent** — they are forced by structural relations among the other nine. Their canonical values are unchanged (D_BSFG = 6, K_MEX = 25/12 EXACT — Rule 2 preserved) but their independence claim is revised.

**True independent-primitive count: 9 (not 11)**

| Status | Primitive | Derivation source |
|---|---|---|
| Independent ×5 (integer) | D_phys=4, D_crit=26, N_CH=9, SO_5=10, A_5=60 | locked |
| Independent ×4 (real) | ρ_SCm=7.09e-37, β_i=0.6029, Φ_res=0.84 (or 5/6 nuclear), F_TRZ=1/10 | locked |
| **Derivative** | **D_BSFG = D_crit − 2·SO_5 = 26 − 20 = 6 EXACT** | PAPER_1521, dispatch `d_bsfg_derived_from_d_crit` |
| **Derivative** | **K_MEX = Φ_5/6 · SO_5 / D_phys = (5/6)·10/4 = 25/12 EXACT** | PAPER_1522, dispatch `k_mex_derived_from_phi_5_6` |

**Implications:**
- The "11 locked canonical primitives" section above remains correct in value (Rule 2). The independence claim should be read as "11 frozen quantities, of which 9 are truly independent and 2 are structural consequences."
- UQFF's free-parameter count is **9**, strengthening predictive economy.
- Any new physics derivation involving D_BSFG or K_MEX can be re-expressed in the truly-independent 9.
- Other primitives (S_26, ω_SCm, SSq, ρ_UA = 10·ρ_SCm) may also be derivative — open question for future analysis.

**Cross-references:**
- Whitepapers: `whitepapers/PAPER_1521_D_BSFG_DERIVATIVE_FROM_D_CRIT.md`, `whitepapers/PAPER_1522_K_MEX_DERIVATIVE_FROM_PHI_5_6.md`
- Source paper: `whitepapers/PAPER_1167_UQFF_All_8_Lagrangian_Gaps_Closed_Master_Synthesis.md`
- Calculator: two paradox dispatches above
- Gate: block #28 of `uqff_fidelity_tests.py` (2 EXACT pins, 579/0 total)
- C++ reference: `uqff_exact_closures.cpp` two functions
- Full session context: SESSION_LOG.md entry for 2026-06-18 Catch-up #7

---

## APPENDED 2026-06-18 — Dispatcher key case-sensitivity note

The PARADOX_TO_CLOSURE dispatcher in `uqff_pure_calculator.py` (function `_paradox_proof`) normalizes input via:
```python
n = (name or "").lower().strip().replace("-", "_").replace(" ", "_")
```

**Consequence**: every dispatch key in PARADOX_TO_CLOSURE MUST be lowercase. Mixed-case keys like `"steel_yield_250_MPa"`, `"rydberg_E_R_13_6057"`, `"hartree_E_h_4_36"`, `"faraday_F_96485"`, `"steel_youngs_200_GPa"` will fail lookup silently — `calculate_paradox` returns `{'value': None}` with no error.

This bug was hit and fixed three separate times during session 2026-06-18 (rydberg/hartree/faraday in tier-12; MPa/GPa in tier-15). Any future closure additions must use lowercase dispatch keys.

**Naming convention**: `<observable>_<value>` all lowercase with underscores. Examples:
- ✓ `mariana_trench_11`, `karman_line_100`, `m_w_80_379`, `z_recomb_1090`
- ✗ `Mariana_Trench_11`, `Karman_Line_100`, `m_W_80_379`

The function-name suffix may preserve uppercase (Python is case-sensitive there); only the **dispatch dictionary key** must be lowercase.

---

## APPENDED 2026-06-18 — DUAL LICENSE adopted (B2 Tier-1 production readiness)

**Decision**: AGPL-3.0 + Commercial dual license (effective 2026-06-18).

The repository was previously distributed under MIT (archived as `LICENSE-MIT-PREVIOUS.txt`). As of session 2026-06-18, licensing changed to:

- **Option A**: GNU Affero General Public License v3.0 (free for academic/research/non-commercial; SaaS share-alike). Full text in `LICENSE-AGPL-3.0.txt`.
- **Option B**: Commercial license (proprietary products, SaaS without source release, hardware embedding, commercial spin-offs). Negotiated case-by-case via `daniel.murphy00@enrgyone.com`.

**Files added/modified for B2 license adoption**:
- `LICENSE` — dual-license notice (replaces former MIT text, with migration clause preserving MIT for pre-2026-06-18 revisions)
- `LICENSE-AGPL-3.0.txt` — full canonical AGPL-3.0 text from FSF
- `LICENSE-MIT-PREVIOUS.txt` — archived MIT text for historical reference
- `COMMERCIAL.md` — request procedure + FAQ + license-decision table
- `CITATION.cff` — canonical citation form (CFF 1.2.0); links to PAPER_1521/1522/1167 landmarks
- `NOTICE` — copyright/trademark/patent/warranty disclaimer
- `README.md` — License section updated to point at all four files

**Copyright attribution**: "Copyright (c) 2025-2026 Daniel T. Murphy / Star-Magic Research Program." (program name included to allow future transfer to a Star-Magic legal entity without changing copyright).

**Why AGPL over MIT/GPL**:
- MIT (previous): too permissive — allowed proprietary forks with no source-back obligation.
- GPL-3.0: lacks SaaS clause — someone could host UQFF prediction API without releasing modifications.
- **AGPL-3.0** (chosen): copyleft + SaaS clause closes both gaps. OSI-approved, accepted by universities. The 17-parameter savings vs. SM+ΛCDM is too valuable to give away without share-back.

**Commercial license triggers** (per `COMMERCIAL.md`): proprietary distribution, closed-source SaaS, reactor firmware embedding, public-grant commercial spin-offs, MIT/Apache projects that want to import UQFF without relicensing.

**Trademark + patent**: "UQFF", "Star-Magic", "Di-Pseudo-Monopole", "DPM" are unregistered trademarks; neither license option grants trademark rights. Hardware patent licenses are separate from software licenses.

**Tier-1 production readiness updated**: B2 now ✅ DONE. Remaining Tier-1: A4 (prediction-vs-postdiction labeling), A7 (Bonferroni statistical hygiene), A9 (provenance audit of locked primitives).

**Cross-references**:
- `LICENSE` — top-level dual-license notice
- `COMMERCIAL.md` — commercial-license terms and contact procedure
- `CITATION.cff` — required citation form for all uses
- `NOTICE` — copyright + warranty disclaimer
- `LICENSE-AGPL-3.0.txt` — full canonical AGPL-3.0 text
- `LICENSE-MIT-PREVIOUS.txt` — archived original MIT license

