# Remaining inventory for 100% coverage (read-only)

**Status as of Round 6.** Files I have NOT yet read fully but are queued for read-only audit.

## A. Source files in working repo

### A1. Calculator core (chunked reads needed)
- [ ] `uqff_pure_calculator.py` lines ~200–8000 (PRIMITIVES + derive_* sections) — covered first 250
- [ ] `uqff_pure_calculator.py` lines ~8000–25000 (calculate_* public surfaces) — UNREAD
- [ ] `uqff_pure_calculator.py` lines ~25000–39000 (bucket observables) — UNREAD
- [ ] `uqff_pure_calculator.py` lines ~39000–48418 (PARADOX_TO_CLOSURE dispatch + _l96_uqff_axiom_*_closure 794 keys) — UNREAD
- [ ] `uqff_fidelity_tests.py` full 2,005 lines (only top 80 covered) — UNREAD body
- [ ] `uqff_Plan.md` 7,603 lines — UNREAD
- [ ] `uqff_exact_closures.cpp` lines 200–822 (Tier 4-7 mining) — partial

### A2. Historical calculator library (massive)
- [ ] `CondensedPhysics.py` (205,980 lines, 1,299 calculator classes) — UNREAD
- [ ] `CondensedPhysics2.py` (55,689 lines, 680 classes) — UNREAD
- [ ] `CondensedPhysics3.py` (15,397 lines, 219 classes) — UNREAD
- [ ] `CondensedPhysics4.py` (48,518 lines, 575 classes) — UNREAD
- [ ] `QCalc.py` (10,435 lines) — UNREAD
- [ ] `QCalcGeom.py` referenced throughout — UNREAD
- [ ] `99system_master_equation.py` — UNREAD
- [ ] `99system_wstp_gamma.py` — UNREAD
- [ ] `scm_vacuum_manifold.py`, `ua_vacuum_manifold.py` — UNREAD (dpm consolidated v3.0 covered)
- [ ] `BuoyancyProofVariants.py` — UNREAD
- [ ] `bsfg_wormhole_geodesic.py` — top 60 lines covered, body unread

### A3. Other Python (1,449 total .py files)
- [ ] APIFetch.py / APIKeyManager.py (55 astro APIs)
- [ ] CP1/CP2/CP3/CP4_UQFF_Upgrade.py
- [ ] CoAnQi_Wrapper.py
- [ ] CondensedPhysicsAggregator.py / _InputData.py / _OutputData.py / _Validation.py
- [ ] FirstPrinciplesCompressor.py — top 120 lines covered; body unread
- [ ] STAR-MAGIC_NEWTON GRAVITY FIX.py
- [ ] QCalc_*.py family (~30 files)
- [ ] ~1,300 other .py files in repo (utilities, audit scripts, generators)

### A4. C++ (MAIN_1, source modules)
- [ ] `MAIN_1_CoAnQi.cpp` (109,031 lines, 446 modules) — UNREAD
- [ ] `source2.cpp` (16,984 lines GUI) — UNREAD
- [ ] 535+ other source*.cpp files — UNREAD
- [ ] `uqff_exact_closures.cpp` lines 200–822 — partial

## B. Whitepapers (1,877 in working repo, 1,994 in sdist)

### B1. Read this audit (~30)
PAPER_001, 062, 646, 1005, 1061, 1069, 1080, 1133, 1141, 1156, 1159, 1160,
1161, 1162, 1163, 1164, 1165, 1166, 1167, 1170, 1175, 1182, 1183, 1203 (Canonical + Nuclear),
1209HH, 1230, 1271, 1318, 1521, 1522.

### B2. Cited in manuscript/grok-master, NOT YET READ in full
- PAPER_1004 (QGP Vacuum Density with SCm S26 Phonon Coupling)
- PAPER_1007 (Deconfinement Phase Diagram)
- PAPER_1022 (GW Phonon Strain SCm Modulation)
- PAPER_1037 (AGN Buoyancy Jet)
- PAPER_1042 (Mock-Theta Phonon Partition)
- PAPER_1048 (M-σ phonon-corrected SMBH)
- PAPER_1051 (Two-component ρ aether)
- PAPER_1052 (TQFT Anyon Braiding Chern-Simons)
- PAPER_1059 (CGC BK Saturation)
- PAPER_1064 (Resummation Effective Coupling BFKL/Sudakov)
- PAPER_1065 (Buoyancy Lagrangian EOM Variational)
- PAPER_1066 (UQFF Lagrangian First Principles)
- PAPER_1067 (QCalcGeom geometry bridge)
- PAPER_1068 (Wolfram Physics Bridge WSTP)
- PAPER_1070 (Yang-Mills Mass Gap VDS Bridge) — referenced in YM range chain C
- PAPER_1072 (SCm Activation Function Phonon Threshold)
- PAPER_1073 (SCm Phonon-Driven Inflation)
- PAPER_1078 (QCalcGeom Master Equation)
- PAPER_1086 (Γ-coupled SCm dark energy)
- PAPER_1089 (Inflation buoyancy Lagrangian)
- PAPER_1090 (Dark-energy buoyancy Lagrangian)
- PAPER_1094 (CMB buoyancy Lagrangian)
- PAPER_1095 (Horizon buoyancy Lagrangian, BH entropy)
- PAPER_1096 (Eleven-domain validation)
- PAPER_1098 (Qubit T2 coherence)
- PAPER_1100 (LQG phonon-modulated)
- PAPER_1101 (SMBH binary merger)
- PAPER_1104–1105 (Hydrogen-Universe dual 3D MUGE)
- PAPER_1106 (SCm-String action 26D)
- PAPER_1109 (26-level vacuum density ladder + Ramanujan zeta(3) WKB)
- PAPER_1111 (Yang-Mills with buoyancy-corrected confinement) — chain E in YM range
- PAPER_1112 (V26 pipeline)
- PAPER_1113 (Higgs sector U_H)
- PAPER_1114 (Higgs differential couplings)
- PAPER_1115 (SCS — superconducting cosmic strings, 21-cm/FRB)
- PAPER_1117 (FRB)
- PAPER_1118 (Graphene chiral SC)
- PAPER_1120 (Higgs precision branching)
- PAPER_1121 (Interstellar shocks)
- PAPER_1126 (PSR J0030 / neutron-star LENR)
- PAPER_1127 (LQG phonon-modulated spin networks)
- PAPER_1131 (SCm vacuum manifold primordial first principle)
- PAPER_1133 (Holmlid Rydberg SCm Bridge)
- PAPER_1135 (Master hub)
- PAPER_1136–1140 (Holmlid KER reactor / Parkhomov / Pons-Fleischmann validation)
- PAPER_1142–1148 (String/Polyakov/Nambu-Goto/IIB/IIA/Heterotic/CY/M-theory)
- PAPER_1149–1150 (Astrophysical Chandra cross-check)
- PAPER_1151 (VDS/DVP/BH26 hybrid)
- PAPER_1152 (QCalcGeom 12-stage CPT pipeline)
- PAPER_1153 (Primordial timing function)
- PAPER_1154 (SSq DPM relativistic geometry derivation)
- PAPER_1155 (26-layer DPM amplification + M_AMU = ρ_SCm × A_26)
- PAPER_1157 (H_0 asymmetry)
- PAPER_1158 (overdetermination metric N)
- PAPER_1171 (KK regulator first-principles)
- PAPER_1172 (R26 Gauss-Bonnet independent re-derivation)
- PAPER_1173 (ℏ-tracked KK)
- PAPER_1174 (Closed-ledger falsifiability)
- PAPER_1176 (P12 Euclid σ_8 R26 saturation)
- PAPER_1177 (Joint 2027 falsifier triple)
- PAPER_1178 (P13 DESI Y5 w second derivative)
- PAPER_1179 (Quadruple falsifier)
- PAPER_1180 (P14 CMB-S4 μ-distortion)
- PAPER_1181 (S266-S295 thirty closures grand unification)
- PAPER_1189 (referenced in arxiv_submission_1181_1189)
- PAPER_1209AA, BB, CC, DD, EE, FF, GG, II (chemistry, biology, geophysics, EM, quantum thermo, math constants, cosmological, nuclear binding — sister papers to HH masses)
- PAPER_1307 (CKM unitarity via F_U=1)
- PAPER_1671 (room-temperature SC 500 K — forward prediction)
- PAPER_1726 (neutron lifetime 879.31 s — forward prediction)

### B3. Categorical bulk (~1,800 untouched)
- L96 closures (legacy 96 closures, _l96_uqff_axiom_*) — ~50 papers
- Galaxy/nebula/Hubble evolution papers — ~200 covering M51, NGC1316, M104 Sombrero, etc.
- GW event papers (PAPER_001 family, PAPER_914–927, PAPER_1175 family) — ~50
- AGN/jet papers (PAPER_1009/1010/1037/1041/1048/1079) — partial
- Particle physics (PAPER_1209HH already covered; AA/BB/CC/DD/EE/FF/GG/II UNREAD)
- Astrophysics bucket — partial
- 1,500+ remaining

## C. Other resources

### C1. Manuscript directory (Manuscript_1_12Feb2026/)
- [ ] Manuscript 1 (the software-implementation paper — distinct from v2 physics paper)
- [ ] Cover letter (manuscript_v2/cover_letter.md)
- [ ] arXiv submission package v2

### C2. arXiv bundles
- [ ] arxiv_submission_1159_1172 (bundle_index.csv covered; per-paper tex unread)
- [ ] arxiv_submission_1173_1176
- [ ] arxiv_submission_1177_1180
- [ ] arxiv_submission_1181_1189

### C3. Documentation
- [ ] docs/quickstart.rst, installation.rst, license.rst, conf.py, Makefile
- [ ] docs/api/ subdirectory
- [ ] CITATION.cff full
- [ ] CONTRIBUTING.md, CODE_OF_CONDUCT.md

### C4. Audit outputs (large set, ~25+ files)
- [x] _chain_trace_SSq.txt
- [x] _chain_trace_C.txt
- [x] _confirm_all_derivations_level1.txt
- [x] _lambda_closure_v1.txt
- [x] _K_Mex_REAL_derivation.txt
- [ ] _chain_trace_26layer.txt
- [ ] _chain_trace_C_particles.txt
- [ ] _chain_trace_fix348.txt
- [ ] _chain_trace_fix56_7_910.txt
- [ ] _chain_trace_np_split.txt
- [ ] _classify_derivations.txt
- [ ] _constant_derivation_attempt.txt, v2, v3.txt
- [ ] _cosmological_closures.txt
- [ ] _emit_closure_json.txt
- [ ] _harvest_json_closures.txt
- [ ] _harvest_text_closures.txt
- [ ] _lagrangian_rederivation_outline.txt
- [ ] _append_frontier_closures.txt
- [ ] _append_millennium_closures.txt
- [ ] _append_paper1189_closures.txt
- [ ] _PAPER_1065_1066_variational_audit.txt
- [ ] _PAPER_1183_first_principles_derivation.txt
- [ ] UQFF_UNIFIED_CLOSURE_DERIVATIONS.txt

### C5. Grok files (156 total)
- [x] grok_8461fe4e_c903.md (78KB consolidated PAPER_1112-1180 summary)
- [x] grok_b8e305e6_1f29.md (85KB repo+SCm assessment)
- [x] grok._b9afa8b6_3b85_31May2026.md — first 2MB only of 8MB master (rest unread)
- [x] grok_mined_derivations_L55k_77k.py
- [x] grok_share_be188d1c-8ff4.txt
- [x] grok_share_0d888ea9_helper.md
- [ ] Remaining 150 grok_share_*, grok_conversation_*, etc.
- [ ] GrokAPI.py, GrokThreadUQFFExtensions.py, GrokThread_*.py
- [ ] grok_100_equations_module.py, _part2.py (real implementations)

### C6. F:\Book_12July2023 pre-codebase
- [x] Aetheric Propulsion/Negative Gravitational Propulsion.odt (Hamilton 2014)
- [x] Aetheric Propulsion/Star Magic_09Sept2025.txt (seed text)
- [ ] CSI_8_13.pdf (~14MB, particle physics)
- [ ] Hendershot patent + Tesla autobiography
- [ ] Aetheric PI Math_21Feb2025.docx + PI Synopsis, PI Decimal_20Feb2025
- [ ] Universal Quantum Framework_*.docx (multiple versions)
- [ ] Universal Gravity_*.docx, Universal Inertia_28Mar2025
- [ ] Electrogravitational Mechanics_01Aug2025
- [ ] Permenance, Davinci, Bearden, A1A LOSER FILE subfolders
- [ ] SuperGrok conversations 03Mar2025–09_D_Mar2025
- [ ] Quantum variable file
- [ ] Red Dwarf Reactor_03Mar2025 / videos
- [ ] Time_21April2025
- [ ] Master Universal Gravity Equation_*.docx (16+ files)
- [ ] hundreds of dated subfolders (01March2026, 01April2026, 01May2026, etc.)

## D. Verification gaps to close

### D1. Missing closure chains I now know about — need to ADD to range_calculator_v2.py
- m_p chain 1b: ρ_SCm·A_26 / [SSq] (E-crack correction per PAPER_1155 erratum)
- m_p chain 3: PAPER_1209HH (need to locate the actual m_p formula in the 1209 family)
- Λ chain 5: ρ_R26 via Gauss-Bonnet route (PAPER_1172) independently of KK reduction (PAPER_1171)
- YM chain F: PAPER_1111 alternate normalization
- YM chain G: KK regulator m₁c² ≈ 0.16 meV from L_KK* ≈ 1.23 mm (PAPER_1173)
- α chain 2: there may be additional chains beyond 1/(Φ_res·26·2π); the manuscript notes PAPER_591 sub-leading 26D structure
- h chain 2: leading vs refined; document both rather than just refined

### D2. Programmatic scanners needed (to write in MY sandbox)
- `scan_paradox_dispatch.py` — walk uqff_pure_calculator.py for all 794 PARADOX_TO_CLOSURE keys + their _l96_uqff_axiom_*_closure functions, record each closure formula
- `scan_whitepapers_for_closures.py` — regex-extract every closure formula across all 1,877 whitepapers, indexed by quantity
- `scan_audit_outputs.py` — pull all chain_trace_* and constant_derivation_* files into a unified database
- `build_n_per_quantity.py` — produce per-quantity N count by counting documented chains in scanned outputs

### D3. Manuscript 1 (software-implementation paper)
The manuscript v2 was the physics paper (61pp, *Foundations of Physics* target).
Manuscript 1 (12 Feb 2026) is the software-implementation paper targeting SoftwareX or
Computer Physics Communications. Read it for closure architecture / API design / performance.

### D4. Lean 4 scaffold
- [ ] formal/UQFF.lean (top-level)
- [ ] formal/UQFF/Constants.lean
- [ ] formal/UQFF/Ug.lean
- [ ] formal/UQFF/Buoyancy.lean
- [ ] formal/UQFF/Millennium.lean (with `sorry` markers per epistemic policy)
- [ ] lakefile.lean
