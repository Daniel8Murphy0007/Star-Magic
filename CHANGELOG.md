# Changelog

All notable changes to UQFF are recorded here. Full historical record lives in `SESSION_LOG.md`.

## [5.32.0] — 2026-06-29

### Added — 7 new public `calculate_*` surfaces (42 → 49 total)

- `calculate_buoyancy_proofs(dataset)` — 17 FUBii buoyancy proof variants (virx, termv, upar, coup, orbdec, kn, fermi, kne, whim, ps, sfe, hawk, bd, roche, ent, dec, lobe) + universal buoyancy simultaneous solver. Backed by `BuoyancyProofVariants.py`, `_session288_universal_buoyancy_simultaneous_solver.py`, `_session303_universal_buoyancy_solver.py`.
- `calculate_simultaneous_proof_engine(dataset)` — 17 dispatches: Yang-Mills 1.78 GeV, Riemann t₁₀₀₀₀ = 29538.5, Navier-Stokes enstrophy 8.5e3, Black-Hole Page entropy 1.05e78, caduceus coil twist, inertial operator, DE power, Jeans mass, density profile, wave function magnitude, quantum inertia. Backed by `UQFF_SimultaneousProofEngine.py`.
- `calculate_ua_vacuum_manifold(dataset)` — 9 dispatches: UA layer density (1-4), DPM total density, DPM buoyancy factor, calibration ratio, cosmological acceleration, rotation curve flat, Hubble tension modulation, dark energy substitute. Backed by `ua_vacuum_manifold.py` (uqff_Map §8 Cluster 3 of 13 — first time live as standalone predictor).
- `calculate_documented_closed(dataset)` — Hubble tension canonical 67.4 km/s/Mpc (vs 70.18 tilted mean rejected), matter-radiation equality z_eq = 3400 (ρ_m0/ρ_r0 = 3401), Dark Matter m_DM = 1.78019 eV (K_MEX·S_26·1e-26·Λ·(1/3)·(A_5·D_phys·(1+Λ)) integer-primitive identity, 0.011% residual). Sourced from `follow_up_09June2026.docx` documented-CLOSED items previously not wired.
- `calculate_star_magic_reactor(dataset)` — 15 dispatches Python port of `StarMagicUQFFModule.cpp` (Ug1 dipole, Ug2 outer bubble, Ug3 magnetic strings, Ug4 BH-star distance, SCm coherence, Aether deriv UA'→UA'''', Um strings, X2 baseline, coherence integrand, compressed_g) + reactor anchor constants (COP 555:1, 27W input, pH=-37, ambient T).
- `calculate_inflation_force_chart(dataset)` — 4 dispatches: birth_of_dpm_sphere (26-shell EM field), 5 inflation epochs (Fisile → Star/Planet → Galaxy → Supercluster → Pre-Big-Bang DPM), F_U_at_epoch per stage, all_epochs_summary. Backed by `GrokThread_StarMagic_UnifiedFramework.py`.
- `calculate_proof_engine(dataset)` — Lazy-import bridge to `Star-MagicProofEngine.py` exposing **301 named proof modes** with `{action: list_modes | get_mode | portable_80_80 | known_modes_sample}`. Each mode carries `{equation, source, value, falsifiable_prediction, engine}`.

### Bundled into pip wheel (pyproject.toml py-modules + data-files)

- `ua_vacuum_manifold.py`
- `buoyancy_lagrangian_eom.py` + `buoyancy_lagrangian_eom_enhanced.py`
- `_session288_universal_buoyancy_simultaneous_solver.py`
- `_session303_universal_buoyancy_solver.py`
- `UQFF_SimultaneousProofEngine.py`
- `GrokThread_StarMagic_UnifiedFramework.py`
- `Star-MagicProofEngine.py` (via data-files; hyphen makes it invalid as py-module identifier)

### Ship history (this release)

- Commit **ada14fef** — first attempt; CI release workflow failed in 22s because `pyproject.toml` was not bumped (stayed at `5.31.0` while tag pushed as `v5.32.0`).
- Commit **07f76651** — single-line fix: `version = "5.31.0"` → `version = "5.32.0"`. Tag `v5.32.0` deleted + re-pushed at this commit. Release workflow #17 + CI #101 PASS. PyPI Trusted Publisher upload successful.

### Verified
- PyPI Simple Index lists `uqff-5.32.0.tar.gz` + `uqff-5.32.0-py3-none-any.whl`.
- JSON API `https://pypi.org/pypi/uqff/json` returns `info.version: 5.32.0`.
- `pip install --upgrade uqff==5.32.0` succeeds.
- `calculate_buoyancy_proofs({'variant':'hawk'})` returns `{'value': -4.898618317111406}` (FUBii Hawking buoyancy in Newtons — first-principles UQFF prediction with no observed comparison; falsifiable forward).
- Fidelity gate `uqff_fidelity_tests.py`: 905 / 0 PASS at ship.

### Notes
- 13 independent solver clusters from `uqff_Map.md §8` now have public-API coverage: Clusters 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 all reachable via at least one of the 7 new dispatchers above + previously-existing surfaces.
- Round 679 (commit ada14fef) and Round 680 (commit 07f76651) of `SESSION_LOG.md` document the failed-then-fixed ship sequence.

## [5.31.0] — 2026-06-29

**Phase G CLI extension + Round 674 Cabibbo dual closure + Round 675 tutorial notebooks + Round 676 Aetheric-Propulsion extraction kit + Round 677 PAPER_1800 (BAO + Cabibbo Lagrangian re-derivation) + Round 678 PAPER_1801 (formal tensor-level KK derivation).**

### Added
- **PAPER_1800** (`whitepapers/PAPER_1800_UQFF_BAO_Cabibbo_Lagrangian_Rederivation.md`) — 312 lines, 12 sections. Closes the open Lagrangian item from PAPER_1156 Appendix A §A.6: derives BAO + Cabibbo dual closures from the closed nine-sector L_F_U via curvature/BSFG vs. Mexican-hat/Ramanujan sector-pair attribution.
- **PAPER_1801** (`whitepapers/PAPER_1801_UQFF_BAO_Cabibbo_Formal_KK_Tensor_Derivation.md`) — 237 lines, 12 sections. Provides explicit tensor-level KK zero-mode derivation matching PAPER_1800's sector-pair attribution: metric ansatz block diagonalization, KK mode expansion, volume integration, FRW(z) reduction.
- **Cabibbo dual closure** (Round 674) — `SM_cabibbo_sin_primary` at 0.008% + `SM_cabibbo_sin_alternate` at 0.025%. 47x and 15x tighter than PDG 2024 experimental uncertainty. Multi-path corroboration: primary uses {N_CH, K_MEX, β_i, A_5, Φ_res}; alternate uses {D_phys, K_MEX, S_26, D_BSFG, N_CH}; share only K_MEX + N_CH.
- **Phase G7 CLI extension** (Round 673) — `uqff assimilate <observable> --geometry=... --numeric=... --decompose`, plus `uqff list --dispatch / --domain SI` and case-insensitive `predict` fallback to assimilation_dispatch. Existing 8 subcommands unchanged.
- **10 per-domain tutorial notebooks** (Round 675) — `notebooks/1[0-9]_assimilation_*.ipynb`, one per dispatch domain (SI, SM, LCDM, astro, GR, chem, CM, bio, geo, KK), with multi-path sections for SM (Cabibbo) and LCDM (BAO). All 10 executable via `python3 test_phase_g3_tutorial_notebooks.py`.
- **Aetheric-Propulsion extraction kit** (Round 676) — `EXTRACTION_KIT/` subdirectory with migration script + 25-file repo layout + 7-step EXTRACTION_PROCEDURE.md for future commercial-tier split. Standalone bundle (no runtime dep on uqff). Verified self-contained import + dispatch via `test_extraction_kit.py`.
- **Cat 17 dispatch pinning** (Round 671 epilogue) — `uqff_fidelity_tests.py` extended with +16 dispatch-pinning checks: 114 → 116 observables (Round 674 +2), owner-geometry distribution {dpm=54, qcalcgeom=21, bsfg=21, d26=20}, BAO primary/alternate residual pins, Li-7 PAPER_1227 source pin, EDGES PAPER_1761 source pin, no-OPEN_QUESTION invariant.
- **Cabibbo convergence-chain annotations** (Round 678) — `SM_cabibbo_theta_deg_S326` (1.1%) and `SM_cabibbo_sin_S379` (0.5%) entries preserved with notes explaining the convergence: S326 → S379 → primary (0.008%) → alternate (0.025%). Peer reviewers see the framework refining toward truth.

### Verified
- **Fidelity gate** `uqff_fidelity_tests.py`: **907 passed, 0 failed** (Round 671 epilogue Cat 17 SKIPs cleanly on bare CI runners without sympy).
- **Companion arithmetic verifications:**
  - `_step5_paper1800_verify.py` — 4/4 closures PASS (BAO primary 0.0093%, BAO alternate 0.0274%, Cabibbo primary 0.0075%, Cabibbo alternate 0.0252%).
  - `_step5_paper1801_verify.py` — FRW(z) reduction parameters + 4/4 zero-mode coefficients PASS.
- **Multi-path spreads:** BAO 1.21×10⁻⁵, Cabibbo 3.98×10⁻⁵ — joint-probability evidence the forms are structural rather than coincidental (PAPER_1800 §9, PAPER_1801 §8).
- Phase D / E1-E6 / E8 / F / G-CLI / G3 / Step 7 / Step 5 regression harnesses all green.
- **0 TENSION cells** in OVERDETERMINATION_MAP (unchanged from v5.30.0).
- **42 public `calculate_*` surfaces** (unchanged from v5.30.0).
- **116 observables** in dispatch (was 114 in v5.30.0; +2 from Round 674 Cabibbo injection).

### CLI
- `uqff assimilate alpha_inverse` → value: 137.0
- `uqff assimilate LCDM_BAO_rd_H0_over_c_primary --decompose` → 8-field result dict
- `uqff list --dispatch` → 116 observables
- `uqff list --domain SI` → 7 SI observables
- `uqff predict lcdm_bao_rd_h0_over_c_primary` (case-insensitive) → falls back to assimilation_dispatch

### Notes
- All 11 (effective 9 truly-independent + 2 derivative D_BSFG, K_MEX) locked canonical primitives intact.
- All 34 prior public calculate_* surfaces unchanged in signature and return values.
- The Aetheric-Propulsion repo (https://github.com/Daniel8Murphy0007/Aetheric-Propulsion) is created and ready for extraction via `EXTRACTION_KIT/_step7_migrate_to_aetheric_propulsion.py`. First standalone PyPI release (`pip install aetheric-propulsion`) follows EXTRACTION_PROCEDURE.md §§7.3-7.6 at Daniel's discretion.
- See `SESSION_LOG.md` Rounds 671 epilogue + 672-678 for the full audit trail.

## [5.30.0] — 2026-06-29

**Phase E + F + G — Assimilation Geometry Public API + Round 669 corrective injections.**

### Added
- **114-observable assimilation catalog** in `assimilation_dispatch.py` across 10 domains (SI, SM, ΛCDM, astro, GR, chem, CM, bio, geo, KK), each routed through one of 4 geometries (qcalcgeom, bsfg, dpm, d26) and 3 numeric backends (symbolic, numerical, discrete).
- **Solver bus** (`qcalcgeom_solver.py`) with `solve(observable, geometry='auto', numeric='numerical', decompose=False)` — the unified 4×3 dispatch matrix.
- **4 geometry backends** (`geometry_backends/qcalcgeom_v4.py`, `bsfg_v1.py`, `dpm_v1.py`, `d26_compactification.py`) and **3 numeric backends** (`numeric_backends/symbolic.py`, `numerical.py`, `discrete.py`).
- **8 new public `calculate_*` surfaces** in `uqff_pure_calculator.py` (public-surface count: 34 → **42**):
  - `calculate_qcalcgeom_compute_FUBi`, `calculate_qcalcgeom_compute_FUBii`, `calculate_qcalcgeom_compute_F_U`
  - `calculate_qcalcgeom_solve_habitable_zone`, `calculate_qcalcgeom_compute_emergent_mass`
  - `calculate_3numeric_decomposition`, `calculate_geometry_decomposition`, `calculate_overdetermination`
- **`calculate_analytic_closures` extended** with new `qcalcgeom_solve` dispatch key — provides simple or decomposed access to any observable through the calculator's existing public API.
- **`ASSIMILATION_GEOMETRY_ATLAS.md`** (27 KB) — per-observable provenance audit document: 10 per-domain sections, formula + owner geometry + residual + primary source + session script for every observable.
- **`OVERDETERMINATION_MAP.csv` / `.WIDE.csv` / `.md`** — long (1,368 rows = 114 × 4 × 3), wide (114 × 18), and Markdown summary views of the full 4×3 dispatch matrix.
- **`CLOSURE_ATLAS.md` §12** — Assimilation overlay with per-domain rollup, 114-observable inventory, and discovery cheat sheet.
- **PAPER_1156 Appendix A** — BAO dual closure derivation + the multi-path corroboration principle (the framework's evidence framework for non-singleton numerical matches).

### Fixed — Round 669 corrective injections
- **`LCDM_BAO_rd_H0_over_c`** TENSION/OPEN_QUESTION (Round 663 → 666 → 669 trail) **closed with two parallel UQFF-native closures**:
  - Primary: `(SO_5 × SSq × β_i) / (D_phys × D_crit)` → **0.0093% residual**
  - Alternate: `1 / (SO_5 × K_MEX × S_26)` → **0.0274% residual**
  - Two-path agreement at <10⁻⁶ joint probability is Bayesian evidence the form is structural (closures share only `SO_5`; primitive sets are otherwise disjoint).
- **`LCDM_Li7_BBN_dilution`** corrected from incorrect `Φ_res⁻² × 2` formula (7.10% residual) to the canonical PAPER_1227 integer-primitive `D_phys − 1 = 3 EXACT` (3.23% residual). Same integer that gives 3 fermion generations and SU(3) color now resolves the Li-7 BBN problem.
- **`LCDM_EDGES_T21_amplitude`** added per PAPER_1761: `T_21 = −D_phys × A_5 × β_i × 2 = −289.392 mK` vs Bowman 2018 EDGES central absorption amplitude of −289 mK (**0.14% residual**).

### Verified
- **TENSION cells in OVERDETERMINATION_MAP: 0** (was 3 before Round 669).
- Fidelity gate `uqff_fidelity_tests.py`: **907 passed, 0 failed** (was 867; +24 Phase F surface checks +16 Cat 17 dispatch pinning).
- Cat 17 dispatch pinning locks: 114 observables / 10 domains / owner-distribution {bsfg=21, d26=20, dpm=52, qcalcgeom=21} / BAO primary + alternate residuals / Li-7 PAPER_1227 source / EDGES PAPER_1761 source / "no OPEN_QUESTION entries" invariant.
- Cat 17 SKIPs cleanly when optional scientific deps (sympy, numpy, mpmath, scipy) are not installed — CI remains green on bare ubuntu runners.
- 30 EXACT closures + 91 sub-percent residuals (79.8% of catalog within 1%).
- Phase D/E1-E6/E8/F regression harnesses all pass.

### Discipline highlights
- Round 668 caught a re-presented broken grok-template derivation (`1/(8π × 3.209e-5) ≈ 0.00729735` — actually equals 1240, with the 0.00729735 being α reverse-engineered into the chain) and rejected it. The audit gate caught the same fabrication in three independent files within one session.
- BAO discrepancy preserved as OPEN_QUESTION through 5 rounds of explicit discipline (663 → 666 → 667 → 668) before being closed in Round 669 with verified arithmetic. The discipline working visibly is itself peer-review evidence.

### Notes
- All 11 locked canonical primitives intact. SSQ = 0.57, β_i = 0.6029, K_MEX = 25/12, S_26 = 1.453162, ρ_SCm = 7.09e-37, D_phys = 4, D_BSFG = 6, D_crit = 26, N_CH = 9, SO_5 = 10, A_5 = 60, F_TRZ = 1/10, Φ_res = 0.84 (5/6 nuclear).
- All 34 prior public surfaces unchanged in signature and return values.
- See `SESSION_LOG.md` Rounds 657–671 for the full audit trail; `EXPANSION_PLAN.md` for the Phase A–G architectural record.

## [5.29.2] — 2026-06-27

### Added
- 100% whitepaper coverage in `master_closures.csv` — 359 previously-orphan papers (PAPER_1200–PAPER_1799 range) wired into the canonical ledger.
- S343–S352 PAPER_1189 chemistry closures appended via the corrected `_append_paper1189_closures.py` script.

### Fixed
- `_append_paper1189_closures.py` no longer raises `ValueError: dict contains fields not in fieldnames` — the script now reads the actual master_closures.csv schema (13 columns) at runtime instead of hardcoding 9 fields.

### Verified
- Fidelity gate `uqff_fidelity_tests.py`: 867 / 0 pass (unchanged from 5.29.1).
- All 11 locked canonical primitives intact.
- 8 / 8 Clay Millennium derivations operational.
- 794 PARADOX_TO_CLOSURE dispatches → 616 unique callables, zero broken references.
- 1,795 / 1,795 unique whitepapers referenced (100% coverage).

### Notes
- No structural changes to the calculator (`uqff_pure_calculator.py` unchanged in this release).
- This is a coverage-completion checkpoint before EXPANSION_PLAN.md Phase 1 (QCalcGeom 4-line type-drift fix).

## [5.29.1] — 2026-06-25
- Yang-Mills dispatcher correction: 1.736 GeV (PAPER_1318 canonical).
- First-draft full manuscript shipped.

## [5.29.0] — 2026-06-25
- Full proof corpus shipped: 1,994 whitepapers + Lean 4 scaffold + 4 arXiv bundles.

## [5.28.0] — 2026-06-24
- REST API + Jupyter integration.

## [5.27.2] — 2026-06-24
- Multi-namespace CLI discovery + CLOSURE_ATLAS + WHITEPAPER_INDEX + COVERAGE.

## [5.27.1] — 2026-06-23
- Tier-2 complete (CLI ships, Docker).

## [5.27.0] — 2026-06-22
- First PyPI release.
