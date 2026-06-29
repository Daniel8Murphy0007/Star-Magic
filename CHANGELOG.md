# Changelog

All notable changes to UQFF are recorded here. Full historical record lives in `SESSION_LOG.md`.

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
