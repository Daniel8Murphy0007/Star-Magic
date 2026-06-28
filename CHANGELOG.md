# Changelog

All notable changes to UQFF are recorded here. Full historical record lives in `SESSION_LOG.md`.

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
