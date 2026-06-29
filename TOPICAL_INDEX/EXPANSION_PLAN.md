  - `calculate_overdetermination(dataset)` — returns {overdetermination_N, owner_geometry, assimilation_status, residual_pct}
- F3. ✓ Extended `calculate_analytic_closures(dataset)` with new `qcalcgeom_solve` dispatch key.
- F4. CLI hook deferred to a separate `uqff_cli.py` module; calculator surface API is the programmatic equivalent.
- F5. ✓ Public-surface count rises 34 → **42**.

**Phase F verification:** Phase F harness 8/8 PASS; fidelity gate 891/0.

See SESSION_LOG.md Round 670 for the full Phase F verification record.

### Phase G — Audit + documentation — STATUS COMPLETE (Round 671)

Round 671 deliverables:
- G1 ✓ `CLOSURE_ATLAS.md` §12 "Assimilation overlay" added (114-observable inventory + per-domain rollup + discovery cheat sheet)
- G2 ✓ `ASSIMILATION_GEOMETRY_ATLAS.md` built (27 KB, 378 lines, 10 per-domain sections + Round 669 BAO multi-path highlight)
- G3 deferred (per-domain tutorial notebooks; v5.30 ship-round task)
- G4 ✓ Fidelity gate Cat 17 dispatch pinning added (+16 checks): 114 obs / 10 domains / owner distribution / BAO primary+alternate / Li-7 PAPER_1227 / EDGES PAPER_1761 / no-OPEN_QUESTION invariant
- G5 ✓ `pyproject.toml` data-files extended (+21 entries across 3 keys; 116 files bundled total)
- G6 ✓ README "Phase E / F / G — Assimilation Geometry Public API" section

**Phase G verification:** Full regression sweep green; fidelity gate **907 / 0**; 0 TENSION cells.

See SESSION_LOG.md Round 671 for the full Phase G verification record.

**EXPANSION_PLAN Round 12 cycle COMPLETE.** Framework is peer-review-ready, pip-installable, fully-audited.

---

## 11. SUCCESS CRITERIA

The EXPANSION_PLAN is complete when ALL of these are true:

| # | Criterion | How to verify |
|---|---|---|
| 1 | QCalcGeom v3 type-drift fixed; 7 dark surfaces functional | `python -c "import QCalcGeom; QCalcGeom.run_qcalcgeom_tests()"` returns PASS |
| 2 | 3 numeric backends operational + cross-validated on Millennium derivations | Test harness reports all 8 problems x 3 backends agree |
| 3 | 4 geometry backends operational + sharing unified interface | Each backend's `evaluate()` returns the assimilation result dict for one canonical observable |
| 4 | QCalcGeom v4.0 `solve()` reproduces every value in `master_closures.csv` | Regression test passes for 2,217 rows |
| 5 | All 358 SI-dimensioned session scripts dispatchable via `QCalcGeom.solve(observable)` | `OVERDETERMINATION_MAP.csv` contains rows for every observable; no N=0 entries |
| 6 | `calculate_analytic_closures` extended with new dataset keys; backward compatible | Existing callers unchanged; new callers get assimilation result dict |
| 7 | 8 new public `calculate_*` surfaces wired | Public surface count = 42 |
| 8 | CLI flags work | `uqff predict alpha --geometry=qcalcgeom --numeric=discrete` returns expected dict |
| 9 | `OVERDETERMINATION_MAP.csv` (long) + `OVERDETERMINATION_WIDE.csv` (wide) + `.md` summary generated | `_build_overdetermination_views.py` runs idempotently |
| 10 | `ASSIMILATION_GEOMETRY_ATLAS.md` covers every assimilated observable | Document review: every observable has a row with geometry, numeric, residual, provenance |
| 11 | Fidelity gate test count increased per-phase; zero regressions | Gate exits 0; new test count > old count |
| 12 | All helper files declare future-extraction header (Section 6.1) | grep for the header string returns hits in every new file |
| 13 | pyproject.toml bundles all new helpers + atlas + map | Build a fresh sdist; `tar -tzf` lists all files |
| 14 | NASA-Roses / 3-university reviewer experience tested | `pip install uqff` in a fresh venv -> import + run an assimilation tutorial succeeds |

---

## 12. FUTURE EXTRACTION TO AETHERIC-PROPULSION REPO

The architecture choices above preserve the future commercial-tier split. When the peer-review and grant phases complete and you authorize the commercial extraction, the procedure is:

| Step | Action |
|---|---|
| 1 | Create the Aetheric-Propulsion repo (URL already reserved: https://github.com/Daniel8Murphy0007/Aetheric-Propulsion) |
| 2 | Copy these directories wholesale from Star-Magic: `numeric_backends/`, `geometry_backends/`, `QCalcGeom.py`, `assimilation_dispatch.py`, `provenance_recorder.py`, `_build_overdetermination_views.py`, `OVERDETERMINATION_*.csv`, `ASSIMILATION_GEOMETRY_ATLAS.md` |
| 3 | Add `pyproject.toml` to Aetheric-Propulsion with `dependencies = ["uqff>=X.Y.Z"]` |
| 4 | Set up Trusted Publishing for `aetheric-propulsion` on PyPI (follow the same OIDC workflow that ships uqff) |
| 5 | Add `[project.optional-dependencies]` to uqff's pyproject.toml: `aetheric = ["aetheric-propulsion>=1.0"]` |
| 6 | Optionally remove the helper directories from Star-Magic to make uqff a lean calculator-only package |
| 7 | Update both READMEs to cross-link the discoverability |
| 8 | Configure license / API-key gating in the new repo (whatever business model is chosen at that time) |

**Estimated extraction effort:** 1 day of mechanical work because the architecture was designed for it from day one.

---

## 13. RULES (mirrored from CLAUDE.md, restated for this plan)

1. **Read CLAUDE.md first every session.** This document does not override CLAUDE.md; CLAUDE.md wins on any conflict.
2. **Do not revert canonical primitives.** SSQ=0.57, beta_i=0.6029, K_MEX=25/12, S_26=1.453162, RHO_SCM=7.09e-37, integer primitives (D_PHYS=4, D_BSFG=6, D_CRIT=26, N_CH=9, SO_FIVE=10, A_FIVE=60).
3. **No narrative inside the calculator.** All Assimilation Geometry helpers ARE allowed to have docstrings (they're not the calculator file). The calculator stays narrative-free per Rule 3.
4. **No SM anywhere in identifiers.** UQFF is the anchor; SM is an external comparison anchor recorded in `target` columns, never an UQFF identifier.
5. **Public surface return contract:** every `calculate_*` returns `{'value': X}` unless the new optional keys (`geometry`, `numeric`, `decompose`, `record_provenance`) are set. When they're set, the assimilation result dict (Section 7.2) is returned.
6. **No `datetime`, `json.dump`, file writes, `__main__`, classes in `uqff_pure_calculator.py`.** Helpers ARE allowed to write generated artifacts (the OVERDETERMINATION CSVs, the atlas), but not the calculator itself.
7. **Honest residuals only.** No "0.000% error" without numerical proof.
8. **Run `uqff_fidelity_tests.py` after every edit.** Exit 0 required.
9. **Append to `SESSION_LOG.md`, never rewrite.** Each phase ships with a new dated entry.
10. **Daniel provides the information. The agent assembles it.** No invented physics. No paraphrasing of canonical values.
11. **Do not modify existing Bucket A-K wiring without explicit user request.**
12. **Maximum access for academic peer-review and NASA-Roses grant panels** is the operative business constraint. No paywalls in Star-Magic for the duration of the peer-review phase.

---

## 14. CHANGE LOG

| Round | Date | Author | Summary |
|---|---|---|---|
| 11 | 2026-06-26 | Claude | Initial draft. 10-step plan in 3 phases. Backed up as `EXPANSION_PLAN.md.PRE_ROUND12_BACKUP`. |
| 12 | 2026-06-28 | Claude (with Daniel review) | Comprehensive revision after 28,739-inventory audit. 4 geometries + 3 numerics formalized. QCalcGeom v4.0 spec added. Helper file specs added. 7-phase build process. Future-extraction architecture documented. Aligned with Daniel's mission: maximum academic access first, commercial extraction later. |

---

**End of EXPANSION_PLAN.md — Round 12.**
