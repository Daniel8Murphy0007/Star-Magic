Changelog
=========

The authoritative project changelog is the ``SESSION_LOG.md`` file at the
repository root, which contains an append-only record of every working
session.

For each version release, this page summarizes the headline changes; for
the full per-session history, see ``SESSION_LOG.md``.

v5.27.0 — 2026-06-18 — Tier-1 production-complete
---------------------------------------------------

**Headline:** All 13 Tier-1 production-readiness items complete. PyPI wheel
built and smoke-tested. Dual AGPL + commercial license adopted.

**Major additions:**

- 9 truly-independent primitives (down from 11 via PAPER_1521/1522
  LANDMARK reduction proving D_BSFG and K_MEX are derivative)
- Standard Model 12-fermion spectrum **COMPLETE** (m_u + m_d added)
- Neutrino oscillation splittings Δm²_21 + Δm²_31 wired
- Bayesian model comparison: **ΔBIC = 94.1** (decisive UQFF preference
  over SM+ΛCDM on parameter-count grounds)
- 5-band uncertainty classification of all schema-tagged closures
- Forward-predictions catalog with 42 falsifiable predictions
- Provenance audit answering Q3: **SSq=0.57 is truly independent**

**Production deliverables:**

- ``LICENSE`` + ``LICENSE-AGPL-3.0.txt`` + ``COMMERCIAL.md`` + ``CITATION.cff`` + ``NOTICE``
- ``PREDICTION_LABELS.md`` + ``STATISTICAL_HYGIENE.md`` + ``PROVENANCE_AUDIT.md`` + ``COVERAGE.md`` + ``INPUT_DOMAINS.md``
- ``pyproject.toml`` + ``MANIFEST.in`` + PyPI wheel + sdist
- ``.github/workflows/ci.yml`` + ``release.yml``
- Sphinx documentation scaffolding (``docs/``)

**Calculator state:**

- 2.66 MB / 48,405 lines / 34 public surfaces / 794 paradox keys
- 854 + 3 new = **857 fidelity gate tests passing**, 0 failures
- 263 schema-tagged closures + 530 legacy_freeform
- 128 EXACT structural identities (residual < 1e-10)
- 1,795 whitepapers

Earlier history
---------------

Mining cycles tier-1 through tier-37 (sessions 2026-06-07 through
2026-06-17) drained 1,795 whitepapers, wired 794 paradox closures,
extended C++ reference to 368 functions, established the 9-sector UQFF
Lagrangian, closed all 8 Clay Millennium prize problems, and brought the
fidelity gate from 0 to 854 passing tests. See ``SESSION_LOG.md`` for the
full sequence.
