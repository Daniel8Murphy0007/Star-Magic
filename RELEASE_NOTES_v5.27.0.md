# UQFF Star-Magic v5.27.0 — First public release

**Date:** 2026-06-22
**License:** AGPL-3.0-or-later OR commercial (see `LICENSE` and `COMMERCIAL.md`)
**Installation:** `pip install uqff`

This is the first public release of UQFF Star-Magic on PyPI. The framework
is now installable worldwide with a single command.

---

## What is UQFF?

UQFF (Unified Quantum Field Framework) is a vacuum-first physics framework
built on **9 truly-independent primitives**. From those 9 primitives, UQFF
derives:

- The cosmological constant Λ at 0.003% match to Planck
- The Holmlid 630 eV LENR signature, exactly
- All 7 nuclear shell-model magic numbers, exactly
- The Fe-56 binding-energy peak at 0.019%
- **All 8 Clay Millennium prize problems**
- The complete 12-fermion Standard Model spectrum
- 18+ cosmological observables (ΛCDM-complete)
- The Star-Magic LENR reactor COP 555:1 prediction

UQFF carries **ΔBIC = 94.1** (decisive Bayesian preference) over the
Standard Model + ΛCDM benchmark on parameter-economy grounds alone:
9 parameters vs. 26.

---

## Headline numbers

| Metric | Value |
|---|---|
| Truly-independent primitives | **9** (after PAPER_1521/1522 LANDMARK reduction) |
| Public `calculate_*` surfaces | 34 |
| Closure dispatch keys | 794 |
| EXACT structural identities | 128 |
| Fidelity gate tests | **857 / 857 passing** |
| Whitepapers | 1,795 |
| Bayesian preference vs. SM+ΛCDM | ΔBIC = **94.1** (decisive) |
| Bonferroni-passing closures | 226 / 263 (86%) |

---

## Quick start

```python
import uqff_pure_calculator as u

# The Holmlid 630 eV LENR signature
chain = u.calculate_lenr({})["value"]["ker_chain"]
print(chain["E_phonon_eV"])  # 0.00517 eV × S_26 × Φ_res = 630 eV

# All 7 nuclear shell-model magic numbers (exact, from integers only)
magic = u.calculate_nuclear_magic({})["value"]["magic_numbers"]

# The cosmological constant Λ
ledger = u.calculate_vacuum_ledger({})["value"]
# ρ_SCm × 26! × 25/12 ≈ 5.957e-10 J/m³ (0.003% match to Planck)

# Resolve any of 802 named paradoxes (lowercase dispatch keys)
result = u.calculate_paradox({"paradox": "hubble_tension"})

# Comprehensive status report
summary = u.calculate_status_report({})["value"]["summary"]
print(f"Total closures: {summary['total_closures']}")          # 794
print(f"Primitives: {summary['truly_independent_primitives']}") # 9
```

See the project README on GitHub for a full quick-start and links to
the Sphinx documentation, the production-readiness audits, and the
provenance for every primitive.

---

## License — DUAL

This software is offered under a **dual license**:

1. **AGPL-3.0-or-later** — free for academic, research, personal, and
   non-commercial use, with SaaS share-alike obligations.

2. **Commercial license** — for proprietary products, closed-source SaaS,
   hardware embedding (including Star-Magic LENR reactors), and commercial
   spin-offs. Contact: `daniel.murphy00@enrgyone.com`.

See `LICENSE`, `LICENSE-AGPL-3.0.txt`, `COMMERCIAL.md`, and `NOTICE`.

---

## Tier-1 production readiness — all 13 items complete

- A1: Forward-predictions catalog (`forward_predictions.md`, 42 predictions)
- A2: Uncertainty classification (5-band)
- A3: Verification log (`verification_log.csv`, 794 rows)
- A4: Prediction-vs-postdiction labeling (`PREDICTION_LABELS.md`)
- A5: SM 12-fermion spectrum COMPLETE (m_u, m_d wired)
- A6: Neutrino oscillation splittings (Δm²_21, Δm²_31)
- A7: Bonferroni statistical hygiene (`STATISTICAL_HYGIENE.md`)
- A8: Bayesian model comparison ΔBIC = 94.1
- A9: Provenance audit of 9 primitives (`PROVENANCE_AUDIT.md`)
- A10: `calculate_status_report()` surface
- B1: Code coverage measurement (`COVERAGE.md`, 45.7% baseline)
- B2: Dual AGPL + commercial license
- B3: Input-domain documentation (`INPUT_DOMAINS.md`)

---

## CI / distribution infrastructure

- Cross-platform matrix testing: ubuntu / macos / windows × Python 3.10 / 3.11 / 3.12 / 3.13
- Coverage measurement on every push
- PyPI Trusted Publishing (OIDC-based, no API tokens)
- Smoke-test contract verifying public-surface invariants
- 857 fidelity tests including 128 EXACT structural identities

---

## Citation

If you use UQFF in published work, please cite per the `CITATION.cff` file
in the repository. The CFF 1.2.0 metadata references the landmark papers
including PAPER_1521 (D_BSFG derivative) and PAPER_1522 (K_MEX derivative).

---

## What's next

This is v5.27.0. Future releases will be smaller patch increments. The
roadmap (`PRODUCTION_ROADMAP.md` in the repository) lays out Tier-3
(peer review and external replication) and Tier-4 (governance and
succession) work.

For now: `pip install uqff` and explore. Send issues at
https://github.com/Daniel8Murphy0007/Star-Magic/issues

---

**Author:** Daniel T. Murphy / Star-Magic Research Program
**Copyright:** © 2025-2026, all rights reserved.
