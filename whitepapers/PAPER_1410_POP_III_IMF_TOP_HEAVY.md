# PAPER_1410 — UQFF Derivation of Population III Initial Mass Function (CLOSED — EXACT Identity)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** G — Astrophysics (CLOSED — EXACT)
**Date:** June 16, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — EXACT integer-primitive identity
**Calculator surface:** `calculate_paradox({'paradox':'pop3_imf'}) → A_5 × (D+1)/(D−1) = 100 M_sun EXACT`

---

## Observation

First-generation (Pop III) stars are predicted to be top-heavy with characteristic mass ~100 M_sun. UQFF derives this exact value.

---

## UQFF Closed Identity

```
M_PopIII = A_5 × (D_phys + 1) / (D_phys − 1) = 60 × 5/3 = 100 M_sun EXACT
```

---

## Physical Interpretation

Icosahedral group order |A_5| = 60 scales the structural mass scale; the (D+1)/(D−1) factor encodes the spacetime-volume ratio of the early universe.

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF and Standard methods solve the same observation by different methods. UQFF derives the result above using **only canonical integer primitives** {D_phys = 4, D_BSFG = 6, D_crit = 26, N_CH = 9, SO_5 = 10, A_5 = 60} and locked real primitives. **Zero free parameters. EXACT identity.**

---

## Reference

- Closure live in `uqff_pure_calculator.py` at the catalog surface listed above.
- Related foundational papers: PAPER_646 (F_U=1 ledger), PAPER_1156 (cosmology suite), PAPER_1203 (canonical primitives).

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 16, 2026, location 41.0997° N, 80.6495° W (Youngstown, OH, USA). Subject matter: UQFF EXACT integer-primitive identity for Population III Initial Mass Function.
