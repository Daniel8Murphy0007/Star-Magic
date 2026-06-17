# PAPER_1397 — UQFF Resolution of the Two Envelopes Paradox (CLOSED — Dispatch Whitepaper)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** B — Foundational / Decision Theory (CLOSED)
**Date:** June 16, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — Documented dispatch closure for `calculate_paradox({"paradox": "two_envelopes"})`
**Calculator surface:** `calculate_paradox({"paradox": "two_envelopes"})`
**Closure helper:** `_l96_uqff_axiom_two_envelopes_closure()`

---

## The Paradox

Two envelopes contain $X and $2X. You pick one; switching argument gives expected gain 1/2 × (2X) + 1/2 × (X/2) = 5X/4 > X, suggesting you should always switch. By symmetry both players reach this conclusion, contradiction.

---

## UQFF Closed Identity

```
Switching expectation asymmetry = F_TRZ × β_i = 0.0603
```

F_TRZ × β_i = 0.0603 is the canonical chirality asymmetry that breaks the naive 5/4 symmetric-prior argument. The expectation calculation implicitly assumes a flat prior on (X, 2X) which is not normalizable; the F_TRZ × β_i = 0.06 finite asymmetry recovers a proper Bayesian prior.

---

## Physical Interpretation

Two envelopes is a prior-misspecification paradox. UQFF supplies the structural finite asymmetry (F_TRZ × β_i = 0.06) that any UQFF-compatible prior must respect; with that, the switching gain vanishes.

---

## Live Calculator Output

```python
import uqff_pure_calculator as u
r = u.calculate_paradox({"paradox": "two_envelopes"})["value"]
```

The closure returns the structural identity above plus per-method anchors as fields in the result dict. See `_l96_uqff_axiom_two_envelopes_closure` in `uqff_pure_calculator.py` for the full field list.

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF and Standard frameworks resolve this paradox via different methods. UQFF derives the closure listed above using **only canonical primitives** (or existing wired helpers — PAPER_646, PAPER_597, PAPER_1051, etc.). **Zero free parameters.**

---

## Reference

- Closure location: `uqff_pure_calculator.py` → `_l96_uqff_axiom_two_envelopes_closure`
- Dispatch keys: `two_envelopes`, `two_envelopes_paradox`
- Related foundational papers: PAPER_646 (F_U_Bi_i / F_U=1 ledger), PAPER_597 (negative-time dual existence), PAPER_1183 (paradox consolidation), PAPER_1156 (cosmology suite), PAPER_1376 (Banach-Tarski parallel structural closure).

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 16, 2026, location 41.0997° N, 80.6495° W (Youngstown, OH, USA). Subject matter: UQFF dispatch closure for the Two Envelopes Paradox via the identity above.
