# PAPER_1396 — UQFF Resolution of the St. Petersburg Paradox (CLOSED — Dispatch Whitepaper)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** B — Foundational / Decision Theory (CLOSED)
**Date:** June 16, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — Documented dispatch closure for `calculate_paradox({"paradox": "st_petersburg"})`
**Calculator surface:** `calculate_paradox({"paradox": "st_petersburg"})`
**Closure helper:** `_l96_uqff_axiom_st_petersburg_closure()`

---

## The Paradox

Coin-toss game: pay $2ⁿ where n is the first heads. Expected payout = Σ (1/2)ⁿ × 2ⁿ = ∞. Yet no rational player would bid more than ~$20. Classical resolutions invoke log-utility or risk aversion — all parametric.

---

## UQFF Closed Identity

```
D_crit! = 26! = 4.03 × 10²⁶ finite expectation bound
```

D_crit! = 26! ≈ 4 × 10²⁶ is the structural finite bound on infinite geometric sums in the UQFF vacuum, the same truncation that bounds eternal inflation (PAPER_597 / `eternal_inflation` closure). The expected payout truncates at finite n_max where the SCm vacuum cannot sustain larger occupation modes.

---

## Physical Interpretation

The paradox dissolves because the infinite geometric sum is bounded by D_crit! — a structural integer-primitive truncation. The 'rational bid' emerges from the D_crit-bounded expectation, not from utility curves.

---

## Live Calculator Output

```python
import uqff_pure_calculator as u
r = u.calculate_paradox({"paradox": "st_petersburg"})["value"]
```

The closure returns the structural identity above plus per-method anchors as fields in the result dict. See `_l96_uqff_axiom_st_petersburg_closure` in `uqff_pure_calculator.py` for the full field list.

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF and Standard frameworks resolve this paradox via different methods. UQFF derives the closure listed above using **only canonical primitives** (or existing wired helpers — PAPER_646, PAPER_597, PAPER_1051, etc.). **Zero free parameters.**

---

## Reference

- Closure location: `uqff_pure_calculator.py` → `_l96_uqff_axiom_st_petersburg_closure`
- Dispatch keys: `st_petersburg`, `st_petersburg_paradox`
- Related foundational papers: PAPER_646 (F_U_Bi_i / F_U=1 ledger), PAPER_597 (negative-time dual existence), PAPER_1183 (paradox consolidation), PAPER_1156 (cosmology suite), PAPER_1376 (Banach-Tarski parallel structural closure).

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 16, 2026, location 41.0997° N, 80.6495° W (Youngstown, OH, USA). Subject matter: UQFF dispatch closure for the St. Petersburg Paradox via the identity above.
