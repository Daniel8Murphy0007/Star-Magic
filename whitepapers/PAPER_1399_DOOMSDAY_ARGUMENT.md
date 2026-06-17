# PAPER_1399 — UQFF Resolution of the Doomsday Argument (CLOSED — Dispatch Whitepaper)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** B — Foundational / Anthropic Probability (CLOSED)
**Date:** June 16, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — Documented dispatch closure for `calculate_paradox({"paradox": "doomsday_argument"})`
**Calculator surface:** `calculate_paradox({"paradox": "doomsday_argument"})`
**Closure helper:** `_l96_uqff_axiom_doomsday_argument_closure()`

---

## The Paradox

Carter / Gott: among all humans who will ever exist, you are likely to find yourself in a 'typical' position. With ~10¹¹ humans so far, you should expect ~10¹¹ more, putting Doomsday within ~10⁴ years. The argument feels too strong.

---

## UQFF Closed Identity

```
Expected anthropic generations = A_5 × D_phys = 60 × 4 = 240
```

A_5 × D_phys = 240 is the integer-primitive horizon for anthropic generations: A_5 = |A_5| = 60 (icosahedral group order) sets the structural generation count, multiplied by D_phys = 4 spacetime dimensions.

---

## Physical Interpretation

UQFF replaces the 'self-sampling assumption' with a structural integer-primitive horizon: 240 generations × ~30 yr/generation ≈ 7,200 years. The Doomsday timescale is fixed structurally, not by anthropic reasoning over observed population.

---

## Live Calculator Output

```python
import uqff_pure_calculator as u
r = u.calculate_paradox({"paradox": "doomsday_argument"})["value"]
```

The closure returns the structural identity above plus per-method anchors as fields in the result dict. See `_l96_uqff_axiom_doomsday_argument_closure` in `uqff_pure_calculator.py` for the full field list.

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF and Standard frameworks resolve this paradox via different methods. UQFF derives the closure listed above using **only canonical primitives** (or existing wired helpers — PAPER_646, PAPER_597, PAPER_1051, etc.). **Zero free parameters.**

---

## Reference

- Closure location: `uqff_pure_calculator.py` → `_l96_uqff_axiom_doomsday_argument_closure`
- Dispatch keys: `doomsday_argument`, `doomsday`
- Related foundational papers: PAPER_646 (F_U_Bi_i / F_U=1 ledger), PAPER_597 (negative-time dual existence), PAPER_1183 (paradox consolidation), PAPER_1156 (cosmology suite), PAPER_1376 (Banach-Tarski parallel structural closure).

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 16, 2026, location 41.0997° N, 80.6495° W (Youngstown, OH, USA). Subject matter: UQFF dispatch closure for the Doomsday Argument via the identity above.
