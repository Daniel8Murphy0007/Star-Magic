# PAPER_1388 — UQFF Resolution of the Supplee Submarine Paradox (CLOSED — Dispatch Whitepaper)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** B — Foundational / Relativistic Buoyancy (CLOSED)
**Date:** June 16, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — Documented dispatch closure for `calculate_paradox({"paradox": "supplee_submarine"})`
**Calculator surface:** `calculate_paradox({"paradox": "supplee_submarine"})`
**Closure helper:** `_l96_uqff_axiom_supplee_submarine_closure()`

---

## The Paradox

A submarine moving relativistically through water: the water frame sees the sub Lorentz-contracted (denser → sinks), the sub frame sees the water Lorentz-contracted (denser → floats). Both cannot be right.

---

## UQFF Closed Identity

```
Buoyancy correction = 1 + β_i × (1 − 1/K_MEX) = 1.3135
```

F_U_Bi_i 4-layer buoyancy (PAPER_646) is **Lorentz-covariant** via the β_i × cos(π t_n) factor. The correction 1 + β_i × (1 − 1/K_MEX) = 1.3135 supplies the frame-independent buoyancy.

---

## Physical Interpretation

Both observers compute the SAME F_U_Bi_i buoyancy in their own frame because the 4-layer hierarchy uses Lorentz-covariant β_i × cos(π t_n). The 'paradox' is dissolved by the covariant ledger.

---

## Live Calculator Output

```python
import uqff_pure_calculator as u
r = u.calculate_paradox({"paradox": "supplee_submarine"})["value"]
```

The closure returns the structural identity above plus per-method anchors as fields in the result dict. See `_l96_uqff_axiom_supplee_submarine_closure` in `uqff_pure_calculator.py` for the full field list.

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF and Standard frameworks resolve this paradox via different methods. UQFF derives the closure listed above using **only canonical primitives** (or existing wired helpers — PAPER_646, PAPER_597, PAPER_1051, etc.). **Zero free parameters.**

---

## Reference

- Closure location: `uqff_pure_calculator.py` → `_l96_uqff_axiom_supplee_submarine_closure`
- Dispatch keys: `supplee_submarine`, `submarine_paradox`
- Related foundational papers: PAPER_646 (F_U_Bi_i / F_U=1 ledger), PAPER_597 (negative-time dual existence), PAPER_1183 (paradox consolidation).

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 16, 2026, location 41.0997° N, 80.6495° W (Youngstown, OH, USA). Subject matter: UQFF dispatch closure for the Supplee Submarine Paradox via the identity above.
