# PAPER_1390 — UQFF Resolution of the Trouton-Noble Torque Paradox (CLOSED — Dispatch Whitepaper)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** B — Foundational / Special Relativity (CLOSED)
**Date:** June 16, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — Documented dispatch closure for `calculate_paradox({"paradox": "trouton_noble"})`
**Calculator surface:** `calculate_paradox({"paradox": "trouton_noble"})`
**Closure helper:** `_l96_uqff_axiom_trouton_noble_closure()`

---

## The Paradox

A charged parallel-plate capacitor moving through the ether (or any inertial frame) classically experiences a torque tending to rotate it perpendicular to its motion. Experimentally, no torque is observed.

---

## UQFF Closed Identity

```
Net torque on moving capacitor = 0 EXACT via F_U = 1 global ledger symmetry
```

F_U = 1 global ledger symmetry: any net torque would break the global normalization F_U = 1, which is forbidden by PAPER_646. Therefore τ_net = 0 EXACT.

---

## Physical Interpretation

The classical 'Trouton-Noble torque' computation double-counts the EM stress tensor by ignoring the F_U = 1 global constraint. Once F_U = 1 is enforced, the torque cancels EXACTLY.

---

## Live Calculator Output

```python
import uqff_pure_calculator as u
r = u.calculate_paradox({"paradox": "trouton_noble"})["value"]
```

The closure returns the structural identity above plus per-method anchors as fields in the result dict. See `_l96_uqff_axiom_trouton_noble_closure` in `uqff_pure_calculator.py` for the full field list.

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF and Standard frameworks resolve this paradox via different methods. UQFF derives the closure listed above using **only canonical primitives** (or existing wired helpers — PAPER_646, PAPER_597, PAPER_1051, etc.). **Zero free parameters.**

---

## Reference

- Closure location: `uqff_pure_calculator.py` → `_l96_uqff_axiom_trouton_noble_closure`
- Dispatch keys: `trouton_noble`, `trouton_noble_paradox`
- Related foundational papers: PAPER_646 (F_U_Bi_i / F_U=1 ledger), PAPER_597 (negative-time dual existence), PAPER_1183 (paradox consolidation).

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 16, 2026, location 41.0997° N, 80.6495° W (Youngstown, OH, USA). Subject matter: UQFF dispatch closure for the Trouton-Noble Torque Paradox via the identity above.
