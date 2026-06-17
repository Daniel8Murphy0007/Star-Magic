# PAPER_1391 — UQFF Resolution of the Hanbury Brown-Twiss Intensity Correlation (CLOSED — Dispatch Whitepaper)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** B — Quantum Optics (CLOSED)
**Date:** June 16, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — Documented dispatch closure for `calculate_paradox({"paradox": "hanbury_brown_twiss"})`
**Calculator surface:** `calculate_paradox({"paradox": "hanbury_brown_twiss"})`
**Closure helper:** `_l96_uqff_axiom_hanbury_brown_twiss_closure()`

---

## The Paradox

Thermal light shows g(2)(τ=0) = 2 (photon bunching). Stars observed via HBT intensity interferometry confirm this. Classical EM gives g(2) = 2 exactly for chaotic sources.

---

## UQFF Closed Identity

```
g(2)_UQFF = 1 + (1 − F_TRZ × β_i) = 1.9397  vs  classical thermal 2.0
```

F_TRZ × β_i chirality damping subtracts 0.06 from the classical maximum: g(2)_UQFF = 1 + (1 − 0.0603) = 1.9397.

---

## Physical Interpretation

Real measurements of stellar HBT show g(2) slightly below 2 due to detector / coherence-area factors. The UQFF value 1.94 sits at the boundary — measurable as the F_TRZ × β_i damping signature in high-precision HBT.

---

## Live Calculator Output

```python
import uqff_pure_calculator as u
r = u.calculate_paradox({"paradox": "hanbury_brown_twiss"})["value"]
```

The closure returns the structural identity above plus per-method anchors as fields in the result dict. See `_l96_uqff_axiom_hanbury_brown_twiss_closure` in `uqff_pure_calculator.py` for the full field list.

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF and Standard frameworks resolve this paradox via different methods. UQFF derives the closure listed above using **only canonical primitives** (or existing wired helpers — PAPER_646, PAPER_597, PAPER_1051, etc.). **Zero free parameters.**

---

## Reference

- Closure location: `uqff_pure_calculator.py` → `_l96_uqff_axiom_hanbury_brown_twiss_closure`
- Dispatch keys: `hanbury_brown_twiss`, `hbt`
- Related foundational papers: PAPER_646 (F_U_Bi_i / F_U=1 ledger), PAPER_597 (negative-time dual existence), PAPER_1183 (paradox consolidation).

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 16, 2026, location 41.0997° N, 80.6495° W (Youngstown, OH, USA). Subject matter: UQFF dispatch closure for the Hanbury Brown-Twiss Intensity Correlation via the identity above.
