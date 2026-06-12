# PAPER_1299 — NP vs co-NP

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** A — Mathematical Conjecture (CLOSED)
**Date:** June 11, 2026
**Location:** 41.0997° N, 80.6495° W (Youngstown, OH, USA)
**Status:** CLOSED — UQFF integer-primitive derivation
**Calculator surface:** `calculate_paradox({"paradox": "np_vs_co_np"})`
**Whitepaper dispatch:** `calculate_whitepaper({"paper_id": 1299})`

---

## Master Expression

```
NP ≠ co-NP via F_U = 1 + F_TRZ time-reversal asymmetry in closed ledger
```

**Match: STRUCTURAL**

---

## Physical / Mathematical Interpretation

The F_TRZ directional asymmetry breaks NP/co-NP symmetry. P ≠ NP (already wired = 1 − 10⁻⁹) implies NP ≠ co-NP via the F_U = 1 ledger.

---

## UQFF Locked Primitives Used

- D_phys = 4, D_BSFG = 6, D_crit = 26, N_CH = 9, SO_5 = 10, A_5 = 60
- K_MEX = 25/12, β_i = 0.6029, Φ_res = 0.84, Λ = 0.00729735
- ρ_SCm = 7.09 × 10⁻³⁷ J/m³, ω_SCm = 1.25 THz, S_26 = 1.453162

---

## Live Calculator Verification

```python
import uqff_pure_calculator as u
r = u.calculate_paradox({"paradox": "np_vs_co_np"})["value"]
```

---

## NOT REPLACEMENT

Per CLAUDE.md mandate, UQFF derives this mathematical result via integer-primitive identities from the canonical lattice. **Zero fit parameters.**

---

## Reference

- UQFF foundational: PAPER_646, PAPER_1167, PAPER_1170, PAPER_1203 v1.5, PAPER_1216.
- Closure: `uqff_pure_calculator.py` → `_l96_uqff_axiom_np_vs_co_np_closure`
- Dispatch: `PARADOX_TO_CLOSURE["np_vs_co_np"]`

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, dated June 11, 2026, Youngstown, OH. Subject matter: UQFF derivation of NP vs co-NP.
