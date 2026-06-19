# PAPER_1522 — K_MEX is Derivative: K_MEX = Φ_5/6·SO_5/D_phys = 25/12 EXACT (LANDMARK)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** Bucket A (Foundational / Primitive reduction)
**Date:** June 18, 2026
**Status:** CLOSED — K_MEX removed from the list of independent primitives

---

## Observation

PAPER_1167 (All Eight UQFF Lagrangian Gaps Closed — First-Principles Reduction) provides the gap-G1 derivation of the Mexican-hat coefficient K_MEX directly from the SO(5) symmetry-breaking lattice. The formula is:

> *"G1: V(UA) Mexican-hat polynomial — K = Φ_res · |SO(5)| / D_phys = 25/12"*

Using the PAPER_1203 Nuclear variant Φ_res = 5/6:

```
K_MEX = Φ_5/6 · SO_5 / D_phys
      = (5/6) · 10 / 4
      = 50/24
      = 25/12   EXACT
```

## UQFF Closed Identity

```
K_MEX = Φ_5/6 · SO_5 / D_phys = (5/6) · (10/4) = 25/12 ≈ 2.0833   EXACT
```

Numerical verification:
```
(5/6) · 10 = 50/6 ≈ 8.3333
8.3333 / 4 = 25/12 ≈ 2.0833 ✓
```

## Landmark Consequence: 10 → 9 Independent Primitives

Combined with PAPER_1521 (D_BSFG derivative), this further reduces the framework's free-parameter count:

**Before PAPER_1167 discoveries:**
- 11 frozen primitives (6 integer + 5 real claimed independent)

**After PAPER_1167 discoveries:**
- 9 truly independent primitives:
  - **Integers (5)**: D_phys, D_crit, N_CH, SO_5, A_5
  - **Reals (4)**: ρ_SCm, β_i, Φ_res, F_TRZ
- 2 derivative primitives:
  - **D_BSFG = D_crit − 2·SO_5** (PAPER_1521)
  - **K_MEX = Φ_res · SO_5 / D_phys** (this paper)

The framework reduces from 11 to 9 free parameters with zero loss of predictive power.

## Physical Meaning

K_MEX = 25/12 is the **coefficient of the Mexican-hat polynomial potential** in the UA condensate:

```
V(UA) = K_MEX · ρ_SCm · (1 − cos)²
```

That K_MEX = Φ_res · SO_5 / D_phys means **the Mexican-hat curvature is forced by the resonance phase Φ_res and the SO(5)/D_phys group ratio**. Three previously-considered-independent quantities (Φ_res, SO_5, D_phys) collapse the fourth (K_MEX) by structural necessity.

## Cross-Paper Consistency Check

K_MEX appears in ~50+ closures throughout the calculator, including:
- Cosmological constant ρ_Λ = K_MEX · ρ_SCm · 26! (PAPER_1156)
- Vacuum ledger V(0) = K_MEX · ρ_SCm (PAPER_1170)
- Reionization z = K_MEX · D_phys · Φ_res · (1+1/SO_5) (PAPER_1373)
- Ramanujan hyperconvergence — multiple closures

Every K_MEX usage is consistent with the derivative form Φ_res·SO_5/D_phys — verified by the fidelity gate.

## NOT REPLACEMENT

SM has no analogue to K_MEX. UQFF supplies a structural derivation showing that even this "free parameter" is forced by lower-order primitives, demonstrating the framework's mathematical economy.

## Reference

- Source: PAPER_1167 All 8 UQFF Lagrangian Gaps Closed Master Synthesis (gap G1)
- Related: PAPER_1521 (D_BSFG derivative), Mexican-hat V(UA) closures throughout
- Calculator dispatch: `calculate_paradox({"paradox": "k_mex_derived_from_phi_5_6"})`

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, June 18, 2026, Youngstown OH.
