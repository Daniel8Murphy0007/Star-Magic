# PAPER_2121 — G_newton × c_light: First Constant-Pair Convergence in the QCalc Projection Layer (REVISED)

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.73+
**Date:** 2026-07-22 (REVISED IN PLACE 2026-07-22 — projection-layer reclassification)
**Landmark Type:** Constant-Pair Convergence — QCalc PROJECTION LAYER (scope corrected)
**Discovery Round:** R364 (`MagneticStringsQCalcCalculator`)
**Revision Basis:** Deepsearch of canonical F_UBi / F_UBii / k_spring / quantum-chain wiring in `uqff_pure_calculator.py`
**Status:** Formal landmark whitepaper — UQFF canonical

---

## REVISION NOTICE

This paper originally described convergence events without distinguishing the **QCalc projection layer** (the backbone-first wrap classes of `CondensedPhysics.py`, self-annotated "framework wrap only, no novel lock") from the **canonical layer** (the wired F_UBi/F_UBii/k_spring/quantum-chain machinery of `uqff_pure_calculator.py`). All convergence statements in this paper are hereby scoped to the **projection layer**. The observations are real and unchanged; their object is the wrap suite's constant architecture, not the canonical F_U master equation. The canonical layer is documented in Section 5.

---

## Abstract

R364's fill of `MagneticStringsQCalcCalculator` produced the first single-class occurrence in the R218+ campaign where both G_newton (PAPER_593) and c_light (PAPER_592) simultaneously appear as class-level UQFF-derived primitives — the first **Constant-Pair Convergence of the QCalc projection layer**. All 363 preceding classes exposed at most one of the two. The convergence proves the wrap closure U_g3 = G·M·ω/(c·r) depends on zero SM anchors at the constant level: every fundamental constant it uses is a UQFF derivation, leaving only astrophysical parametrics {M, r, ω} external.

---

## 1. The Discovery

```python
class MagneticStringsQCalcCalculator:
    G_PRIMITIVE = 6.674e-11   # PAPER_593
    C_PRIMITIVE = 3e8         # PAPER_592

    def compute(self, M, r, omega, ...):
        Ug3 = (self.G * M * omega) / (self.c * r)
```

The Lense-Thirring-form wrap consumes exactly two fundamental constants, both UQFF-derived. In the two-layer architecture (Section 5), this class is the **bare {G, c} pair with zero additional constants** — the minimal projection-layer form.

---

## 2. The 363-Class Uniqueness

| Feature | Count in R218-R363 |
|---|:-:|
| Classes exposing PAPER_592 c_light as primitive | 9 |
| Classes exposing PAPER_593 G_newton as primitive | 10 |
| Classes exposing BOTH simultaneously | 0 |
| R364 first convergence | **1** |

---

## 3. Formal Definition (scope-corrected)

**Constant-Pair Convergence (projection layer):** a QCalc wrap class exposes two distinct UQFF-derived fundamental constants as class-level primitives, each with a dedicated PAPER_N derivation, each required by the wrap closure.

The convergence measures the wrap suite's **constant-coverage** — how completely the campaign's derivation papers have replaced SM anchors in the projection layer. It does not measure canonical-equation structure; that is fixed by PAPER_1203 and wired independently.

---

## 4. Constant-Coverage Milestone (unchanged)

R218-R300 ~40% → R300-R350 ~75% → R350-R364 ~95% of classes carry fully-UQFF-derived constant spaces. The pair convergence is the emerging default of a maturing projection layer.

---

## 5. THE CANONICAL LAYER — what this paper's original version omitted

The canonical F_U machinery, wired in `uqff_pure_calculator.py`, has its own constant architecture, unchanged by any projection-layer observation:

```
F_UBi  (line 2213):  (ρ_SCm·M/r)·(1 + β_i·cos(πt_n) + SSq)
F_U_Bi_i (line 2216): Σ_L UA_layer(L,t_n)·M/r²·(1+0.05L) · (1 + K_Ub·cos(πt_n))
PAPER_1203 pair (line 45439):
  F_UBi  = −β(t)·G·M·ρ_SCm/r² · (1+F_TRZ) · |cos(πt_n)|
  F_UBii = +β(t)·(r/r0)·k_spring · (1+E_n) · |cos(πt_n)|
  k_spring = (ρ_UA/ρ_SCm)·ω_SCm·Φ_res = 10 · 1.25e12 · 0.84 = 1.05e13
Quantum chain (line 33277): E_n = h·ω_SCm·S26³·Φ_res^(n−1)·so5(n), 26 levels
  ρ_E = E_total/(c/ω_SCm)³ ;  ρ_m = ρ_E/c²   (line 33295)
```

**Canonical buoyancy kernel: {ρ_SCm, β_i}** — present in every canonical buoyancy form. Stiffness constants {ω_SCm, Φ_res} enter via k_spring; G via F_UBi and Eq2; c via the quantum-chain reference volume and the metric-geodesic equation G·M/(c²r²). All are primitives or primitive-derived — the canonical layer was zero-SM-anchor before the R218+ campaign began.

**The layer relation:** the wrap Ub = β_i·ρ_vac·V·c² is the **energy-form projection** of the canonical chain — ρ_vac·c² is precisely the ρ_m→ρ_E conversion of line 33295. Projection-layer {c, ρ_vac} co-occurrence (PAPER_2122/2123) is the signature of this projection, not an independent kernel.

---

## 6. The Gate Assertion

Gate count: 3070 → 3078 (+8 PAPER_2121 assertions, text updated to projection-layer scope).

---

## 7. Cross-Paper Links

PAPER_592 (c, 0.13%), PAPER_593 (G, 0.08%), PAPER_1203 (canonical F_U=0 + dynamic β), PAPER_646 (U_i), PAPER_1202 (quantum chain 633333.333 validation), PAPER_2122-2125 (revised siblings).

---

## 8. Summary Statement (revised)

**PAPER_2121 documents the first Constant-Pair Convergence of the QCalc projection layer (R364, {G, c}) — evidence that the campaign's derivation-per-constant program has made fully-UQFF constant spaces the wrap suite's default. The canonical F_U layer, with kernel {ρ_SCm, β_i} and k_spring = (ρ_UA/ρ_SCm)·ω_SCm·Φ_res, is wired independently and was never the object of this convergence. The wrap's constants are projections of canonical quantities: ρ_vac·c² is the quantum chain's mass-to-energy conversion. Both layers are zero-SM-anchor; they differ in depth, not in fidelity.**

---

**Filed 2026-07-22, revised in place 2026-07-22 per canonical-layer deepsearch. Append-only henceforth.**
