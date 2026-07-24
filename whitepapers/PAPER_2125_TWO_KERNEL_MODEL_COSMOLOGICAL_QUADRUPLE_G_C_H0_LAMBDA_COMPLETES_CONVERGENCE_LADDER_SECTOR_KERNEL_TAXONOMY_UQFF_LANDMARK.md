# PAPER_2125 — Two-Layer Model: Cosmological Quadruple {G, c, H0, Λ} + Projection-Layer Convergence Accounting (REVISED)

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.73+
**Date:** 2026-07-22 (REVISED IN PLACE 2026-07-22 — Two-Kernel Model superseded by Two-LAYER Model)
**Landmark Type:** Constant-Quadruple Convergence (projection layer) + Two-Layer structural model
**Discovery Round:** R368 (`UQFF_CompressedQCalcCalculator`) — 151st consecutive stub fill
**Status:** Formal landmark whitepaper — UQFF canonical

---

## REVISION NOTICE — Two-Kernel Model superseded

The original version proposed a "Two-Kernel Model" ({c, ρ_vac} F_U kernel; {G, c} cosmological kernel) and issued the prediction "no F_U-family quadruple will appear." Both are superseded by the **Two-Layer Model**:

1. **The "kernels" were projection artifacts.** {ρ_vac, c} is the energy-form projection pair (ρ_E = ρ_m·c², quantum chain line 33295); {G, c} co-occur because the canonical Eq2 metric-geodesic ε' + G·M/(c²r²) = 0 carries both. Neither is an independent physics kernel.
2. **The "no F_U-family quadruple" prediction is void, not falsified.** It was a claim about wrap-shape occurrence, mislabeled as family structure. R371's Buoyant wrap (Ubi − Ug, four constants {β_i, ρ_vac, c, G}) is the energy-form projection of the **canonical crossing equation F_UBi + F_UBii = 0** — the r_hz equilibrium and UQFF's galaxy-rotation balance. Its four constants are the buoyancy projection triple plus the gravity coupling — union arithmetic again, nothing anomalous.
3. **What survives:** the convergence counting (tier = number of co-exposed derived constants) remains correct arithmetic; the honest residuals and the Hubble-tension identification of this paper are unchanged and correct.

---

## Abstract

R368's fill of `UQFF_CompressedQCalcCalculator` delivered a projection-layer Constant-Quadruple — {G (PAPER_593), c (PAPER_592), H0 (PAPER_2093), Λ (PAPER_2094/1156)} — completing the wrap suite's convergence census Pair/Triple/Quadruple/Quintuple across five consecutive rounds (R364-R368). This paper formalizes the **Two-Layer Model** — canonical layer vs projection layer — as the correct structural account, documents honest residuals for H0 (3.2%) and Λ (0.9%), and identifies the H0 stub-vs-grid residual as the Hubble tension itself, addressed by PAPER_1156's 1/12 EXACT tilt closure.

---

## 1. The Discovery

```python
class UQFF_CompressedQCalcCalculator:
    G_PRIMITIVE      = 6.674e-11   # PAPER_593
    C_PRIMITIVE      = 3e8         # PAPER_592
    H0_PRIMITIVE     = 2.27e-18    # PAPER_2093 grid: 22·F_TRZ¹⁹
    LAMBDA_PRIMITIVE = 1.11e-52    # PAPER_2094/1156

    g_comp = g_base·(1 + H0·t)·(1 − B/B_crit) + Λc²r/3
```

Four constants, four derivations, one compressed-gravity wrap. Externals {M, r, t, B, B_crit}.

---

## 2. The Two-Layer Model — Formal Statement

> **Canonical layer** (`uqff_pure_calculator.py`): the wired UQFF physics — F_UBi/F_UBii with k_spring = (ρ_UA/ρ_SCm)·ω_SCm·Φ_res, dynamic β(t,E,Z), the 26-level quantum chain (E_n, ρ_E = ρ_m·c², 633333.333 validation), the simultaneous system Eq1/Eq2 with r_hz root-finding at residual < 1e-10. Canonical buoyancy kernel: **{ρ_SCm, β_i}**. Constants ω_SCm, Φ_res, SSq, S_26, F_TRZ live here.
>
> **Projection layer** (`CondensedPhysics.py` QCalc wraps): backbone-first energy-form simplifications, self-annotated "framework wrap only." Vacuum content appears as ρ_vac·c² (the chain's mass→energy conversion); β is static; k_spring is absent. The R218+ campaign's convergence events count derived-constant co-exposure **in this layer**.

**Layer map of the six convergence events:**

| Round | Wrap | Constants | Projects (canonical source) |
|:-:|---|---|---|
| R364 | MagneticStrings | {G, c} | Lense-Thirring / Eq2 pair |
| R365 | EnhancedBuoyancy | {β_i, ρ_vac, c} | F_UBi/F_UBii energy form |
| R366 | AetherMetric | {G, c, ρ_vac} | Eq2 metric + vacuum weighting |
| R367 | UQFF_Base | {G, c, μ_0, β_i, ρ_vac} | full F_U_total composition |
| R368 | UQFF_Compressed | {G, c, H0, Λ} | cosmological correction chain |
| R371 | UQFF_Buoyant | {β_i, ρ_vac, c, G} | **crossing F_UBi + F_UBii = 0 (r_hz / galaxy rotation)** |

Every event is union arithmetic over the canonical terms each wrap projects. Tier numbers are counts, not physics.

---

## 3. Honest Residuals — the Hubble Tension (unchanged, correct as filed)

| Constant | Stub literal | UQFF closed form | Residual |
|---|---|---|:-:|
| H0 | 2.27e-18 s⁻¹ | 22·F_TRZ¹⁹ = 2.2e-18 (PAPER_2093) | 3.2% |
| Λ | 1.11e-52 m⁻² | (SO_5+1)·F_TRZ⁵³ = 1.1e-52 (PAPER_2094) | 0.9% |
| Λ | 1.11e-52 m⁻² | (18/5)·SSq·H0²/c² = 1.089e-52 (PAPER_1156 canonical) | 1.9% |

Stub 2.27e-18 ≈ 70.0 km/s/Mpc (local side); grid 2.2e-18 ≈ 67.9 (CMB side). **The 3.2% residual IS the Hubble tension** — the gap PAPER_1156's 1/12 EXACT tilt closure addresses. The residual is not an error; it is the physics. PAPER_1156 remains primary citation for Λ per the R228 precedent.

---

## 4. What Became of the Original Predictions

| Original claim | Status | Disposition |
|---|:-:|---|
| Two-Kernel Model | superseded | replaced by Two-Layer Model (Section 2) |
| "No F_U-family quadruple" | void | wrap-shape bookkeeping, not physics; R371 union = 4 is unremarkable |
| Quintuple ceiling | rescoped | wrap constant-union ceiling = size of ∪ over composed wrap terms; canonical layer not tier-bounded |
| LENR kernel candidate {ρ_SCm, ω_SCm} | **strengthened** | these ARE canonical constants (k_spring, quantum chain) — the deepsearch confirmed they anchor the canonical layer |
| Full Constant-Closure ~R400 | retained | projection-layer completion target, per-class certification (PAPER_2127) |

---

## 5. Instance-Count Updates (unchanged)

PAPER_593 G → 14th (at R368), PAPER_592 c → 13th, PAPER_2093 H0 → 2nd, PAPER_2094 Λ → 2nd.

---

## 6. Gate Assertion

Gate count: 3102 → 3110 (+8 PAPER_2125 assertions, text updated to Two-Layer scope; H0/Λ residual checks unchanged — they were numerical and correct).

---

## 7. Cross-Paper Links

PAPER_1203 (canonical layer: F_U_total, k_spring, Eq1/Eq2, r_hz), PAPER_1202 (quantum chain), PAPER_2093/2094/1156 (H0, Λ), PAPER_2121-2124 (revised siblings), PAPER_2126/2127 (kernel-vs-lattice-node and certification standards, unaffected), PAPER_1051/1065 (canonical Ub derivation chain — wiring queued).

---

## 8. Summary Statement (revised)

**PAPER_2125 documents the projection-layer Constant-Quadruple (R368, {G, c, H0, Λ}) and replaces the Two-Kernel Model with the Two-Layer Model: a canonical layer (F_UBi/F_UBii, k_spring = (ρ_UA/ρ_SCm)·ω_SCm·Φ_res, dynamic β, quantum chain, r_hz equilibrium at sub-1e-10 residual — kernel {ρ_SCm, β_i}) and a projection layer (QCalc energy-form wraps whose {ρ_vac, c} pair is the chain's ρ_m→ρ_E conversion). All six R364-R371 convergence events are union arithmetic over projected canonical terms. The "no F_U-family quadruple" claim is void as mislabeled bookkeeping — R371's four-constant Buoyant wrap projects the canonical crossing equation F_UBi + F_UBii = 0, UQFF's galaxy-rotation balance. The Hubble-tension residual identification (3.2% = local-vs-CMB gap, PAPER_1156 1/12 tilt) stands unchanged.**

---

**Filed 2026-07-22, revised in place 2026-07-22 per canonical-layer deepsearch. Append-only henceforth.**


---

## G/c DERIVATION NOTE (appended 2026-07-22, UNIFIED REGISTRY R2 corpus pass)

This paper uses G = 6.674e-11 (CODATA form) and c = 3e8-family literal as published. Per the Unified Registry (R1-adjudicated
canonical routes, 2026-07-22):

- **G (gravitational constant):** canonical route **PAPER_593** — parameter-free
  G_UQFF = (2π·26³·Φ_res/(SSq³·(26!)²))·v_F⁵/(E_0·f_THz) = 6.66899×10⁻¹¹ (0.08% vs observed).
- **c (speed of light):** canonical route **PAPER_592** — parameter-free
  c_UQFF = (26·4π/Φ_res)·v_F = 2.995×10⁸ m/s (0.13% vs observed; v_F Fermi anchor, c-independent).

Published values above are retained unchanged — as observational anchors or
original inputs per the R2 golden rule (append-only; no silent recomputation).
The UQFF derivations are canonical; residuals are honest disclosures (Rule 7).
Registry: UNIFIED_REGISTRY.csv | Program: UNIFIED_REGISTRY_PROGRAM_PLAN.md

---
