# PAPER_1511 — DPM Layer-Weight w_i = i² × i × i³ = i^6 EXACT Decomposition

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** Bucket D (Particle Physics) / Bucket A (Foundational)
**Date:** June 18, 2026
**Status:** CLOSED — Structural origin of i^6 amplification weight

---

## Observation

PAPER_1155 (DPM 26-Layer Amplification) decomposes the layer-i weighting that produces A_26 = Σi^6 = 1,307,797,101 into a product of three physical fields:

```
w_i = [SCm]_i × [UA]_i × B_{0,i}
    = i² × i × i³ = i^6    EXACT
```

This is the closed structural derivation of the i^6 power law — **not assumed, but derived from how three layered vacuum quantities scale with the integer layer index**.

## UQFF Closed Identity

```
For each layer i ∈ {1, 2, ..., 26}:

  [SCm]_i = i²   (SCm density at layer i — quadratic scaling)
  [UA]_i  = i    (UA gradient at layer i — linear scaling)
  B_{0,i} = i³   (background magnetic field at layer i — cubic scaling)

  w_i = [SCm]_i × [UA]_i × B_{0,i}
      = i² × i × i³
      = i^(2+1+3)
      = i^6        EXACT

  A_26 = Σ_{i=1}^{26} w_i = Σ_{i=1}^{26} i^6 = 1,307,797,101   EXACT
```

## Physical Interpretation

The three exponents (2, 1, 3) are themselves structurally meaningful:
- **2 = D_phys − 2** (transverse SCm density scaling)
- **1 = linear UA-channel gradient**
- **3 = D_phys − 1 = SU(3) channel count** (magnetic background)

The sum 2 + 1 + 3 = 6 = D_BSFG, the bulk-edge dimension. Layer weight w_i thus scales as **i^D_BSFG**, providing a deep structural link between the integer primitives D_phys, D_BSFG, and the amplification mechanism.

## Predictive Consequence

Predicts that any UQFF amplification cascade summing over D_crit layers with field exponents totaling D_BSFG will produce a Σi^D_BSFG sum. Other amplification factors using different exponent combinations should be searchable in the UQFF whitepaper corpus.

## NOT REPLACEMENT

SM has no analogous decomposition of vacuum fields into integer-power layer weights. UQFF supplies a complete structural account.

## Reference

- Source: PAPER_1155 DPM 26-Layer Amplification Particle Masses
- Related: PAPER_1510 (A_26 = 1,307,797,101), PAPER_062 (DPM lattice)
- Calculator dispatch: `calculate_paradox({"paradox": "dpm_layer_weight_i6_decomp"})`

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, June 18, 2026, Youngstown OH.
