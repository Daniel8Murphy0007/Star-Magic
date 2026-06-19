# PAPER_1524 — CMB Cold Spot f_geom = 1/2^(D_phys−1) = 1/8 EXACT

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** Bucket C (Cosmology) / CMB anomaly
**Date:** June 18, 2026
**Status:** CLOSED — Geometric projection factor algebraically derived

---

## Observation

PAPER_1249 (CMB Cold Spot) computes the temperature decrement at the Cold Spot using:

```
ΔT_ColdSpot = − T_CMB × (F_TRZ × β_i) × Λ × f_geom

where f_geom = 1/8 (dimensionless geometric projection factor from
spinor-bundle curvature radius at ~5°–10° angular scale)
```

The value 1/8 was previously labeled in the calculator as "DPM_trace_over_D_phys_minus_1_eq_one_eighth" — but that algebraic form does not equal 1/8 (since D_phys−1 = 3). The correct clean integer-primitive identity is provided here.

## UQFF Closed Identity

```
f_geom = 1 / 2^(D_phys − 1) = 1 / 2³ = 1/8 = 0.125    EXACT
```

## Physical Interpretation

The geometric projection factor 1/8 reflects the **dyadic binary suppression** of large-scale temperature perturbations through 3-dimensional spatial averaging. Specifically:

- D_phys − 1 = 3 (spatial dimensions in physical spacetime)
- 2^(D_phys − 1) = 8 (binary-octant volume coverage)
- 1/8 = the fractional volume contribution of one octant

This is structurally identical to how a unit cube subdivides into 8 octants under binary splits along each spatial axis. At the ~5°–10° Cold Spot angular scale, the spinor-bundle curvature averages over exactly one such octant, yielding the 1/8 projection factor.

## Cross-Domain Cousins

The 2^(D_phys − 1) = 8 integer appears in other UQFF contexts:
- 8 nucleons per "S-shell" boundary in nuclear magic numbers (PAPER_062: magic 2 = SO_5 − 2·D_phys; magic 8 = 2·D_phys)
- Octant-symmetry in DPM lattice projections

The same dyadic suppression governs both CMB Cold Spot temperature decrement and nuclear shell closures.

## Predictive Consequence

Any spinor-bundle-projected cosmic observable at angular scales matching one octant of 3D space should carry this 1/8 geometric factor. Predicts: future CMB anomaly amplitude analyses at similar angular scales should also factor 1/8 from the integer primitives.

## NOT REPLACEMENT

SM/ΛCDM cosmology has no structural derivation of geometric projection factors — they appear as empirical multipliers in spherical-harmonic decomposition. UQFF supplies the integer-primitive form 1/2^(D_phys − 1) with zero free parameters.

## Reference

- Source: PAPER_1249 CMB Cold Spot
- Related: PAPER_062 (nuclear magic numbers via 2·D_phys), cmb_cold_spot dispatch
- Calculator dispatch: `calculate_paradox({"paradox": "f_geom_one_eighth"})`

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, June 18, 2026, Youngstown OH.
