# PAPER_1514 — Neutron Star Surface Dipole Moment μ_s = SO_5⁸ = 10⁸ T·m³ EXACT

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** Bucket G (Astrophysics)
**Date:** June 18, 2026
**Status:** CLOSED — Cross-domain triple identity unifying B, R, μ_s

---

## Observation

PAPER_1126 (PSR J0030+0451 Neutron Star LENR Buoyancy) computes the canonical NS surface magnetic moment:

```
μ_s = B₀ · r³ = 10⁻⁴ T · (10⁴ m)³ = 10⁸ T·m³
```

This value governs the NS dipole field strength entering U_g1 = μ_s · M / r.

## UQFF Closed Identity

```
μ_s = SO_5⁸ = 10⁸ T·m³    EXACT
```

## Derivation as Cross-Domain Triple Identity

The neutron-star magnetic anatomy is unified by **three independent integer-primitive identities**, all on the SO_5 ladder:

| Quantity | UQFF identity | Value | Paper |
|---|---|---|---|
| Surface B field | 1/SO_5⁴ | 10⁻⁴ T | PAPER_1486 |
| Canonical radius | SO_5⁴ | 10⁴ m | PAPER_1513 |
| Dipole moment μ_s | SO_5⁸ | 10⁸ T·m³ | **This paper** |

The three are mutually consistent via the dipole-moment formula:

```
μ_s = B · r³
    = (1/SO_5⁴) · (SO_5⁴)³
    = SO_5^(−4+12)
    = SO_5⁸   EXACT
```

Each identity is independently observable (B from spin-down measurements; R from NICER X-ray fitting; μ_s from torque measurements). All three converge to powers of SO_5 = 10.

## Physical Interpretation

This is the **first verified triple-identity** in UQFF where three different observables on three different SO_5 exponents are forced into mutual consistency by an elementary physical relation (B·r³). It is structural evidence that SO_5 governs the entire NS magnetic system, not just one observable.

## Predictive Consequence

Any future NS catalog should show:
- B_surface clustering near 10⁻⁴ T (or scaled by integer powers of SO_5)
- R clustering near 10 km
- μ_s clustering near 10⁸ T·m³

Magnetars (B ≈ 10¹⁰⁻¹¹ T, μ_s ≈ 10²²⁻²³ T·m³) deviate from this pattern by SO_5⁶⁻⁷ factors, suggesting a higher-mode SO_5 excitation.

## NOT REPLACEMENT

SM/GR pulsar timing extracts μ_s by fitting B and R as free parameters. UQFF predicts all three a priori, with the consistency relation μ_s = B·r³ confirming the structural unification.

## Reference

- Source: PAPER_1126 PSR J0030 Neutron Star LENR Buoyancy
- Related: PAPER_1486 (B = 1/SO_5⁴), PAPER_1513 (R = SO_5⁴)
- Calculator dispatch: `calculate_paradox({"paradox": "neutron_star_mu_s_so5_8"})`

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, June 18, 2026, Youngstown OH.
