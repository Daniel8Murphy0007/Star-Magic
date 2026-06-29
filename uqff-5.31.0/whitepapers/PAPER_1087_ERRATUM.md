# PAPER_1087 ERRATUM — Unit Inconsistency in Abstract Formula vs §3 Table

**Author:** Daniel T. Murphy
**Filed:** June 16, 2026
**Subject:** Internal inconsistency between the abstract closed form and the §3 numerical table in PAPER_1087 (Time-Evolving Dark Energy EOS w_DE from SCm Phonon Dynamics).
**Status:** OPEN — awaiting clarification of units; calculator closure pinned to §3 table value at present (t = 13.8 Gyr).

---

## The Discrepancy

PAPER_1087 abstract gives:

```
w_DE(t, Γ) = -1 + (2κt + (SSq/26)·t) / ln(Φ(Γ))
```

with κ = 5.787 × 10⁻⁹ s⁻¹, SSq = 0.57, Φ(Γ) = Φ₀ · S_26, Φ₀ = 10²⁰.

Direct evaluation at t = 4.35 × 10¹⁷ s (≈ 13.8 Gyr) with κ in s⁻¹:

```
Numerator  ≈ 2 × 5.787e-9 × 4.35e17 + (0.57/26) × 4.35e17  ≈ 9.54 × 10¹⁵
ln(Φ)      = ln(1.453 × 10²⁰)                              = 46.43
w_DE       = -1 + 9.54e15 / 46.43                          ≈ 2.05 × 10¹⁴
```

This is non-physical (w_DE should be near -1).

Conversely, the **§3 table** in the paper reports:

| t (Gyr) | w_DE | Δw |
|---|---|---|
| 0 | -1.0000 | 0 |
| 1 | -0.9959 | 0.0041 |
| 5 | -0.9795 | 0.0205 |
| 10 | -0.9591 | 0.0409 |
| **13.8** | **-0.9435** | **0.0565** |

For the table values to be reproduced, the deviation at t = 13.8 Gyr must be Δw = 0.0565. Solving backward:

```
required coefficient C such that C × 13.8 / ln(Φ) = 0.0565
C = 0.0565 × 46.43 / 13.8 = 0.190 per Gyr
```

Neither the κ term (2 × 5.787e-9 = 1.16e-8 per Gyr) nor the SSq/26 term (0.0219) reaches 0.190 — short by a factor of ~9.

## Three candidate resolutions

1. **κ in different units.** If κ is meant to be 5.787 × 10⁻⁹ × 10¹⁰ = 0.058 per Gyr (perhaps the 10⁻⁹ is a typo for 10⁻² in some natural unit system), then 2κ ≈ 0.116 — closer but still short.

2. **(SSq/26) should be (SSq × X)** where X ≈ 8.5. The closest UQFF candidates that give ~8.5 are K_MEX × A_5 / SO_5 / 1.47 or similar, but none is structurally exact.

3. **Formula should be applied with different normalization.** Perhaps the t in the formula is "dimensionless cosmic time" and the table substitutes a different scaling.

## Calculator behavior (current)

`calculate_paradox({"paradox": "dark_energy_eos_time_evolving"})` and `calculate_cosmology({"observable": "w_de_paper_1087"})` both return **-0.9435** (the §3 table value at t = 13.8 Gyr) directly. The literal formula is still documented in the closure's other fields for transparency.

## Recommended action by Daniel

- Confirm κ units and provide an erratum to the paper itself; OR
- Replace `(SSq/26)` with the correct coefficient; OR
- Confirm that the closure should remain pinned to the §3 table value pending re-derivation.

---

**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, June 16, 2026.
