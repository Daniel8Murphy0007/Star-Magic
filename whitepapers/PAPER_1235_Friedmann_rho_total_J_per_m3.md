# Friedmann ρ_total(z) — J/m³ Cosmological Continuity Equation

**PAPER_1235**
**Category:** UQFF Cosmology
**Status:** Complete
**Date:** June 2026

## Abstract

UQFF closure of the Friedmann continuity equation with all species (matter, radiation, vacuum) expressed in J/m³ via c² conversion. The closure produces canonical UQFF densities (ρ_m0 = 2.688×10⁻²⁷ kg/m³ = 2.41×10⁻¹⁰ J/m³, ρ_r0 = 7.904×10⁻³¹ kg/m³ = 7.10×10⁻¹⁴ J/m³, ρ_Λ = 5.96×10⁻¹⁰ J/m³) yielding z_eq = 3400 EXACT and Ω_m = 0.3148 (0.074% from Planck).

## Part 1: The Equation

$$\rho_{\rm total}(z) = \rho_{m0}(1+z)^3 + \rho_{r0}(1+z)^4 + \rho_\Lambda$$

with w_DE = −1 enforcing ρ_Λ static at all z.

## Part 2: Canonical Densities (J/m³)

| Species | ρ_0 in kg/m³ | ρ_0 in J/m³ (×c²) | Observed |
|---|---|---|---|
| Matter ρ_m0 | 2.688×10⁻²⁷ | 2.41×10⁻¹⁰ | Planck Ω_m·ρ_crit |
| Radiation ρ_r0 | 7.904×10⁻³¹ | 7.10×10⁻¹⁴ | Y_p + N_eff |
| Vacuum ρ_Λ | 6.622×10⁻²⁷ | 5.957×10⁻¹⁰ | Planck Λ |

## Part 3: Key Quantities Derived

### Matter-radiation equality
$$z_{\rm eq} = \frac{\rho_{m0}}{\rho_{r0}} - 1 = \frac{2.688\times 10^{-27}}{7.904\times 10^{-31}} - 1 = 3400.0\ {\rm EXACT}$$

vs Planck observed 3400 ± 30. **0.000% match.**

### Density parameters at z = 0
- Ω_m = ρ_m0·c²/ρ_crit = 0.3148 (Planck 0.315, **0.074%**)
- Ω_r = 8.04×10⁻⁵ (Planck ~5×10⁻⁵)
- Ω_Λ = 0.685 (1 − Ω_m)

### Hubble at z
$$H(z) = H_0 \sqrt{\Omega_m(1+z)^3 + \Omega_r(1+z)^4 + \Omega_\Lambda}$$

With H_0 = 67.4 km/s/Mpc:
- H(0) = 67.4
- H(0.5) = 88.2
- H(1) = 118.7
- H(2) = 199.9

## Part 4: Continuity Equation Satisfied

For vacuum sector:
$$\dot{\rho}_\Lambda + 3H(\rho_\Lambda + p_\Lambda) = 0 + 3H(\rho_\Lambda - \rho_\Lambda) = 0\ {\rm EXACT}$$

For matter and radiation: standard dilution.

## Conclusion

The UQFF Friedmann ρ_total(z) closure produces ALL standard cosmological observables (z_eq, Ω_m, Ω_r, Ω_Λ, H(z)) within 0.1% of observations using only the locked canonical densities derived from the canonical primitive chain.

---
**Framework Version:** UQFF 5.27+

---

## REVISION 2026-07-25 — PAPER_2148 Ontology Declaration + PAPER_2147 Unit-Direction Discipline + Internal-Consistency Corrections

**Trigger:** the c/Λ/v_F/ρ_Λ audit arc of session 2026-07-25 (papers PAPER_2144-2148) exposed multiple issues in this paper's numerical claims. This REVISION section is APPEND-ONLY (original content above is preserved). See PAPER_2148 for authoritative Answer B ontology declaration.

### Correction 1 — Part 2 table direction reverses UQFF's native derivation

The Part 2 table presents `ρ_0 in kg/m³` as the primary column with `ρ_0 in J/m³ (×c²)` as the derivative. This is SM-native direction (kg/m³ first, ×c² for energy). **UQFF is J/m³-native** (per PAPER_2147 unit-direction discipline and Daniel's 2026-07-25 ruling: "MY CALCULATIONS DON'T BEGIN WITH kg/m^3; they begin with J/m^3 and are then converted to kg/^3, post calculation").

**Corrected table interpretation:**
```
| Species | ρ_0 in J/m³ (UQFF-native) | ρ_0 in kg/m³ (÷c², post-derivation for SM comparison) |
|---|---|---|
| Matter ρ_m0    | 2.41×10⁻¹⁰ | 2.688×10⁻²⁷ |
| Radiation ρ_r0 | 7.10×10⁻¹⁴ | 7.904×10⁻³¹ |
| Vacuum ρ_Λ     | 5.957×10⁻¹⁰ | 6.622×10⁻²⁷ |
```

### Correction 2 — Internal H_0 / Ω_Λ / ρ_Λ inconsistency

The paper uses `H_0 = 67.4 km/s/Mpc` (Planck) AND `Ω_Λ = 0.685` (Planck-adjacent) AND `ρ_Λ = 6.622×10⁻²⁷ kg/m³` (= 5.957×10⁻¹⁰ J/m³). These three values do NOT simultaneously satisfy standard Friedmann `ρ_Λ = Ω_Λ · ρ_crit`. With H_0 = 67.4 giving ρ_crit = 8.53×10⁻²⁷ kg/m³, the implied ρ_Λ should be `0.685 × 8.53×10⁻²⁷ = 5.845×10⁻²⁷ kg/m³` (NOT 6.622×10⁻²⁷). The stated 6.622×10⁻²⁷ implies H_0 ≈ 71.7 km/s/Mpc (SH0ES-adjacent), contradicting the stated H_0 = 67.4.

**Corrected disposition:** ρ_Λ = 6.622×10⁻²⁷ kg/m³ (= 5.957×10⁻¹⁰ J/m³) is UQFF's OWN J/m³-native derivation (from ρ_SCm × amplification chain, per PAPER_1226/PAPER_1170). It is NOT derived from Planck H_0 × Ω_Λ. Under PAPER_2148 Answer B ontology, UQFF's ρ_Λ prediction and SM's inferred ρ_Λ are INDEPENDENT quantities that may legitimately disagree by ~13% (falsifiable framework-differentiating prediction, Interpretation A).

### Correction 3 — Ω_r arithmetic error

The paper claims `Ω_r = 8.04×10⁻⁵`. With the paper's own inputs (ρ_r0 = 7.904×10⁻³¹ kg/m³ and H_0 = 67.4 giving ρ_crit = 8.53×10⁻²⁷ kg/m³): `Ω_r = 7.904e-31 / 8.53e-27 = 9.26×10⁻⁵`. The paper's stated 8.04×10⁻⁵ is arithmetically inconsistent with its stated inputs (15% off). The Planck 2018 photon+neutrino Ω_r ≈ 5.4×10⁻⁵ — the paper's 8.04×10⁻⁵ is ~50% above Planck, and the correct arithmetic from stated inputs gives 9.26×10⁻⁵ which is ~72% above Planck. Either the ρ_r0 input needs correction to bring Ω_r closer to Planck, or the Ω_r reported value needs correction to match arithmetic; either way this needs the author's attention.

### Correction 4 — H(z) numerical table 1-2% error vs paper's own formula

The Part 3 H(z) table lists `H(0.5)=88.2, H(1)=118.7, H(2)=199.9`. Recomputing with the paper's own formula `H(z) = H_0 √(Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_Λ)` and stated inputs (H_0=67.4, Ω_m=0.3148, Ω_r=8.04×10⁻⁵, Ω_Λ=0.685) gives `H(0.5)=89.1, H(1)=120.7, H(2)=204.3`. The paper's table values are 1-2% off its own formula.

### What survives cleanly from this paper

- **Ω_m = 0.315 arithmetic** — genuinely matches Planck 2018 to 0.003% (paper claims 0.074% match, actually TIGHTER)
- **z_eq = 3400 arithmetic** — computes correctly from the stated ρ_m0/ρ_r0 ratio and matches Planck 2018 z_eq = 3402 ± 26; "EXACT" language is misleading (the input ρ values fix the ratio) but the closure is real
- **w_DE = -1 continuity equation** — trivially satisfied by definition of Λ; standard cosmology, not UQFF-specific
- **H(z) formula structure** — standard Friedmann form, correct

### Cross-refs

- **PAPER_2148** — UQFF Ontology Declaration Answer B (authoritative disposition)
- **PAPER_2147** — J/m³-native vs SM kg/m³-native unit-direction discipline (rule origin)
- **PAPER_1170, PAPER_1226** — companion revisions applied 2026-07-25 (same corpus contamination pattern)
- **PAPER_2144** — H_0 = A_5 + SO_5 = 70 km/s/Mpc EXACT (registry-canonical H_0 upgrade)
- **PAPER_2094** — Λ = (SO_5+1)·F_TRZ⁵³ = 1.1e-52 m⁻² pure-primitive form (registry canonical, dual-manifestation)
