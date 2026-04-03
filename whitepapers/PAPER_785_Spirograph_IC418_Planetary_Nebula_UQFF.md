# PAPER_785: Spirograph Nebula IC 418 — UQFF Planetary Nebula Fast Wind

**Author:** Daniel T. Murphy  
**Framework:** UQFF (Universal Quantum Field Superconductive Framework)  
**Session:** 181 | v5.42  
**Date:** 2026  
**CP4 Class:** #369 — SpirographNebulaIC418UQFFCalculator  

---

## Abstract

IC 418, the "Spirograph Nebula," is one of the most intricately structured planetary nebulae, located ~2,000 light-years away (z ≈ 0.0007) in the constellation Lepus. Hubble WFPC2 imaging reveals a complex pattern of nested shells and radial spokes resembling a spirograph drawing. The central white dwarf (T_eff ~36,000 K) drives a fast stellar wind at ~1,500 km/s (v = 1.5×10⁶ m/s), ionizing the ejected AGB envelope. Under UQFF, the fast wind velocity (v = 1.5×10⁶ m/s) with standard PN B-field (B = 10⁻⁵ T) yields g_IC418 ≈ 1.580×10⁻² m/s².

---

## 1. Introduction

IC 418 represents a late-stage AGB star (~1–3 M☉ progenitor) that ejected its outer envelope ~3,000 years ago, creating the current planetary nebula shell. The Hubble image reveals the "spirograph" interference pattern from multiple overlapping AGB pulsation shells that were ejected at slightly different angles. The central star's fast wind at ~1,500 km/s (confirmed by UV spectroscopy) is the highest drive velocity in the slow/fast wind PN interaction model. UQFF encodes this through v = 1.5×10⁶ m/s, which is 1.5× the LBV eruptive velocity (and thus 15× the standard HII velocity), yielding g = 1.580×10⁻² m/s².

---

## 2. Master UQFF Gravity Equation

```
g_IC418(r, t) = (G × M) / r² × (1 + H(z)×t) × (1 - E_rad) × (1 + f_TRZ) + a_EM
```

### 2.1 Parameters

| Parameter | Symbol | Value | Source |
|-----------|--------|-------|--------|
| Nebula mass (envelope) | M | ~0.6 M☉ = 1.193×10³⁰ kg | Hubble |
| Nebula radius | r | 1×10¹⁶ m (~0.1 pc) | Hubble |
| Age | t | ~3,000 yr = 9.468×10¹⁰ s | Expansion age |
| E_rad | — | 0.20 | EUV photoionization loss |
| Redshift | z | 0.0007 | Distance |
| v_EM | v | 1.5×10⁶ m/s | Central star fast wind |
| B_EM | B | 10⁻⁵ T | PN field |
| f_TRZ | — | 0.05 | UQFF |

---

## 3. Long-Form Derivation

### Step 1: Base Gravitational Term
```
g_grav = 6.6743e-11 × 1.193e30 / (1e16)² = 7.962e-13 m/s²
```

### Step 2: Cosmic Expansion Factor
```
H(z)×t ≈ negligible (2.268e-18 × 9.468e10 ≈ 2.15e-7); factor ≈ 1.0000002
```

### Step 3: Radiation Energy Loss (EUV Photoionization)
```
E_rad = 0.20; 1 - E_rad = 0.80
```

### Step 4: Time-Reversal Correction
```
f_TRZ = 0.05; 1 + f_TRZ = 1.05
```

### Step 5: Gravitational Total
```
g_grav_total = 7.962e-13 × 1.0 × 0.80 × 1.05 = 6.689e-13 m/s²
```

### Step 6: Aether EM Correction (Fast Wind)
```
v = 1.5×10⁶ m/s, B = 10⁻⁵ T
a_EM = (1.602e-19 × 1.5e6 × 1e-5 / 1.673e-27) × 11 × 1e-12 = 1.580e-2 m/s²
```

### Step 7: Final Solution
```
g_IC418 = 6.689e-13 + 1.580e-2 ≈ 1.580e-2 m/s²
```

---

## 4. Physical Interpretation

IC 418's result (g = 1.580×10⁻²) falls precisely between the standard HII (1.053×10⁻³) and LBV eruptive (1.053×10⁻²) regimes. The 1.5× velocity multiplier (1.5×10⁶ vs. 1.0×10⁶ m/s) directly yields a 1.5× larger result. IC 418 establishes the "fast wind PN" UQFF subcategory at v = 1.5×10⁶ m/s and will be shared with NGC 6307+7027 (PAPER_788) and M57 (PAPER_791) as the canonical planetary nebula fast-wind value.

---

## 5. Conclusions

UQFF applied to IC 418 Spirograph Nebula yields g_IC418 ≈ 1.580×10⁻² m/s², establishing the planetary nebula fast-wind UQFF parameter class. The 1.5×10⁶ m/s central star wind velocity is directly observed in UV spectroscopy, providing an observational anchor for the PN fast-wind UQFF subcategory.

*PAPER_785, CP4 class #369. v5.42.*
