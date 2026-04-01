# PAPER_661: UQFF Primordial Black Hole Dark Matter
**Subtitle:** Demonstrates UQFF elevates PBH lifetimes by ~30x, reopening the dark matter window for masses M ~ 10¹⁰–10¹⁵ g.
**Module:** UQFFPBHDarkMatter  
**Session:** Session 172  
**Date:** April 2, 2026  
**Version:** v5.29  
**Status:** Complete — CP4 #245 | UQFF Session 172

---

## Abstract
Primordial black holes (PBHs) with masses M < 10¹² kg evaporate via Hawking radiation before the present epoch, eliminating them as dark matter candidates. We show that UQFF elevates their evaporation timescale by a factor of ~11–30, converting many previously evaporating PBHs into viable dark matter candidates.

## 1. Standard Hawking Evaporation
$$\tau_{std} = \frac{5120\pi G^2 M^3}{\hbar c^4}$$

For M = 10¹² kg: τ_std ≈ 10¹⁰ yr (comparable to Hubble time).

## 2. UQFF Lifetime Enhancement

### Step 1 — f_TRZ Negentropic Extension
τ' = τ_std / (1 − f_TRZ) ≈ 1.11 × τ_std

### Step 2 — [UA]/[SCm] Aether Suppression
τ'' = τ' × (ρ_UA/ρ_SCm) ≈ 10 × τ'

### Step 3 — U_m String Anchoring
τ_UQFF = τ'' × exp(U_m / k_B T_H)

$$\tau_{UQFF} = \frac{\tau_{std}}{1-f_{TRZ}} \cdot \frac{\rho_{UA}}{\rho_{SCm}} \cdot \exp\!\left(\frac{U_m}{k_B T_H}\right)$$

Total factor: ~30× for typical U_m/k_BT_H ≈ 1.

## 3. Dark Matter Window
| Mass (g) | τ_std | τ_UQFF | DM candidate (UQFF) |
|----------|-------|--------|---------------------|
| 10¹⁰ | < t_H | > 0.1 t_H | marginal |
| 10¹² | ~ t_H | ~30 t_H | stable |
| 10¹⁵ | >> t_H | >> t_H | stable |

UQFF reopens M ~ 10¹⁰–10¹⁵ g as a dark matter window.

## 4. C++ Module
`UQFFPBHDarkMatter.h / .cpp` — Session 172
CP4 #245 — `UQFFPBHDarkMatterCalculator`

## References
- Carr, B. & Hawking, S.W. (1974). *MNRAS* 168, 399.
- UQFF calibrated constants: PAPER_631 (Millennium Prize context).


---
*PAPER_661 | Session 172 | Star-Magic UQFF Framework v5.29 | Daniel Murphy*
