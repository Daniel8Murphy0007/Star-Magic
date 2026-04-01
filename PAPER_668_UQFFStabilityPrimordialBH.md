# PAPER_668: UQFF Primordial Black Hole Stability Analysis
**Subtitle:** Mass classification of PBHs under UQFF: stable_DM / marginal / evaporating. UQFF minimum DM-viable mass derived.
**Module:** UQFFStabilityPrimordialBH  
**Session:** Session 172  
**Date:** April 2, 2026  
**Version:** v5.29  
**Status:** Complete — CP4 #252 | UQFF Session 172

---

## Abstract
We systematically classify primordial black holes under UQFF by their evaporation prospects: stable dark matter, marginal, or evaporating. The minimum DM-viable mass shifts significantly compared to standard Hawking theory.

## 1. Standard Evaporation Timescale
$$\tau_{std} = \frac{5120\pi G^2 M^3}{\hbar c^4}$$

Standard boundary: τ_std = t_H (Hubble time) → M_boundary ≈ 5×10¹¹ kg.

## 2. UQFF Timescale
$$\tau_{UQFF} = \frac{\tau_{std}}{1-f_{TRZ}} \cdot \frac{\rho_{UA}}{\rho_{SCm}} \cdot \exp\!\left(\frac{U_m}{k_B T_H}\right)$$

Factor ≈ 30×.

## 3. Mass Classifications
| Category | Condition | Mass range (standard) | Mass range (UQFF) |
|----------|-----------|-----------------------|-------------------|
| stable_DM | τ_UQFF ≥ t_H | M ≥ 5×10¹¹ kg | M ≥ 1.6×10¹¹ kg |
| marginal | 0.1 t_H ≤ τ < t_H | narrow window | wider window |
| evaporating | τ < 0.1 t_H | M < 4×10¹¹ kg | M < 1.3×10¹¹ kg |

## 4. Numerical (M = 10¹² kg)
- τ_std ≈ 10¹⁰ yr 
- τ_UQFF ≈ 3×10¹¹ yr **> Hubble time** → viable DM

## 5. C++ Module
`UQFFStabilityPrimordialBH.h / .cpp` — Session 172
CP4 #252 — `UQFFStabilityPrimordialBHCalculator`


---
*PAPER_668 | Session 172 | Star-Magic UQFF Framework v5.29 | Daniel Murphy*
