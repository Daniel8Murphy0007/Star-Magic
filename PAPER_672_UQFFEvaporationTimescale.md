# PAPER_672: UQFF Evaporation Timescale
**Subtitle:** Full UQFF evaporation timescale formula, mass boundary calculations for universe-age crossings, and U_m sensitivity sweep.
**Module:** UQFFEvaporationTimescale  
**Session:** Session 172  
**Date:** April 2, 2026  
**Version:** v5.29  
**Status:** Complete — CP4 #256 | UQFF Session 172

---

## Abstract
We derive the complete UQFF black hole evaporation timescale and calculate the mass at which τ_UQFF equals the Hubble time, demarcating the boundary between stable and evaporating black holes.

## 1. Standard Timescale
$$\tau_{Hawking} = \frac{5120\pi G^2 M^3}{\hbar c^4}$$

## 2. UQFF Timescale
$$\tau_{UQFF} = \frac{\tau_{Hawking}}{1-f_{TRZ}} \cdot \frac{\rho_{UA}}{\rho_{SCm}} \cdot \exp\!\left(\frac{U_m}{k_B T_H}\right)$$

Factor ≈ **30×** (1.11 × 10 × 2.718).

## 3. Universe-Age Boundary Masses
Standard: $M_{cross,std} = \left(\frac{\hbar c^4 t_H}{5120\pi G^2}\right)^{1/3} \approx 5.5\times10^{11}$ kg

UQFF: $M_{cross,UQFF} = M_{cross,std} / (30)^{1/3} \approx 1.8\times10^{11}$ kg

UQFF shifts the evaporation boundary 3.1× lower in mass.

## 4. Sensitivity to U_m
| U_m/k_BT_H | τ factor over τ_std |
|-------------|---------------------|
| 0 | 11.1 |
| 1 | 30.2 |
| 2 | 82 |
| 3 | 222 |
| 5 | 1660 |

## 5. Physical Implications
PBHs in the previously evaporating mass range 1.8×10¹¹–5.5×10¹¹ kg become stable under UQFF.

## 6. C++ Module
`UQFFEvaporationTimescale.h / .cpp` — Session 172
CP4 #256 — `UQFFEvaporationTimescaleCalculator`


---
*PAPER_672 | Session 172 | Star-Magic UQFF Framework v5.29 | Daniel Murphy*
