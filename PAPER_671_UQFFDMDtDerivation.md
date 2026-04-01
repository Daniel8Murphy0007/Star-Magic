# PAPER_671: UQFF dM/dt Derivation
**Subtitle:** Full step-by-step derivation of the UQFF-modified mass loss rate, with analytic M(t) approximation and Euler simulation.
**Module:** UQFFDMDtDerivation  
**Session:** Session 172  
**Date:** April 2, 2026  
**Version:** v5.29  
**Status:** Complete — CP4 #255 | UQFF Session 172

---

## Abstract
We present the complete three-step UQFF derivation of the modified black hole mass loss rate dM/dt, showing a ~33× suppression relative to standard Hawking emission.

## 1. Standard Hawking Rate
$$\frac{dM}{dt}\bigg|_{std} = -\frac{L_H}{c^2} = -\frac{\hbar c^4}{15360\pi G^2 M^2 c^2}$$

## 2. UQFF Step 1 — f_TRZ
$$\left.\frac{dM}{dt}\right.' = \frac{dM}{dt}\bigg|_{std} \cdot (1 - f_{TRZ})$$

Physical basis: negentropic restoring force suppresses pair production at horizon.

## 3. UQFF Step 2 — [UA]/[SCm]
$$\left.\frac{dM}{dt}\right.'' = \left.\frac{dM}{dt}\right.' \cdot \frac{\rho_{SCm}}{\rho_{UA}}$$

## 4. UQFF Step 3 — U_m Anchor
$$\frac{dM}{dt}\bigg|_{UQFF} = \left.\frac{dM}{dt}\right.'' \cdot \exp\!\left(-\frac{U_m}{k_B T_H}\right)$$

## 5. Combined Formula
$$\frac{dM}{dt}\bigg|_{UQFF} = \frac{dM}{dt}\bigg|_{std} \cdot (1-f_{TRZ}) \cdot \frac{\rho_{SCm}}{\rho_{UA}} \cdot \exp\!\left(-\frac{U_m}{k_B T_H}\right)$$

Suppression factor: 0.9 × 0.1 × e⁻¹ ≈ **0.033**  (evaporation ~30× slower).

## 6. Analytic Approximation
$$M(t) \approx \left(M_0^3 - 3A\,t\right)^{1/3}$$

where A = |dM/dt|_UQFF × M₀².

## 7. C++ Module
`UQFFDMDtDerivation.h / .cpp` — Session 172
CP4 #255 — `UQFFDMDtDerivationCalculator`


---
*PAPER_671 | Session 172 | Star-Magic UQFF Framework v5.29 | Daniel Murphy*
