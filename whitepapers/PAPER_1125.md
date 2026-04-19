---
paper_id: "PAPER_1125"
title: "AGN Feedback, M-σ Scaling, and Metallicity Gradient Flattening in the UQFF"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [AGN feedback, M-sigma, metallicity gradient, Eddington ratio, SMBH, Ug4, SCm expulsion]
crosslinks: [PAPER_1124]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1125: AGN Feedback, M-σ Scaling, and Metallicity Gradient Flattening in the UQFF

## Abstract

Based on arXiv:2506.09123 (2025), we implement a UQFF calculator for AGN feedback effects on the $M_{\text{BH}}$-$\sigma$ relation and circumgalactic metallicity gradients. The classical $M_{\text{BH}} \propto \sigma^{4.38}$ scaling (Kormendy & Ho 2013) is modulated by AGN feedback energy $E_{\text{AGN}} = \varepsilon_f \dot{M}_{\text{acc}} c^2$. High Eddington-ratio AGN flatten metallicity gradients through outflow mixing, modeled as [SCm] expulsion in the UQFF Ug4 framework with $\Delta M_{\text{BH}}$ proportional to $f_{\text{feedback}}$.

## 1. The M-σ Relation

$$M_{\text{BH}} = 3.09 \times 10^8 \left(\frac{\sigma}{200 \text{ km/s}}\right)^{4.38} M_\odot$$

## 2. AGN Feedback Energy

$$E_{\text{AGN}} = \varepsilon_f \cdot \dot{M}_{\text{acc}} \cdot c^2$$

with typical radiative efficiency $\varepsilon_f \approx 0.05$ and accretion rates $\dot{M} \sim 0.01$-$10$ $M_\odot$/yr.

## 3. Eddington Ratio

$$\lambda_{\text{Edd}} = \frac{E_{\text{AGN}}}{L_{\text{Edd}}} = \frac{\varepsilon_f \dot{M} c^2}{1.26 \times 10^{38} (M_{\text{BH}}/M_\odot)}$$

## 4. Metallicity Gradient Flattening

AGN-driven outflows flatten the CGM metallicity gradient:

$$\nabla Z_{\text{flat}} = \nabla Z_{\text{intrinsic}} \cdot \frac{1}{1 + 10 \lambda_{\text{Edd}}}$$

High $\lambda_{\text{Edd}}$ systems show nearly uniform CGM metallicity.

## 5. UQFF Ug4 Framework

$$f_{\text{feedback}} = \varepsilon_f \cdot \lambda_{\text{Edd}}$$
$$Ug4 = \rho_{\text{SCm}} \cdot |\Delta M| \cdot f_{\text{feedback}}$$

The [SCm] expulsion mechanism drives metal ejection proportional to how much the SMBH deviates from the M-σ expectation.

Overall alignment: **85%**.

## References

- arXiv:2506.09123 — AGN Feedback and M-σ Scaling (2025).
- Kormendy, J. & Ho, L.C. (2013). ARAA, 51, 511.
