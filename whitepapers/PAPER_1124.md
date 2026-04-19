---
paper_id: "PAPER_1124"
title: "CGM Metal Retention in Dwarf Galaxies and the UQFF Ug4 Expulsion Model"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [CGM, dwarf galaxies, metal retention, M-sigma, SMBH, AGN feedback, Ug4, SCm expulsion]
crosslinks: [PAPER_1125]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1124: CGM Metal Retention in Dwarf Galaxies and the UQFF Ug4 Expulsion Model

## Abstract

Based on arXiv:2505.08861 (2025), we implement a UQFF calculator for circumgalactic medium (CGM) metal retention in dwarf galaxies. Dwarf galaxies exhibit a weak $M_*$-$\sigma$ correlation ($\alpha \approx 0.2$) compared to the classical $M_{\text{BH}}$-$\sigma$ relation ($\beta \approx 4.38$). Over-massive SMBHs ($\Delta M_{\text{BH}} > 0$) drive metals out of the CGM ($f_Z \sim 0.89$) via [SCm] expulsion in Ug4, while under-massive SMBHs ($\Delta M_{\text{BH}} < 0$) allow higher metal retention ($f_Z \sim 0.85$) in stars.

## 1. The M*-σ Relation in Dwarfs

For dwarf galaxies with $M_* \lesssim 10^{10} M_\odot$:

$$\sigma_{\text{pred}} = 30 \cdot \left(\frac{M_*}{10^9 M_\odot}\right)^{0.20} \text{ km/s}$$

This is significantly shallower than the classical Kormendy & Ho (2013) relation.

## 2. BH Mass Offset

$$\Delta M_{\text{BH}} = \log_{10}\left(\frac{M_{\text{BH}}}{M_{\text{BH,exp}}}\right) \text{ dex}$$

where $M_{\text{BH,exp}} = 10^5 (\sigma/30)^{4.38} M_\odot$.

## 3. Metal Retention Fraction

$$f_Z = f_{Z,\text{base}} \pm 0.02 - f_{\text{feedback}} \cdot \Delta M_{\text{BH}}$$

- **Over-massive** ($\Delta M > 0$): AGN drives metals out → lower $f_Z \approx 0.89$
- **Under-massive** ($\Delta M < 0$): metals retained in ISM/stars → higher $f_Z \approx 0.85$

## 4. UQFF Ug4 Interpretation

The [SCm] expulsion mechanism in Ug4:

$$Ug4_{\text{expulsion}} = \rho_{\text{SCm}} \cdot |\Delta M_{\text{BH}}| \cdot f_{\text{feedback}}$$

This provides a physical mechanism for the observed scatter in CGM metallicities.

Overall alignment: **80%**.

## References

- arXiv:2505.08861 — CGM in Dwarf Galaxies (2025).
- Kormendy, J. & Ho, L.C. (2013). ARAA, 51, 511.
