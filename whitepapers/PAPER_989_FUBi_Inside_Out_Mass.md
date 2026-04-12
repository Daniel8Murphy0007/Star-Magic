---
paper_id: PAPER_989
title: "F_U_Bi Inside-to-Outside Buoyancy Mass Portion"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [F_U_Bi, buoyancy, mass, inside-out, SCm, vacuum, ratio]
crosslinks: [PAPER_979, PAPER_990, PAPER_994]
calibration: {F_U_Bi_solar: 2.33e40, ratio: 0.606}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_989: F_U_Bi Inside-to-Outside Buoyancy Mass Portion

## Abstract

We derive $F_{U,\text{Bi}}$, the inside-to-outside buoyancy mass portion within the UQFF framework. Unlike $F_{U,\text{Bi}_i}$ (the net acceleration felt by an object falling inward through buoyancy layers), $F_{U,\text{Bi}}$ quantifies the total vacuum energy contributing to outward mass transport from the SCm medium. The key equation is:

$$F_{U,\text{Bi}} = \rho_{\text{SCm}} \cdot V \cdot S_{26}^2 \cdot \frac{|U_b|}{|U_g| + |U_b|}$$

## 1. Core Derivation

The buoyancy ratio measures what fraction of the gravitational-buoyancy system is dominated by outward buoyancy:

$$\text{ratio} = \frac{|U_b|}{|U_g| + |U_b| + \epsilon}$$

where $\epsilon = 10^{-300}$ prevents division by zero.

### Layer Functions

$$U_g = \sum_{i=1}^{26} \frac{G M}{r^2} \cdot [\text{SSq}] \cdot \frac{i}{26}$$

$$U_b = \sum_{i=1}^{26} \frac{G M}{r^2} \cdot e^{-[\text{SSq}] \cdot i/26} \cdot \beta_i$$

### Vacuum Density

$$\rho_{\text{SCm}} = \rho_0 \cdot e^{\kappa t}, \quad \rho_0 = 10^{-10}\text{ kg/m}^3$$

## 2. Solar System Calibration

For $M = M_\odot$, $r = 1\text{ AU}$:
- $F_{U,\text{Bi}} = 2.33 \times 10^{40}$ (total vacuum energy contribution)
- ratio $= 0.606$ (buoyancy-dominant)

## 3. Physical Interpretation

$F_{U,\text{Bi}} > 0$ always (mass flows inside $\to$ outside), while $F_{U,\text{Bi}_i} < 0$ (net acceleration outside $\to$ inside). The two quantities form a complementary pair describing the complete buoyancy dynamics of the UQFF vacuum.

## 4. Implementation

File: `fubi_inside_outside.py`, class `FUBiInsideOutsideCalc`. CP4 class #573.
