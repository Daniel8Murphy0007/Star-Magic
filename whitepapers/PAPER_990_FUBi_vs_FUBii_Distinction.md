---
paper_id: PAPER_990
title: "F_U_Bi vs F_U_Bi_i Distinction — Direction, Magnitude, Dimensionality"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [F_U_Bi, F_U_Bi_i, distinction, direction, buoyancy, sign]
crosslinks: [PAPER_989, PAPER_979, PAPER_991]
calibration: {F_U_Bi: 2.33e40, F_U_Bi_i: -2.41e-02}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_990: F_U_Bi vs F_U_Bi_i Distinction

## Abstract

We formalize the distinction between two complementary UQFF buoyancy quantities that are frequently confused:

| Property | $F_{U,\text{Bi}}$ | $F_{U,\text{Bi}_i}$ |
|----------|--------------------|----------------------|
| Direction | Inside → Outside | Outside → Inside (net) |
| Sign | Always positive | Typically negative |
| Units | Energy (J) | Acceleration (m/s²) |
| Magnitude | $\sim 10^{40}$ | $\sim 10^{-2}$ |
| Physical meaning | Vacuum mass transport | Net gravitational acceleration |

## 1. F_U_Bi (Inside-to-Outside)

$$F_{U,\text{Bi}} = \rho_{\text{SCm}} \cdot V \cdot S_{26}^2 \cdot \frac{|U_b|}{|U_g| + |U_b|}$$

This is the total vacuum energy flowing outward through the 26-layer buoyancy structure. It is always positive and cosmologically large because it includes $V_{\text{region}} \sim 10^{48}\text{ m}^3$.

## 2. F_U_Bi_i (Outside-to-Inside)

$$F_{U,\text{Bi}_i} = U_g + U_m + U_A - U_b + F_n \cdot S_{26} \cdot \Phi \cdot E_{\text{net}}$$

This is the per-particle net acceleration through the 6-layer canonical structure. Negative values indicate net inward pull (gravity dominates buoyancy at the acceleration level).

## 3. Complementarity

The ratio $|F_{U,\text{Bi}_i}| / F_{U,\text{Bi}}$ gives the fractional acceleration per unit vacuum energy — the "efficiency" of gravitational focusing through the buoyancy medium.

## 4. Implementation

File: `fubi_inside_outside.py`, class `FUBiDistinctionCalc`. CP4 class #574.
