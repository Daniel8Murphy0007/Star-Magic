---
paper_id: PAPER_991
title: "Centaurus A AGN F_U_Bi_i Curves — 7-Point Gamma Sweep with Jet Modulation"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [Centaurus_A, AGN, jet, F_U_Bi_i, Gamma, sweep, CenA, SMBH]
crosslinks: [PAPER_989, PAPER_993, PAPER_930]
calibration: {M_BH: "5.5e7 Msun", a_spin: 0.70, B_gauss: 3000, gamma_points: 7}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_991: Centaurus A AGN F_U_Bi_i Curves

## Abstract

We compute numerical $F_{U,\text{Bi}_i}$ curves for Centaurus A ($M_{\text{BH}} = 5.5 \times 10^7\,M_\odot$, $a = 0.70$, $B = 3000\text{ G}$) across 7 linewidth values $\Gamma \in \{0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0\}\text{ THz}$.

## 1. Jet Modulation Factor

$$M_{\text{jet}}(\Gamma) = 1 + A_{\text{jet}} \cdot \exp\!\left(-\frac{(\Gamma - \Gamma_0)^2}{2\sigma_\Gamma^2}\right)$$

with $A_{\text{jet}} = 1.5$, $\Gamma_0 = 2\pi \times 0.1\text{ THz}$, $\sigma_\Gamma = 0.08 \times 2\pi\text{ THz}$.

## 2. Jet Power

$$P_{\text{jet}}(\Gamma) = \frac{B^2}{8\pi} \left(\frac{r_H}{c}\right)^2 a^2 c \cdot M_{\text{jet}}(\Gamma)$$

## 3. Combined F_U_Bi_i at Horizon

$$F_{U,\text{Bi}_i}^{\text{CenA}}(\Gamma) = U_g(r_H) - U_b(r_H) + P_{\text{jet}}(\Gamma) \times 10^{-45}$$

The jet power coupling at $10^{-45}$ ensures dimensional consistency with the acceleration-scale F_U_Bi_i.

## 4. Implementation

File: `fubi_inside_outside.py`, class `CentaurusAFUBiCurves`. CP4 class #575.
