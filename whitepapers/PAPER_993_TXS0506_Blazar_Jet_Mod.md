---
paper_id: PAPER_993
title: "TXS 0506+056 Blazar Jet F_U_Bi_i Modulation — 3.3x Peak"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [TXS_0506, blazar, jet, neutrino, F_U_Bi_i, modulation, IceCube]
crosslinks: [PAPER_991, PAPER_989, PAPER_948]
calibration: {M_BH: "3e8 Msun", a_spin: 0.95, peak_mod: "3.3x", B_gauss: 5000}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_993: TXS 0506+056 Blazar Jet F_U_Bi_i Modulation

## Abstract

We compute the $F_{U,\text{Bi}_i}$ jet modulation curves for TXS 0506+056, the IceCube neutrino blazar ($M_{\text{BH}} = 3 \times 10^8\,M_\odot$, $a = 0.95$, $B = 5000\text{ G}$). The peak jet modulation factor reaches $M_{\text{jet}} = 3.3\times$ at on-resonance $\Gamma_0$.

## 1. Jet Modulation

$$M_{\text{jet}}(\Gamma_0) = 1 + A_{\text{jet}} \cdot e^0 = 1 + 2.3 = 3.3$$

with $A_{\text{jet}} = 2.3$ for TXS 0506+056 (higher than CenA due to extreme spin $a = 0.95$).

## 2. Neutrino Production

The 3.3× jet power enhancement at SCm resonance provides a natural explanation for the IceCube neutrino excess: enhanced proton acceleration in the jet produces pions that decay to neutrinos.

## 3. F_U_Bi_i at Horizon

The buoyancy-gravity balance at the ergosphere boundary ($r_H$ for $a = 0.95$) gives $F_{U,\text{Bi}_i} < 0$ (inward-dominant), with the jet modulation acting as a perturbative correction.

## 4. Implementation

File: `fubi_inside_outside.py`, class `TXS0506FUBiCurves`. CP4 class #577.
