---
paper_id: PAPER_992
title: "GW190425 NS Merger F_U_Bi_i Curves — 47% Peak Strain Suppression"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [GW190425, NS_merger, strain, phonon, suppression, LIGO]
crosslinks: [PAPER_916, PAPER_989, PAPER_993]
calibration: {M_total: "3.4 Msun", d_Mpc: 159, peak_suppression: "47%"}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_992: GW190425 NS Merger F_U_Bi_i Curves

## Abstract

We compute $F_{U,\text{Bi}_i}$ numerical curves for the GW190425 neutron star merger ($M_{\text{total}} = 3.4\,M_\odot$, $d = 159\text{ Mpc}$) and derive the phonon-suppressed gravitational wave strain with 47% peak reduction at on-resonance $\Gamma = \Gamma_0$.

## 1. Strain Suppression

$$h_{\text{UQFF}} = h_{\text{GR}} \cdot S_{26} \cdot (1 - 0.47 \cdot e^0) = h_{\text{GR}} \cdot 0.530 \cdot S_{26}$$

at $\Gamma = \Gamma_0$ (on-resonance). The 47% factor comes from the buoyancy-to-gravity ratio at the merger radius.

## 2. Gamma Dependence

At 7 linewidth values, the strain reduction ranges from near-zero (far off-resonance, $\Gamma = 10\text{ THz}$) to the full 47% (on-resonance, $\Gamma = 0.1\text{ THz}$).

## 3. Mass-Gap Classification

With $m_1 = 2.52\,M_\odot$, the heavier component sits in the NS/BH mass gap. The phonon-modulated F_U_Bi_i provides additional classification power: $P(\text{NS}) = 49\%$, $P(\text{BH}) = 51\%$.

## 4. Implementation

File: `fubi_inside_outside.py`, class `GW190425FUBiCurves`. CP4 class #576.
