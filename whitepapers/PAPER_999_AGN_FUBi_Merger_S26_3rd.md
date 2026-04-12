---
paper_id: PAPER_999
title: "AGN F_U_Bi Merger with 3rd-Order Ramanujan S₂₆⁽³⁾"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [AGN, F_U_Bi, merger, Ramanujan, S26, 3rd-order, SMBH]
crosslinks: [PAPER_991, PAPER_979, PAPER_989]
calibration: {S26_3rd: 0.095, M_bh: 5.5e7, a_spin: 0.70}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_999: AGN F_U_Bi Merger with S₂₆⁽³⁾

## Abstract

We derive the buoyancy force F_U_Bi at the SMBH horizon of AGN merger systems using the 3rd-order Ramanujan summation S₂₆⁽³⁾ = 0.095. The recursive coefficient R_n^{(d,k)} = Σ_{j=0}^{k-1} (-1)^j C(k-1,j)/(n+j)! provides alternating-sign corrections that reduce S₂₆ from 0.57 (1st-order) to 0.095 (3rd-order).

## 1. Core Equation

F_{U,Bi}^{AGN} = ρ_SCm · V · S₂₆⁽³⁾² · |Ub|/(|Ug| + |Ub|)

For Centaurus A (M_BH = 5.5×10⁷ M☉, a = 0.70): F_U_Bi = 5.47×10³⁵ m/s².

## 2. Ramanujan S₂₆⁽³⁾

S₂₆⁽³⁾ = Σ_{n=1}^{26} (z^n/n²⁶) · R_n^{(26,3)} = 9.500×10⁻²

Distinct from S₂₆⁽¹⁾ = 0.57 due to higher-order alternating cancellation.

## 3. Implementation

File: `fubi_agn_ns_mergers.py`, class `AGNFUBiMergerCalc`. CP4 class #583.
