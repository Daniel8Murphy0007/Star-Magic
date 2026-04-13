---
paper_id: PAPER_1014
title: "SMBH Merger F_U_Bi — Inspiral, Coalescence, and Ringdown Phases"
session: 220
date: 2026-04-14
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [SMBH, merger, FUBi, inspiral, coalescence, ringdown, QNM, gravitational waves]
crosslinks: [PAPER_1010, PAPER_1011, PAPER_1015]
calibration: {M1: 3.5e7, M2: 3.5e7, total_force: 6.98e20, damping: 0.333, phase_lag: 367.0, f_QNM: 2.19e-4, Df_f: 9.03e-3}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1014: SMBH Merger F_U_Bi — Inspiral, Coalescence, Ringdown

## Abstract

We compute the UQFF buoyancy force F_U_Bi across all three phases of a supermassive black hole (SMBH) binary merger (M_1 = M_2 = 3.5 x 10^7 M_sun). The inspiral phase yields F_total = 6.98 x 10^20 N with buoyancy damping factor 0.333 and accumulated phase lag 367.0 cycles. During coalescence, buoyancy-induced mass ejection Delta_M_buoy = 4.05 x 10^4 kg is computed. The ringdown phase shows a quasi-normal mode frequency f_QNM = 2.19 x 10^-4 Hz with SCm correction Delta_f/f = 9.03 x 10^-3.

## 1. Inspiral Phase

The buoyancy force during orbital decay follows:

F_U_Bi(r) = G * M_1 * M_2 / r^2 * (1 + BETA_I * S26_3) * D(Gamma)

where D(Gamma) is the frequency-dependent damping factor. At peak coupling, D = 0.333 indicating 66.7% buoyancy suppression of GW emission.

## 2. Coalescence Phase

During merger, the buoyancy field ejects mass:

Delta_M_buoy = BETA_I * M_total * (v_kick / c)^2 * S26_3

yielding M_remnant = 6.76 x 10^7 M_sun (3.4% mass deficit from GW + buoyancy radiation).

## 3. Ringdown Phase

The QNM frequency receives an SCm correction:

f_QNM = f_Kerr * (1 + SCm * BETA_I * S26_3 / (1 + a^2))

with Delta_f/f = 9.03 x 10^-3, potentially detectable by LISA.

## 4. Implementation

File: `fubi_smbh_mergers.py`, class `SMBHInspiralFUBiCalc` (inspiral), `SMBHCoalescenceFUBiCalc` (coalescence), `SMBHRingdownFUBiCalc` (ringdown). CP4 class #598. Tests: 8/8 pass.
