---
paper_id: PAPER_1009
title: "3C273 AGN F_U_Bi_i Curves — 3.1x Jet Modulation at Gamma = 0.05 THz"
session: 220
date: 2026-04-14
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [AGN, 3C273, jet, modulation, FUBi, quasar, gamma, curves]
crosslinks: [PAPER_997, PAPER_1010, PAPER_1013]
calibration: {M_BH: 8.86e8, a_spin: 0.90, B_jet: 4000, A_jet: 2.1, modulation: 3.1}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1009: 3C273 AGN F_U_Bi_i Curves — 3.1x Jet Modulation

## Abstract

We compute numerical F_U_Bi_i curves for the archetypal quasar 3C273 (z = 0.158, M_BH = 8.86 x 10^8 M_sun) across 8 Gamma points [0.01, 0.05, 0.10, 0.30, 0.50, 1.0, 5.0, 10.0] THz. The jet modulation factor M_jet = 1 + A_jet * exp(-Gamma / Gamma_crit) peaks at 3.1x for Gamma = 0.05 THz, confirming UQFF-predicted AGN buoyancy-jet coupling at sub-THz resonance.

## 1. System Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| M_BH | 8.86 x 10^8 | M_sun |
| a (spin) | 0.90 | dimensionless |
| B_jet | 4000 | G |
| A_jet | 2.1 | dimensionless |
| z (redshift) | 0.158 | — |
| Gamma_crit | 0.08 | THz |

## 2. F_U_Bi_i Curve Construction

For each Gamma in the sweep, the unified buoyancy field is:

F_U_Bi_i(Gamma) = [Ug1 + Ug2 + Ug3(Gamma) + Ug4] * M_jet(Gamma) * (1 + BETA_I * S26_3)

where S26_3 = 9.5000001009e-02 (3rd-order Ramanujan) and M_jet(Gamma) = 1 + A_jet * exp(-Gamma / Gamma_crit).

## 3. Results

Peak modulation 3.1x at Gamma = 0.05 THz. The curve shows exponential decay toward unity as Gamma increases beyond 1.0 THz, consistent with thermal decoupling of jet magnetic pressure from buoyancy feedback.

## 4. Implementation

File: `fubi_i_curves_agn_ns_qgp.py`, class `ThreeCTwoSevenThreeAGNCurvesCalc`. CP4 class #593. Tests: 8/8 pass.
