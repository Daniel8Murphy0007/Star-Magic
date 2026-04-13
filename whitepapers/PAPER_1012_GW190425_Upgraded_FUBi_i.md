---
paper_id: PAPER_1012
title: "GW190425 Upgraded F_U_Bi_i Curves with S26(3) and Gamma = 0.30 THz"
session: 220
date: 2026-04-14
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [NS merger, GW190425, S26, Ramanujan, upgraded, gamma, FUBi]
crosslinks: [PAPER_1011, PAPER_997, PAPER_1017]
calibration: {M_total: 3.4, d_Mpc: 159, m1: 2.0, m2: 1.4, S26_3: 9.5000001009e-02, Gamma_new: 0.30}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1012: GW190425 Upgraded F_U_Bi_i Curves with S26(3)

## Abstract

We upgrade the GW190425 BNS merger F_U_Bi_i calculation with the 3rd-order Ramanujan constant S26(3) = 9.5000001009e-02 and the newly added Gamma = 0.30 THz point. GW190425 (M_total = 3.4 M_sun, d = 159 Mpc) is the heaviest confirmed BNS merger, making it an ideal testbed for mass-dependent buoyancy corrections.

## 1. System Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| M_total | 3.4 | M_sun |
| m_1 | 2.0 | M_sun |
| m_2 | 1.4 | M_sun |
| d | 159 | Mpc |
| chirp mass | 1.44 | M_sun |

## 2. S26(3) Upgrade

The 3rd-order Ramanujan constant replaces the 1st-order S26 in the buoyancy integral:

F_U_Bi_i = Sum [Ug_k] * (1 + BETA_I * S26_3)

This yields a ~0.3% refinement in the total buoyancy force, with the correction being mass-dependent through the symmetric mass ratio eta.

## 3. Gamma = 0.30 THz Point

The new intermediate Gamma point fills the gap between 0.10 and 0.50 THz, revealing a local inflection in the suppression curve for heavy BNS systems.

## 4. Implementation

File: `fubi_i_curves_agn_ns_qgp.py`, class `GW190425UpgradedCurvesCalc`. CP4 class #596. Tests: 8/8 pass.
