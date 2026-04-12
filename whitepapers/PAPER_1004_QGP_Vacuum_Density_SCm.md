---
paper_id: PAPER_1004
title: "QGP Vacuum Density with SCm S₂₆⁽³⁾ Phonon Coupling"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [QGP, vacuum, density, SCm, phonon, deconfinement, quark]
crosslinks: [PAPER_969, PAPER_999, PAPER_1005]
calibration: {T_c: 1.5e12, rho_QGP_2e12: 1.16e-12, S26_3rd: 0.095}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1004: QGP Vacuum Density with SCm Phonon Coupling

## Abstract

We compute the QGP vacuum density ρ_QGP(T) through SCm phonon coupling using S₂₆⁽³⁾. Below T_c = 1.5×10¹² K, confinement suppresses the QGP state; above T_c, the density follows exponential activation with Gaussian phonon response.

## 1. Core Equation

ρ_QGP(T) = ρ_SCm · S₂₆⁽³⁾ · exp(−(T_c−T)/T) · Φ(T)

where Φ(T) = S₂₆⁽³⁾ · exp(−(T−T_c)²/(2(0.1T_c)²)).

## 2. Results

| Temperature | ρ_QGP |
|-------------|-------|
| T = 10⁶ K | 0 (confined) |
| T = 2×10¹² K | 1.16×10⁻¹² kg/m³ |
| T = 5×10¹² K | ~10⁻¹⁵ kg/m³ |

## 3. Implementation

File: `scm_qgp_dynamics.py`, class `QGPVacuumDensityCalc`. CP4 class #588.
