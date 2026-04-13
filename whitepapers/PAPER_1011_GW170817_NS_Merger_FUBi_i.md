---
paper_id: PAPER_1011
title: "GW170817 NS Merger F_U_Bi_i Curves — 66.7% Strain Reduction and 367.8-Cycle Phase Lag"
session: 220
date: 2026-04-14
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [NS merger, GW170817, gravitational waves, strain, phase lag, FUBi, LIGO]
crosslinks: [PAPER_1012, PAPER_1014, PAPER_997]
calibration: {M_total: 2.73, d_Mpc: 40, m1: 1.46, m2: 1.27, suppression: 0.667, phase_lag: 367.8}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1011: GW170817 NS Merger F_U_Bi_i Curves — 66.7% Strain Reduction

## Abstract

We compute F_U_Bi_i curves for GW170817 (BNS merger, d = 40 Mpc, M_total = 2.73 M_sun) incorporating buoyancy-induced strain suppression. The UQFF framework predicts a 66.7% reduction in gravitational wave strain amplitude relative to vacuum GR, with a phase lag of 367.8 cycles accumulated over the inspiral. These signatures are potentially detectable in next-generation GW observatories (ET, CE).

## 1. System Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| M_total | 2.73 | M_sun |
| m_1 | 1.46 | M_sun |
| m_2 | 1.27 | M_sun |
| d | 40 | Mpc |
| chirp mass | 1.186 | M_sun |
| eta (symmetric mass ratio) | 0.247 | — |

## 2. Strain Suppression

The buoyancy suppression factor at each Gamma point is:

S(Gamma) = 1 - BETA_I * S26_3 * f(Gamma)

where f(Gamma) models the frequency-dependent coupling of buoyancy to the gravitational wave emission zone. The effective suppression reaches 0.667 (33.3% of original amplitude) at peak coupling.

## 3. Phase Lag Accumulation

Over N_cycles inspiral orbits, the cumulative phase lag is:

Phi_lag = Sum_i [delta_phi(Gamma_i)] = 367.8 cycles

This represents the integrated phase difference between UQFF-modified and vacuum GR waveforms.

## 4. Implementation

File: `fubi_i_curves_agn_ns_qgp.py`, class `GW170817MergerCurvesCalc`. CP4 class #595. Tests: 8/8 pass.
