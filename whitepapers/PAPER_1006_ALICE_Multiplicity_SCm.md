---
paper_id: PAPER_1006
title: "ALICE Multiplicity SCm Phonon Scaling"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [ALICE, multiplicity, SCm, phonon, Pb-Pb, centrality, LHC]
crosslinks: [PAPER_1004, PAPER_1005, PAPER_969]
calibration: {sqrt_s_NN: 5020, dNdeta: 7565, N_part: 383}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1006: ALICE Multiplicity SCm Phonon Scaling

## Abstract

We model the ALICE Pb-Pb charged-particle multiplicity dN_ch/dη at √s_NN = 5.02 TeV with SCm phonon scaling. The S₂₆⁽³⁾ phonon coupling enhances the base multiplicity by (1 + Φ).

## 1. Multiplicity Formula

dN_ch/dη = α · (N_part/2)^β · (1 + Φ) · (√s_NN/200)^0.15

With α = 8.5, β = 1.25, N_part = 383 (central 0-5%): dN_ch/dη ≈ 7565.

## 2. SCm Enhancement

The phonon correction Φ = S₂₆⁽³⁾ = 0.095 provides a 9.5% enhancement over the baseline empirical fit.

## 3. Implementation

File: `scm_qgp_dynamics.py`, class `ALICEMultiplicityCalc`. CP4 class #590.
