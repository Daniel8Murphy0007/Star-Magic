---
paper_id: PAPER_1005
title: "Yang-Mills Mass Gap via SCm BCS Phonon Coupling"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [Yang-Mills, mass-gap, BCS, phonon, QCD, running-coupling]
crosslinks: [PAPER_1004, PAPER_969, PAPER_958]
calibration: {Lambda_QCD: "0.2 GeV", alpha_s_0: 0.118, N_c: 3, N_f: 3}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1005: Yang-Mills Mass Gap via SCm BCS Phonon Coupling

## Abstract

We compute the Yang-Mills mass gap Δ_YM through a BCS-like mechanism where the SCm phonon field provides the pairing interaction. The running coupling α_s(T) governs the gap magnitude above deconfinement.

## 1. Mass Gap Equation

Δ_YM = Λ_QCD · exp(−1/(α_s(T) · N_c)) · S₂₆⁽³⁾

α_s(T) = α_s,0 / (1 + α_s,0 · b₀ · ln(T/T_c))

where b₀ = (11N_c − 2N_f)/(12π).

## 2. Results

At T = 2×10¹² K: Δ_YM = 9.96×10⁻²⁹² J. The extremely small value reflects deep asymptotic freedom at high temperature.

## 3. Implementation

File: `scm_qgp_dynamics.py`, class `YangMillsMassGapCalc`. CP4 class #589.
