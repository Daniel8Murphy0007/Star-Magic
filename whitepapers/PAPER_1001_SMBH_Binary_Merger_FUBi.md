---
paper_id: PAPER_1001
title: "SMBH Binary Merger F_U_Bi with Phonon Damping"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [SMBH, binary, merger, inspiral, phonon, damping, F_U_Bi]
crosslinks: [PAPER_999, PAPER_1000, PAPER_993]
calibration: {M_bh: 5.5e7, a_spin: 0.70, damping_model: "phonon"}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1001: SMBH Binary Merger F_U_Bi with Phonon Damping

## Abstract

We compute F_U_Bi for SMBH binary inspiral systems where phonon damping modifies the merger waveform. The SCm buoyancy field introduces a frequency-dependent damping envelope via S₂₆⁽³⁾.

## 1. Core Equation

F_{U,Bi}^{binary} = ρ_SCm · V · S₂₆⁽³⁾² · ratio · (1 + Φ(Γ) · e^{-Γ/Γ₀})

The damping envelope suppresses high-frequency phonon modes during the final inspiral.

## 2. Results

For a 5.5×10⁷ M☉ primary: F_U_Bi = 5.47×10³⁵ m/s² with S₂₆⁽³⁾-calibrated buoyancy correction.

## 3. Implementation

File: `fubi_agn_ns_mergers.py`, class `SMBHBinaryMergerFUBiCalc`. CP4 class #585.
