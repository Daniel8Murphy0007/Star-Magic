---
paper_id: PAPER_1007
title: "Deconfinement Phase Diagram with SCm Phonon Boundary"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [deconfinement, phase-diagram, QGP, hadron, SCm, phonon, boundary]
crosslinks: [PAPER_1004, PAPER_1005, PAPER_1006]
calibration: {T_c: 1.5e12, phases: ["hadron", "mixed", "QGP"]}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1007: Deconfinement Phase Diagram with SCm Boundary

## Abstract

We construct the T vs μ_B phase diagram for the hadron-QGP transition with an SCm phonon boundary. The phase transition at T_c = 1.5×10¹² K is modified by the running coupling α_s(T) and the S₂₆⁽³⁾ phonon response.

## 1. Phase Classification

| Temperature | Phase |
|-------------|-------|
| T < 0.9 T_c | Hadron (confined) |
| 0.9 T_c ≤ T ≤ 1.1 T_c | Mixed (crossover) |
| T > 1.1 T_c | QGP (deconfined) |

## 2. SCm Phonon Boundary

The on-resonance phonon factor Φ(Γ₀) = S₂₆⁽³⁾ = 0.095 modulates the crossover width.

## 3. Implementation

File: `scm_qgp_dynamics.py`, class `DeconfinementPhaseDiagramCalc`. CP4 class #591.
