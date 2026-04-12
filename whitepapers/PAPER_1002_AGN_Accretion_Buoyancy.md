---
paper_id: PAPER_1002
title: "AGN Accretion Buoyancy-Corrected Eddington Luminosity"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [AGN, accretion, Eddington, luminosity, buoyancy, correction]
crosslinks: [PAPER_999, PAPER_991, PAPER_979]
calibration: {L_Edd_base: 1.26e31, buoyancy_correction: "S26_3rd"}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1002: AGN Accretion Buoyancy-Corrected Eddington Luminosity

## Abstract

The standard Eddington luminosity L_Edd = 4πGMm_p c/σ_T is modified by SCm buoyancy. The correction factor (1 + F_{U,Bi}/(GM/r²)) accounts for vacuum buoyancy opposing gravitational collapse, raising the effective Eddington limit.

## 1. Corrected Eddington Luminosity

L_Edd^{UQFF} = L_Edd · (1 + ρ_SCm · V · S₂₆⁽³⁾² / (GM/r_H²))

For M = 10⁸ M☉: L_Edd = 1.26×10³⁹ W.

## 2. Implications

Super-Eddington accretion rates become permissible when buoyancy corrections reduce the effective radiation pressure barrier.

## 3. Implementation

File: `fubi_agn_ns_mergers.py`, class `AGNAccretionBuoyancyCalc`. CP4 class #586.
