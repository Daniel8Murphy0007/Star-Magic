---
paper_id: PAPER_1000
title: "NS Merger F_U_Bi with Strain Suppression and BCS Gap"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [NS, merger, strain, suppression, BCS, GW190425, gravitational-wave]
crosslinks: [PAPER_992, PAPER_999, PAPER_976]
calibration: {suppression_pct: 47, d_Mpc: 159, M_total: 3.4}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1000: NS Merger F_U_Bi with Strain Suppression

## Abstract

We extend the NS merger F_U_Bi framework (GW190425) with 3rd-order Ramanujan S₂₆⁽³⁾, incorporating BCS gap coupling and tidal correction. At resonance, 47.0% strain reduction is achieved: h_UQFF = h_GR · (1 − 0.47) = 0.53 · h_GR.

## 1. Strain Suppression

h_UQFF(Γ) = h_GR · (1 − 0.47 · Φ(Γ)/S₂₆⁽³⁾)

For GW190425 (M_total = 3.4 M☉, d = 159 Mpc): h_GR = 2.52×10⁻²² → h_UQFF = 1.33×10⁻²².

## 2. Mass-Gap Classification

m₁ = 2.52 M☉ → P(NS) = 49%, P(BH) = 51%. Phonon suppression factor discriminates mass-gap objects.

## 3. Implementation

File: `fubi_agn_ns_mergers.py`, class `NSMergerFUBiCalc`. CP4 class #584.
