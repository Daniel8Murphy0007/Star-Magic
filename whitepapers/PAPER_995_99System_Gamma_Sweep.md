---
paper_id: PAPER_995
title: "99-System Gamma Sweep Computation — Aggregate F_U_Bi_i at 7 Linewidths"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [99-system, Gamma, sweep, linewidth, aggregate, F_U_Bi_i, catalogue]
crosslinks: [PAPER_984, PAPER_996, PAPER_989]
calibration: {systems: 99, gamma_values: 7, aggregate_0p1: -6.11e13}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_995: 99-System Gamma Sweep Computation

## Abstract

We compute the aggregate $F_{U,\text{Bi}_i}$ across all 99 astrophysical systems in the UQFF catalogue at 7 linewidth values $\Gamma \in \{0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0\}\text{ THz}$.

## 1. System Catalogue

The 99 systems span 6 categories:
- **Stars** (20): Main sequence, giants, supergiants ($0.1$–$100\,M_\odot$)
- **Galaxies** (20): Spirals through ellipticals ($10^9$–$10^{13}\,M_\odot$)
- **Nebulae** (15): Planetary through H II regions
- **Compact Objects** (15): Neutron stars ($1.4$–$2.5\,M_\odot$) and stellar BHs ($3$–$100\,M_\odot$)
- **Clusters** (15): Globular through galaxy superclusters ($10^{13}$–$10^{16}\,M_\odot$)
- **Cosmological** (14): Large-scale structure ($10^{15}$–$10^{17}\,M_\odot$)

## 2. Aggregate at Reference Linewidth

$$F_U^{(99)}(\Gamma_0 = 0.1\text{ THz}) = -6.11 \times 10^{13}\text{ m/s}^2$$

The negative sign confirms global inward dominance across all mass scales.

## 3. Stability

100% of systems return finite, stable $F_{U,\text{Bi}_i}$ values across all 7 $\Gamma$ values — no divergences or instabilities detected.

## 4. Implementation

File: `99system_wstp_gamma.py`, class `NinetyNineSystemGammaSweep`. CP4 class #579.
