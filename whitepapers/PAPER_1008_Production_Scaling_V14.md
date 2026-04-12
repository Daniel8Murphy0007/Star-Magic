---
paper_id: PAPER_1008
title: "Production Scaling v14 — 600k calc/s with 24 Benchmark Kernels"
session: 219
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [production, scaling, benchmark, throughput, kernels, v14, 600k]
crosslinks: [PAPER_997, PAPER_999, PAPER_1004]
calibration: {target: 600000, kernels: 24, v13_carry: 20, v14_new: 4}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1008: Production Scaling v14 — 600k calc/s

## Abstract

We upgrade the production benchmark from v13 (550k calc/s, 20 kernels) to v14 (600k calc/s, 24 kernels). Four new kernels using S₂₆⁽³⁾:

| Kernel | Description |
|--------|-------------|
| `kernel_agn_merger_fubi` | AGN merger F_U_Bi with S₂₆⁽³⁾ |
| `kernel_qgp_vacuum_density` | QGP vacuum density ρ_QGP(T) |
| `kernel_alice_multiplicity` | ALICE dN_ch/dη SCm scaling |
| `kernel_ym_mass_gap` | Yang-Mills Δ_YM via BCS phonon |

## 1. Scaling History

| Version | Target (calc/s) | Kernels |
|---------|-----------------|---------|
| v11 | 500,000 | 16 |
| v13 | 550,000 | 20 |
| **v14** | **600,000** | **24** |

## 2. REST API

20 total routes including new endpoints:
- `POST /api/fubi/agn-merger` — AGN merger F_U_Bi with S₂₆⁽³⁾
- `POST /api/qgp/scm-dynamics` — QGP SCm phonon dynamics

## 3. Implementation

File: `production_scaling_v14.py`, class `ProductionScalingV14`. CP4 class #592.
