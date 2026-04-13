---
paper_id: PAPER_1018
title: "Production Scaling v15 — 650k calc/s with 30 Benchmark Kernels"
session: 220
date: 2026-04-14
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [production, scaling, benchmark, throughput, kernels, v15, 650k]
crosslinks: [PAPER_1008, PAPER_1009, PAPER_1014, PAPER_1015]
calibration: {target: 650000, kernels: 30, v14_carry: 24, v15_new: 6}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1018: Production Scaling v15 — 650k calc/s

## Abstract

We upgrade the production benchmark from v14 (600k calc/s, 24 kernels) to v15 (650k calc/s, 30 kernels). Six new kernels cover the Session 220 physics domains:

| Kernel | Description |
|--------|-------------|
| `kernel_3c273_agn_fubi` | 3C273 AGN F_U_Bi_i jet modulation |
| `kernel_ton618_agn_fubi` | TON618 ultramassive AGN curves |
| `kernel_gw170817_merger` | GW170817 BNS strain suppression |
| `kernel_smbh_merger_fubi` | SMBH merger inspiral F_U_Bi |
| `kernel_dm_halo_nfw` | DM halo NFW with SCm coupling |
| `kernel_txs0506_3gamma` | TXS 0506+056 3-Gamma-point profile |

## 1. Scaling History

| Version | Target (calc/s) | Kernels |
|---------|-----------------|---------|
| v11 | 500,000 | 16 |
| v13 | 550,000 | 20 |
| v14 | 600,000 | 24 |
| **v15** | **650,000** | **30** |

## 2. Validation Checks

- TON618 (2.19e-1) > 3C273 (1.78e-1): Mass hierarchy confirmed
- GW170817 damped (2.83e-23) < undamped (8.51e-23): Buoyancy suppression confirmed
- All 30 kernels finite and non-negative
- Total throughput exceeds 650k target

## 3. REST API

22 total routes including new endpoints:
- `POST /api/fubi/smbh-merger` — SMBH merger F_U_Bi with S26(3)
- `POST /api/dm/halo-nfw` — DM halo NFW profile with SCm coupling

## 4. Implementation

File: `production_scaling_v15.py`, class `ProductionScalingV15`. CP4 class #602. Tests: 8/8 pass.
