---
paper_id: PAPER_997
title: "Production Scaling v13 — 550k calc/s with 20 Benchmark Kernels"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [production, scaling, benchmark, throughput, kernels, v13]
crosslinks: [PAPER_968, PAPER_978, PAPER_989]
calibration: {target: 550000, kernels: 20, v11_carry: 16, v13_new: 4}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_997: Production Scaling v13 — 550k calc/s

## Abstract

We upgrade the production benchmark from v11 (500k calc/s, 16 kernels) to v13 (550k calc/s, 20 kernels). Four new kernels are added:

| Kernel | Description |
|--------|-------------|
| `kernel_fubi_inside_out` | F_U_Bi inside-to-outside mass portion |
| `kernel_99sys_gamma_sweep` | 99-system aggregate at given Γ |
| `kernel_agn_cena_fubi` | Centaurus A AGN F_U_Bi_i with jet modulation |
| `kernel_ns_merger_gw190425` | GW190425 NS merger strain with phonon suppression |

## 1. Scaling History

| Version | Target (calc/s) | Kernels |
|---------|-----------------|---------|
| v4 | 100,000 | 4 |
| v7 | 300,000 | 8 |
| v11 | 500,000 | 16 |
| **v13** | **550,000** | **20** |

## 2. Benchmark Architecture

```python
class ProductionScalingV13:
    TARGET = 550_000
    def single_pass(self) -> float:    # Execute all 20 kernels
    def simulate(self, duration_s)     # Timed benchmark loop
    def compute(self, dataset)         # Full pipeline with metadata
```

## 3. Implementation

File: `production_scaling_v13.py`, class `ProductionScalingV13`. CP4 class #581.
