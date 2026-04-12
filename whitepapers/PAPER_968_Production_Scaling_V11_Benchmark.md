# PAPER_968: Production Scaling v11 Benchmark (500k calc/s)

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 215
**Source:** production_scaling_v11.py (ProductionScalingV11)
**Calculator:** ProductionScalingV11BenchmarkCalc (CP4 #552)
**CVW:** v2.0.0 compliant

---

## Abstract

We benchmark the UQFF production pipeline at 500,000 calculations per second, an 11.1% improvement over v10 (450k). Two new kernels — 26D Ramanujan summation and Triadic solver — bring the total to 16 simultaneously benchmarked kernels.

---

## 1. Scaling History

| Version | Target (calc/s) | Kernels | Session |
|---------|-----------------|---------|---------|
| v4 | 100,000 | 4 | 198 |
| v5 | 150,000 | 6 | 200 |
| v6 | 200,000 | 8 | 204 |
| v7 | 300,000 | 10 | 208 |
| v8 | 350,000 | 11 | 210 |
| v9 | 400,000 | 12 | 213 |
| v10 | 450,000 | 14 | 214 |
| **v11** | **500,000** | **16** | **215** |

## 2. New v11 Kernels

### kernel_ramanujan_26d_sum
$S_{26}(z) = \Sigma z^n / n^{26} \cdot R_n^{(26)}$ — accelerated polylog for hot-loop evaluation.

### kernel_triadic_solver
Combined Compressed + Resonant + Buoyancy Triadic evaluation in a single kernel call.

## 3. Throughput

$$\text{rate} = \frac{N_\text{iter} \times 16}{t_\text{elapsed}} \geq 500{,}000 \text{ calc/s}$$

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
