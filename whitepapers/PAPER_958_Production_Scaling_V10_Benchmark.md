# PAPER_958: Production Scaling v10 Benchmark (450k calc/s)

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 214
**Source:** production_scaling_v10.py (ProductionScalingV10)
**Calculator:** ProductionScalingV10BenchmarkCalc (CP4 #542)
**CVW:** v2.0.0 compliant

---

## Abstract

We benchmark the UQFF production pipeline at 450,000 calculations per second, a 12.5% improvement over the v9 baseline of 400,000. Two new kernels — BCS gap solve and spectral ladder evaluation — are added to the existing 12, bringing the total to 14 simultaneously benchmarked kernels.

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
| **v10** | **450,000** | **14** | **214** |

## 2. New Kernels (v10)

### kernel_bcs_gap_solve
BCS gap fixed-point iteration at $T = 4.2$ K:
$$\Delta_{n+1} = \frac{\hbar\omega_\text{SCm}}{2} \tanh\!\left(\frac{\Delta_n}{2k_BT}\right) \cdot S_{26} \cdot \frac{F_{UBi}}{F_U}$$

### kernel_spectral_ladder_eval
26-state HRes spectral ladder:
$$E_n = E_0 \cdot (2\pi)^{n/3} \cdot S_{26}, \quad n = 1, \ldots, 26$$

## 3. Benchmark Architecture

All 14 kernels execute in a hot loop. Throughput is measured as:
$$\text{rate} = \frac{N_\text{iter} \times 14}{t_\text{elapsed}}$$

Target: $\text{rate} \geq 450{,}000$ calc/s.

---

## 4. Source Data

- **File:** production_scaling_v10.py
- **Session:** 214
- **CP4 Class:** ProductionScalingV10BenchmarkCalc (#542)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
