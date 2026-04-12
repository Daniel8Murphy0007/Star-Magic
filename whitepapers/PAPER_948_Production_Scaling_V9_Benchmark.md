# PAPER_948: Production Scaling V9 Benchmark (400k calc/s)

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** production_scaling_v9.py (ProductionScalingV9)
**Calculator:** ProductionScalingV9BenchmarkCalc (CP4 #532)
**CVW:** v2.0.0 compliant

---

## Abstract

We present the v9 production scaling benchmark targeting 400,000 calculations per second, a 14% increase over the v8 target of 350,000. The benchmark suite comprises 12 kernels: the 10 v8 kernels plus two new additions -- Centaurus A jet power evaluation and SMBH merger strain damping computation. Performance is measured in calculations per second with automated pass/fail against the 400k target.

---

## 1. Benchmark Kernels

| Kernel | Operation | Source |
|--------|-----------|--------|
| 1--10 | v8 base kernels | production_scaling_v8.py |
| 11 | CenA Jet P_BZ | blazar_jet_power_curves_extended.py |
| 12 | SMBH Merger D_total | smbh_binary_mergers.py |

### Target

$$\bar{R} = \frac{1}{K} \sum_{k=1}^{K} R_k \geq 400{,}000 \text{ calc/s}$$

where $K = 12$ is the total kernel count.

---

## 2. Scaling History

| Version | Target (calc/s) | Kernels | Date |
|---------|-----------------|---------|------|
| v4 | 100,000 | 4 | Dec 2025 |
| v5 | 150,000 | 5 | Jan 2026 |
| v6 | 200,000 | 6 | Feb 2026 |
| v7 | 300,000 | 8 | Mar 2026 |
| v8 | 350,000 | 10 | Apr 2026 |
| **v9** | **400,000** | **12** | **Apr 2026** |

---

## 3. New Kernels

**Kernel 11 -- CenA Jet:** Evaluates $P_\text{BZ}$ for Centaurus A at three linewidths. Single-pass complexity O(1).

**Kernel 12 -- SMBH Merger Strain:** Computes $D_\text{total}(q) = 0.333 + 0.197(1 - q)$ and applies to fiducial strain. Single-pass complexity O(1).

---

## 4. Source Data

- **File:** production_scaling_v9.py
- **Session:** 213
- **CP4 Class:** ProductionScalingV9BenchmarkCalc (#532)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Production benchmark history: v4 (Session 191) through v9 (Session 213)
