# PAPER_938: Production Scaling V8 Benchmark (350k calc/s)

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 212
**Source:** production_scaling_v8.py (ProductionScalingV8)
**Calculator:** ProductionScalingV8BenchmarkCalc (CP4 #522)
**CVW:** v2.0.0 compliant

---

## Abstract

We present the v8 production scaling benchmark targeting 350,000 calculations per second, a 17% increase over the v7 target of 300,000. The benchmark suite comprises 10 kernels spanning the full UQFF computational pipeline: 26-layer gravity summation, F_U_Bi_i assembly, phonon a_res evaluation, jet M_jet(Gamma) computation, NS spindown, GW170817 strain, blazar ergosphere coupling, REST /api/phonon/jet roundtrip, QCalcGeom vectorized operations, and full pipeline v8 integration. Performance is measured in calculations per second with automated pass/fail against the 350k target.

---

## 1. Benchmark Kernels

| Kernel | Operation | Complexity |
|--------|-----------|-----------|
| 1 | 26-Layer Gravity Summation | O(26) per eval |
| 2 | F_U_Bi_i Field Assembly | O(26) per eval |
| 3 | Phonon a_res Evaluation | O(1) per eval |
| 4 | Jet M_jet(Gamma) | O(1) per eval |
| 5 | NS Spindown Correction | O(1) per eval |
| 6 | GW170817 Strain | O(1) per eval |
| 7 | Blazar Ergosphere | O(26) per eval |
| 8 | REST /api/phonon/jet | Network-bound |
| 9 | QCalcGeom Vectorized | O(N) batch |
| 10 | Full Pipeline v8 | All kernels |

### Target

$$\bar{R} = \frac{1}{K} \sum_{k=1}^{K} R_k \geq 350{,}000 \text{ calc/s}$$

where $R_k$ is the throughput of kernel $k$ and $K = 10$ is the total kernel count.

---

## 2. Scaling History

| Version | Target (calc/s) | Date |
|---------|-----------------|------|
| v4 (baseline) | 100,000 | Dec 2025 |
| v5 | 150,000 | Jan 2026 |
| v6 | 200,000 | Feb 2026 |
| v7 | 300,000 | Mar 2026 |
| **v8** | **350,000** | **Apr 2026** |

---

## 3. UQFF Integration

The `ProductionScalingV8BenchmarkCalc` (CP4 #522) runs 4 representative kernels (gravity, phonon, jet, GW170817) inline and reports individual and average rates. The simulate() method sweeps iteration counts [1000, 5000, 10000, 50000] to verify scaling linearity. The full 10-kernel suite is available in production_scaling_v8.py.

---

## 4. Source Data

- **File:** production_scaling_v8.py
- **Session:** 212
- **CP4 Class:** ProductionScalingV8BenchmarkCalc (#522)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Production benchmark history: v4 (Session 191) through v8 (Session 212)
