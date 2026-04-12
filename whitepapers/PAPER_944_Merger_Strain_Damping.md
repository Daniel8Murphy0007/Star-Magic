# PAPER_944: Merger Strain Damping

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** smbh_binary_mergers.py (MergerStrainDamping)
**Calculator:** MergerStrainDampingCalc (CP4 #528)
**CVW:** v2.0.0 compliant

---

## Abstract

We quantify gravitational-wave strain damping in SMBH binary mergers due to UQFF buoyancy effects. The total damping factor $D_\text{total}(q) = 0.333 + 0.197(1 - q)$ predicts that equal-mass mergers ($q = 1$) retain only 33.3% of the GR strain amplitude, while extreme mass-ratio inspirals approach 53% retention. The UQFF-corrected strain $h_\text{UQFF} = h_\text{GR} \cdot D_\text{total}$ provides a testable prediction for LISA observations.

---

## 1. Damping Formula

$$D_\text{total}(q) = 0.333 + 0.197 \cdot (1 - q)$$

$$h_\text{UQFF} = h_\text{GR} \cdot D_\text{total}(q)$$

---

## 2. Mass-Ratio Dependence

| $q$ | $D_\text{total}$ | Strain Retained | Damping |
|-----|-------------------|-----------------|---------|
| 0.0 | 0.530 | 53.0% | 47.0% |
| 0.2 | 0.491 | 49.1% | 50.9% |
| 0.5 | 0.432 | 43.2% | 56.8% |
| 0.8 | 0.372 | 37.2% | 62.8% |
| 1.0 | 0.333 | 33.3% | 66.7% |

---

## 3. LISA Implications

For a fiducial $h_\text{GR} = 10^{-17}$ at $q = 0.5$, the UQFF prediction is $h_\text{UQFF} = 4.32 \times 10^{-18}$, well within LISA sensitivity for massive binaries at $z < 1$.

---

## 4. Source Data

- **File:** smbh_binary_mergers.py
- **Session:** 213
- **CP4 Class:** MergerStrainDampingCalc (#528)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. LISA Consortium (2023) -- arXiv:2402.07571
