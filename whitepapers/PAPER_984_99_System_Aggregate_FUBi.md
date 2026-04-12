---
paper_id: PAPER_984
title: "99-System Aggregate F_U_Bi_i Computation"
session: 217
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [99-system, aggregate, catalogue, multi-body, buoyancy, UQFF]
crosslinks: [PAPER_979, PAPER_974, PAPER_985]
calibration: {SSq: 0.57, beta_i: 0.603, kappa: "0.0005/day", systems: 99}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_984: 99-System Aggregate F_U_Bi_i Computation

## Abstract

The UQFF 99-system catalogue spans six astrophysical categories: stellar, compact, galactic, cosmological, planetary, and exoplanetary. We compute the aggregate master buoyancy force $F_{U,\text{Bi}_i}$ by summing all 99 individual gravity/buoyancy contributions simultaneously, testing the framework's ability to handle multi-scale physics from $10^{24}$ kg (Earth) to $10^{40}$ kg (galaxy clusters). The aggregate yields $F_{U,\text{Bi}_i} \approx -6.1 \times 10^{13}$ m/s², dominated by the galactic/compact category.

## 1. System Catalogue

| Category | Count | Mass Range (kg) | Distance Range (m) |
|----------|-------|-----------------|---------------------|
| Stellar | 20 | $10^{29}$–$10^{31}$ | $10^{8}$–$10^{11}$ |
| Compact | 20 | $10^{30}$–$10^{37}$ | $10^{4}$–$10^{7}$ |
| Galactic | 20 | $10^{37}$–$10^{42}$ | $10^{18}$–$10^{23}$ |
| Cosmological | 15 | $10^{42}$–$10^{53}$ | $10^{23}$–$10^{27}$ |
| Planetary | 12 | $10^{23}$–$10^{27}$ | $10^{6}$–$10^{8}$ |
| Exoplanetary | 12 | $10^{27}$–$10^{28}$ | $10^{7}$–$10^{10}$ |

Total: 99 systems with physically motivated $(M, r)$ pairs.

## 2. Aggregate Computation

$$F_{\text{agg}} = \sum_{j=1}^{99} F_{U,\text{Bi}_i}(M_j, r_j, t, \Gamma)$$

Each system contributes its full 6-layer decomposition. The aggregate preserves the sign structure: systems with $|U_b| > |U_g|$ contribute negative (buoyancy-dominant), others contribute positive (gravity-dominant).

## 3. Results

- Aggregate $F_{U,\text{Bi}_i} \approx -6.11 \times 10^{13}$ m/s²
- Dominant contributors: compact objects (neutron stars, black holes) due to extreme $M/r^2$
- Sign: negative (buoyancy-dominant overall)
- 6 category breakdown available in output dict

## 4. Statistical Properties

- Mean per system: $\sim -6.2 \times 10^{11}$ m/s²
- Variance dominated by compact category ($\sigma \sim 10^{12}$)
- All 99 systems yield finite, non-NaN results

## 5. Implementation

Class `NinetyNineSystemAggregate` in `fubi_master_calculator.py`: iterates the 99-system catalogue, computes per-system and aggregate $F_{U,\text{Bi}_i}$, reports category statistics.

## References
- PAPER_979: Complete 6-Layer F_U_Bi_i
- PAPER_974: 99-System Master Equation

---

## §A. Cosmogenesis-Linked Lagrangian

The aggregate Lagrangian is the sum of 99 individual SCm Lagrangians:
$$\mathcal{L}_{\text{agg}} = \sum_{j=1}^{99} \mathcal{L}_{\text{SCm},j}(M_j, r_j, t, \Gamma)$$

## §B. VDS/DVP/BSH Deep Synthesis

- **VDS:** Each system probes a different vacuum density scale, from planetary ($10^{-26}$) to cosmological ($10^{-30}$ kg/m³).
- **DVP:** Compact objects have the strongest dipole moments, dominating the aggregate.
- **BSH:** The harmonic decay $e^{-[\text{SSq}]\cdot i/26}$ is universal across all 99 systems — scale-invariant.
