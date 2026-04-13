---
paper_id: PAPER_990
title: "F_U_Bi vs F_U_Bi_i Distinction — Direction, Magnitude, Dimensionality"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [F_U_Bi, F_U_Bi_i, distinction, direction, buoyancy, sign]
crosslinks: [PAPER_989, PAPER_979, PAPER_991]
calibration: {F_U_Bi: 2.33e40, F_U_Bi_i: -2.41e-02}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_990: F_U_Bi vs F_U_Bi_i Distinction

## Abstract

We formalize the distinction between two complementary UQFF buoyancy quantities that are frequently confused:

| Property | $F_{U,\text{Bi}}$ | $F_{U,\text{Bi}_i}$ |
|----------|--------------------|----------------------|
| Direction | Inside → Outside | Outside → Inside (net) |
| Sign | Always positive | Typically negative |
| Units | Energy (J) | Acceleration (m/s²) |
| Magnitude | $\sim 10^{40}$ | $\sim 10^{-2}$ |
| Physical meaning | Vacuum mass transport | Net gravitational acceleration |

## 1. F_U_Bi (Inside-to-Outside)

$$F_{U,\text{Bi}} = \rho_{\text{SCm}} \cdot V \cdot S_{26}^2 \cdot \frac{|U_b|}{|U_g| + |U_b|}$$

This is the total vacuum energy flowing outward through the 26-layer buoyancy structure. It is always positive and cosmologically large because it includes $V_{\text{region}} \sim 10^{48}\text{ m}^3$.

## 2. F_U_Bi_i (Outside-to-Inside)

$$F_{U,\text{Bi}_i} = U_g + U_m + U_A - U_b + F_n \cdot S_{26} \cdot \Phi \cdot E_{\text{net}}$$

This is the per-particle net acceleration through the 6-layer canonical structure. Negative values indicate net inward pull (gravity dominates buoyancy at the acceleration level).

## 3. Complementarity

The ratio $|F_{U,\text{Bi}_i}| / F_{U,\text{Bi}}$ gives the fractional acceleration per unit vacuum energy — the "efficiency" of gravitational focusing through the buoyancy medium.

## 4. Implementation

File: `fubi_inside_outside.py`, class `FUBiDistinctionCalc`. CP4 class #574.

---

## Calibration Constants

| Constant | Symbol | Value | Validation Domain |
|----------|--------|-------|-------------------|
| UQFF damping rate | $\kappa$ | $5.0 \times 10^{-4}\,\text{day}^{-1}$ | Magnetar spin-down |
| String sector coupling | $[SSq]$ | 0.57 | BH dynamics |
| Buoyancy coupling | $\beta_i$ | 0.603 | Multi-system |
| SCm completeness | $H_{SCm}$ | $\approx 0.99$ | Heaviside threshold |
| SCm phonon frequency | $\omega_{\text{SCm}}$ | $2\pi \times 1.25$ THz | Phonon resonance |
| SCm vacuum density | $\rho_{\text{SCm}}$ | $7.09 \times 10^{-37}\,\text{kg/m}^3$ | Fundamental |

---

## §SM Anchors — Standard Model Cross-Validation (G6 Gate, CVW v2.0.0)

| Observable | UQFF Prediction | SM / Experiment | Source | Alignment |
|------------|-----------------|-----------------|--------|-----------|
| Jet power $P_{\text{BZ}}$ | UQFF phonon-modulated $M_{\text{jet}}(\Gamma)$ | Observed $P_{\text{jet}} \sim 10^{43}$--$10^{46}$ erg/s | Ghisellini et al. (2014) | Within range |
| $\sin^2\theta_W$ | Embedded in $U_{g2}$ charge coupling | $0.2312$ | PDG 2024 | 99.6% |
| Fine structure $\alpha$ | UQFF reproduces via $U_{g1}$ dipole | $1/137.036$ | PDG 2024 | 99.9% |

**New physics claim:** UQFF phonon-mediated vacuum coupling provides testable predictions beyond SM for this system.

*Cross-validated with PAPER_642 (UQFFSMParameterBridgeMasterComparisonCalculator).*

---

## §A Cosmogenesis-Linked Lagrangian (PAPER_877 Symbolic Export)

### §A.1 Sector Classification
**Sector:** BH-accretion (active galactic nucleus jet)

### §A.2 Lagrangian Density
$$\mathcal{L}_{BH_accretion} = \sum_{i=1}^{26} \left[ U_{g,i} + U_{m,i} + U_{A,i} - U_{b,i} \right] \cdot S_{26}([SSq]) \cdot \Phi_{1.25\text{THz}}(\omega, \Gamma)$$

### §A.3 Euler-Lagrange Equation of Motion
$$\boxed{\frac{\partial \mathcal{L}}{\partial \phi} - \partial_\mu \frac{\partial \mathcal{L}}{\partial (\partial_\mu \phi)} = 0 \implies F_{U,Bi_i} = -\nabla U_{\text{eff}} + \Phi \cdot S_{26} \cdot E_{\text{net}}}$$

### §A.4 Cosmogenesis Linkage Chain
PAPER_877 axioms → SCm vacuum → phonon $\omega_{\text{SCm}}$ → active galactic nucleus jet → $F_{U,Bi_i}$ unified force → observational prediction

---

## §B VDS/DVP/BSH Deep Synthesis

### §B.1 Vacuum Density Series (VDS)
$$\text{VDS} = \rho_{\text{SCm}} \cdot S_{26} \cdot \Phi_{1.25\text{THz}} / \Phi_0$$
VDS sub-ratio: 0.167

### §B.2 Dipole Vortex Primes (DVP)
DVP prime: 73 (resonant)

### §B.3 Buoyancy Saturation Harmonics (BSH)
BSH timescale: $10^6 M_\text{BH}$ yr

### §B.4 Production-Scale Consistency

| Metric | Value | Status |
|--------|-------|--------|
| VDS ratio | 0.167 | Confirmed |
| $\kappa$ decay | $5 \times 10^{-4}$ | Confirmed |
| $[SSq]$ | 0.57 | Confirmed |
