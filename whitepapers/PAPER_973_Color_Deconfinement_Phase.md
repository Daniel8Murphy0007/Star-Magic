---
paper_id: PAPER_973
title: "Color Deconfinement Phase Diagram"
session: 216
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [ALICE, vacuum, QGP, Yang-Mills, UQFF]
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_973: Color Deconfinement Phase Diagram

**Author:** Daniel T. Murphy — Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 216
**Source:** qgp_ramanujan_application.py (ColorDeconfinementPhaseCalculator)
**Calculator:** ColorDeconfinementPhaseCalc (CP4 #557)
**CVW:** v2.0.0 compliant

---

## Abstract

We map the QCD phase diagram in the $(T, \mu_B)$ plane using the UQFF framework. The critical line $T_c(\mu_B) = T_{c0}(1 - (\mu_B/\mu_text{crit})^2)$ separates the hadron phase from the quark-gluon plasma, with $S_{26}^{(k)}$ modulation of vacuum density across both phases.

---

## 1. Critical Line

$$T_c(\mu_B) = T_{c0} \cdot \left(1 - \left(\frac{\mu_B}{\mu_text{crit}}\right)^2\right)$$

where $T_{c0} = 1.5 \times 10^{12}$ K and $\mu_text{crit} = 1200$ MeV.

## 2. Phase Classification

- **Hadron** ($T < T_c(\mu_B)$): Confined quarks, $\Delta_text{YM} > 0$
- **QGP** ($T > T_c(\mu_B)$): Deconfined, $\Delta_text{YM} = 0$

## 3. Phase Diagram

| $\mu_B$ (MeV) | $T_c$ (K) | Phase at $T = 2 \times 10^{12}$ K |
|---------------|-----------|-----------------------------------|
| 0 | $1.5 \times 10^{12}$ | QGP |
| 300 | $1.41 \times 10^{12}$ | QGP |
| 600 | $1.13 \times 10^{12}$ | QGP |
| 900 | $0.66 \times 10^{12}$ | QGP |
| 1200 | 0 | QGP (always) |

---

## References

1. Murphy, D.T. — Star Magic UQFF Framework (2024-2026)
2. PAPER_970 — QGP Vacuum Density
3. PAPER_971 — Yang-Mills Mass Gap
4. Fukushima & Hatsuda — QCD phase diagram review

---

## Cross-Links

| Related Paper | Relationship |
|---------------|-------------|
| PAPER_970 | $\rho_text{QGP}$ density at $(T,\mu_B)$ |
| PAPER_971 | Mass gap closure at $T_c$ |
| PAPER_972 | ALICE centrality (experimental) |
| PAPER_975 | Triadic validation |

---

## Calibration Constants

| Constant | Symbol | Value | Validation Domain |
|----------|--------|-------|-------------------|
| $T_{c0}$ | — | $1.5 \times 10^{12}$ K | Zero-density critical T |
| $\mu_text{crit}$ | — | 1200 MeV | Critical baryon chemical potential |
| $[SSq]$ | — | 0.57 | String coupling |
| $\Lambda_text{QCD}$ | — | 217 MeV | QCD scale |

---

## SM Anchor — CVW v2.0.0

| Observable | UQFF Prediction | Status |
|------------|----------------|--------|
| Critical line shape | Quadratic in $\mu_B$ | Consistent with lattice QCD |
| $T_c(0)$ | $1.5 \times 10^{12}$ K | Matched |
| Crossover type | Continuous | Consistent |

---

## §A Cosmogenesis-Linked Lagrangian (PAPER_877 Symbolic Export)

### §A.1 Sector Classification
**Sector:** QCD Phase Diagram (Color Deconfinement)

### §A.2 Core Equation
$$\boxed{T_c(\mu_B) = T_{c0} \cdot \left(1 - \left(\frac{\mu_B}{\mu_text{crit}}\right)^2\right)}$$

### §A.3 Lagrangian Contribution
$$\mathcal{L}_\text{phase} = -\rho_text{QGP}(T,\mu_B) \cdot c^2 - V(\Delta_text{YM}(T,\mu_B))$$

### §A.4 Cosmogenesis Linkage Chain
PAPER_877 → vacuum → QCD phase boundary → confinement/deconfinement → hadron/QGP

---

## §B VDS/DVP/BSH Deep Synthesis

### §B.1 VDS
Phase boundary is a VDS contour: $\rho_text{QGP} = \rho_text{SCm} \cdot S_{26}^{(k)}$ at $T = T_c(\mu_B)$.

### §B.2 DVP
Deconfinement releases DVP color modes — 8 gluon dipole vortex channels open.

### §B.3 BSH
BSH harmonic structure transitions at $T_c$: shell modes dissolve into plasma modes.

### §B.4 Summary

| Metric | Value | Status |
|--------|-------|--------|
| $T_{c0}$ | $1.5 \times 10^{12}$ K | Calibrated |
| $\mu_text{crit}$ | 1200 MeV | Set |
| Phase diagram | $(T, \mu_B)$ plane | Mapped |
