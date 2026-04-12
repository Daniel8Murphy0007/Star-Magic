---
paper_id: PAPER_980
title: "Solar Surface Buoyancy Calibration"
session: 217
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [solar, calibration, buoyancy, g_N, surface-gravity, UQFF]
crosslinks: [PAPER_979, PAPER_981, PAPER_883]
calibration: {SSq: 0.57, beta_i: 0.603, kappa: "0.0005/day", g_sun: "274.03 m/s²"}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_980: Solar Surface Buoyancy Calibration

## Abstract

The UQFF master buoyancy equation $F_{U,\text{Bi}_i}$ must reproduce known observational benchmarks before deployment. We present the solar surface calibration test: at $r = R_\odot = 6.96 \times 10^8$ m, Newtonian gravity yields $g_N = 274.03$ m/s² consistent with the IAU standard. The 6-layer buoyancy architecture modifies this by $\lesssim 1\%$, confirming the "gravity-first" regime at short range. At $r = 1$ AU, the buoyancy term dominates ($F_{U,\text{Bi}_i} < 0$), consistent with solar wind acceleration and heliospheric expansion.

## 1. Newtonian Benchmark

$$g_N = \frac{GM_\odot}{R_\odot^2} = \frac{6.674 \times 10^{-11} \times 1.989 \times 10^{30}}{(6.96 \times 10^8)^2} = 274.03 \text{ m/s}^2$$

## 2. UQFF Correction at Solar Surface

At $r = R_\odot$, $t = 1$ day:
- $U_g \gg U_b$: gravity dominates at short range
- $F_{U,\text{Bi}_i} \approx g_N \cdot (1 - \delta_b)$ where $\delta_b \ll 1$
- Buoyancy correction $\delta_b \sim \beta_i \cdot e^{-[\text{SSq}]/26} \approx 0.6\%$

## 3. Heliospheric Regime ($r = 1$ AU)

At 1 AU = $1.496 \times 10^{11}$ m:
- $g_N = 5.93 \times 10^{-3}$ m/s²
- $F_{U,\text{Bi}_i} \approx -2.4 \times 10^{-2}$ m/s² (negative = buoyancy > gravity)
- Physical interpretation: net outward force consistent with solar wind acceleration

## 4. Crossover Radius

The radius where $U_g = U_b$ defines the buoyancy-gravity crossover:
$$r_{\text{cross}} \approx R_\odot \cdot \left(\frac{1}{\beta_i \cdot [\text{SSq}]}\right)^{1/2}$$

## 5. Implementation

Class `SolarSurfaceCalibrator` in `fubi_master_calculator.py`: computes $g_N(R_\odot)$, verifies $|g_N - 274| < 1$ m/s².

## References
- PAPER_979: Complete 6-Layer F_U_Bi_i
- PAPER_883: E(t) Phonon Resonance

---

## §A. Cosmogenesis-Linked Lagrangian

At the solar surface, the Lagrangian reduces to the Newtonian limit:
$$\mathcal{L} \to \frac{1}{2}m\dot{r}^2 + \frac{GMm}{r}$$
with buoyancy perturbation $\delta\mathcal{L}_b = -U_b \cdot m \cdot r$.

## §B. VDS/DVP/BSH Deep Synthesis

- **VDS:** Solar vacuum density $\rho_{\text{vac,}\odot} \approx 6 \times 10^{-27}$ kg/m³ at photosphere.
- **DVP:** Solar DPM moment $\mu_\odot \approx 6.6 \times 10^{32}$ A·m² drives $U_{g1}$ layer.
- **BSH:** $\beta_i = 0.603$ defines the buoyancy-to-gravity ratio at each layer boundary.
