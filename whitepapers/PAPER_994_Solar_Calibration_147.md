---
paper_id: PAPER_994
title: "Solar Calibration g_eff Convergence via 26-Layer Buoyancy Correction"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [solar, calibration, g_eff, buoyancy, correction, Sun, surface_gravity]
crosslinks: [PAPER_989, PAPER_990, PAPER_980]
calibration: {g_N: 274, g_eff: 108.05, S26_ratio: 1.529}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_994: Solar Calibration g_eff Convergence

## Abstract

We derive the effective solar surface gravity including the 26-layer buoyancy correction:

$$g_{\text{eff}} = \frac{g_N}{1 + \frac{\beta_i \cdot S_{26}}{[\text{SSq}] \cdot 13.5}}$$

With $g_N = GM_\odot / R_\odot^2 = 274\text{ m/s}^2$, $\beta_i = 0.603$, $S_{26} \approx 19.5$, $[\text{SSq}] = 0.57$:

$$g_{\text{eff}} = \frac{274}{1 + \frac{0.603 \times 19.5}{0.57 \times 13.5}} = \frac{274}{1 + 1.529} = \frac{274}{2.529} \approx 108.05\text{ m/s}^2$$

## 1. Physical Interpretation

The buoyancy correction reduces the effective surface gravity by a factor of $\sim 2.53$. This represents the SCm vacuum medium partially supporting mass against gravitational collapse — the "buoyancy floor" of the UQFF framework.

## 2. Convergence

The 99-system $\Gamma$ sweep shows that this solar calibration value is stable across all 7 linewidth values tested ($\Gamma \in \{0.01, ..., 10.0\}\text{ THz}$), confirming that the solar calibration is Γ-independent (as expected for a static surface gravity measurement).

## 3. Implementation

File: `fubi_inside_outside.py`, class `SolarCalibration147Calc`. CP4 class #578.
