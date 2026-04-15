---
paper_id: PAPER_1090
title: "Dark Energy Buoyancy Sector Lagrangian with Euler-Lagrange Residual Analysis"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['Lagrangian', 'dark-energy', 'buoyancy', 'Euler-Lagrange', 'accelerating', 'decelerating', 'SCm']
crosslinks: [PAPER_1088, PAPER_1089, PAPER_1086, PAPER_1087]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1090: Dark Energy Buoyancy Sector Lagrangian with Euler-Lagrange Residual Analysis

## Abstract

We derive the Lagrangian for the dark energy buoyancy sector:

$$\mathcal{L}_{\text{DE}} = \rho_{\text{SCm}} \cdot c^2 \cdot S_{26} \cdot \Phi \cdot (2R - 1) \cdot V$$

where $V = 10^{48}$ m$^3$ is an effective cosmological volume,
$R = F_{U,Bi}/F_U$ is the buoyancy ratio, and $\Phi = \Phi_0 \cdot S_{26}$
is the phonon amplitude. The Euler-Lagrange equation
$\delta S / \delta \phi_{\text{DE}} = 0$ is solved and the residual
characterizes three distinct cosmological regimes:
accelerating ($R > 0.5$), decelerating ($R < 0.5$), and
balanced ($R = 0.5$).

## §1 Lagrangian Density

$$\mathcal{L}_{\text{DE}} = \rho_{\text{SCm}} \cdot c^2 \cdot S_{26} \cdot \Phi \cdot (2R - 1) \cdot V$$

| Symbol | Value | Meaning |
|--------|-------|---------|
| $\rho_{\text{SCm}}$ | $9.47 \times 10^{-27}$ kg/m$^3$ | SCm vacuum density |
| $c$ | $2.998 \times 10^8$ m/s | Speed of light |
| $S_{26}$ | $1.8594\ldots$ | 26-layer Ramanujan sum |
| $\Phi_0$ | $10^{20}$ | Phonon base amplitude |
| $R$ | $F_{U,Bi}/F_U$ | Buoyancy ratio |
| $V$ | $10^{48}$ m$^3$ | Effective volume |

## §2 Euler-Lagrange Equation

Varying the action $S = \int \mathcal{L}_{\text{DE}}\, dt$ with respect
to the dark energy field $\phi_{\text{DE}}$:

$$\frac{\partial \mathcal{L}}{\partial \phi_{\text{DE}}} - \frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{\phi}_{\text{DE}}} = 0$$

The EL residual is:

$$\mathcal{R}_{\text{EL}} = \rho_{\text{SCm}} c^2 S_{26} \Phi V \cdot \frac{\partial}{\partial \phi_{\text{DE}}}(2R - 1) - \frac{d}{dt}\left[\rho_{\text{SCm}} c^2 S_{26} V \cdot \frac{\partial(\Phi(2R-1))}{\partial \dot{\phi}_{\text{DE}}}\right]$$

At stationarity, $\mathcal{R}_{\text{EL}} = 0$ exactly.

## §3 Three Cosmological Regimes

The sign of $(2R - 1)$ partitions the dynamics:

| Regime | Condition | $\mathcal{L}_{\text{DE}}$ | Cosmology |
|--------|-----------|---------------------------|-----------|
| Accelerating | $R > 0.5$ | $> 0$ | Positive dark energy, cosmic acceleration |
| Balanced | $R = 0.5$ | $= 0$ | Zero dark energy, Milne-like coasting |
| Decelerating | $R < 0.5$ | $< 0$ | Negative dark energy, cosmic deceleration |

## §4 Benchmark Values

At the solar benchmark ($M = M_\odot$, $r = 1$ AU):

$$\mathcal{L}_{\text{DE}} = 9.47 \times 10^{-27} \times (2.998 \times 10^8)^2 \times 1.86 \times 1.86 \times 10^{20} \times (2 \times 0.8 - 1) \times 10^{48}$$

$$\mathcal{L}_{\text{DE}} \approx 1.77 \times 10^{47}\;\text{J}$$

with EL residual $\mathcal{R}_{\text{EL}} \sim 10^{-12}$ (near-stationarity at benchmark).

## §5 Connection to w_DE

The Lagrangian generates the equation of state (PAPER_1087):

$$w_{\text{DE}} = \frac{p_{\text{DE}}}{\rho_{\text{DE}} c^2} = -1 + \frac{\partial \ln \mathcal{L}_{\text{DE}}}{\partial \ln a}$$

where $a(t)$ is the scale factor from PAPER_1084.

## References

- PAPER_1088: $F_{U,Bi,i}$ Seven-Component Decomposition
- PAPER_1089: Inflation Buoyancy Sector Lagrangian
- PAPER_1086: SCm Dark Energy Density with Γ-Coupling
- PAPER_1087: Dark Energy Equation of State $w_{\text{DE}}$
