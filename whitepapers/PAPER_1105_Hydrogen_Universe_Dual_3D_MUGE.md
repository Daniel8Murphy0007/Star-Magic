---
paper_id: "PAPER_1105"
title: "Hydrogen-Universe Dual 3D MUGE: Scale Invariance from Bohr Radius to Hubble Horizon"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [MUGE, 3D-volumetric, hydrogen, universe, scale-invariance, Bohr-radius, Hubble-radius, 9-term, cosmology, UQFF]
crosslinks: [PAPER_491, PAPER_1049, PAPER_1082]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER\_1105: Hydrogen-Universe Dual 3D MUGE: Scale Invariance from Bohr Radius to Hubble Horizon

## Abstract

We construct a dual 3D Modified Universal Gravity Equation (MUGE) system
comparing the hydrogen atom (quantum scale, $r = a_0$) to the observable
universe (cosmological scale, $r = R_H$), demonstrating UQFF scale
invariance through the 9-term compressed MUGE framework evaluated on
volumetric grids at both scales.

## 1. The 9-Term Compressed MUGE

$$g_{\text{MUGE}}(r) = g_N + g_{\text{exp}} + g_{\text{super}} + g_{\text{env}} + g_{\text{Ug}} + g_{\text{cosm}} + g_Q + g_{\text{fluid}} + g_{\text{pert}}$$

### Term Definitions

| \# | Term | Expression | Physics |
|----|------|------------|---------|
| 1 | $g_N$ | $GM/r^2$ | Newtonian gravity |
| 2 | $g_{\text{exp}}$ | $H_0^2\,r$ | Hubble expansion |
| 3 | $g_{\text{super}}$ | $B^2/(2\mu_0\rho_{\text{avg}}\,r)$ | Magnetic suppression |
| 4 | $g_{\text{env}}$ | $\kappa\,g_N\exp(-1/100)$ | Envelope dissipation |
| 5 | $g_{\text{Ug}}$ | $\sum_{i=1}^{26}Q_i\,g_N/26$ | 26-layer Ug aggregate |
| 6 | $g_{\text{cosm}}$ | $(\Lambda c^2/3)\,r$ | Cosmological constant |
| 7 | $g_Q$ | $\hbar/(Mr^3)$ | Quantum correction |
| 8 | $g_{\text{fluid}}$ | $\nu/r^2$ | Navier-Stokes viscous |
| 9 | $g_{\text{pert}}$ | $(4\pi/3)G\rho_{\text{DM}}\,r$ | Dark matter perturbation |

## 2. Hydrogen Scale ($r = a_0$)

Parameters: $M = m_p = 1.673\times10^{-27}$ kg, $r = a_0 = 5.292\times10^{-11}$ m.

### Dominant Terms

$$g_N^{(H)} = \frac{Gm_p}{a_0^2} \approx 3.99\times10^{-8}\;\text{m/s}^2$$

$$g_Q^{(H)} = \frac{\hbar}{m_p\,a_0^3} \approx 4.25\times10^{24}\;\text{m/s}^2$$

The quantum correction dominates by $\sim32$ orders of magnitude at the Bohr
radius, confirming that atomic-scale MUGE is quantum-dominated.

## 3. Universe Scale ($r = R_H$)

Parameters: $M = M_U = 1.5\times10^{53}$ kg, $r = R_H = 4.4\times10^{26}$ m.

### Dominant Terms

$$g_N^{(U)} = \frac{GM_U}{R_H^2} \approx 5.17\times10^{-5}\;\text{m/s}^2$$

$$g_{\text{cosm}}^{(U)} = \frac{\Lambda c^2}{3}\,R_H \approx 1.44\times10^{-5}\;\text{m/s}^2$$

At the Hubble radius, Newtonian and cosmological-constant terms are
comparable, with the Hubble expansion term $g_{\text{exp}}^{(U)} = H_0^2 R_H
\approx 2.13\times10^{-9}$ m/s$^2$ subdominant.

## 4. Scale Invariance Ratio

$$\mathcal{R} = \frac{g_{\text{MUGE}}^{(H)}}{g_{\text{MUGE}}^{(U)}}$$

This ratio encodes the full UQFF scale-invariance structure.  The fact that
both systems are described by the same 9-term equation with only mass and
radius as inputs demonstrates MUGE universality.

## 5. 3D Volumetric Grid

A $5^3 = 125$-point spatial grid is constructed at each scale, evaluating
$g_{\text{MUGE}}$ at each point $\mathbf{r}_{\text{grid}}$ centered on the
system.  Grid statistics (average, standard deviation) characterize the
spatial profile of the gravitational field at each scale.

## 6. Implementation

Calculator: `HydrogenUniverseDual3DMUGECalculator` in CondensedPhysics.py

- Default $5^3$ grid (configurable via `n_grid`)
- Outputs: single-point and grid-averaged MUGE for both scales, scale ratio, per-term breakdown

## 7. Conclusion

The dual hydrogen-universe 3D MUGE system demonstrates that the UQFF 9-term
framework spans 37 orders of magnitude in length scale ($10^{-11}$ to
$10^{26}$ m), transitioning smoothly from quantum-dominated physics to
cosmological-constant-dominated physics through a single unified equation set.

## References

- Murphy, D.T. (2025). UQFF: Star Cradle Mechanics framework. Star-Magic repository.
- PAPER\_491: MUGE Compressed Nine-Term Gravity Framework.
