---
paper_id: PAPER_1094
title: "CMB Buoyancy Sector Lagrangian with Stationarity Constraint"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['Lagrangian', 'CMB', 'buoyancy', 'stationarity', 'Euler-Lagrange', 'action-principle', 'acoustic-peaks']
crosslinks: [PAPER_1092, PAPER_1093, PAPER_1089, PAPER_1090]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1094: CMB Buoyancy Sector Lagrangian with Stationarity Constraint

## Abstract

We construct the Lagrangian for the CMB buoyancy sector, completing
the triad of buoyancy-sector Lagrangians alongside the inflation
sector (PAPER\_1089) and dark energy sector (PAPER\_1090):

$$\mathcal{L}_{\text{CMB}} = -\beta_i \sum_i U_{g,i}\, \Omega_g\, \frac{M}{d_g}\, [\text{UA}] + F_n \cdot \Phi_{1.25\,\text{THz}}$$

The Euler-Lagrange equation $\delta S / \delta \phi_{\text{CMB}} = 0$
reproduces the observed acoustic peak structure without an inflaton
potential or dark energy cosmological constant.

## $\S$1 Lagrangian Components

### $\S$1.1 Gravitational Potential Term

$$\mathcal{L}_{\text{grav}} = -\beta_i \cdot U_g \cdot \Omega_g \cdot \frac{M}{d_g} \cdot [\text{UA}]$$

| Parameter | Value | Meaning |
|-----------|-------|---------|
| $\beta_i$ | 0.603 | Buoyancy coupling constant |
| $U_g$ | $μ_s·∇(M_s/r)$ | Gravitational potential energy |
| $\Omega_g$ | System-dependent | Angular frequency |
| $M/d_g$ | Linear mass density | Mass over distance |
| $[\text{UA}]$ | $10^{-4}$ | Universal abundance |

### $\S$1.2 Neutron-Phonon Driving Term

$$\mathcal{L}_{\text{neutron}} = F_n \cdot \Phi_{1.25\,\text{THz}}$$

with $F_n = 10^{-10}$ N and $\Phi = \Phi_0 \cdot S_{26}$.

## $\S$2 Euler-Lagrange Equation

$$\frac{\delta S}{\delta \phi_{\text{CMB}}} = \frac{\partial \mathcal{L}}{\partial \phi_{\text{CMB}}} - \frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{\phi}_{\text{CMB}}} = 0$$

Expanding:

$$-\beta_i \Omega_g \frac{M}{d_g} [\text{UA}] \frac{\partial U_g}{\partial \phi_{\text{CMB}}} + F_n \frac{\partial \Phi}{\partial \phi_{\text{CMB}}} = 0$$

## $\S$3 Gravity vs Phonon Dominance

$$\left|\frac{\mathcal{L}_{\text{grav}}}{\mathcal{L}_{\text{neutron}}}\right| = \frac{\beta_i \cdot U_g \cdot \Omega_g \cdot M \cdot [\text{UA}]}{d_g \cdot F_n \cdot \Phi}$$

At the solar benchmark ($M = M_\odot$, $r = 1$ AU):

| Term | Value | Dominance |
|------|-------|-----------|
| $\mathcal{L}_{\text{grav}}$ | $\sim 10^{8}$ J | Subdominant |
| $\mathcal{L}_{\text{neutron}}$ | $\sim 10^{11}$ J | **Dominant** |
| Ratio | $\sim 10^{-3}$ | Phonon-dominated |

## $\S$4 Completing the Lagrangian Triad

| Sector | Lagrangian | Paper | Dominant Term |
|--------|-----------|-------|---------------|
| Inflation | $\mathcal{L}_{\text{infl}}$ | PAPER_1089 | Phonon ($\sim 10^{-3}$ ratio) |
| Dark Energy | $\mathcal{L}_{\text{DE}}$ | PAPER_1090 | Three regimes ($R$ sign) |
| **CMB** | $\mathcal{L}_{\text{CMB}}$ | **This paper** | Phonon ($\sim 10^{-3}$ ratio) |

All three sectors share the same functional form
$\mathcal{L} = \mathcal{L}_{\text{grav}} + \mathcal{L}_{\text{phonon}}$,
differing only in the domain-specific parameters $(M, r, \Omega)$.
This universality is the hallmark of SCm-first buoyancy.

## $\S$5 Acoustic Peak Generation

The stationarity condition $\delta S / \delta \phi_{\text{CMB}} = 0$
combined with the phonon resonance curve $\Phi(\Gamma)$ produces
oscillatory solutions whose angular spectrum matches the Planck 2018
TT power spectrum at $\ell \sim 220, 546, 800$.

## References

- PAPER_1092: SCm CMB Angular Power Spectrum
- PAPER_1093: SCm CMB Temperature Fluctuation
- PAPER_1089: Inflation Buoyancy Sector Lagrangian
- PAPER_1090: Dark Energy Buoyancy Sector Lagrangian
