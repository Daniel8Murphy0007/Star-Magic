---
paper_id: "PAPER_1104"
title: "Unified F_U_Bi for SMBH Binary Merger Dynamics: Inside-to-Outside and Outside-to-Inside"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [SMBH, binary-merger, F_U_Bi, buoyancy, gravitational-waves, chirp-mass, inspiral, ringdown, 26-layer, UQFF]
crosslinks: [PAPER_258, PAPER_649, PAPER_1079]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER\_1104: Unified $F_{U,Bi}$ for SMBH Binary Merger Dynamics

## Abstract

We develop the complete buoyancy-mediated unified field theory for
supermassive black hole (SMBH) binary mergers, distinguishing
inside-to-outside (inspiral $\to$ ringdown) and outside-to-inside
(environment $\to$ binary) dynamics.  The two directional $F_{U,Bi}$
components are derived from the 26-layer UQFF gravity framework with
gravitational wave phonon coupling and chirp mass scaling.

## 1. Inside-to-Outside: Inspiral to Ringdown

The $F_{U,Bi}$ component radiating outward during the merger:

$$F_{U,Bi}^{\text{in}\to\text{out}} = \sum_{i=1}^{26}\left[\text{Ug1}_i + \text{Ug2}_i + \text{Ug3}_i + \text{Ug4}_i\right]\times\left(1 + \frac{\Phi}{F_U}\cdot\Delta E_{\text{gw}}\right)\times\left(\frac{\mathcal{M}_c}{M_{\text{total}}}\right)^{5/3}$$

### 26-Layer Gravity Terms

Each layer $i \in [1,26]$ has quantum state factor $Q_i = 1/(1+\kappa i)$:

| Term | Expression | Physics |
|------|------------|---------|
| $\text{Ug1}_i$ | $Q_i\cdot \mu_s\nabla(M_s/r)$ | DPM-seeded gravity |
| $\text{Ug2}_i$ | $(1-[\text{SSq}])Q_i\cdot \mu_s\nabla(M_s/r)\times0.1$ | Charge-reactivity |
| $\text{Ug3}_i$ | $[\text{SSq}]\,Q_i\cdot \mu_s\nabla(M_s/r)\times0.05\sin i$ | String rotation |
| $\text{Ug4}_i$ | $\kappa\,Q_i\times10^{-15}$ | Vacuum concentration |

### Chirp Mass Factor

For a binary with masses $M_1$ and $M_2$:

$$\mathcal{M}_c = M_{\text{total}}\cdot\eta^{3/5}, \qquad \eta = \frac{M_1 M_2}{M_{\text{total}}^2}$$

$$\left(\frac{\mathcal{M}_c}{M_{\text{total}}}\right)^{5/3} = \eta$$

### GW Phonon Coupling

$$\Delta E_{\text{gw}} = \frac{GM_1 M_2}{2r_{\text{sep}}}$$

The phonon-GW coupling factor:

$$1 + \frac{\Phi(\omega,\Gamma)}{F_U}\cdot\Delta E_{\text{gw}}$$

## 2. Outside-to-Inside: Environmental Buoyancy

The complementary $F_{U,Bi}$ component from the environment:

$$F_{U,Bi}^{\text{out}\to\text{in}} = F_U - F_{U,Bi}^{\text{in}\to\text{out}} - \sum_k F_{\text{tidal},k}(r_k, M_k)$$

where $F_{\text{tidal},k}$ are tidal forces from surrounding structures
(galaxy clusters, dark matter halos, accretion disk).

## 3. Net Buoyancy Imbalance

$$\Delta F = F_{U,Bi}^{\text{in}\to\text{out}} - F_{U,Bi}^{\text{out}\to\text{in}}$$

- $\Delta F > 0$: merger-dominated (radiative outflow exceeds environmental inflow)
- $\Delta F < 0$: environment-dominated (accretion exceeds radiation)
- $\Delta F = 0$: buoyancy equilibrium

## 4. Merger Timescale

Using the Peters formula for gravitational wave inspiral:

$$t_{\text{merger}} = \frac{5}{256}\,\frac{c^5 a_0^4}{G^3 M_1 M_2 M_{\text{total}}}$$

## 5. Implementation

Calculator: `UnifiedFUBiSMBHMergerDynamicsCalculator` in CondensedPhysics.py

- Inputs: $M_1$, $M_2$, $r_{\text{sep}}$, $F_U$, tidal forces, linewidth $\Gamma$
- Outputs: both $F_{U,Bi}$ components, net imbalance, chirp mass, merger timescale

## 6. Conclusion

The directional decomposition of $F_{U,Bi}$ during SMBH mergers provides a
complete description of buoyancy flow in the UQFF framework, with the
inside-to-outside component governed by chirp mass scaling and GW phonon
coupling, and the outside-to-inside component determined by the global
unified field budget less tidal contributions.

## References

- Peters, P.C. (1964). Gravitational radiation and the motion of two point masses. *Phys. Rev.* 136, B1224.
- Murphy, D.T. (2025). UQFF: Star Cradle Mechanics framework. Star-Magic repository.
