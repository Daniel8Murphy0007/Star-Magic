---
paper_id: PAPER_1088
title: "F_U_Bi_i Seven-Component Force Decomposition: Phonon, Inflation, BCS, VDS, DVP, BSH, and QCalcGeom Sectors"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['F_U_Bi_i', 'decomposition', 'seven-component', 'BCS', 'VDS', 'DVP', 'BSH', 'QCalcGeom', 'phonon', 'inflation']
crosslinks: [PAPER_1089, PAPER_1090, PAPER_892, PAPER_893]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1088: $F_{U,Bi,i}$ Seven-Component Force Decomposition

## Abstract

The unified buoyancy force $F_{U,Bi,i}$ is decomposed into
seven orthogonal sectors that account for all fundamental
buoyancy-mediated interactions in the UQFF framework:

$$F_{U,Bi,i} = F_{\text{phonon}} + F_{\text{inflation}} + F_{\text{BCS}} + F_{\text{VDS}} + F_{\text{DVP}} + F_{\text{BSH}} + F_{\text{QCalcGeom}}$$

Each sector contributes a fractional share to the total buoyancy
budget measured against the base gravitational acceleration
$g_{	ext{base}} = \mu_s\nabla(M_s/r)$.

## §1 Base Gravity and Component Definitions

$$g_{\text{base}} = \frac{G M}{r^2}, \quad G = 6.674 \times 10^{-11}\;\text{m}^3/(\text{kg}\cdot\text{s}^2)$$

The seven force components are:

### §1.1 Phonon Sector

$$F_{\text{phonon}} = \Phi_{1.25\,\text{THz}} \cdot g_{\text{base}} \cdot M$$

Vacuum phonon resonance at 1.25 THz, coupling buoyancy to
the SCm lattice oscillation spectrum.

### §1.2 Inflation Sector

$$F_{\text{inflation}} = \beta_i \cdot U_g \cdot \frac{M}{d} \cdot [\text{UA}]$$

Buoyancy force driving inflationary expansion,
with $\beta_i = 0.603$ and $[\text{UA}] = 10^{-4}$.

### §1.3 BCS (Buoyancy-Condensate-Superfluid) Sector

$$F_{\text{BCS}} = \Delta_{\text{BCS}}^2 \cdot g_{\text{base}} \cdot M$$

Cooper-pair condensate contribution from the BCS
gap parameter $\Delta_{\text{BCS}}$.

### §1.4 VDS (Vacuum Density Series) Sector

$$F_{\text{VDS}} = \rho_{\text{vac}} \cdot V_{\text{eff}} \cdot g_{\text{base}}$$

Vacuum density series contribution encoding
$k_\eta = 10^{-113}$.

### §1.5 DVP (Dipole Vortex Primes) Sector

$$F_{\text{DVP}} = \sum_{p \in \mathcal{P}} \frac{\mu_p}{r_p^3} \cdot g_{\text{base}} \cdot M$$

Vortex-line contributions indexed by
prime-number dipole configurations.

### §1.6 BSH (Buoyancy Shell Harmonics) Sector

$$F_{\text{BSH}} = \sum_{\ell=0}^{26} Y_\ell^0(\theta) \cdot g_{\text{base}} \cdot M$$

Spherical-harmonic shell decomposition across 26 layers.

### §1.7 QCalcGeom (Quantum Calculator Geometry) Sector

$$F_{\text{QCalcGeom}} = \alpha_{\text{QCG}} \cdot g_{\text{base}} \cdot M$$

Geometric correction from the quantum calculator
manifold structure.

## §2 Fractional Contribution Table

| Sector | Force Component | Fraction of $F_{U,Bi,i}$ |
|--------|----------------|--------------------------|
| Phonon | $F_{\text{phonon}}$ | Dominant at resonance |
| Inflation | $F_{\text{inflation}}$ | $\beta_i / \sum$ |
| BCS | $F_{\text{BCS}}$ | $\sim 10^{-3}$ |
| VDS | $F_{\text{VDS}}$ | $\sim 10^{-5}$ |
| DVP | $F_{\text{DVP}}$ | $\sim 10^{-4}$ |
| BSH | $F_{\text{BSH}}$ | $\sim 10^{-2}$ |
| QCalcGeom | $F_{\text{QCalcGeom}}$ | $\sim 10^{-3}$ |

## §3 Budget Closure

$$\sum_{k=1}^{7} f_k = 1, \quad f_k \equiv \frac{F_k}{F_{U,Bi,i}}$$

The seven fractions are computed dynamically per system and
must satisfy exact budget closure to machine precision.

## References

- PAPER_1089: Inflation Buoyancy Sector Lagrangian
- PAPER_1090: Dark Energy Buoyancy Sector Lagrangian
- PAPER_892: FUBiUnifiedFieldBuoyancyForceCalc
- PAPER_893: GravBuoyancyForceCalculator
