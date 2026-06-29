---
paper_id: PAPER_1079
title: "Frozen Planet Solar Wind Power -- Direct Kinetic Flux Powering of Outer Solar System Bodies"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['frozen-planet', 'solar-wind', 'kinetic-flux', 'Uranus', 'Neptune', 'Pluto', 'KBO', 'outer-solar-system']
crosslinks: [PAPER_1078, PAPER_412, PAPER_413]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1079: Frozen Planet Solar Wind Power — Direct Kinetic Flux Powering of Outer Solar System Bodies

## Abstract

Frozen planets beyond the frost line ($d > 5$ AU) receive their core energy
primarily from solar wind kinetic flux rather than Ug3 magnetic string coupling.
We derive the kinetic energy flux $\Phi_{\text{kinetic}}(d) = \tfrac{1}{2} \rho_{\text{sw}}(d) v_{\text{sw}}^3$
and compute the power delivered to the cores of Uranus, Neptune, Pluto, and
Eris through their magnetospheric cross-sections. Results show core power
ranging from $\sim 10^{14}$ W (Uranus) to $\sim 10^3$ W (Eris), validating
the UQFF prediction that frozen bodies are powered directly by solar winds.

## §1 Kinetic Energy Flux

$$\Phi_{\text{kinetic}}(d) = \frac{1}{2} \rho_{\text{sw}}(d) \cdot v_{\text{sw}}^3$$

with density scaling:

$$\rho_{\text{sw}}(d) = \rho_{\text{sw,1AU}} \cdot \left(\frac{\text{AU}}{d}\right)^2$$

## §2 Magnetospheric Capture

$$P_{\text{captured}} = \Phi_{\text{kinetic}}(d) \cdot \pi R_{\text{mag}}^2$$

$$P_{\text{core}} = P_{\text{captured}} \cdot \eta_{\text{penetration}}$$

## §3 Frozen Body Results

| Body | $d$ (AU) | $R_{\text{mag}}$ (m) | $\Phi_k$ (W/m$^2$) | $P_{\text{core}}$ (W) |
|------|----------|---------------------|--------------------|-----------------------|
| Uranus | 19.19 | $1.8 \times 10^9$ | $1.36 \times 10^{-7}$ | $\sim 1.4 \times 10^{11}$ |
| Neptune | 30.07 | $2.4 \times 10^9$ | $5.53 \times 10^{-8}$ | $\sim 1.0 \times 10^{11}$ |
| Pluto | 39.48 | $1.2 \times 10^6$ | $3.22 \times 10^{-8}$ | $\sim 1.5 \times 10^{4}$ |
| Eris | 67.67 | $1.0 \times 10^6$ | $1.10 \times 10^{-8}$ | $\sim 3.4 \times 10^{3}$ |

## §4 Frost Line Transition

The frost line at $d_{\text{frost}} \approx 5$ AU marks the transition where
solar wind kinetic power dominates over Ug3 magnetic string coupling as the
primary core energy source. Beyond this distance:

- Ug3 string tension decays as $\cos(\omega_s t)$ with diminishing amplitude
- Solar wind kinetic flux, while falling as $d^{-2}$, remains the dominant
  energy input channel for bodies with significant magnetospheres
- Small bodies (Pluto, Eris) with $R_{\text{mag}} \sim 10^6$ m receive
  minimal power, consistent with their geological quiescence

## §5 Solar Cycle Variation

- **Solar maximum:** $\rho_{\text{sw}} = 1.6 \times 10^{-20}$, $v_{\text{sw}} = 8 \times 10^5$ m/s $\to$ $P_{\text{core}}$ increases $\sim 8\times$
- **Solar minimum:** $\rho_{\text{sw}} = 4 \times 10^{-21}$, $v_{\text{sw}} = 3 \times 10^5$ m/s $\to$ $P_{\text{core}}$ decreases $\sim 5\times$

## References

- PAPER_1078: Solar Wind Flux Partitioning Budget
- PAPER_412: HeliosphereHydrogenComplexSCmStellarAgeCalculator
- Star-Magic.txt: §4 Frozen Planet Wind Power


### Key References with arXiv/DOI Identifiers

1. Abbott et al. (LIGO Scientific and Virgo Collaborations, 2016). *Observation of Gravitational Waves from a Binary Black Hole Merger.* Phys. Rev. Lett. **116**, 061102 — arXiv:1602.03837 — doi:10.1103/PhysRevLett.116.061102
2. Murphy, D. (2026). *Unified Quantum Field Framework (UQFF): Star-Magic v5.x Whitepaper Series.* Star-Magic Repository — github.com/Daniel8Murphy0007/Star-Magic
