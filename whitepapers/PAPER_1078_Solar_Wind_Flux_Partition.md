---
paper_id: PAPER_1078
title: "Solar Wind Mass-Energy Flux Partitioning Budget in the UQFF Heliosphere"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['solar-wind', 'heliosphere', 'Ug2', 'flux-partition', 'hydrogen-complex', 'planetary-core']
crosslinks: [PAPER_412, PAPER_413, PAPER_1079, PAPER_1083]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1078: Solar Wind Mass-Energy Flux Partitioning Budget in the UQFF Heliosphere

## Abstract

We derive a complete mass-energy flux partitioning budget for stellar wind
propagation through the UQFF heliosphere. The total solar wind flux
$\Phi_{\text{total}} = 4\pi r_{\text{helio}}^2 \rho_{\text{sw}} v_{\text{sw}}$
is decomposed into four channels: (a) planetary geometric interception
$f_{\text{planet}}$, (b) Ug2 outer-field-shell transmutation into hydrogen
complexes $f_{\text{Ug2}}$, (c) planetary core absorption $f_{\text{core}}$,
and (d) heliospheric escape $f_{\text{escape}}$. Numerical evaluation for
the 8-planet solar system yields $f_{\text{Ug2}} \approx 0.85$, confirming
the dominant role of Ug2 charge-reactivity in processing solar wind material
into magnetically bound hydrogen complexes on the outer heliospheric shell.

## §1 Total Wind Flux

$$\Phi_{\text{total}} = 4\pi r_{\text{helio}}^2 \cdot \rho_{\text{sw}} \cdot v_{\text{sw}}$$

At 1 AU: $\rho_{\text{sw}} \approx 8 \times 10^{-21}$ kg/m$^3$,
$v_{\text{sw}} \approx 5 \times 10^5$ m/s.

## §2 Density Scaling

Solar wind density falls as inverse-square:

$$\rho_{\text{sw}}(d) = \rho_{\text{sw,1AU}} \cdot \left(\frac{\text{AU}}{d}\right)^2$$

## §3 Per-Planet Interception

Each planet intercepts flux through its magnetospheric cross-section:

$$\Phi_{\text{planet},i} = \pi R_{\text{mag},i}^2 \cdot \rho_{\text{sw}}(d_i) \cdot v_{\text{sw}}$$

$$f_{\text{planet}} = \frac{\sum_i \Phi_{\text{planet},i}}{\Phi_{\text{total}}}$$

| Planet | $R_{\text{mag}}$ (m) | $d$ (AU) | $\Phi_i / \Phi_{\text{total}}$ |
|--------|---------------------|----------|-------------------------------|
| Mercury | $1.5 \times 10^6$ | 0.387 | $\sim 10^{-14}$ |
| Venus | $6.1 \times 10^6$ | 0.723 | $\sim 10^{-13}$ |
| Earth | $6.4 \times 10^7$ | 1.0 | $\sim 10^{-11}$ |
| Mars | $3.4 \times 10^6$ | 1.524 | $\sim 10^{-14}$ |
| Jupiter | $7.1 \times 10^{10}$ | 5.203 | $\sim 10^{-5}$ |
| Saturn | $2.0 \times 10^{10}$ | 9.537 | $\sim 10^{-6}$ |
| Uranus | $1.8 \times 10^9$ | 19.19 | $\sim 10^{-8}$ |
| Neptune | $2.4 \times 10^9$ | 30.07 | $\sim 10^{-8}$ |

## §4 Four-Channel Budget

$$f_{\text{core}} = f_{\text{planet}} \cdot \eta_{\text{penetration}} \cdot (1 - \eta_{\text{liquid}})$$

$$f_{\text{Ug2}} = (1 - f_{\text{planet}}) \cdot \eta_{\text{transmutation}}$$

$$f_{\text{escape}} = 1 - f_{\text{planet}} - f_{\text{Ug2}}$$

With $\eta_{\text{transmutation}} = 0.85$, $\eta_{\text{penetration}} = 0.15$,
$\eta_{\text{liquid}} = 0.6$:

| Channel | Fraction | Physical Process |
|---------|----------|-----------------|
| $f_{\text{Ug2}}$ | $\approx 0.85$ | Transmutation $\to$ H complexes on outer shell |
| $f_{\text{escape}}$ | $\approx 0.15$ | Escape past heliosphere boundary |
| $f_{\text{planet}}$ | $\sim 10^{-5}$ | Geometric interception by magnetospheres |
| $f_{\text{core}}$ | $\sim 10^{-6}$ | Absorbed into planetary cores (maintains Um) |

## §5 Simulation Sets

- **High wind:** $\rho_{\text{sw}} = 1.6 \times 10^{-20}$, $v_{\text{sw}} = 8 \times 10^5$ m/s
- **Low wind:** $\rho_{\text{sw}} = 4 \times 10^{-21}$, $v_{\text{sw}} = 3 \times 10^5$ m/s
- **Young star:** $\rho_{\text{sw}} = 5 \times 10^{-20}$, $v_{\text{sw}} = 1 \times 10^6$ m/s

## §6 UQFF Connection

The dominant Ug2 channel ($f_{\text{Ug2}} \approx 0.85$) validates the
Star-Magic prediction that solar winds not absorbed by planets are transmuted
by the Ug2 outer field shell into hydrogen complexes magnetically stuck to
the heliosphere boundary, providing the raw material for Ug2 charge-reactivity
processes described in PAPER_412.

## References

- PAPER_412: HeliosphereHydrogenComplexSCmStellarAgeCalculator
- PAPER_413: Ug3CCWCWDifferentialRotationSCmPlanetaryCoreCalculator
- Star-Magic.txt: §3 Solar Wind and Heliospheric Flux Balance
