---
paper_id: PAPER_1083
title: "Planetary Core Energy Maintenance via Solar Wind Absorption"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['planetary-core', 'solar-wind', 'Um', 'Ug3', 'energy-budget', 'core-maintenance', 'depletion']
crosslinks: [PAPER_1078, PAPER_1079, PAPER_413]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1083: Planetary Core Energy Maintenance via Solar Wind Absorption

## Abstract

We derive the energy budget for planetary core maintenance through solar
wind absorption. The core receives power $P_{\text{wind}} = \eta_{\text{pen}} \cdot \Phi_{\text{sw}}(d) \cdot \pi R_{\text{mag}}^2$
from wind kinetic flux, which must balance Um maintenance cost
$P_{U_m} = \mu_{\text{core}}^2 \omega_{\text{core}} / \tau_{U_m}$ and
radiative losses $P_{\text{rad}} = \sigma T_{\text{core}}^4 \cdot 4\pi R_{\text{core}}^2$.
The net energy rate $dE/dt = P_{\text{wind}} - P_{U_m} - P_{\text{rad}}$
determines whether a planet's core is in thermal equilibrium or depleting.
The core SCm+UA is exclusively interactive with Ug3.

## §1 Wind Power Absorbed by Core

$$P_{\text{wind}} = \eta_{\text{pen}} \cdot (1 - \eta_{\text{liquid}}) \cdot \frac{1}{2} \rho_{\text{sw}}(d) \cdot v_{\text{sw}}^3 \cdot \pi R_{\text{mag}}^2$$

## §2 Um Maintenance Cost

$$P_{U_m} = \frac{\mu_{\text{core}}^2 \cdot \omega_{\text{core}}}{\tau_{U_m}}$$

where $\mu_{\text{core}}$ is the planetary dipole moment, $\omega_{\text{core}}$ the
rotation rate, and $\tau_{U_m} \sim 10^{15}$ s the magnetic diffusion timescale.

## §3 Radiative Loss

$$P_{\text{rad}} = \sigma T_{\text{core}}^4 \cdot 4\pi R_{\text{core}}^2$$

## §4 Core Energy Budget

$$\frac{dE_{\text{core}}}{dt} = P_{\text{wind}} - P_{U_m} - P_{\text{rad}}$$

$$\tau_{\text{depletion}} = \frac{E_{\text{core}}}{|dE/dt|}$$

Equilibrium condition: $dE/dt \geq 0 \iff P_{\text{wind}} \geq P_{U_m} + P_{\text{rad}}$.

## §5 Earth Example

| Parameter | Value | Source |
|-----------|-------|--------|
| $d$ | 1.0 AU | Orbital distance |
| $R_{\text{mag}}$ | $6.4 \times 10^7$ m | Earth magnetosphere |
| $R_{\text{core}}$ | $3.486 \times 10^6$ m | Inner+outer core |
| $T_{\text{core}}$ | 5700 K | Core temperature |
| $\mu_{\text{core}}$ | $8 \times 10^{22}$ A$\cdot$m$^2$ | Geomagnetic dipole |
| $P_{\text{wind}}$ | $\sim 10^{10}$ W | Wind power to core |
| $P_{U_m}$ | $\sim 10^{15}$ W | Magnetic maintenance |
| $P_{\text{rad}}$ | $\sim 10^{19}$ W | Radiative loss |
| $dE/dt$ | $< 0$ | Net depletion |

## §6 Ug3 Exclusivity

Core SCm+UA interacts exclusively with Ug3 (magnetic string rotation), not
with Ug1 (dipole) or Ug2 (charge-reactivity). This means:
- Core energy maintenance depends on Ug3 string coupling distance
- Beyond the frost line, Ug3 coupling weakens and wind power dominates (PAPER_1079)
- Core dynamo longevity correlates with Ug3 string tension retention

## §7 Planetary Comparison

| Planet | $P_{\text{wind}}$ | $P_{U_m}$ | $P_{\text{rad}}$ | Regime |
|--------|-------------------|-----------|-----------------|--------|
| Earth | $\sim 10^{10}$ W | $\sim 10^{15}$ W | $\sim 10^{19}$ W | Depleting (slowly) |
| Jupiter | $\sim 10^{14}$ W | $\sim 10^{18}$ W | $\sim 10^{21}$ W | Depleting |
| Neptune | $\sim 10^{11}$ W | $\sim 10^{13}$ W | $\sim 10^{17}$ W | Depleting |
| Mars | $\sim 10^{6}$ W | $\sim 0$ W | $\sim 10^{16}$ W | Dead dynamo |

## References

- PAPER_1078: Solar Wind Flux Partitioning Budget
- PAPER_1079: Frozen Planet Solar Wind Power
- PAPER_413: Ug3CCWCWDifferentialRotationSCmPlanetaryCoreCalculator
- Star-Magic.txt: §5 Core Energy Maintenance
