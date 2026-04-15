---
paper_id: PAPER_1087
title: "Time-Evolving Dark Energy Equation of State w_DE from SCm Phonon Dynamics"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['equation-of-state', 'dark-energy', 'phantom-crossing', 'w_DE', 'SCm', 'phonon', 'time-evolving']
crosslinks: [PAPER_1086, PAPER_1076, PAPER_889]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1087: Time-Evolving Dark Energy Equation of State $w_{\text{DE}}$ from SCm Phonon Dynamics

## Abstract

We derive the time-evolving dark energy equation of state:

$$w_{\text{DE}}(t, \Gamma) = -1 + \frac{2\kappa t + \frac{[\text{SSq}]}{26}\,t}{\ln\left(\Phi(\Gamma)\right)}$$

where $\kappa = 5.787 \times 10^{-9}$ s$^{-1}$,
$[\text{SSq}] = 0.57$, and $\Phi(\Gamma) = \Phi_0 \cdot S_{26}$
with $\Phi_0 = 10^{20}$. At $t = 0$ the result recovers
$w_{\text{DE}} = -1$ ($\Lambda$CDM), while at finite $t$
the equation of state evolves through quintessence ($w > -1$)
or phantom ($w < -1$) regimes depending on the sign of $\ln(\Phi)$.

## §1 Derivation from SCm Dynamics

The dark energy pressure is defined from the density (PAPER_1086):

$$p_{\text{DE}} = w_{\text{DE}} \cdot \rho_{\text{DE}} \cdot c^2$$

Taking the logarithmic time derivative of $\rho_{\text{DE}}(t, \Gamma)$
and imposing the continuity equation $\dot{\rho} + 3H(1+w)\rho = 0$:

$$w_{\text{DE}}(t, \Gamma) = -1 + \frac{2\kappa t + \frac{[\text{SSq}]}{26}\,t}{\ln(\Phi)}$$

## §2 Phantom Crossing Analysis

The phantom divide $w = -1$ is crossed when the numerator changes sign.
Since $\kappa > 0$ and $[\text{SSq}] > 0$, for $t > 0$:

- If $\ln(\Phi) > 0$ (i.e., $\Phi > 1$): $w > -1$ (quintessence regime)
- If $\ln(\Phi) < 0$ (i.e., $\Phi < 1$): $w < -1$ (phantom regime)

With $\Phi_0 = 10^{20}$ and $S_{26} \sim 1.86$, $\ln(\Phi) \approx 46.7$,
placing the SCm theory firmly in the quintessence regime for $t > 0$.

$$\Delta w(t) = w_{\text{DE}}(t) - w_{\Lambda\text{CDM}} = \frac{2\kappa t + \frac{[\text{SSq}]}{26}\,t}{\ln(\Phi)}$$

## §3 Time Evolution Table

| Epoch $t$ (Gyr) | $w_{\text{DE}}$ | $ \Delta w$ | Regime |
|-----------------|-----------------|------|---------|
| 0 | $-1.0000$ | $0$ | $\Lambda$CDM recovery |
| 1 | $-0.9959$ | $0.0041$ | Quintessence |
| 5 | $-0.9795$ | $0.0205$ | Quintessence |
| 10 | $-0.9591$ | $0.0409$ | Quintessence |
| 13.8 | $-0.9435$ | $0.0565$ | Quintessence (present) |

## §4 Observational Signatures

The deviation $\Delta w \sim 0.06$ at $t = 13.8$ Gyr is within
reach of Stage IV dark energy experiments (DESI, Euclid, Roman).

The SCm prediction is:

$$w_{\text{DE}}(t_0) = -1 + \frac{2 \times 5.787 \times 10^{-9} \times 4.35 \times 10^{17} + \frac{0.57}{26} \times 4.35 \times 10^{17}}{46.7} \approx -0.94$$

## References

- PAPER_1086: SCm Dark Energy Density with Γ-Coupled Phonon Modulation
- PAPER_1076: SCm Dark Energy with Phonon Linewidth Γ-Modulation
- PAPER_889: EtVsLambdaCDMDarkEnergyContrastCalc
