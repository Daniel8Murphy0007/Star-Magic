---
title: "UQFF Prediction P13: DESI Y5 Dark Energy w(z) Second-Derivative Strict-Static Test"
session: 259
date: 2026-05-11
status: CLOSED
cvw: v2.0.0
paper_id: PAPER_1178
---

# UQFF Prediction P13: DESI Y5 Dark Energy $w(z)$ Second-Derivative Strict-Static Test

## 1. Setup

DESI 2024 BAO + SN Ia + CMB combined analysis reports a possible drift in the dark-energy
equation of state, parametrized as

$$
w(z) = w_0 + w_a \cdot \frac{z}{1 + z}, \qquad (w_0, w_a) = (-0.827 \pm 0.063,\; -0.75 \pm 0.29) .
$$

P7 (PAPER_1170) already constrains the *first* derivative: UQFF predicts $w_0 = -1$, $w_a = 0$
because the R26 vacuum energy is strictly static.

P13 extends to the *second* derivative.

## 2. UQFF Closed Form (P13)

Because $\rho_\Lambda = (V(0) + \rho_{R26} + \rho_{KK} + \rho_{BSFG})$ is *strictly* time-independent
in the closed Lagrangian (PAPER_1167, PAPER_1170), all derivatives of $w(z)$ vanish identically:

$$
\boxed{\frac{d^n w}{dz^n}\bigg|_{\rm UQFF} = 0 \quad \text{for all } n \ge 1.}
$$

In particular, the second derivative

$$
\frac{d^2 w}{dz^2} = 0 \pm 0
$$

is a *zero-parameter* prediction — UQFF has no freedom to evade it.

## 3. Numerical Comparison

Standard $w_0 w_a$CDM has $d^2 w/dz^2 = 2 w_a/(1+z)^3$. At $z = 0.5$:

| Model | $d^2 w/dz^2$ at $z=0.5$ |
|---|---|
| $\Lambda$CDM | 0 |
| DESI 2024 best fit ($w_a = -0.75$) | $-0.444$ |
| **UQFF (P13)** | **0** |

## 4. Falsifiable Prediction (P13)

DESI Y5 (2028) is forecast to constrain $w_a$ at the $\pm 0.05$ level after BAO + SN + CMB
combination. The implied second-derivative bound at $z = 0.5$ is

$$
|d^2 w/dz^2| < 0.030 \quad \text{(DESI Y5 1-}\sigma\text{)} .
$$

**Falsifier:** if DESI Y5 confirms $w_a$ deviating from 0 at the $\ge 3\sigma$ level
(equivalently $|d^2 w/dz^2| > 0.030$ at $z = 0.5$ at $\ge 3\sigma$), UQFF static-vacuum is excluded.

## 5. Cross-Check with CP4 #262

CP4 calculator `UQFFDarkEnergySecondDerivativeCalculator` returns
$d^2 w/dz^2 = 0.0$ at all $z$ for default inputs, with the DESI Y5 bound included for comparison.

## 6. Status

CLOSED. P13 added to predictions table alongside P6--P12.
P13 is unique among the predictions for being a **strict zero with no UQFF freedom**.
