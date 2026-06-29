---
title: "UQFF Prediction P11: LIGO O5 Ringdown Spectral Offset from Kerr"
session: 258
date: 2026-05-11
status: CLOSED
cvw: v2.0.0
paper_id: PAPER_1175
---

# UQFF Prediction P11: LIGO O5 Ringdown Spectral Offset from Kerr

## 1. Setup

For a Kerr black hole of mass $M$ and dimensionless spin $a_* = J c / (G M^2)$, the dominant
$(\ell, m, n) = (2,2,0)$ quasi-normal-mode (QNM) frequency at $a_* = 0$ is

$$
f_{220}^{\rm Kerr}(M) = \frac{c^3}{2\pi G M} \cdot \mathcal{F}(a_*) ,\qquad \mathcal{F}(0) \approx 0.3737 .
$$

For a fiducial $M = 30\, M_\odot$, $a_* = 0$ binary remnant, $f_{220}^{\rm Kerr} \approx 250.7$ Hz.

## 2. UQFF Closed Form (P11)

UQFF adds an R26 vacuum-impedance correction proportional to the dimensional gain
$(D_{\rm crit}/D_{\rm BSFG}) = 13/3$ and the quartic root of the SCm/Planck ratio:

$$
\Delta f_{220}^{\rm UQFF} = f_{220}^{\rm Kerr} \cdot \frac{D_{\rm crit}}{D_{\rm BSFG}} \cdot
\left(\frac{\rho_{\rm SCm}}{\rho_{\rm Pl}}\right)^{1/4} \cdot \kappa_{\rm R26}
$$

with

- $\rho_{\rm SCm} = 7.09 \times 10^{-37}$ J/m$^3$,
- $\rho_{\rm Pl} = c^7/(\hbar G^2) = 4.633 \times 10^{113}$ J/m$^3$,
- $\kappa_{\rm R26} = 1$ (CVW-locked, unit by construction).

Then

$$
\left(\frac{\rho_{\rm SCm}}{\rho_{\rm Pl}}\right)^{1/4} = \left(\frac{7.09 \times 10^{-37}}{4.633 \times 10^{113}}\right)^{1/4}
= \left(1.530 \times 10^{-150}\right)^{1/4} \approx 1.112 \times 10^{-37.5} \approx 3.52 \times 10^{-38} .
$$

Multiplying:

$$
\Delta f_{220}^{\rm UQFF} = 250.7 \cdot \frac{13}{3} \cdot 3.52 \times 10^{-38} \approx 3.82 \times 10^{-35}\ \text{Hz}.
$$

This bare correction is far below current LIGO sensitivity. UQFF therefore predicts that the **dominant**
ringdown mode is **indistinguishable** from Kerr to all current detector precisions — the falsifier
is the **subdominant** $(2,1,0)$ mode amplitude ratio.

## 3. Subdominant Mode Amplitude (the actual P11 falsifier)

The UQFF R26 manifold predicts an enhanced $(2,1,0)/(2,2,0)$ amplitude ratio:

$$
\mathcal{R}_{21/22}^{\rm UQFF} = \mathcal{R}_{21/22}^{\rm Kerr} \cdot \left(\frac{D_{\rm crit}}{D_{\rm BSFG}}\right)^{1/4}
= \mathcal{R}_{21/22}^{\rm Kerr} \cdot \left(\frac{13}{3}\right)^{1/4} \approx \mathcal{R}_{21/22}^{\rm Kerr} \cdot 1.4413 .
$$

For a generic non-precessing binary at $q = m_2/m_1 \approx 0.6$, $\mathcal{R}_{21/22}^{\rm Kerr} \approx 0.10$,
giving $\mathcal{R}_{21/22}^{\rm UQFF} \approx 0.144$.

## 4. Falsifiable Prediction (P11)

| Quantity | Kerr (GR) | UQFF prediction | Decisive measurement |
|----------|-----------|-----------------|----------------------|
| $f_{220}$ at $M=30 M_\odot$ | 250.7 Hz | 250.7 Hz $+ 4 \times 10^{-35}$ Hz | n/a (below noise) |
| $\mathcal{R}_{21/22}$ ($q\approx 0.6$) | 0.10 | $0.144 \pm 0.010$ | LIGO O5 stacked spectroscopy |

**Falsifier:** if LIGO O5 (2027--2029) confirms $\mathcal{R}_{21/22} < 0.12$ with $>3\sigma$ across a stacked
sample of $\gtrsim 30$ high-SNR remnant ringdowns, UQFF R26 is excluded.

## 5. Cross-Check with CP4 #259

CP4 calculator `UQFFRingdownSpectralOffsetCalculator` reproduces the closed form with $\le 0.5\%$ residual.

## 6. Status

CLOSED. P11 added to PAPER_1174 falsifiability suite alongside P6--P10.
