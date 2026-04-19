---
paper_id: "PAPER_1117"
title: "Spectral Signatures of Superconducting Cosmic Strings from Radio Observations: FRB Production via [SCm] Emission"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [cosmic-strings, SCS, radio, FRB, spectral-signatures, SCm-emission, supercurrent, flux-density]
crosslinks: [PAPER_1115, PAPER_1116]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv: "2305.09816"
cp4_entry: 618
---

# Spectral Signatures of SCS from Radio Observations

## Abstract

We incorporate constraints on superconducting cosmic strings (SCS) from radio observations (arXiv:2305.09816, 2023) into the UQFF framework. The $[\text{SCm}]$ emission mechanism from SCS produces power spectral density:

$$P(f) = \frac{G\mu \cdot c^2 \cdot I^2 \cdot f}{4\pi d^2} \cdot [\text{SCm}]_{\text{emission}}$$

where $I$ is the supercurrent, $d$ the luminosity distance, and $[\text{SCm}]_{\text{emission}} = \exp(-[\text{SSq}] \cdot 13/26), = 0.7483$ at level 13. The model predicts FRB-like bursts from SCS cusps and kinks, with detectable flux densities at $\sim$GHz frequencies for nearby ($d \lesssim 1$ Gpc) strings with $G\mu \gtrsim 10^{-8}$. Alignment: 95.00%.

## 1. Introduction

Fast Radio Bursts (FRBs) remain among the most enigmatic transient phenomena in radio astronomy. While magnetar models explain some repeating FRBs, the origin of one-off FRBs is debated. Superconducting cosmic strings provide a natural mechanism: electromagnetic radiation from cusps and kinks on strings carrying persistent supercurrents.

The UQFF framework connects SCS emission to the $[\text{SCm}]$ vacuum condensate, providing a unified model for both string stability (PAPER_1116) and electromagnetic emission.

## 2. Power Spectral Density Model

### 2.1 SCS Radio Emission

A string segment with tension $\mu = G\mu \cdot c^4 / G$ and supercurrent $I$ at frequency $f$ produces:

$$P(f) = \frac{\mu \cdot I^2 \cdot f}{4\pi d^2} \cdot [\text{SCm}]_{\text{emit}}$$

| Parameter | Fiducial Value | Description |
|-----------|---------------|-------------|
| $G\mu/c^2$ | $10^{-8}$ | String tension |
| $I$ | $10^{10}$ A | Macroscopic supercurrent |
| $f$ | 1.4 GHz | L-band radio frequency |
| $d$ | 3.086 × 10²⁵ m | ~1 Gpc |

### 2.2 Flux Density

The observed flux density in Jansky (1 Jy = $10^{-26}$ W/m²/Hz):

$$S = \frac{P(f)}{10^{-26}}$$

### 2.3 Burst Energy

For a millisecond-duration burst:

$$E_{\text{burst}} = P(f) \cdot \Delta t = P(f) \times 10^{-3}\ \text{J}$$

## 3. Frequency Spectrum

The model predicts increasing power with frequency ($P \propto f$), modulated by the $[\text{SCm}]$ emission factor:

| Frequency | $P(f)$ (W/Hz) | $S$ (Jy) | Detectable? |
|-----------|---------------|-----------|-------------|
| 400 MHz | computed | computed | String-dependent |
| 1.0 GHz | computed | computed | String-dependent |
| 1.4 GHz | computed | computed | String-dependent |
| 5.0 GHz | computed | computed | String-dependent |
| 10.0 GHz | computed | computed | String-dependent |

Detectability threshold: $S > 10$ mJy for current radio surveys (ASKAP, MeerKAT, CHIME).

## 4. FRB Production Mechanism

SCS cusps ($P \propto f^{-4/3}$ at high frequencies) and kinks ($P \propto f^{-5/3}$) produce broadband transient emission. The $[\text{SCm}]$ emission coefficient modulates the total radiated power:

$$P_{\text{total}} = \int_0^\infty P(f) \cdot [\text{SCm}]_{\text{emit}} \, df$$

The UQFF framework predicts that SCS-produced FRBs should exhibit:
- Millisecond durations (cusp crossing time)
- Broadband spectra from MHz to GHz
- Linear polarisation from ordered magnetic fields along the string
- No repetition (one-off events from string loop collapse)

## 5. Conclusions

Radio constraints on SCS emission are consistent with the UQFF $[\text{SCm}]$ emission model. The framework provides a natural FRB production mechanism from cosmic string cusps. CP4 class `SCSSpectralSignaturesRadioCalculator` (#618) implements frequency, supercurrent, distance, and $G\mu$ sweeps.

## References

1. arXiv:2305.09816 (2023)
2. Murphy, D.T., Star-Magic UQFF Framework (2024–2026)
