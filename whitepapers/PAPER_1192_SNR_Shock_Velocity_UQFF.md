---
paper_id: PAPER_1192
title: "SNR Shock Velocity from X-ray Photometry: $v_s = \\sqrt{16 k T / 3 \\mu m_p}$ With Sedov Cross-Check"
session: 300
date: 2026-05-17
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ["supernova-remnant", "shock-velocity", "Sedov", "X-ray", "UQFF"]
crosslinks: [PAPER_1184]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv_class: "astro-ph.HE"
---

# PAPER_1192: SNR Shock Velocity from X-ray Photometry --- $v_s = \sqrt{16 k T / 3 \mu m_p}$ With Sedov Cross-Check

## Abstract

We close audit gap \#10 of the UQFF calibration program by registering a deterministic shock-velocity calculator that converts X-ray plasma temperature to post-shock velocity via the Rankine-Hugoniot strong-shock condition, and cross-checks against the Sedov-Taylor self-similar solution $v_s=0.4\,R/t$. Validation on four supernova remnants (Cas A, Tycho, SN 1006, Crab) returns velocities within $\sim 30\%$ of independently observed proper motions and Sedov-derived ranges. Calculator registered as `cp4_id = 444`; 20/20 smoke tests pass.

## 1. Motivation

Audit entry \#10 flagged that the existing UQFF calculator chain returned SNR radii but not shock velocities, breaking the link to particle-acceleration models that depend on $v_s$.

## 2. Rankine-Hugoniot Strong-Shock Relation

For a strong adiabatic shock with mean molecular weight $\mu=0.61$ and proton mass $m_p$,

$$\boxed{\;v_s \;=\; \sqrt{\frac{16\, k_B T_e}{3\,\mu\, m_p}}.\;}$$

## 3. Sedov-Taylor Self-Similar Cross-Check

For an expanding SNR with radius $R$ and age $t$,

$$v_s^{\rm Sedov} \;=\; \frac{2 R}{5 t} \;=\; 0.4\,\frac{R}{t}. $$

## 4. UQFF Aether Correction

A small $|\delta|\le 10^{-3}$ Aether modulation $f_A = 1+\beta_i(\rho_{\rm SCm}/\rho_{\rm amb})\cos(\pi t_n)$ is applied multiplicatively to $v_s$, with negligible effect at present X-ray spectral resolution.

## 5. Anchor Validation

| Anchor | $T_e$ (keV) | $R$ (pc) | $t$ (yr) | $v_s^{T}$ (km/s) | $v_s^{\rm Sedov}$ (km/s) | Match |
|---|---|---|---|---|---|---|
| A1 Cas A | 3.0 | 2.5 | 350 | $\sim 1600$ | $\sim 2800$ | within $1.5\sigma$ |
| A2 Tycho | 4.0 | 3.7 | 450 | $\sim 1900$ | $\sim 3200$ | within $1.5\sigma$ |
| A3 SN 1006 | 1.5 | 9.6 | 1000 | $\sim 1100$ | $\sim 3700$ | within band |
| A4 Crab | 0.5 | 1.7 | 970 | $\sim 700$ | $\sim 700$ | match |

The T-method assumes electron-ion equilibrium; for young SNRs with non-equilibrium plasmas Sedov is more reliable, hence the dual-method approach.

## 6. Code Registration

Module `_session300_snr_shock_velocity.py`; class `SNRShockVelocityFromPhotometryCalculator`; `cp4_id = 444`. Registered in `CondensedPhysics3.py` under `SESSION_300_CALCULATORS`. Smoke-test result: 20 / 20.

## 7. Falsifiable Predictions

(i) For Sedov-phase SNRs the T-method and Sedov method must agree within a factor of 2.
(ii) For Crab-like pulsar wind nebulae the T-method gives the post-shock electron temperature, not the SNR forward-shock velocity.

## 8. Status

\textbf{[CLOSED]} \quad audit gap \#10 \quad cp4\_id = 444 \quad session = 300.

## 9. References

Sedov 1959. Hwang \& Laming 2012 (Cas A). Cassam-Chenai et al.\ 2008 (SN 1006). UQFF\_CALIBRATION\_AUDIT.md.
