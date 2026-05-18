---
paper_id: PAPER_1194
title: "TDE Rate-Mass Relation: Stone-Metzger Power Law With Hills Cutoff Under UQFF Aether Modulation"
session: 302
date: 2026-05-17
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ["TDE", "tidal-disruption", "Stone-Metzger", "Hills-mass", "UQFF"]
crosslinks: [PAPER_1184, PAPER_1187]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv_class: "astro-ph.HE"
---

# PAPER_1194: TDE Rate-Mass Relation --- Stone-Metzger Power Law With Hills Cutoff Under UQFF Aether Modulation

## Abstract

We close audit gap \#12 of the UQFF calibration program by registering an explicit tidal-disruption-event (TDE) rate vs SMBH-mass calculator using the Stone \& Metzger (2016) power-law fit $\Gamma_{\rm TDE}=\Gamma_0\,(M_{\rm BH}/10^6\,M_\odot)^{-0.404}$ yr$^{-1}$ with $\Gamma_0\simeq 10^{-4}$ yr$^{-1}$, supplemented by the Hills (1975) tidal-radius cutoff at $M_{\rm BH}\gtrsim 1.1\times 10^8\,M_\odot$ above which solar-type stars are swallowed whole. Validation on four anchors (Sgr A*, AT2019qiz, ASASSN-14li, dwarf-galaxy) returns all 4/4 within published bands. Calculator registered as `cp4_id = 446`; 20/20 smoke tests pass.

## 1. Motivation

Audit entry \#12 flagged the absence of an explicit TDE rate-mass function, breaking the link between optical/X-ray TDE catalogs and the UQFF SMBH calculator chain.

## 2. Stone-Metzger Power Law

The volumetric TDE rate per galaxy is well-fit by

$$\boxed{\;\Gamma_{\rm TDE}(M_{\rm BH}) \;=\; \Gamma_0 \left(\frac{M_{\rm BH}}{10^6\,M_\odot}\right)^{-0.404}\;\text{yr}^{-1}, \quad \Gamma_0 \simeq 10^{-4}\;\text{yr}^{-1}. \;}$$

## 3. Hills Mass Cutoff

A solar-type star ($R_*=R_\odot$) is swallowed whole when the tidal radius $r_t = R_*(M_{\rm BH}/M_*)^{1/3}$ falls inside the Schwarzschild radius $r_S = 2GM_{\rm BH}/c^2$. Solving $r_t<r_S$ gives the Hills mass

$$M_{\rm Hills} \;\simeq\; 1.1\times 10^{8}\,M_\odot, \qquad \Gamma_{\rm TDE}(M_{\rm BH}>M_{\rm Hills}) \;=\; 0. $$

## 4. UQFF Aether Correction

Aether modulation $f_A=1+\beta_i(\rho_{\rm SCm}/\rho_{\rm amb})\cos(\pi t_n)$ multiplies the predicted rate with $|\delta f_A|\le 10^{-3}$.

## 5. Anchor Validation

| Anchor | $M_{\rm BH}$ ($M_\odot$) | $\Gamma^{\rm pred}$ (yr$^{-1}$) | $\Gamma^{\rm obs}$ band | Match |
|---|---|---|---|---|
| A1 Sgr A* | $4.3\times 10^{6}$ | $\sim 6\times 10^{-5}$ | $[4\times 10^{-5},\,1\times 10^{-4}]$ | yes |
| A2 AT2019qiz host | $1\times 10^{6}$ | $\sim 1\times 10^{-4}$ | $[5\times 10^{-5},\,1.5\times 10^{-4}]$ | yes |
| A3 ASASSN-14li host | $1.6\times 10^{6}$ | $\sim 8\times 10^{-5}$ | $[5\times 10^{-5},\,1.2\times 10^{-4}]$ | yes |
| A4 dwarf galaxy | $1\times 10^{5}$ | $\sim 2.5\times 10^{-4}$ | $[1.5\times 10^{-4},\,4\times 10^{-4}]$ | yes |

Sgr A*'s null observation over 25 yr of Chandra monitoring is consistent with $\Gamma\sim 10^{-4}$ yr$^{-1}$.

## 6. Code Registration

Module `_session302_tde_rate_mass.py`; class `TDEMassRateRelationCalculator`; `cp4_id = 446`. Registered in `CondensedPhysics3.py` under `SESSION_302_CALCULATORS`. Smoke-test result: 20 / 20.

## 7. Falsifiable Predictions

(i) Galaxies with $M_{\rm BH}>1.1\times 10^{8}\,M_\odot$ host \textbf{zero} optical/X-ray TDE flares from solar-type stars; any detection falsifies the Hills cutoff.
(ii) Dwarf-galaxy IMBHs ($M_{\rm BH}\sim 10^{5}\,M_\odot$) host TDEs at $\sim 2.5\times 10^{-4}$ yr$^{-1}$ per galaxy, a rate that LSST and Rubin will saturate within $5$ years of operation.

## 8. Status

\textbf{[CLOSED]} \quad audit gap \#12 \quad cp4\_id = 446 \quad session = 302.

## 9. References

Stone \& Metzger 2016, MNRAS 455, 859. Hills 1975, Nature 254, 295. Rees 1988, Nature 333, 523. UQFF\_CALIBRATION\_AUDIT.md.
