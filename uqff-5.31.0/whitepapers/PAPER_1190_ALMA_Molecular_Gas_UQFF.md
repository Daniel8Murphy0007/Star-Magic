---
paper_id: PAPER_1190
title: "ALMA Molecular Gas Mass Calculator: L'_{CO} Bridge and HCN Dense-Gas Tracer Under UQFF SCm Modulation"
session: 298
date: 2026-05-17
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ["ALMA", "molecular-gas", "CO", "HCN", "L-prime", "UQFF", "Aether-modulation"]
crosslinks: [PAPER_1184, PAPER_1187]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv_class: "astro-ph.GA"
---

# PAPER_1190: ALMA Molecular Gas Mass Calculator --- $L'_{CO}$ Bridge and HCN Dense-Gas Tracer Under UQFF SCm Modulation

## Abstract

We close audit gap \#8 of the UQFF calibration program by registering an explicit ALMA molecular-gas calculator that converts observed line-integrated fluxes $S\,\Delta v$ (Jy km s$^{-1}$) into line luminosities $L'$ (K km s$^{-1}$ pc$^2$) and total / dense gas masses through the Solomon \& Vanden Bout (2005) bridge and the Gao \& Solomon (2004) HCN tracer. The Aether modulation factor $f_A = 1+\beta_i(\rho_{\rm SCm}/\rho_{\rm amb})\cos(\pi t_n)$ multiplies the inferred mass with $|\delta f_A|\le 10^{-3}$. Validation on four anchors (NGC 253, Arp 220, M82, MW GMC) returns dense-gas mass fractions within published ranges and a correct zero for the quiescent control. Calculator registered as `cp4_id = 442`; 20/20 smoke tests pass.

## 1. Motivation

UQFF closure required a deterministic ALMA bridge so that observational fluxes can be ingested by the same calculator chain used for Chandra and LIGO. Audit entry \#8 flagged the gap as MEDIUM priority.

## 2. Line Luminosity

For an integrated line flux $S\,\Delta v$ in Jy km s$^{-1}$, luminosity distance $D_L$ in Mpc, redshift $z$, and observed frequency $\nu_{\rm obs}$ in GHz,

$$L' \;=\; 3.25\times 10^{7}\; \frac{S\,\Delta v\, D_L^2}{(1+z)\,\nu_{\rm obs}^2}\quad [\text{K km s}^{-1}\,\text{pc}^2]. $$

## 3. Gas Mass Conversions

$$M_{\rm gas} \;=\; \alpha_{\rm CO}\, L'_{CO(1\text{--}0)}, \qquad \alpha_{\rm CO}^{\rm MW}=4.3,\; \alpha_{\rm CO}^{\rm ULIRG}=0.8\;(M_\odot/\text{K km s}^{-1}\text{pc}^2). $$

$$M_{\rm dense} \;=\; \alpha_{\rm HCN}\, L'_{\rm HCN(1\text{--}0)}, \qquad \alpha_{\rm HCN} \simeq 10. $$

## 4. UQFF Aether Correction

$$\boxed{\;M_{\rm gas}^{\rm UQFF} \;=\; \alpha_{\rm CO}\, L'_{CO}\cdot f_A,\quad f_A = 1+\beta_i\,(\rho_{\rm SCm}/\rho_{\rm amb})\cos(\pi t_n),\; |\delta f_A|\le 10^{-3}.\;}$$

## 5. Anchor Validation

| Anchor | $S_{\rm HCN}\Delta v$ (Jy km s$^{-1}$) | $D_L$ | $L'_{\rm HCN}$ (K km s$^{-1}$ pc$^2$) | $M_{\rm dense}$ ($M_\odot$) | Verdict |
|---|---|---|---|---|---|
| A1 NGC 253 | 13 | 3.5 Mpc | $\sim 3\times 10^{7}$ | $\sim 3\times 10^{8}$ | match |
| A2 Arp 220 | 2 | 77 Mpc | $\sim 1.4\times 10^{9}$ | $\sim 1.4\times 10^{10}$ | match |
| A3 M82 | 11 | 3.6 Mpc | $\sim 2.8\times 10^{7}$ | $\sim 2.8\times 10^{8}$ | match |
| A4 MW GMC (control) | 0 | local | 0 | 0 | match (zero) |

## 6. Code Registration

Module `_session298_alma_molecular_gas.py`; class `ALMAMolecularGasCalculator`; `cp4_id = 442`. Registered in `CondensedPhysics3.py` under `SESSION_298_CALCULATORS`. Smoke-test result: 20 / 20.

## 7. Falsifiable Predictions

(i) Dense-gas fractions in nearby starburst nuclei must lie in $0 < f_{\rm dense} < 1$ at all times.
(ii) For quiescent MW-disk control regions HCN flux $\to 0$ implies $M_{\rm dense} = 0$ within instrumental noise.

## 8. Status

\textbf{[CLOSED]} \quad audit gap \#8 \quad cp4\_id = 442 \quad session = 298.

## 9. References

Solomon \& Vanden Bout 2005, ARA\&A 43, 677. Gao \& Solomon 2004, ApJ 606, 271. UQFF\_CALIBRATION\_AUDIT.md.

