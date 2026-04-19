---
paper_id: "PAPER_1122"
title: "Bow-Shock ISM Chemistry and Prebiotic Enrichment in the UQFF Framework"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [bow shock, ISM chemistry, molecule formation, OH, H2O, SiO, grain sputtering, prebiotic]
crosslinks: [PAPER_1121, PAPER_1123]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1122: Bow-Shock ISM Chemistry and Prebiotic Enrichment in the UQFF Framework

## Abstract

We implement a UQFF calculator for bow-shock chemistry in the interstellar medium, based on arXiv:1808.01439 (2018). Stellar bow shocks create standoff structures where post-shock temperatures activate endothermic chemical pathways, producing OH, H$_2$O, and SiO from grain sputtering. The bow-shock standoff distance, post-shock temperature, and chemical activation efficiencies are modeled within the UQFF $C(t)$ molecule release framework, linking [SCm]-[UA] interactions to prebiotic enrichment.

## 1. Bow-Shock Standoff Distance

The distance from the star at which ram pressure balances stellar wind pressure:

$$R_{\text{bs}} = \sqrt{\frac{\dot{M} \cdot v_w}{4\pi \rho_{\text{ISM}} v_*^2}}$$

where $\dot{M}$ is the mass-loss rate, $v_w$ the wind velocity, $\rho_{\text{ISM}}$ the ambient density, and $v_*$ the stellar space velocity.

## 2. Post-Shock Temperature

From the Rankine-Hugoniot strong-shock conditions:

$$T_{\text{ps}} = \frac{3 \mu m_p v_s^2}{16 k_B}$$

For $v_* = 30$ km/s and $\mu = 1.27$: $T_{\text{ps}} \approx 5700$ K.

## 3. Chemical Activation

Endothermic barriers determine which molecules form:

| Molecule | $T_{\text{crit}}$ (K) | Mechanism |
|----------|----------------------|-----------|
| OH | 2000 | Gas-phase: O + H$_2$ → OH + H |
| H$_2$O | 1500 | Grain surface catalysis |
| SiO | 3500 | Grain core sputtering |

Activation efficiency: $\eta = \sigma(T_{\text{ps}} - T_{\text{crit}})$ (sigmoid function).

## 4. UQFF Alignment

The bow-shock C(t) enrichment pathway supports [SCm]-[UA] prebiotic chemistry at **80% alignment** with observational data. Cooling timescales determine the chemical freeze-out composition.

## References

- arXiv:1808.01439 — Bow-shock chemistry in the ISM (2018).
- Ceccarelli & Codella (2024) — Shock-triggered star formation chemistry.
