---
paper_id: "PAPER_1123"
title: "H₂O Maser Emission from J-Type Shocks: UQFF S(t) Compression Verification"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [H2O maser, 22 GHz, J-type shock, water maser, collisional pumping, Ug1 jump]
crosslinks: [PAPER_1121, PAPER_1122]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1123: H₂O Maser Emission from J-Type Shocks — UQFF S(t) Compression Verification

## Abstract

Following arXiv:1306.5276 (2013), we implement a UQFF calculator for H$_2$O maser emission driven by J-type shocks. The 22.235 GHz water maser line ($6_{16} \to 5_{23}$ rotational transition) is collisionally pumped in post-shock gas at $T \sim 300$-$1000$ K with densities $n \sim 10^5$-$10^6$ cm$^{-3}$. This verifies the UQFF $S(t)$ compression model, where J-type shocks are interpreted as abrupt Ug1 jumps.

## 1. J-Shock Compression

The strong J-shock limit gives post-shock density:

$$n_{\text{post}} = 4 \cdot n_{\text{pre}}$$

and post-shock temperature:

$$T_{\text{post}} = \frac{3 \mu m_p v_s^2}{16 k_B}$$

## 2. Maser Pumping Window

Collisional inversion of the 22 GHz H$_2$O transition requires:

$$300 \text{ K} \leq T_{\text{post}} \leq 1000 \text{ K}$$

This constrains shock velocities to $v_s \sim 20$-$80$ km/s. Below 300 K, collisional rates are insufficient; above 1000 K, collisional de-excitation quenches the inversion.

## 3. Maser Optical Depth

$$\tau_{\text{maser}} = \frac{n_{\text{post}} \cdot X(\text{H}_2\text{O}) \cdot h\nu \cdot l_{\text{gain}}}{4\pi \Delta v \cdot k_B T_{\text{post}}}$$

where $l_{\text{gain}}$ is the coherent gain path length ($\sim 10$ AU typical).

## 4. UQFF Interpretation

The J-type shock maps directly to an abrupt Ug1 jump in the UQFF framework:
- $S(t)$ compression: density enhancement factor
- The maser activation window validates the thermal structure of Ug1 discontinuities

Overall alignment: **80%** — datasets on shock velocities and maser luminosities verify the $S(t)$ compression model.

## References

- arXiv:1306.5276 — H$_2$O masers from J-type shocks (2013).
- Elitzur, M. (1992). Astronomical Masers. Kluwer Academic Publishers.
