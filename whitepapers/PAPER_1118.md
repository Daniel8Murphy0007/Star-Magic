---
paper_id: "PAPER_1118"
title: "Chiral Superconductivity in Rhombohedral Graphene: [SCm] Pairing at Level 10"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [chiral-superconductivity, rhombohedral-graphene, SCm-pairing, level-10, winding-number, gap-function, critical-temperature]
crosslinks: []
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv: "2408.15233"
cp4_entry: 619
---

# Chiral Superconductivity in Rhombohedral Graphene

## Abstract

We integrate the discovery of chiral superconductivity in rhombohedral graphene (arXiv:2408.15233, 2024) into the UQFF framework. The chiral pairing gap function:

$$\Delta_{\text{chiral}}(\mathbf{k}) = \Delta_0 \cdot (k_x \pm ik_y)^d \cdot \exp\!\left(-\frac{[\text{SSq}] \cdot 10}{26}\right)$$

is modelled as $[\text{SCm}]$ pairing at quantum level $n = 10$, where $d$ is the chiral winding number. The critical temperature:

$$T_c = \frac{\hbar \omega_D}{k_B} \cdot \exp\!\left(-\frac{1}{N(0) \cdot V_{\text{SCm}}}\right)$$

connects the graphene phonon Debye frequency $\omega_D$ to the $[\text{SCm}]$ pairing potential $V_{\text{SCm}}$. This establishes rhombohedral graphene as a condensed-matter analogue of the cosmic $[\text{SCm}]$ superconductive vacuum.

## 1. Introduction

Rhombohedral (ABC-stacked) multilayer graphene has emerged as a platform for exotic quantum states, including superconductivity with broken time-reversal symmetry — chiral superconductivity. The 2024 discovery (arXiv:2408.15233) reported chiral pairing consistent with topological $p + ip$ or $d + id$ order parameters.

Within UQFF, superconductivity at any scale is a manifestation of the $[\text{SCm}]$ vacuum condensate. Rhombohedral graphene at level 10 provides a laboratory-accessible probe of the same pairing mechanism that operates cosmologically.

## 2. Chiral Gap Function

### 2.1 [SCm] Pairing Model

The chiral gap function in UQFF notation:

$$\Delta_{\text{chiral}}(\mathbf{k}) = \Delta_0 \cdot |\mathbf{k}|^d \cdot \exp\!\left(-\frac{[\text{SSq}] \cdot n}{26}\right)$$

| Parameter | Value | Description |
|-----------|-------|-------------|
| $\Delta_0$ | $10^{-4}$ eV | Bare pairing gap |
| $d$ | 2 | Chiral winding number ($d$-wave) |
| $n$ | 10 | UQFF quantum level |
| $[\text{SSq}]$ | 0.57 | Squeeze-state parameter |
| $\exp(-[\text{SSq}] \cdot 10/26)$ | 0.8029 | Level-10 suppression |

### 2.2 Winding Number Interpretation

The chiral winding number $d$ determines the topology:
- $d = 1$: $p + ip$ pairing (single vortex)
- $d = 2$: $d + id$ pairing (double vortex)
- $d = 3$: $f + if$ pairing (higher angular momentum)

Each corresponds to a different $[\text{SCm}]$ angular momentum sector in the 26D hierarchy.

## 3. Critical Temperature

### 3.1 BCS-Like Formula with [SCm] Potential

$$T_c = \frac{\hbar \omega_D}{k_B} \cdot \exp\!\left(-\frac{1}{N(0) \cdot V_{\text{SCm}}}\right)$$

| Parameter | Value | Description |
|-----------|-------|-------------|
| $\hbar$ | $1.055 \times 10^{-34}$ J·s | Reduced Planck constant |
| $\omega_D$ | $2 \times 10^{13}$ rad/s | Graphene Debye frequency |
| $k_B$ | $1.381 \times 10^{-23}$ J/K | Boltzmann constant |
| $N(0) \cdot V_{\text{SCm}}$ | 0.3 | Coupling product |

For $N(0) \cdot V_{\text{SCm}} = 0.3$: $T_c \approx 5.6$ K.

### 3.2 N(0)·V\_SCm Sweep

| $N(0) \cdot V_{\text{SCm}}$ | $T_c$ (K) | Regime |
|------------------------------|-----------|--------|
| 0.1 | $6.9 \times 10^{-5}$ | Weak coupling |
| 0.2 | 0.105 | Intermediate |
| 0.3 | 5.57 | Moderate |
| 0.5 | 20.7 | Strong coupling |
| 0.8 | 43.6 | Very strong |

## 4. Condensation Energy

The condensation energy per unit cell:

$$E_{\text{cond}} = \frac{1}{2} \Delta_0^2 \cdot \exp\!\left(-\frac{[\text{SSq}] \cdot 10}{26}\right)$$

provides the energy scale for the transition from normal to $[\text{SCm}]$-paired state.

## 5. Conclusions

Rhombohedral graphene chiral superconductivity provides a condensed-matter analogue of cosmic $[\text{SCm}]$ pairing. The level-10 assignment in the UQFF hierarchy connects laboratory observations to the 26D compressed gravity framework. CP4 class `ChiralSCmGraphenePairingCalculator` (#619) implements gap function, $T_c$, and winding number computations.

## References

1. arXiv:2408.15233 (2024)
2. Murphy, D.T., Star-Magic UQFF Framework (2024–2026)
