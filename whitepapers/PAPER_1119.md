---
paper_id: "PAPER_1119"
title: "Lorentz Regauging and Vacuum Energy Extraction: Heaviside Component in UQFF"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [Lorentz-regauging, Heaviside, vacuum-energy, COP, Poynting, negentropy, TRZ, quasi-longitudinal, Bearden]
crosslinks: []
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
cp4_entry: 620
---

# Lorentz Regauging and Vacuum Energy Extraction

## Abstract

We integrate the Bearden (2000) Lorentz regauging formalism into the UQFF framework. The total electromagnetic energy flow comprises both the observed Poynting vector and the far larger Heaviside component:

$$S_{\text{total}} = S_{\text{Poynting}} + S_{\text{Heaviside}}$$

where $S_{\text{Heaviside}} = f_H \cdot 10^{13} \cdot S_{\text{Poynting}} \cdot (\rho_{\text{UA}} / \rho_{\text{SCm}})$. Breaking the Lorentz symmetric regauging condition (3-symmetry → 4-symmetry flow) enables extraction of vacuum energy with coefficient of performance $\text{COP} > 1.0$. The extracted power:

$$P_{\text{extracted}} = S_{\text{Heaviside}} \cdot A \cdot \eta_{\text{TRZ}}$$

is modulated by the TRZ extraction efficiency and the $\rho_{\text{UA}}/\rho_{\text{SCm}} = 10$ vacuum density ratio. Quasi-longitudinal waves in the Um (magnetism) field enable negentropy — thermodynamically open energy extraction from the structured vacuum.

## 1. Introduction

In standard electrodynamics, the Poynting vector $\mathbf{S} = \mathbf{E} \times \mathbf{B} / \mu_0$ represents the measurable energy flow. However, Heaviside (1893) identified an additional divergence-free component that is traditionally discarded by the Lorentz gauge condition. Bearden (2000) argued that this Heaviside component carries $\sim 10^{13}$ times more energy than the Poynting flow, but is rendered inaccessible by the symmetric regauging.

In UQFF, this energy resides in the $[\text{UA}]$ vacuum and is partially accessible through the $[\text{SCm}]$ superconductive condensate via TRZ (Triadic Resonant Zone) coupling.

## 2. Energy Flow Decomposition

### 2.1 Poynting Component

$$S_{\text{Poynting}} = \frac{E \times B}{\mu_0}$$

For $E = 10^6$ V/m and $B = 1.0$ T:

$$S_{\text{Poynting}} = \frac{10^6 \times 1.0}{1.2566 \times 10^{-6}} = 7.96 \times 10^{11}\ \text{W/m}^2$$

### 2.2 Heaviside Component

$$S_{\text{Heaviside}} = f_H \cdot 10^{13} \cdot S_P \cdot \frac{\rho_{\text{UA}}}{\rho_{\text{SCm}}}$$

| Parameter | Value | Description |
|-----------|-------|-------------|
| $f_H$ | 0.01 | Heaviside coupling fraction |
| $10^{13}$ | — | Heaviside/Poynting ratio (Bearden) |
| $\rho_{\text{UA}} / \rho_{\text{SCm}}$ | 10 | Vacuum density ratio |

$$S_{\text{Heaviside}} = 0.01 \times 10^{13} \times 7.96 \times 10^{11} \times 10 = 7.96 \times 10^{23}\ \text{W/m}^2$$

### 2.3 Coefficient of Performance

$$\text{COP} = \frac{S_{\text{Poynting}} + P_{\text{extracted}} / A}{S_{\text{Poynting}}}$$

COP > 1.0 indicates thermodynamically open operation — energy extracted from the structured vacuum exceeds the input.

## 3. Extraction Mechanism

### 3.1 TRZ Coupling

The Triadic Resonant Zone provides the extraction channel:

$$P_{\text{extracted}} = S_{\text{Heaviside}} \cdot A \cdot \eta_{\text{TRZ}}$$

For $A = 10^{-4}$ m² and $\eta_{\text{TRZ}} = 0.1$:

$$P_{\text{extracted}} = 7.96 \times 10^{23} \times 10^{-4} \times 0.1 = 7.96 \times 10^{18}\ \text{W}$$

### 3.2 Quasi-Longitudinal Waves

The Um (magnetism) field supports quasi-longitudinal wave modes:

$$P_{\text{quasi}} = P_{\text{extracted}} \cdot f_{\text{quasi}} \cdot \exp\!\left(-\frac{[\text{SSq}] \cdot 18}{26}\right)$$

These modes carry energy along field lines rather than transversely, enabling directed energy flow from the vacuum.

## 4. Lorentz Regauging Interpretation

### 4.1 Symmetry Breaking

Standard EM uses Lorentz gauge: $\nabla \cdot \mathbf{A} + \mu_0 \epsilon_0 \partial\Phi/\partial t = 0$. This enforces 3-symmetry (equal and opposite energy flows cancel). Breaking to 4-symmetry flow allows net extraction from the Heaviside component.

### 4.2 UQFF Connection

In UQFF, the vacuum is not empty but structured with $[\text{UA}]$ and $[\text{SCm}]$ densities. The Heaviside component represents the $[\text{UA}]$ vacuum energy that flows around but does not interact with conventional detectors. The $[\text{SCm}]$ condensate provides the symmetry-breaking mechanism.

## 5. Results

| Observable | Value | Description |
|-----------|-------|-------------|
| $S_{\text{Poynting}}$ | $7.96 \times 10^{11}$ W/m² | Measurable flow |
| $S_{\text{Heaviside}}$ | $7.96 \times 10^{23}$ W/m² | Vacuum flow |
| $\text{COP}$ | $\sim 10^{12}$ | Extraction ratio |
| $P_{\text{quasi}}$ | computed | Quasi-longitudinal power |

## 6. Conclusions

The Lorentz regauging formalism provides a pathway to vacuum energy extraction within the UQFF framework. The Heaviside component, modulated by the $\rho_{\text{UA}}/\rho_{\text{SCm}}$ ratio, represents the dominant energy flow in any electromagnetic system. CP4 class `LorentzRegaugingVacuumEnergyCalculator` (#620) implements E-field, B-field, efficiency, and area sweeps.

## References

1. Bearden, T.E., "Energy from the Vacuum" (2000)
2. Heaviside, O., "Electromagnetic Theory" (1893)
3. Murphy, D.T., Star-Magic UQFF Framework (2024–2026)
