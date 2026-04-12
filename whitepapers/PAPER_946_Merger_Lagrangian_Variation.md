# PAPER_946: Merger Lagrangian Variation

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** smbh_binary_mergers.py (MergerLagrangianVariation)
**Calculator:** MergerLagrangianVariationCalc (CP4 #530)
**CVW:** v2.0.0 compliant

---

## Abstract

We derive the UQFF Euler-Lagrange equation for SMBH binary merger phonon fields. The stationarity condition $\delta S / \delta\varphi_\text{merger} = 0$ yields a critical radius $r_\text{crit}$ at which the gravitational Lagrangian $\mathcal{L}_\text{grav} = -\beta_i \sum_i U_{g,i} \cdot M/d_g \cdot [\text{UA}]$ balances the phonon Lagrangian $\mathcal{L}_\text{phonon} = F_n \cdot \Phi$. The 26-channel buoyancy sum governs the transition.

---

## 1. Merger Lagrangian

$$\mathcal{L}_\text{merger} = \mathcal{L}_\text{grav} + \mathcal{L}_\text{phonon}$$

$$\mathcal{L}_\text{grav} = -\beta_i \sum_{i=1}^{26} U_{g,i} \cdot \frac{M}{d_g} \cdot [\text{UA}]_i$$

$$\mathcal{L}_\text{phonon} = F_n \cdot \Phi$$

---

## 2. Stationarity Condition

$$\frac{\delta S}{\delta \varphi_\text{merger}} = \frac{\partial}{\partial E_\text{net}}\left(-\beta_i \sum U_{g,i} \cdot \Omega_g \cdot \frac{M}{d_g} \cdot [\text{UA}] + F_n \cdot \Phi\right) = 0$$

This yields the critical radius:

$$r_\text{crit} = \frac{2\beta_i \sum U_{g,i} \cdot M}{|F_n \cdot \Phi|}$$

---

## 3. Physical Interpretation

Below $r_\text{crit}$, gravitational effects dominate, and the merger phonon field acts as a passive tracer. Above $r_\text{crit}$, phonon back-reaction becomes significant and modifies the inspiral dynamics through buoyancy feedback.

---

## 4. Source Data

- **File:** smbh_binary_mergers.py
- **Session:** 213
- **CP4 Class:** MergerLagrangianVariationCalc (#530)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Arnowitt, R., Deser, S. & Misner, C.W. (1962) -- Gravitation (Wiley)
