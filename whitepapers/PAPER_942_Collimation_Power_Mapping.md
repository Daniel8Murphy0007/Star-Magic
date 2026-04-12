# PAPER_942: Collimation-Power Mapping

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** linewidth_jet_modulation.py (CollimationPowerMapping)
**Calculator:** CollimationPowerMappingCalc (CP4 #526)
**CVW:** v2.0.0 compliant

---

## Abstract

We derive the mapping between phonon linewidth $\Gamma$, jet collimation half-angle $\theta_\text{half}$, and brightness contrast for astrophysical jets. Narrow linewidths ($\Gamma \leq 0.07$ THz) produce tightly collimated jets with $\theta_\text{half} \lesssim 3^\circ$, while broad linewidths ($\Gamma > 0.15$ THz) yield wide opening angles $\theta_\text{half} \gtrsim 10^\circ$. The contrast (peak-to-background ratio) scales with $M_\text{jet}$.

---

## 1. Collimation Relation

$$\theta_\text{half} = \max\!\left(0.5^\circ,\; \frac{30^\circ}{Q}\right)$$

where $Q = \omega_\text{SCm} / (2\Gamma)$ is the quality factor.

---

## 2. Mapping Results

| $\Gamma$ (THz) | $Q$ | $\theta_\text{half}$ | Contrast |
|-----------------|-----|---------------------|----------|
| 0.05 | 12.5 | $2.4^\circ$ | High |
| 0.10 | 6.25 | $4.8^\circ$ | Moderate |
| 0.30 | 2.08 | $14.4^\circ$ | Low |

---

## 3. Observational Implications

The $Q$-$\theta_\text{half}$ relation provides a direct VLBI prediction: systems with narrower phonon linewidths should exhibit tighter jet collimation at mas scales, testable with EHT and ngEHT observations.

---

## 4. Source Data

- **File:** linewidth_jet_modulation.py
- **Session:** 213
- **CP4 Class:** CollimationPowerMappingCalc (#526)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Event Horizon Telescope Collaboration (2019) -- ApJL, 875, L1
