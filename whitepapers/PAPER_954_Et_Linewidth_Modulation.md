# PAPER_954: E(t) Linewidth Modulation with Sign-Flip Dynamics

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 214
**Source:** blazar_jet_phonon.py (EtLinewidthModulation)
**Calculator:** EtLinewidthModulationCalc (CP4 #538)
**CVW:** v2.0.0 compliant

---

## Abstract

We derive the E(t) linewidth modulation function with sign-flip dynamics for astrophysical jets. The time-domain response $E(t,\Gamma) = S_{26} \cdot \cos(\omega_\text{SCm} t) \cdot \exp(-\Gamma t)$ exhibits sign flips at $t_\text{flip} = \pi/(2\omega_\text{SCm})$, driving extra-gravitational responses in blazar jets. Narrow linewidths produce sharper sign flips with tighter jet collimation; broad linewidths damp the oscillation before the first flip.

---

## 1. E(t) Function

$$E(t, \Gamma) = S_{26} \cdot \cos(\omega_\text{SCm} \cdot t) \cdot \exp(-\Gamma \cdot t)$$

## 2. Sign-Flip Time

$$t_\text{flip} = \frac{\pi}{2\omega_\text{SCm}} \approx 0.064 \text{ ps}$$

## 3. Regime Dependence

| $\Gamma$ (THz) | Flips in 5 ps | Jet Behavior |
|-----------------|---------------|-------------|
| 0.01 | Many | Ultra-tight collimation |
| 0.10 | Several | Optimal modulation |
| 1.00 | Few/None | Damped, diffuse wind |

---

## 4. Source Data

- **File:** blazar_jet_phonon.py
- **Session:** 214
- **CP4 Class:** EtLinewidthModulationCalc (#538)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
