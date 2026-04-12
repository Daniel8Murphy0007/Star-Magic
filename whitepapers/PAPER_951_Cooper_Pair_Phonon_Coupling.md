# PAPER_951: Cooper Pair Phonon Coupling

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 214
**Source:** bcs_superconductivity_uqff.py (CooperPairPhononCoupling)
**Calculator:** CooperPairPhononCouplingCalc (CP4 #535)
**CVW:** v2.0.0 compliant

---

## Abstract

We derive the effective Cooper-pair coupling strength via SCm phonon exchange at 1.25 THz. The coupling $V_\text{eff}(\omega,\Gamma) = V_\text{SCm} \cdot \Phi_{1.25\,\text{THz}}(\omega, \Gamma)$ is modulated by the phonon resonance profile $\Phi = \exp(-(\omega - \omega_\text{SCm})^2/(2\Gamma^2)) \cdot S_{26}$, yielding pair binding energy $E_\text{pair} = 2\Delta(T)$.

---

## 1. Coupling Strength

$$V_\text{eff}(\omega, \Gamma) = V_\text{SCm} \cdot \exp\!\left(-\frac{(\omega - \omega_\text{SCm})^2}{2\Gamma^2}\right) \cdot S_{26}$$

## 2. Pair Binding Energy

$$E_\text{pair} = 2\Delta(T)$$

At on-resonance ($\omega = \omega_\text{SCm}$), $\Phi = S_{26}$ and coupling is maximized.

---

## 3. Source Data

- **File:** bcs_superconductivity_uqff.py
- **Session:** 214
- **CP4 Class:** CooperPairPhononCouplingCalc (#535)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Cooper, L.N. (1956) -- Phys. Rev. 104, 1189
