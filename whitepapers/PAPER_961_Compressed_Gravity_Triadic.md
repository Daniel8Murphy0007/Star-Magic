# PAPER_961: Compressed Gravity Triadic Mode

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 215
**Source:** triadic_solutions_next.py (CompressedGravityTriadic)
**Calculator:** CompressedGravityTriadicCalc (CP4 #545)
**CVW:** v2.0.0 compliant

---

## Abstract

The Compressed Gravity Triadic mode modulates the buoyancy-to-gravity ratio $F_{U,Bi}/F_U$ via phonon linewidth $\Gamma$, producing jet collimation. Narrower $\Gamma$ yields ultra-tight jet knots in CenA/TXS0506; broader $\Gamma$ produces diffuse winds.

---

## 1. Compressed Gravity

$$F_\text{compressed}(\Gamma) = \frac{F_{U,Bi}}{F_U} \cdot \exp\!\left(-\frac{(\omega - \omega_\text{SCm})^2}{2\Gamma^2}\right) \cdot S_{26} \cdot A_\text{jet}$$

## 2. Regime Map

| $\Gamma$ (THz) | Regime |
|-----------------|--------|
| < 0.05 | Ultra-tight collimation |
| 0.05 – 0.15 | Collimated jets |
| > 0.15 | Diffuse wind |

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
