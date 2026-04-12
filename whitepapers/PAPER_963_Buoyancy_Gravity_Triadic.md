# PAPER_963: Buoyancy Gravity Triadic Mode

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 215
**Source:** triadic_solutions_next.py (BuoyancyGravityTriadic)
**Calculator:** BuoyancyGravityTriadicCalc (CP4 #547)
**CVW:** v2.0.0 compliant

---

## Abstract

The Buoyancy Gravity Triadic mode evaluates $E_\text{net}(t, \Gamma) = S_{26} \cdot \cos(\omega_\text{SCm} t) \cdot \exp(-\Gamma t) - \text{threshold}$. Positive $E_\text{net}$ drives expansion (nebulae, HII regions); negative $E_\text{net}$ drives erosion (filaments, cometary knots).

---

## 1. Net Energy

$$E_\text{net}(t, \Gamma) = S_{26} \cdot \cos(\omega_\text{SCm} \cdot t) \cdot \exp(-\Gamma \cdot t) - \text{threshold}$$

## 2. Regime Classification

| $E_\text{net}$ | Regime | Astrophysical Example |
|-----------------|--------|----------------------|
| > 0.1 | Expansion | Nebulae, HII regions |
| < -0.1 | Erosion | Filaments, pillars |
| [-0.1, 0.1] | Neutral | Transition zones |

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
