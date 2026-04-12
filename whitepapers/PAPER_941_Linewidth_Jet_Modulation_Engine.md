# PAPER_941: Linewidth Jet Modulation Engine

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** linewidth_jet_modulation.py (LinewidthJetModulationSweep)
**Calculator:** LinewidthJetModulationSweepCalc (CP4 #525)
**CVW:** v2.0.0 compliant

---

## Abstract

We present a systematic linewidth-to-jet modulation mapping engine that sweeps $\Gamma$ from 0.01 to 1.0 THz and computes the jet modulation factor $M_\text{jet}$, quality factor $Q$, and operational regime (narrow / optimal / broad) at each point. The engine provides a universal lookup for any astrophysical jet system, decoupling the linewidth physics from specific source parameters.

---

## 1. Core Equations

$$M_\text{jet}(\Gamma) = 1 + A_\text{jet} \exp\!\left(-\frac{(\omega - \omega_\text{SCm})^2}{2\Gamma^2}\right) \cdot S_{26} \cdot \left(\frac{2F_{U\text{Bi}}}{F_U} - 1\right)$$

$$Q = \frac{\omega_\text{SCm}}{2\Gamma}$$

where $S_{26} = \sum_{k=1}^{26} e^{-[\text{SSq}] \cdot k/26}$ and $[\text{SSq}] = 0.57$.

---

## 2. Regime Classification

| Regime | $\Gamma$ Range | Characteristics |
|--------|---------------|-----------------|
| Narrow | $\Gamma \leq 0.07$ THz | High Q, tight collimation |
| Optimal | $0.07 < \Gamma \leq 0.15$ THz | Peak modulation |
| Broad | $\Gamma > 0.15$ THz | Low Q, wide opening angle |

---

## 3. Reference Systems

The `ReferenceSystemMatcher` class compares computed $(M_\text{jet}, Q)$ pairs against five calibrated systems: M87 ($A_\text{jet} = 1.20$), Sgr A* ($0.80$), Centaurus A ($0.95$), TXS 0506+056 ($1.20$), and 3C 273 ($1.05$).

---

## 4. Source Data

- **File:** linewidth_jet_modulation.py
- **Session:** 213
- **CP4 Class:** LinewidthJetModulationSweepCalc (#525)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Blandford, R.D. & Konigl, A. (1979) -- ApJ, 232, 34
