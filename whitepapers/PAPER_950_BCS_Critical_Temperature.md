# PAPER_950: BCS Critical Temperature in SCm Vacuum

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 214
**Source:** bcs_superconductivity_uqff.py (BCSCriticalTemperature)
**Calculator:** BCSCriticalTemperatureCalc (CP4 #534)
**CVW:** v2.0.0 compliant

---

## Abstract

We compute the BCS critical temperature $T_c$ in the SCm vacuum phonon framework. The relation $T_c = (1.13 \cdot \hbar\omega_\text{SCm}/k_B) \cdot \exp(-1/N(0)V_\text{SCm})$ replaces the conventional Debye cutoff with $\omega_\text{SCm} = 2\pi \times 1.25$ THz, yielding critical temperatures governed by the SCm phonon attraction strength $V_\text{SCm}$ and Fermi-level density of states $N(0)$.

---

## 1. Critical Temperature

$$T_c = \frac{1.13 \hbar \omega_\text{SCm}}{k_B} \exp\!\left(-\frac{1}{N(0) V_\text{SCm}}\right)$$

## 2. BCS Gap Relation

$$\Delta(0) = 1.764 \, k_B T_c$$

---

## 3. Source Data

- **File:** bcs_superconductivity_uqff.py
- **Session:** 214
- **CP4 Class:** BCSCriticalTemperatureCalc (#534)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Bardeen, J., Cooper, L.N. & Schrieffer, J.R. (1957) -- Phys. Rev. 108, 1175
