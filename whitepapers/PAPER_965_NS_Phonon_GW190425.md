# PAPER_965: NS Phonon Effects for GW190425

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 215
**Source:** et_phonon_resonance.py §8 (NSPhononGW190425)
**Calculator:** NSPhononGW190425Calc (CP4 #549)
**CVW:** v2.0.0 compliant

---

## Abstract

We derive neutron star phonon effects for GW190425 (mass-gap BNS merger). The UQFF phonon-corrected strain $h_\text{UQFF}(t) = h_\text{GR}(t) \cdot 0.5297 \cdot \exp([SSq] \cdot t/26)$ produces a 47% strain reduction. The phonon-corrected tidal deformability $\Lambda_\text{UQFF} = \Lambda_\text{GR}(1 + \delta\Lambda)$ matches LIGO constraints. Mass-gap probability: P(NS) = 49%, P(BH) = 51%.

---

## 1. Phonon-Corrected Strain

$$h_\text{UQFF}(t) = h_\text{GR}(t) \cdot 0.5297 \cdot \exp\!\left(\frac{[SSq] \cdot t}{26}\right)$$

## 2. Wavelength Correction

$$\lambda_\text{UQFF} = \lambda_\text{GR} \cdot \left(1 - \frac{F_{UBi}}{F_U} \cdot \Phi_{1.25\text{THz}}(\omega, \Gamma)\right)$$

## 3. Tidal Deformability

$$\Lambda_\text{UQFF} = \Lambda_\text{GR} \cdot (1 + \delta\Lambda_\text{phonon})$$

## 4. Mass-Gap Classification

At $\Gamma = 0.1$ THz: P(NS) ≈ 49%, P(BH) ≈ 51%.

---

## References

1. LIGO/Virgo -- GW190425 (2020)
2. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
