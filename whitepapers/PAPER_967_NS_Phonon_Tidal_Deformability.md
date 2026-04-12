# PAPER_967: NS Phonon Tidal Deformability Correction

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 215
**Source:** et_phonon_resonance.py §8 (NSPhononGW190425.tidal_deformability_correction)
**Calculator:** NSPhononTidalDeformabilityCalc (CP4 #551)
**CVW:** v2.0.0 compliant

---

## Abstract

The SCm phonon correction to neutron star tidal deformability $\Lambda$ bridges the gap between GR-only predictions and LIGO/Virgo observations for GW190425. The correction $\delta\Lambda = F_{UBi}/F_U \cdot \Phi_{1.25\text{THz}} \cdot 0.1$ yields a 10% maximal shift in $\Lambda$, consistent with the mass-gap nature of the lighter component.

---

## 1. Tidal Deformability Correction

$$\Lambda_\text{UQFF} = \Lambda_\text{GR} \cdot (1 + \delta\Lambda_\text{phonon})$$

$$\delta\Lambda_\text{phonon} = \frac{F_{UBi}}{F_U} \cdot \Phi_{1.25\text{THz}}(\omega_\text{SCm}, \Gamma) \cdot 0.1$$

## 2. Phonon Occupation

$$\Phi_{1.25\text{THz}}(\omega, \Gamma) = \exp\!\left(-\frac{(\omega - \omega_\text{SCm})^2}{2\Gamma^2}\right) \cdot S_{26}$$

---

## References

1. LIGO/Virgo -- GW190425 (2020)
2. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
