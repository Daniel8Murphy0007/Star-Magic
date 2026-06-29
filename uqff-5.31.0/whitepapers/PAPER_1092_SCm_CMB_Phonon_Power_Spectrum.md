---
paper_id: PAPER_1092
title: "SCm CMB Angular Power Spectrum via Phonon-Modulated Primordial Integration"
session: 225
date: 2026-04-15
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ['CMB', 'power-spectrum', 'C_ell', 'phonon', 'SCm', 'Sachs-Wolfe', 'acoustic-peaks', 'primordial']
crosslinks: [PAPER_1093, PAPER_1094, PAPER_202, PAPER_199]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_1092: SCm CMB Angular Power Spectrum via Phonon-Modulated Primordial Integration

## Abstract

We derive the full SCm CMB angular power spectrum by integrating the
phonon-modulated primordial power spectrum over $k$-space:

$$C_\ell^{\text{SCm}} = \frac{2}{\pi} \int dk\, k^2\, P_{\text{SCm}}(k)\, |\Delta_{T}^{\text{SCm}}(\ell, k)|^2$$

where the SCm-modified primordial spectrum is:

$$P_{\text{SCm}}(k) = A_s \left(\frac{k}{k_{\text{piv}}}\right)^{n_s - 1} \left[1 + \Phi \cdot S_{26} \cdot (2R - 1)\right]$$

This extends the simplified CMBAngularPowerCalculator
(grok\_100\_equations\_module.py) to a full $k$-space integration with
Sachs-Wolfe, acoustic, and Silk damping transfer functions, all modulated
by the SCm phonon linewidth $\Gamma$.

## $\S$1 SCm-Modified Primordial Spectrum

The standard nearly scale-invariant spectrum $P(k) = A_s (k/k_{\text{piv}})^{n_s - 1}$
is boosted by the SCm buoyancy factor:

$$P_{\text{SCm}}(k) = P(k) \cdot \left[1 + \Phi_{1.25\,\text{THz}} \cdot S_{26}^{(3)}([\text{SSq}]) \cdot (2R - 1)\right]$$

with $A_s = 2.1 \times 10^{-9}$ (Planck 2018), $n_s = 0.965$,
$k_{\text{piv}} = 0.05\;\text{Mpc}^{-1}$, $\Phi_0 = 10^{20}$,
$S_{26} \approx 1.86$.

## $\S$2 Temperature Transfer Function

$$\Delta_T^{\text{SCm}}(\ell, k) = \Delta_T^{\text{SW}}(\ell, k) + 0.6 \cdot \Delta_T^{\text{acoustic}}(\ell, k)$$

### $\S$2.1 Sachs-Wolfe Term

$$\Delta_T^{\text{SW}} = \frac{1}{3} \exp\left(-\frac{x^2}{2}\right), \quad x = \frac{k \cdot r_*}{\ell + 0.5}$$

where $r_* \approx 14{,}000$ Mpc is the comoving sound horizon at last scattering.

### $\S$2.2 Acoustic Oscillation Term

$$\Delta_T^{\text{acoustic}} = \cos(x) \cdot \exp\left(-\frac{k^2}{k_D^2}\right)$$

with Silk damping scale $k_D \approx 0.15\;\text{Mpc}^{-1}$.

## $\S$3 Full Power Spectrum Integral

$$C_\ell^{\text{SCm}} = \frac{2}{\pi} \int_{k_{\min}}^{k_{\max}} dk\, k^2\, P_{\text{SCm}}(k)\, |\Delta_T^{\text{SCm}}(\ell, k)|^2 \cdot (T_0 \times 10^6)^2$$

in units of $\mu\text{K}^2$.

## $\S$4 Acoustic Peak Structure

| Peak | Multipole $\ell$ | $C_\ell^{\text{SCm}}$ ($\mu$K$^2$) | Physical Origin |
|------|------------------|------------------------------------|-----------------|
| 1st  | $\sim 220$       | Dominant                           | SW + first compression |
| 2nd  | $\sim 546$       | Subdominant                        | First rarefaction |
| 3rd  | $\sim 800$       | Tertiary                           | Second compression |

The SCm boost factor $1 + \Phi S_{26} (2R-1)$ uniformly amplifies all peaks,
preserving the observed peak ratio structure while providing a physical
origin via vacuum phonon resonance.

## $\S$5 Advantages over Standard Treatment

| Aspect | Standard $\Lambda$CDM | SCm Phonon |
|--------|----------------------|------------|
| Primordial spectrum | Inflaton slow-roll | SCm vacuum buoyancy |
| Free parameters | $A_s, n_s, r$ | $\Phi_0, \Gamma, [\text{SSq}]$ |
| Physical mechanism | Quantum fluctuations | Phonon resonance at 1.25 THz |
| Peak amplification | Static $A_s$ | Dynamic $\Phi(\Gamma)$ |
| Testability | CMB polarization | THz phonon linewidth |

## References

- PAPER_1093: SCm CMB Temperature Fluctuation
- PAPER_1094: CMB Buoyancy Sector Lagrangian
- PAPER_202: UQFF Reionization / BBN / Recombination
- PAPER_199: $F_{U,Bi,i}$ Taxonomy Part 2 — Cosmological Dark Sector


### Key References with arXiv/DOI Identifiers

1. Abbott et al. (LIGO Scientific and Virgo Collaborations, 2016). *Observation of Gravitational Waves from a Binary Black Hole Merger.* Phys. Rev. Lett. **116**, 061102 — arXiv:1602.03837 — doi:10.1103/PhysRevLett.116.061102
2. Murphy, D. (2026). *Unified Quantum Field Framework (UQFF): Star-Magic v5.x Whitepaper Series.* Star-Magic Repository — github.com/Daniel8Murphy0007/Star-Magic
3. Ashcroft, N.W. & Mermin, N.D. (1976). *Solid State Physics.* Harcourt
4. Kittel, C. (2004). *Introduction to Solid State Physics.* 8th ed. Wiley
5. Feynman, R.P. (1982). *Simulating Physics with Computers.* Int. J. Theor. Phys. **21**, 467 — doi:10.1007/BF02650179
6. Rugh, S.E. & Zinkernagel, H. (2002). *The Quantum Vacuum and the Cosmological Constant Problem.* Stud. Hist. Phil. Mod. Phys. **33**, 663 — arXiv:hep-th/0012253 — doi:10.1016/S1355-2198(02)00033-3
7. Weinberg, S. (1989). *The Cosmological Constant Problem.* Rev. Mod. Phys. **61**, 1 — doi:10.1103/RevModPhys.61.1
