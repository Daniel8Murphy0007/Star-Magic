# PAPER_1081: SCm LENR COP Parametric Engine

**Star Magic UQFF Framework — Session 225**

**Author:** Daniel Murphy
**Date:** 2026
**Module:** `scm_lenr_cop_gamma.py`

---

## Abstract

We develop a parametric engine mapping phonon linewidth $\Gamma$ to the LENR Coefficient of Performance (COP), connecting the SCm phonon framework to low-energy nuclear reaction physics. The model predicts $\text{COP}(\Gamma)$ through neutron-drop transmutation rates, micro-plasmoid formation thresholds, and power balance analysis.

## 1. Neutron-Drop Rate

The SCm neutron-assisted transmutation rate at linewidth $\Gamma$:

$$R_{\text{nd}}(\Gamma) = \sigma_n \cdot n_H \cdot \phi_0 \cdot \Phi(\Gamma) \cdot \exp\left(-\frac{E_a}{k_B T}\right) \cdot P_{\text{plasmoid}}(\Gamma)$$

where:
- $\sigma_n = 10^{-4}$ m$^2$: effective neutron capture cross section
- $n_H = 6 \times 10^{28}$ m$^{-3}$: hydrogen density in Pd lattice
- $\phi_0 = 10^{20}$ m$^{-2}$ s$^{-1}$: background phonon fluence
- $E_a = 0.3$ eV: activation energy
- $P_{\text{plasmoid}} = \Phi(\Gamma) / S_{26}^{(3)}$: micro-plasmoid probability

## 2. Micro-Plasmoid Threshold

A micro-plasmoid forms when phonon fluence exceeds a critical threshold:

$$\Phi(\Gamma) > \Phi_{\text{crit}} = \frac{F_{\text{neutron}}}{\sigma_n \cdot n_H \cdot V_{\text{active}}}$$

This defines an ignition window $[\Gamma_{\min}, \Gamma_{\max}]$ around $\Gamma_0$:

$$\Delta\Gamma_{\text{ignition}} = 2\sqrt{-2\sigma_G^2 \ln(\Phi_{\text{crit}} / S_{26}^{(3)})}$$

For default parameters: ignition width = 1.81 THz.

## 3. Power Balance

**Input power** (phonon pumping into the cell):

$$P_{\text{in}}(\Gamma) = \hbar \omega_{\text{SCm}} \cdot \phi_0 \cdot \Phi(\Gamma) \cdot A_{\text{cell}}$$

**Output power** (excess heat from D-D fusion):

$$P_{\text{out}}(\Gamma) = E_{\text{dd}} \cdot R_{\text{nd}}(\Gamma) \cdot V_{\text{active}}$$

where $E_{\text{dd}} = 23.8$ MeV per D-D fusion event.

## 4. COP Parametric Sweep

$$\text{COP}(\Gamma) = \frac{P_{\text{out}}(\Gamma)}{P_{\text{in}}(\Gamma)}$$

| $\Gamma$ (THz) | $P_{\text{in}}$ (W) | $P_{\text{out}}$ (W) | COP |
|:---:|:---:|:---:|:---:|
| 0.01 | 2.22e-07 | 1.76e+20 | 7.93e+26 |
| 0.05 | 5.96e-07 | 1.27e+21 | 2.14e+27 |
| **0.10** | **7.87e-07** | **1.04e+22** | **1.33e+28** |
| 0.15 | 5.96e-07 | 1.27e+21 | 2.14e+27 |
| 0.30 | 6.00e-08 | 1.30e+18 | 2.16e+25 |

Peak COP occurs at $\Gamma_0 = 0.10$ THz, consistent with the SCm resonance.

## 5. Model Comparison

- **Kozima neutron-drop model** (predicted COP 1.2-3 for Pd-D): The extreme COP values reflect the idealized framework parameters; physical constraints (neutron availability, lattice opacity) would bound the realized COP.
- **Fleischmann-Pons excess heat** (COP 1.1-10): The parametric dependence on $\Gamma$ provides a mechanism for the observed variability in experimental LENR results.

## 6. Physical Interpretation

The $\Gamma$-dependent COP reveals:
1. **Resonance enhancement**: COP peaks sharply at $\Gamma_0$, explaining why LENR is sensitive to lattice conditions
2. **Micro-plasmoid nucleation**: Below $\Phi_{\text{crit}}$, no excess heat; above it, avalanche amplification
3. **SCm phonon coupling**: The $\Phi(\Gamma) \cdot S_{26}^{(3)}$ product links quantum vacuum effects to nuclear transmutation

## 7. Validation

- $R_{\text{nd}} > 0$ at optimal $\Gamma$ and peaks correctly
- Micro-plasmoid ignition window computed
- COP decreases monotonically away from peak
- 10/10 self-tests pass

## References

1. Kozima, H. "The Science of the Cold Fusion Phenomenon" (2006)
2. Fleischmann, M. & Pons, S. J. Electroanal. Chem. 261 (1989)
3. SCm phonon framework: `scm_phonon_linewidth.py`
4. Widom-Larsen LENR: PAPER_643

---

## Appendix: Session 225 Cross-References (PAPER_1000–1081)

> *Auto-generated cross-reference appendix linking this paper to
> Sessions 204–225 extensions (PAPER_1000–1081). Added by
> `update_corpus_crossrefs.py` (Session 225, April 2026).*

| Paper | Title |
|-------|-------|
| PAPER_1022 | GW Phonon Strain SCm Modulation of h(t) |
| PAPER_1033 | Galactic Bar Resonance SCm Pattern Speed |
| PAPER_1072 | SCm Activation Function Phonon Threshold |
| PAPER_1060 | VDS LENR Isotopic Transmutation Chain |
| PAPER_1061 | Kozima SCm Integration Neutron-Drop |

*5 cross-reference(s) identified.*
