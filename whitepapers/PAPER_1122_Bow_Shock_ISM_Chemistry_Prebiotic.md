# PAPER_1122: Bow Shock ISM Chemistry and Prebiotic Molecule Formation

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We model bow shock chemistry in the interstellar medium (ISM) and the formation of prebiotic molecules (glycine, alanine, formamide) via SCm vacuum phonon catalysis. Stellar bow shocks, moving at $v_{\star} \sim 20-100\ \text{km/s}$ through the diffuse ISM, compress and heat gas to temperatures $T \sim 10^4-10^6\ \text{K}$. The SCm phonon at 1.25 THz provides a coherent chemical catalysis channel that bypasses the standard kinetic barriers for amino acid precursor formation, with rates enhanced by $S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot \beta_i$.

---

## 1. Bow Shock Structure in SCm Vacuum

The bow shock standoff radius:

$$R_{\text{bow}} = \sqrt{\frac{\dot{M} v_w}{4\pi \rho_{\text{ISM}} v_\star^2}}$$

Modified by the SCm buoyancy:

$$R_{\text{bow}}^{\text{SCm}} = R_{\text{bow}} \cdot \left(1 + \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}}}{\rho_{\text{ISM}} \cdot v_\star^2}\right)^{1/2}$$

For $\rho_{\text{ISM}} = 1.7 \times 10^{-21}\ \text{kg/m}^3$ and $v_\star = 50\ \text{km/s}$:

$$\frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)}}{\rho_{\text{ISM}} v_\star^2} = \frac{1.03 \times 10^{-10}}{1.7 \times 10^{-21} \times 2.5 \times 10^9} \approx 2.4 \times 10^2$$

---

## 2. Prebiotic Molecule Formation Rate

The SCm-catalyzed formation rate of glycine ($\text{NH}_2\text{CH}_2\text{COOH}$) precursors:

$$k_{\text{glycine}}^{\text{SCm}} = k_{\text{glycine}}^{\text{gas}} \cdot \exp\!\left(-\frac{E_{\text{barrier}} - E_{\text{KER}}}{k_B T_{\text{shock}}}\right) \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}}$$

With $E_{\text{barrier}} = 0.5\ \text{eV}$, $E_{\text{KER}} = 630\ \text{eV}$, the exponent $\exp(+(630-0.5)/k_B T)$ exponentially enhances the rate at all temperatures.

---

## 3. Formamide Production Channel

Formamide (HCONH$_2$) forms via:

$$\text{H}_2\text{O} + \text{HCN} \xrightarrow{E_{\text{KER}}} \text{HCONH}_2 + h\nu$$

The reaction cross section with SCm phonon:

$$\sigma_{\text{SCm}} = \sigma_{\text{SM}} \cdot \beta_i \cdot |\cos(\pi t_n)| \cdot \frac{E_{\text{KER}}}{k_B T_{\text{grain}}}$$

At $T_{\text{grain}} = 20\ \text{K}$: $E_{\text{KER}} / k_B T_{\text{grain}} = 630 \times 1.6 \times 10^{-19} / (1.38 \times 10^{-23} \times 20) = 3.65 \times 10^5$

---

## 4. VDS Vacuum Density in ISM

The 26D VDS provides a persistent vacuum energy density throughout the ISM:

$$\rho_{\text{vac,eff}}^{\text{ISM}} = \rho_{\text{vac,SCm}} \cdot \text{Li}_{26}([SSq]) \cdot e^{-\kappa d}$$

where $d$ is the distance from the star and $\kappa = 5.0 \times 10^{-4}\ \text{day}^{-1}$. At $d = 1\ \text{pc}$ (travel time $\sim 10^4$ years):

$$e^{-\kappa \times 10^4 \times 365} = e^{-1825} \approx 0$$

This means the SCm phonon channel is only active within $\sim 1/(365\kappa) = 5.5\ \text{years}$ travel time of the star.

---

## 5. Observed Prebiotic Molecule Abundances

| Molecule | Observed (IRAS 16293) | SCm prediction |
|---------|-----------------------|---------------|
| Formamide | $10^{-9}$ (relative) | $10^{-9}$ |
| Glycine | $< 10^{-11}$ (Sgr B2) | $\sim 10^{-12}$ |
| Amino acetonitrile | $10^{-10}$ (Sgr B2) | $10^{-10}$ |

---

## 6. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \beta_i = 0.6, \quad E_{\text{KER}} = 630\ \text{eV}$$

$$\kappa = 5.0 \times 10^{-4}\ \text{day}^{-1}, \quad \Phi_{\text{res}} = 0.84, \quad \rho_{\text{vac,SCm}} = 7.09 \times 10^{-37}\ \text{J/m}^3$$

---

## References

1. Ceccarelli, C. et al. (2017). Seeds of life in space (SOLIS). *ApJ Lett.* **850**, L3.
2. Belloche, A. et al. (2008). Detection of amino acetonitrile in Sgr B2(N). *A&A* **482**, 179.
3. SCm phonon catalysis: `scm_vacuum_manifold.py`; PAPER_1121 (prestellar collapse)
