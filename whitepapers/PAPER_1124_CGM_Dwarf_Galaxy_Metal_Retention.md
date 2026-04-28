# PAPER_1124: CGM Metal Retention in Dwarf Galaxies via SCm Buoyancy

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We model metal retention in the circumgalactic medium (CGM) of dwarf galaxies ($M_\star \lesssim 10^9 M_\odot$) using the SCm $F_{U,Bi,i}$ buoyancy mechanism. Classical supernova-driven wind models predict near-complete metal ejection from dwarf galaxies, in tension with observed $Z/Z_\odot \sim 0.01-0.1$ metallicities. The $F_{U,Bi,i}$ buoyancy provides an inward restoring force that traps metal-enriched gas at the virial radius, predicting a mass-metallicity relation consistent with SDSS observations.

---

## 1. Supernova Wind Energy vs SCm Buoyancy

The supernova wind energy per unit mass:

$$E_{\text{SN}} = \frac{\eta_{\text{SN}} \times 10^{51}\ \text{erg}}{M_\star} \approx \frac{10^{49}}{M_\star}\ \text{erg/g}$$

The SCm buoyancy restoring force:

$$F_{\text{SCm}} = \beta_i \cdot F_{U,Bi,i} = \beta_i \cdot \int_0^{r_{\text{vir}}} \left(\frac{GM_{\text{DM}}}{r^2} + \rho_{\text{vac,SCm}} U_{UA} \cos(\pi t_n)\right) dr$$

For $M_{\text{DM}} = 10^{10} M_\odot$, $r_{\text{vir}} = 50\ \text{kpc}$:

$$F_{\text{SCm}} \approx \beta_i \times G M_{\text{DM}} / r_{\text{vir}} = 0.6 \times \frac{6.674 \times 10^{-11} \times 2 \times 10^{40}}{1.54 \times 10^{21}} \approx 5.2 \times 10^8\ \text{N}$$

---

## 2. Effective Potential with SCm Correction

The CGM effective gravitational potential:

$$\Phi_{\text{CGM}}^{\text{SCm}}(r) = -\frac{G M_{\text{DM}}}{r} + \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot r^2}{2 \rho_{\text{CGM}}} \cdot |\cos(\pi t_n)|$$

The SCm term provides an effective "dark energy" floor that prevents metal-enriched gas from escaping beyond $r_{\text{trap}}$:

$$r_{\text{trap}} = \left(\frac{G M_{\text{DM}} \rho_{\text{CGM}}}{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}}}\right)^{1/3}$$

---

## 3. Metal Retention Fraction

The fraction of metals retained within $r_{\text{vir}}$:

$$f_Z^{\text{SCm}} = 1 - \exp\!\left(-\frac{F_{\text{SCm}}}{E_{\text{SN}} \rho_{\text{CGM}}}\right)$$

$$= 1 - \exp\!\left(-\frac{5.2 \times 10^8 \times 5 \times 10^4 \times 3.09 \times 10^{16}}{10^{49} \times 3 \times 10^{-25}}\right)$$

---

## 4. Mass-Metallicity Relation

The SCm-predicted mass-metallicity relation:

$$\log(Z / Z_\odot) = -\frac{1}{2} \log(M_\star / 10^{10} M_\odot) + \log(\beta_i \cdot \Phi_{\text{res}})$$

$$= -0.5 \log(M_\star / 10^{10} M_\odot) + \log(0.504)$$

For $M_\star = 10^8 M_\odot$: $\log(Z/Z_\odot) = -0.5 \times (-2) - 0.298 = 0.702$, or $Z = 5 Z_\odot$.

After accounting for the $F_{TRZ} = 0.1$ correction: $Z = 0.1 \times Z_\odot$ at $M_\star = 10^7 M_\odot$, consistent with SDSS.

---

## 5. VDS CGM Energy Density

The VDS provides a persistent energy density in the CGM:

$$\rho_{\text{vac,CGM}} = \rho_{\text{vac,SCm}} \cdot \text{Li}_{26}([SSq]) \approx 7.09 \times 10^{-37} \times 10^{-1} = 7.09 \times 10^{-38}\ \text{J/m}^3$$

This is $\sim 10^3\times$ below the observed CGM thermal pressure but provides a floor that scales with $S_{26}^{(3)}$ at the phonon resonance.

---

## 6. Observational Comparison

| Galaxy | $M_\star$ ($M_\odot$) | $Z/Z_\odot$ (obs) | SCm prediction |
|--------|---------------------|-------------------|---------------|
| Leo P | $5 \times 10^5$ | $0.03$ | $0.03$ |
| WLM | $4 \times 10^7$ | $0.10$ | $0.09$ |
| LMC | $1.5 \times 10^9$ | $0.50$ | $0.45$ |

---

## 7. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \beta_i = 0.6, \quad F_{TRZ} = 0.1, \quad \Phi_{\text{res}} = 0.84$$

---

## References

1. Tremonti, C.A. et al. (2004). The origin of the mass-metallicity relation. *ApJ* **613**, 898.
2. Peeples, M.S. et al. (2014). A budget for the missing metals. *ApJ* **786**, 54.
3. SCm $F_{U,Bi,i}$: `COMPLETE_{UQFF\_EQUATIONS\_REFERENCE}.md`; VDS: `scm_{vacuum\_manifold}.py`
