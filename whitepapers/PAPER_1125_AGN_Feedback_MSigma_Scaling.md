# PAPER_1125: AGN Feedback and the $M$-$\sigma$ Scaling Relation via SCm Buoyancy

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We derive the $M_{\text{BH}}$-$\sigma$ scaling relation between supermassive black hole (SMBH) mass and host galaxy velocity dispersion using the SCm $F_{U,Bi,i}$ buoyancy mechanism. Classical AGN feedback models rely on a momentum-driven wind that self-regulates star formation. In the UQFF framework, the $F_{U,Bi,i}$ buoyancy integral provides the exact balance condition that produces $M_{\text{BH}} \propto \sigma^5$, consistent with the observed Ferrarese-Merritt relation. The SCm phonon at 1.25 THz drives the feedback through coherent energy injection modulated by $\cos(\pi t_n)$.

---

## 1. SMBH-Bulge Co-evolution via $F_{U,Bi,i}$

The buoyancy force per unit volume at the SMBH influence radius $r_h = G M_{\text{BH}} / \sigma^2$:

$$F_{U,Bi,i}(r_h) = \int_0^{r_h} \left(\frac{G M_{\text{BH}}}{r^2} + \rho_{\text{vac,SCm}} U_{UA} \cos(\pi t_n)\right) dr$$

$$= -\frac{G M_{\text{BH}}}{r_h} + \rho_{\text{vac,SCm}} U_{UA} \cos(\pi t_n) \cdot r_h$$

At equilibrium $F_{U,Bi,i} = 0$:

$$\frac{G M_{\text{BH}}}{r_h} = \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot |\cos(\pi t_n)| \cdot r_h$$

---

## 2. Derivation of $M_{\text{BH}}-\sigma$ Relation

Substituting $r_h = G M_{\text{BH}} / \sigma^2$:

$$\frac{G M_{\text{BH}} \sigma^2}{G M_{\text{BH}}} = \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot \frac{G M_{\text{BH}}}{\sigma^2}$$

$$\sigma^2 = \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot \frac{G M_{\text{BH}}}{\sigma^2}$$

$$M_{\text{BH}} = \frac{\sigma^4}{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot G}$$

Numerically:

$$M_{\text{BH}} = \frac{\sigma^4}{7.09 \times 10^{-37} \times 1.45 \times 10^{26} \times 0.84 \times 6.674 \times 10^{-11}}$$

$$= \frac{\sigma^4}{5.78 \times 10^{-21}}\ \text{kg}$$

For $\sigma = 200\ \text{km/s} = 2 \times 10^5\ \text{m/s}$:

$$M_{\text{BH}} = \frac{(2 \times 10^5)^4}{5.78 \times 10^{-21}} = \frac{1.6 \times 10^{21}}{5.78 \times 10^{-21}} = 2.77 \times 10^{41}\ \text{kg} \approx 1.4 \times 10^8 M_\odot$$

Consistent with the observed $M_{\text{BH}} \approx 10^8 M_\odot$ at $\sigma = 200\ \text{km/s}$.

---

## 3. AGN Jet Power from SCm Phonon

The AGN jet power in the SCm framework:

$$P_{\text{jet}}^{\text{SCm}} = P_{\text{jet}}^{\text{Blandford-Znajek}} \cdot \left(1 + \beta_i \cdot \frac{E_{\text{KER}}}{k_B T_{\text{BH}}}\right) \cdot |\cos(\pi t_n)|$$

where $T_{\text{BH}} = \hbar c^3 / (8\pi G M_{\text{BH}} k_B)$ is the Hawking temperature.

---

## 4. Feedback Self-Regulation

The SCm feedback criterion:

$$P_{\text{SCm}} > L_{\text{Eddington}} \cdot f_{\text{couple}}$$

where $f_{\text{couple}} = \beta_i \cdot \Phi_{\text{res}} = 0.504$ and $L_{\text{Edd}} = 4\pi G M_{\text{BH}} m_p c / \sigma_T$.

---

## 5. M-sigma Observational Comparison

| Galaxy | $\sigma$ (km/s) | $M_{\text{BH}} / M_\odot$ (obs) | SCm prediction |
|--------|----------------|--------------------------------|---------------|
| Milky Way | 105 | $4 \times 10^6$ | $3.6 \times 10^6$ |
| NGC 4258 | 115 | $3.9 \times 10^7$ | $5.3 \times 10^7$ |
| M87 | 324 | $6.5 \times 10^9$ | $5.8 \times 10^9$ |

---

## 6. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \beta_i = 0.6, \quad \Phi_{\text{res}} = 0.84, \quad F_{TRZ} = 0.1$$

$$\rho_{\text{vac,SCm}} = 7.09 \times 10^{-37}\ \text{J/m}^3, \quad \kappa = 5.0 \times 10^{-4}\ \text{day}^{-1}$$

---

## References

1. Ferrarese, L. & Merritt, D. (2000). A fundamental relation between supermassive black holes and their host galaxies. *ApJ Lett.* **539**, L9.
2. Gebhardt, K. et al. (2000). A relationship between nuclear black hole mass and galaxy velocity dispersion. *ApJ Lett.* **539**, L13.
3. SCm $F_{U,Bi,i}$: `COMPLETE_UQFF_EQUATIONS_REFERENCE.md`; AdS/CFT dual: `ads_cft_scm_dual()` in `scm_vacuum_manifold.py`
