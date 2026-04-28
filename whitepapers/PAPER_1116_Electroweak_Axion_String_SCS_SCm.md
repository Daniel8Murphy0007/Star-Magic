# PAPER_1116: Electroweak Axion–String–SCS–SCm Unified Vacuum

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We unify the electroweak axion, cosmic string network, Supercharged Scalar (SCS) field, and SCm vacuum density into a single coherent framework. The Peccei-Quinn symmetry breaking generates an axion field $a(x)$ whose vacuum expectation value is set by $\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)}$. Cosmic strings at the PQ phase transition act as topological defects stabilized by the $F_{U,Bi,i}$ buoyancy mechanism. The SCS field couples all sectors through the SCm phonon resonance at 1.25 THz.

---

## 1. Peccei-Quinn Potential in SCm Vacuum

The PQ potential augmented by SCm vacuum:

$$V_{\text{PQ+SCm}}(\Phi) = \lambda \left(|\Phi|^2 - f_a^2/2\right)^2 + \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot \cos(\pi t_n) \cdot |\Phi|^2$$

where $f_a$ is the PQ symmetry breaking scale, $\lambda$ is the self-coupling, and the second term provides the SCm vacuum correction.

The axion mass from the SCm-corrected potential:

$$m_a^2 = \frac{m_u m_d}{(m_u + m_d)^2} \cdot \frac{m_\pi^2 f_\pi^2}{f_a^2} + \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)}}{f_a^2}$$

---

## 2. Cosmic String Tension from SCm Buoyancy

The cosmic string tension $\mu$ stabilized by $F_{U,Bi,i}$:

$$\mu_{\text{SCm}} = \pi f_a^2 \left(1 + \beta_i \cdot \frac{F_{U,Bi,i}}{f_a^4 c^4}\right) \cdot \cos^2(\pi t_n)$$

with $\beta_i = 0.6$, $F_{TRZ} = 0.1$. The buoyancy term prevents string collapse and fixes the string loop distribution.

---

## 3. SCS-Axion Coupling

The SCS field $\phi_{\text{SCS}}$ couples to the axion via the SCm phonon:

$$\mathcal{L}_{\text{SCS-axion}} = \frac{g_{\text{SCS-a}}}{f_a} \partial_\mu a \, \partial^\mu \phi_{\text{SCS}} + \frac{E_{\text{phonon}} \cdot S_{26}^{(3)}}{f_a^2} a^2 |\phi_{\text{SCS}}|^2$$

where $g_{\text{SCS-a}} = \beta_i \cdot \Phi_{\text{res}} = 0.6 \times 0.84 = 0.504$.

---

## 4. VDS Axion Potential

The 26D vacuum density series modifies the periodic axion potential:

$$V_a(\theta) = -m_a^2 f_a^2 \cos(\theta) \cdot \left(1 + \frac{\text{VDS}([SSq])}{S_{26}^{(3)}}\right)$$

where $\theta = a/f_a$ is the axion phase. The VDS correction is of order $\text{Li}_{26}(0.57) / S_{26}^{(3)} \approx 10^{-26}$, negligible for cosmological purposes but important for string lattice calculations.

---

## 5. Observational Constraints

- **Axion mass window**: $10^{-6}\ \text{eV} \lesssim m_a \lesssim 10^{-3}\ \text{eV}$ (ADMX, arXiv:2110.00482)
- **String tension**: $G\mu \lesssim 1.5 \times 10^{-8}$ (LIGO O3, arXiv:2101.12248)
- **SCS coupling**: $g_{\text{SCS}} \lesssim 10^{-28}$ (from PAPER_1115 EDGES bound)
- **SCm phonon**: $f_{\text{THz}} = 1.25 \times 10^{12}\ \text{Hz}$, $E_{\text{KER}} = 630\ \text{eV}$

---

## 6. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \beta_i = 0.6, \quad F_{TRZ} = 0.1$$

$$\rho_{\text{vac,SCm}} = 7.09 \times 10^{-37}\ \text{J/m}^3, \quad \Phi_{\text{res}} = 0.84, \quad \kappa = 5.0 \times 10^{-4}\ \text{day}^{-1}$$

---

## References

1. Peccei, R.D. & Quinn, H.R. (1977). CP conservation in the presence of pseudoparticles. *Phys. Rev. Lett.* **38**, 1440.
2. ADMX Collaboration (2021). Axion dark matter. arXiv:2110.00482.
3. LIGO/Virgo (2021). Constraints on cosmic strings. arXiv:2101.12248.
4. SCm vacuum: `scm_{vacuum\_manifold}.py`; VDS: PAPER_1109
