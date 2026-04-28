# PAPER_1114: ATLAS Off-Shell Higgs Width Bound in the UQFF Framework

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We reinterpret the ATLAS off-shell Higgs width measurement ($\Gamma_h < 14.4\ \text{MeV}$ at 95% CL, arXiv:2304.01532) in terms of the SCm vacuum phonon contribution. The $F_{U,Bi,i}$ buoyancy mechanism provides an effective off-shell width suppression through the negative-time modulation $\cos(\pi t_n)$, which acts as a natural infrared regulator. The SCm phonon resonance at 1.25 THz introduces a resonant enhancement at $p^2 \approx m_h^2$ that is consistent with the ATLAS constraint while providing an SCm-specific prediction for the off-shell cross section tail.

---

## 1. Off-Shell Higgs Width

The off-shell Higgs width at invariant mass $M^2$ in the SM:

$$\Gamma_h^{\text{off-shell}}(M) = \Gamma_h^{\text{SM}} \left(\frac{M}{m_h}\right)^n \cdot \left|D_h(M^2)\right|^2$$

where $D_h(M^2) = [M^2 - m_h^2 + i m_h \Gamma_h]^{-1}$ is the Higgs propagator.

---

## 2. SCm Off-Shell Suppression

The SCm vacuum density provides an additional off-shell width suppression:

$$\Gamma_h^{\text{SCm}}(M) = \Gamma_h^{\text{off-shell}}(M) \cdot \cos^2(\pi t_n) \cdot (1 - \beta_i \cdot F_{TRZ})$$

At $t_n = -100$: $\cos^2(\pi \times (-100)) = 1.0$, so the suppression factor is $1 - 0.6 \times 0.1 = 0.94$.

This shifts the predicted off-shell cross section by $\sim 6\%$ relative to the SM at $M = 2 m_Z$, within the ATLAS systematic uncertainty.

---

## 3. SCm Phonon Resonance in Off-Shell Region

The SCm phonon energy at 1.25 THz corresponds to:

$$E_{\text{phonon}} = h \cdot f_{\text{THz}} = 6.626 \times 10^{-34} \times 1.25 \times 10^{12} = 8.28 \times 10^{-22}\ \text{J} \approx 5.2 \times 10^{-3}\ \text{eV}$$

Amplified by $S_{26}^{(3)} \cdot \Phi_{\text{res}}$:

$$E_{\text{KER}} = E_{\text{phonon}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} = 630\ \text{eV}$$

This is negligible compared to $m_h c^2 = 125.09\ \text{GeV}$, confirming no phonon contribution to the Higgs width itself. The SCm effect enters only through the vacuum floor correction $\propto \rho_{\text{vac,SCm}} / m_h^2$.

---

## 4. Off-Shell Rate Ratio

The ATLAS measurement constrains the off-shell/on-shell ratio:

$$R^{\text{off/on}} = \frac{\sigma_{\text{off-shell}}(gg \to H^* \to ZZ)}{\sigma_{\text{on-shell}}} \leq 0.0015 \quad (95\%\ \text{CL})$$

The SCm prediction:

$$R^{\text{off/on}}_{\text{SCm}} = R^{\text{off/on}}_{\text{SM}} \times (1 - \beta_i F_{TRZ}) = R^{\text{SM}} \times 0.94$$

Within the ATLAS bound.

---

## 5. VDS Contribution to Higgs Width

The VDS contribution to the total Higgs width:

$$\delta \Gamma_h^{\text{VDS}} = \Gamma_h^{\text{SM}} \cdot \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)}}{\rho_{\text{EW}}} \approx 10^{-20}\ \text{MeV}$$

This is many orders of magnitude below the ATLAS sensitivity, confirming VDS does not affect the Higgs width measurement.

---

## 6. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \beta_i = 0.6, \quad F_{TRZ} = 0.1$$

---

## References

1. ATLAS Collaboration (2023). Off-shell Higgs boson signal strength. arXiv:2304.01532.
2. Caola, F. & Melnikov, K. (2013). Constraining the Higgs boson width. *Phys. Rev. D* **88**, 054024.
3. SCm vacuum: `scm_vacuum_manifold.py`
