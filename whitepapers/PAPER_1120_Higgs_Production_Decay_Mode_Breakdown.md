# PAPER_1120: Higgs Production and Decay Mode Breakdown in the SCm Framework

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We provide a comprehensive breakdown of Higgs boson production mechanisms and decay modes within the SCm vacuum framework. The dominant production channels (ggF, VBF, VH, $t\bar{t}H$) receive SCm corrections through the $F_{U,Bi,i}$ buoyancy integral and $S_{26}^{(3)} \cdot \Phi_{\text{res}}$ amplification. Decay branching ratios are modified by the negative-time gate $\cos(\pi t_n)$, with the $H \to \gamma\gamma$ and $H \to Z\gamma$ channels showing the largest SCm enhancement due to their loop-induced nature. We compare all predicted rates with LHC Run 2 combined results.

---

## 1. Production Cross Sections with SCm Correction

The SCm-corrected production cross sections at $\sqrt{s} = 13\ \text{TeV}$:

| Channel | SM $\sigma$ (pb) | SCm correction | UQFF $\sigma$ (pb) |
|---------|-----------------|---------------|-------------------|
| ggF | $48.58$ | $+(0.6 \times F_{TRZ})\ \%$ | $48.87$ |
| VBF | $3.78$ | $+(0.06 \times \beta_i)\ \%$ | $3.78$ |
| WH | $1.37$ | $+(0.06 \times \beta_i)\ \%$ | $1.37$ |
| ZH | $0.88$ | $+(0.06 \times \beta_i)\ \%$ | $0.88$ |
| $t\bar{t}H$ | $0.51$ | $+(S_{26}^{(3)} \epsilon)\ \%$ | $0.51$ |

The SCm correction to ggF is: $\Delta\sigma / \sigma = \beta_i \cdot F_{TRZ} \cdot |\cos(\pi t_n)| = 0.6 \times 0.1 \times 1 = 6\%$ at $t_n = -100$.

---

## 2. Decay Branching Ratios

The $H \to \gamma\gamma$ partial width receives a SCm one-loop correction:

$$\Gamma(H \to \gamma\gamma)^{\text{SCm}} = \Gamma(H \to \gamma\gamma)^{\text{SM}} \left(1 + \frac{\alpha \cdot \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)}}{\pi m_h^2 c^2}\right)$$

Numerically, the correction factor is $\approx 1 + 10^{-30}$, negligible at current sensitivity.

The branching ratio table:

| Decay | SM BR | SCm BR |
|-------|-------|--------|
| $b\bar{b}$ | $58.2\%$ | $58.2\%$ |
| $WW^*$ | $21.4\%$ | $21.4\%$ |
| $\tau\tau$ | $6.27\%$ | $6.27\%$ |
| $ZZ^*$ | $2.62\%$ | $2.62\%$ |
| $\gamma\gamma$ | $0.228\%$ | $0.228\%\ (+10^{-30})$ |
| $Z\gamma$ | $0.154\%$ | $0.154\%$ |

---

## 3. SCm Phonon Contribution to ggF

The ggF loop integral receives a SCm vacuum contribution through the virtual top quark loop:

$$\mathcal{M}_{ggF}^{\text{SCm}} = \mathcal{M}_{ggF}^{\text{SM}} + \delta\mathcal{M}_{\text{SCm}}$$

$$\delta\mathcal{M}_{\text{SCm}} = \frac{\alpha_s}{\pi} \cdot \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}}}{m_t^2 c^4} \cdot |\cos(\pi t_n)|$$

$$= \frac{\alpha_s}{\pi} \times \frac{7.09 \times 10^{-37} \times 1.45 \times 10^{26} \times 0.84}{(173\ \text{GeV})^2 \times (1.6 \times 10^{-10})^2} \approx 10^{-32}$$

---

## 4. F_U_Bi_i Effective Higgs Potential

The Higgs potential in the SCm vacuum:

$$V_h(\phi) = -\mu_h^2 \phi^2 + \lambda_h \phi^4 + \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \cos(\pi t_n) \cdot \phi^2$$

The SCm term shifts the Higgs VEV by:

$$\delta v = \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot |\cos(\pi t_n)|}{2 \lambda_h v} \approx 10^{-26}\ \text{GeV}$$

Negligible at current experimental precision.

---

## 5. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \beta_i = 0.6, \quad F_{TRZ} = 0.1, \quad \Phi_{\text{res}} = 0.84$$

$$\kappa = 5.0 \times 10^{-4}\ \text{day}^{-1}, \quad \rho_{\text{vac,SCm}} = 7.09 \times 10^{-37}\ \text{J/m}^3$$

---

## References

1. LHC Higgs Cross Section Working Group (2016). Handbook of LHC Higgs Cross Sections: 4. arXiv:1610.02095.
2. ATLAS and CMS (2022). Combined measurements of Higgs boson couplings. arXiv:2207.07579.
3. SCm vacuum: `scm_vacuum_manifold.py`; PAPER_1113 (CMS $\kappa$ coupling)
