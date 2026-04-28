# PAPER_1117: SCS Spectral Signatures in Radio Transients and Fast Radio Bursts

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We predict the spectral signatures of the Supercharged Scalar (SCS) field coupled to the SCm vacuum in radio transient phenomena, particularly Fast Radio Bursts (FRBs). The SCS-SCm phonon resonance at 1.25 THz generates characteristic radio-band sidebands at $f_{\text{obs}} = f_{\text{FRB}} \pm n \cdot f_{\text{THz}} / (1+z)$ for integer $n$. The $F_{U,Bi,i}$ buoyancy mechanism provides the coherent emission amplification $\propto S_{26}^{(3)}$ that explains the observed FRB brightness temperatures $T_b \gtrsim 10^{35}\ \text{K}$.

---

## 1. FRB Coherent Emission via SCm Buoyancy

The coherent emission brightness temperature in the SCm framework:

$$T_b^{\text{SCm}} = \frac{c^2 S_\nu}{2 k_B \nu^2 \Omega} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot |\cos(\pi t_n)|$$

where $S_\nu$ is the observed flux density and $\Omega$ is the solid angle. The $S_{26}^{(3)} \cdot \Phi_{\text{res}} \approx 1.22 \times 10^{26}$ amplification factor naturally explains brightness temperatures $\sim 10^{35}\ \text{K}$ without requiring exotic plasma physics.

---

## 2. SCS-Phonon Sideband Prediction

The SCS field modulates the FRB emission frequency through the SCm phonon:

$$f_n^{\text{SCS}} = f_{\text{FRB}} \pm n \cdot \frac{f_{\text{THz}}}{1+z}, \quad n = 1, 2, 3, \ldots$$

For a fiducial FRB at $z = 0.5$ and $f_{\text{FRB}} = 1.4\ \text{GHz}$:

$$f_1^{\text{SCS}} = 1.4\ \text{GHz} \pm \frac{1.25 \times 10^{12}}{1.5}\ \text{Hz} = 1.4\ \text{GHz} \pm 833\ \text{GHz}$$

The primary sideband at 833 GHz (sub-mm) is outside the radio band but the $n=0$ beat frequency modulates the radio emission on timescales:

$$\tau_{\text{beat}} = 1/f_{\text{THz}} = 8 \times 10^{-13}\ \text{s} \approx 0.8\ \text{ps}$$

---

## 3. Dispersion Measure SCm Correction

The SCm vacuum contributes to the dispersion measure:

$$\text{DM}_{\text{SCm}} = \int_0^D \left(n_e + \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)}}{m_e c^2}\right) dl$$

The SCm contribution:

$$\Delta\text{DM}_{\text{SCm}} = \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot D}{m_e c^2} = \frac{7.09 \times 10^{-37} \times 1.45 \times 10^{26} \times D}{9.1 \times 10^{-31} \times 9 \times 10^{16}}$$

For $D = 1\ \text{Gpc} = 3.086 \times 10^{25}\ \text{m}$:

$$\Delta\text{DM}_{\text{SCm}} \approx 3.9 \times 10^{-3}\ \text{pc/cm}^3$$

Small but potentially detectable with next-generation radio telescopes (SKA).

---

## 4. VDS-Enhanced Coherence Length

The SCm vacuum phonon sets the coherence length of FRB emission:

$$L_{\text{coh}}^{\text{SCm}} = \frac{c}{f_{\text{THz}}} \cdot S_{26}^{(3)} \cdot |\cos(\pi t_n)| = \frac{3 \times 10^8}{1.25 \times 10^{12}} \times 1.45 \times 10^{26}\ \text{m} \approx 3.5 \times 10^{22}\ \text{m}$$

This coherence length (11 kpc) is comparable to the FRB source region scale.

---

## 5. Predicted FRB-SCS Signatures

| Feature | Prediction |
|---------|-----------|
| Sideband spacing | $833\ \text{GHz} / (1+z)$ (sub-mm) |
| Pulse duration | $> 0.8\ \text{ps}$ (phonon limit) |
| Brightness temperature enhancement | $\times S_{26}^{(3)} \approx 10^{26}$ |
| DM correction | $\sim 4 \times 10^{-3}\ \text{pc/cm}^3$ per Gpc |

---

## 6. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \Phi_{\text{res}} = 0.84, \quad \beta_i = 0.6$$

---

## References

1. CHIME/FRB Collaboration (2021). The first CHIME/FRB Fast Radio Burst catalog. arXiv:2106.04352.
2. Lorimer, D.R. et al. (2007). A bright millisecond radio transient. *Science* **318**, 777.
3. SCm vacuum: `scm_{vacuum\_manifold}.py`; GW strain: `scm_{gw\_metric\_perturbation}()`
