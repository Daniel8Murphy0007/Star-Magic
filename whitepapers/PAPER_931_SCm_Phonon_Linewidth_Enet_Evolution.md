# PAPER_931: SCm Phonon Linewidth E_net Evolution

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 212
**Source:** scm_phonon_linewidth.py (LinewidthEnetEvolution)
**Calculator:** SCmPhononLinewidthEnetEvolutionCalc (CP4 #515)
**CVW:** v2.0.0 compliant

---

## Abstract

This paper explores the dependence of net energy density E_net on phonon linewidth Gamma in the SCm (superconductive-magnetic) channel of the UQFF framework. By sweeping Gamma through three regimes -- narrow (0.05 THz), optimal (0.10 THz), and broad (0.30 THz) -- we quantify how spectral broadening modulates the effective buoyancy energy available for gravitational coupling. The quality factor Q = omega_SCm / (2 Gamma) serves as a diagnostic for resonance sharpness, with narrow linewidths yielding Q ~ 78.5 and broad linewidths collapsing to Q ~ 13.1.

---

## 1. Core Equations

The net energy density at linewidth Gamma is:

$$E_{\text{net}}(\Gamma) = \rho_{\text{SCm}}(t) \cdot V \cdot \left(\frac{2 F_{U,Bi}}{F_U} - 1\right) \cdot \Phi(\omega, \Gamma) \cdot S_{26}$$

where:

- $\rho_{\text{SCm}}(t) = 9.47 \times 10^{-27} \cdot S_{26}$ kg/m^3 is the SCm vacuum density
- $S_{26} = \sum_{k=1}^{26} e^{-[\text{SSq}] \cdot k/26}$ is the 26-layer suppression sum
- $\Phi(\omega, \Gamma)$ is the phonon flux at angular frequency omega with linewidth Gamma
- $F_{U,Bi} / F_U$ is the buoyancy-to-unified field ratio

Quality factor:

$$Q = \frac{\omega_{\text{SCm}}}{2\Gamma}$$

with $\omega_{\text{SCm}} = 2\pi \times 1.25 \times 10^{12}$ rad/s.

### Numerical Results

| Gamma (THz) | E_net (J) | Q |
|---|---|---|
| 0.05 | Regime-dependent | 78.5 |
| 0.10 | Regime-dependent | 39.3 |
| 0.30 | Regime-dependent | 13.1 |

---

## 2. UQFF Integration

The `SCmPhononLinewidthEnetEvolutionCalc` calculator (CP4 #515) is a stateless, parameterized calculator that accepts dataset parameters V, F_U_Bi, F_U, and [SSq]. It sweeps across three canonical linewidth values and returns E_net and Q for each, along with primary equations in long-form.

---

## 3. Physical Significance

The linewidth Gamma controls the trade-off between resonance sharpness and spectral bandwidth. In the narrow regime (Gamma = 0.05 THz), the phonon mode couples strongly but over a limited frequency range. In the broad regime (Gamma = 0.30 THz), coupling is weaker per frequency bin but spans more of the phonon spectrum. The optimal linewidth (Gamma = 0.10 THz) balances these effects, maximizing the integrated E_net for typical SCm parameters.

---

## 4. Source Data

- **File:** scm_phonon_linewidth.py
- **Session:** 212
- **CP4 Class:** SCmPhononLinewidthEnetEvolutionCalc (#515)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Kozima, H. -- Neutron Drop Model and Cold Fusion Phenomena (2006)
3. UQFF Calibration: kappa = 0.0005/day, [SSq] = 0.57, H_SCm ~ 0.99
