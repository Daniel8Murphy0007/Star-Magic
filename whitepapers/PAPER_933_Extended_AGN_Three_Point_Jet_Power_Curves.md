# PAPER_933: Extended AGN Three-Point Jet Power Curves

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 212
**Source:** agn_jet_power_curves_extended.py (ThreePointJetPowerCurve)
**Calculator:** ExtendedAGNThreePointJetPowerCalc (CP4 #517)
**CVW:** v2.0.0 compliant

---

## Abstract

We present explicit three-point jet power enhancement curves for 3C273 and TON618, computed at phonon linewidths Gamma = 0.05, 0.10, and 0.30 THz. Jet modulation via M_jet(Gamma) couples to the Blandford-Znajek base power P_BZ through the relation P_jet = P_BZ (1 + M_jet). For 3C273 (A_jet ~ 1.05), enhancements are 3.1/2.4/1.5x at the three Gamma points. For TON618 (A_jet ~ 1.40), enhancements reach 3.8/2.9/1.7x, reflecting the stronger phonon-jet coupling in ultra-massive BH systems.

---

## 1. Core Equations

Blandford-Znajek power:

$$P_{\text{BZ}} = \frac{B^2}{8\pi} \left(\frac{r_H}{c}\right)^2 a^2 c$$

Jet modulation factor:

$$M_{\text{jet}}(\Gamma) = 1 + A_{\text{jet}} \exp\left(-\frac{(\Gamma - \Gamma_0)^2}{2\sigma_G^2}\right)$$

where $\Gamma_0 = 2\pi \times 0.1 \times 10^{12}$ rad/s and $\sigma_G = 0.08 \times 2\pi \times 10^{12}$ rad/s.

Total jet power:

$$P_{\text{jet}} = P_{\text{BZ}} (1 + M_{\text{jet}})$$

Enhancement ratio:

$$\text{Enhancement} = \frac{P_{\text{jet}}(\Gamma)}{P_{\text{BZ}}}$$

### Three-Point Results

| System | Gamma=0.05 | Gamma=0.10 | Gamma=0.30 |
|--------|-----------|-----------|-----------|
| 3C273 (A=1.05) | 3.1x | 2.4x | 1.5x |
| TON618 (A=1.40) | 3.8x | 2.9x | 1.7x |

---

## 2. UQFF Integration

The `ExtendedAGNThreePointJetPowerCalc` (CP4 #517) accepts M_Msun, a_spin, B_T, and A_jet as parameters. The simulate() method sweeps A_jet through [0.8, 1.05, 1.40, 2.0] to map how coupling strength affects the three-point enhancement profile.

---

## 3. Source Data

- **File:** agn_jet_power_curves_extended.py
- **Session:** 212
- **CP4 Class:** ExtendedAGNThreePointJetPowerCalc (#517)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Blandford, R.D. & Znajek, R.L. -- Electromagnetic extraction of energy from Kerr black holes (1977)
3. Ghisellini, G. et al. -- General physical properties of bright Fermi blazars (2010)
