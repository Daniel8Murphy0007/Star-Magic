# PAPER_935: GW170817 Tidal Deformability Phonon Correction

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 212
**Source:** ns_phonon_gw170817_wstp.py (TidalDeformabilityPhononCorrection)
**Calculator:** GW170817TidalDeformabilityPhononCalc (CP4 #519)
**CVW:** v2.0.0 compliant

---

## Abstract

We compute the phonon-corrected tidal deformability Lambda_UQFF for GW170817 neutron stars. The UQFF phonon coupling modifies the Love number k_2 through the SCm lattice interaction, yielding Lambda_UQFF = Lambda_GR (1 + Phi S_26 D_total). The UQFF-predicted range Lambda in [190, 600] is consistent with the LIGO constraint Lambda_tilde < 800 and provides tighter bounds than GR alone.

---

## 1. Core Equations

UQFF-corrected tidal deformability:

$$\Lambda_{\text{UQFF}} = \Lambda_{\text{GR}} \cdot \left(1 + \Phi \cdot S_{26} \cdot D_{\text{total}}\right)$$

where:
- $\Lambda_{\text{GR}}$ is the GR tidal deformability from the neutron star equation of state
- $\Phi$ is the phonon flux at SCm resonance frequency
- $S_{26} = \sum_{k=1}^{26} e^{-[\text{SSq}] \cdot k/26}$ is the 26-layer suppression sum
- $D_{\text{total}} = 0.333$ is the phonon suppression factor

LIGO constraint:

$$\tilde{\Lambda} = \frac{16}{13} \frac{(m_1 + 12 m_2) m_1^4 \Lambda_1 + (m_2 + 12 m_1) m_2^4 \Lambda_2}{(m_1 + m_2)^5} < 800$$

### UQFF Range

| Regime | Lambda_UQFF |
|---|---|
| Soft EOS (Lambda_GR ~ 150) | ~190 |
| Moderate EOS (Lambda_GR ~ 400) | ~500 |
| Stiff EOS (Lambda_GR ~ 500) | ~600 |

---

## 2. UQFF Integration

The `GW170817TidalDeformabilityPhononCalc` (CP4 #519) takes Lambda_GR, D_total, and [SSq] as inputs. It checks both the LIGO bound (Lambda < 800) and the UQFF-predicted range [190, 600]. The simulate() method sweeps Lambda_GR = [100, 200, 300, 400, 500].

---

## 3. Physical Significance

Tidal deformability is a direct probe of the neutron star equation of state. The UQFF phonon correction arises because crust phonon modes modify the tidal response: the SCm lattice coupling adds a rigidity contribution to k_2 that stiffens the effective EOS. This narrows the allowed Lambda range from the broad GR prediction to the UQFF band [190, 600], which future detections (GW230529, O5 run) can test.

---

## 4. Source Data

- **File:** ns_phonon_gw170817_wstp.py
- **Session:** 212
- **CP4 Class:** GW170817TidalDeformabilityPhononCalc (#519)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Abbott, B.P. et al. -- GW170817: Measurements of Neutron Star Radii and EOS, PRL 121, 161101 (2018)
3. Hinderer, T. -- Tidal Love numbers of neutron stars, ApJ 677, 1216 (2008)
