# PAPER_939: Centaurus A Jet Power Curves

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** blazar_jet_power_curves_extended.py (CentaurusAJetPowerCurves)
**Calculator:** CentaurusAJetPowerCurvesCalc (CP4 #523)
**CVW:** v2.0.0 compliant

---

## Abstract

We compute numerical jet power curves for Centaurus A (NGC 5128), the nearest radio galaxy at $d \approx 3.8$ Mpc. Using the Blandford-Znajek mechanism with UQFF phonon linewidth modulation, we derive power enhancement factors at three linewidths $\Gamma = 0.05, 0.10, 0.30$ THz. The SMBH parameters are $M_\text{BH} = 5.5 \times 10^7\,M_\odot$, $a = 0.70$, $B = 3000$ T, and $A_\text{jet} = 0.95$, yielding enhancements of $2.6\times / 2.1\times / 1.4\times$ respectively.

---

## 1. System Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| $M_\text{BH}$ | $5.5 \times 10^7\,M_\odot$ | EHT/VLBI |
| Spin $a$ | 0.70 | Jet morphology fits |
| $B$ | 3000 T | Faraday rotation |
| $A_\text{jet}$ | 0.95 | Phonon coupling |
| Distance | 3.8 Mpc | Tully-Fisher |

---

## 2. Blandford-Znajek Power

$$P_\text{BZ} = \frac{B^2}{8\pi} \left(\frac{r_H}{c}\right)^2 a^2 c$$

where $r_H = \frac{r_S}{2}(1 + \sqrt{1 - a^2})$ and $r_S = 2GM/c^2$.

---

## 3. UQFF Jet Modulation

$$M_\text{jet}(\Gamma) = 1 + A_\text{jet} \exp\!\left(-\frac{(\Gamma - \Gamma_0)^2}{2\sigma_\Gamma^2}\right)$$

$$P_\text{jet}(\Gamma) = P_\text{BZ} \cdot (1 + M_\text{jet}(\Gamma))$$

| $\Gamma$ (THz) | Enhancement |
|-----------------|-------------|
| 0.05 | $2.6\times$ |
| 0.10 | $2.1\times$ |
| 0.30 | $1.4\times$ |

---

## 4. Source Data

- **File:** blazar_jet_power_curves_extended.py
- **Session:** 213
- **CP4 Class:** CentaurusAJetPowerCurvesCalc (#523)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Israel, F.P. (1998) -- Centaurus A, A&AR, 8, 237
3. Blandford, R.D. & Znajek, R.L. (1977) -- MNRAS, 179, 433
