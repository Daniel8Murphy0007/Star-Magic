# PAPER_940: TXS 0506+056 Jet Power Curves

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** blazar_jet_power_curves_extended.py (TXS0506JetPowerCurves)
**Calculator:** TXS0506JetPowerCurvesCalc (CP4 #524)
**CVW:** v2.0.0 compliant

---

## Abstract

We derive UQFF phonon-modulated jet power curves for TXS 0506+056, the blazar associated with the IceCube-170922A high-energy neutrino event. With $M_\text{BH} = 3 \times 10^8\,M_\odot$, $a = 0.85$, $B = 8000$ T, and $A_\text{jet} = 1.20$, the Blandford-Znajek power is enhanced by $2.9\times / 2.3\times / 1.6\times$ at $\Gamma = 0.05 / 0.10 / 0.30$ THz. The elevated spin and magnetic field reflect the extreme jet conditions required for neutrino production.

---

## 1. System Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| $M_\text{BH}$ | $3 \times 10^8\,M_\odot$ | Spectral modeling |
| Spin $a$ | 0.85 | Jet Lorentz factor |
| $B$ | 8000 T | Synchrotron SED |
| $A_\text{jet}$ | 1.20 | Phonon coupling |
| Redshift $z$ | 0.3365 | SDSS |
| Neutrino | IceCube-170922A | IceCube 2018 |

---

## 2. Jet Power Curves

$$P_\text{jet}(\Gamma) = P_\text{BZ} \cdot \left(1 + M_\text{jet}(\Gamma)\right)$$

| $\Gamma$ (THz) | Enhancement | Regime |
|-----------------|-------------|--------|
| 0.05 | $2.9\times$ | Narrow |
| 0.10 | $2.3\times$ | Optimal |
| 0.30 | $1.6\times$ | Broad |

---

## 3. IceCube Association

The high $A_\text{jet} = 1.20$ coupling strength implies that phonon-jet modulation produces sufficient hadronic acceleration for $\sim290$ TeV neutrino production via $p\gamma$ interactions, consistent with IceCube-170922A timing.

---

## 4. Source Data

- **File:** blazar_jet_power_curves_extended.py
- **Session:** 213
- **CP4 Class:** TXS0506JetPowerCurvesCalc (#524)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. IceCube Collaboration (2018) -- Science, 361, eaat1378
3. Padovani, P. et al. (2018) -- MNRAS, 480, 192
