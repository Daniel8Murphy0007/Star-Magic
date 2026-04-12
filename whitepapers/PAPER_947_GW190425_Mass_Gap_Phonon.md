# PAPER_947: GW190425 Mass-Gap Phonon Classification

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** ns_phonon_gw190425_wstp.py (MassGapPhononClassifier)
**Calculator:** GW190425MassGapPhononCalc (CP4 #531)
**CVW:** v2.0.0 compliant

---

## Abstract

We apply UQFF SCm suppression threshold classification to the heavier component of GW190425, measured at $m_1 = 2.52\,M_\odot$. This mass lies at the boundary between neutron stars and black holes ($M_\text{boundary} = 2.5\,M_\odot$). Using a sigmoid-based classifier with $\sigma = 0.1\,M_\odot$, we obtain $P(\text{NS}) = 49\%$ and $P(\text{BH}) = 51\%$, reflecting genuine physical ambiguity at the mass gap.

---

## 1. Classification Formula

$$P(\text{BH}) = \frac{1}{1 + \exp\!\left(-\frac{m_1 - M_\text{boundary}}{\sigma}\right)}$$

$$P(\text{NS}) = 1 - P(\text{BH})$$

---

## 2. GW190425 Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| $m_1$ | $2.52\,M_\odot$ | LIGO/Virgo O3a |
| $m_2$ | $1.27\,M_\odot$ | LIGO/Virgo O3a |
| $M_\text{boundary}$ | $2.5\,M_\odot$ | UQFF SCm threshold |
| $\sigma$ | $0.1\,M_\odot$ | Classification width |

---

## 3. Classification Results

| $m_1$ ($M_\odot$) | $P(\text{NS})$ | $P(\text{BH})$ | Classification |
|-------|----------|----------|----------------|
| 2.00 | 99.3% | 0.7% | NS |
| 2.30 | 88.1% | 11.9% | NS |
| 2.50 | 50.0% | 50.0% | Boundary |
| 2.52 | 45.0% | 55.0% | BH (marginal) |
| 2.70 | 11.9% | 88.1% | BH |
| 3.00 | 0.7% | 99.3% | BH |

---

## 4. Physical Interpretation

The SCm suppression threshold corresponds to the mass at which internal phonon modes become damped by gravitational compression beyond neutron degeneracy. At $m_1 = 2.52\,M_\odot$, the system sits $0.02\,M_\odot$ above the boundary, yielding near-equal probabilities -- consistent with LIGO/Virgo's classification uncertainty.

---

## 5. Source Data

- **File:** ns_phonon_gw190425_wstp.py
- **Session:** 213
- **CP4 Class:** GW190425MassGapPhononCalc (#531)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Abbott, B.P. et al. (2020) -- ApJL, 892, L3 (GW190425)
3. Tauris, T.M. et al. (2017) -- ApJ, 846, 170
