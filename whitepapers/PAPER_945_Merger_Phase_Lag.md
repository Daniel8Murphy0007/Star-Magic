# PAPER_945: Merger Phase Lag

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** smbh_binary_mergers.py (MergerPhaseLag)
**Calculator:** MergerPhaseLagCalc (CP4 #529)
**CVW:** v2.0.0 compliant

---

## Abstract

We derive the UQFF phonon-induced gravitational-wave phase lag for SMBH binary mergers in the LISA band (1--100 mHz). The cumulative phase shift $\Delta\Phi = 2\pi(f_\text{max} - f_0) \cdot D_\text{total}(q) \cdot S_{26}$ yields 200--400 cycles depending on mass ratio, providing a distinctive observational signature separable from GR waveform templates.

---

## 1. Phase Lag Formula

$$\Delta\Phi = 2\pi (f_\text{max} - f_0) \cdot D_\text{total}(q) \cdot S_{26}$$

where $S_{26} = \sum_{k=1}^{26} e^{-[\text{SSq}] \cdot k/26}$ with $[\text{SSq}] = 0.57$.

---

## 2. Phase Lag vs Mass Ratio

| $q$ | $D_\text{total}$ | $\Delta\Phi$ (rad) | Cycles |
|-----|-------------------|---------------------|--------|
| 0.2 | 0.491 | $\sim 1900$ | $\sim 302$ |
| 0.5 | 0.432 | $\sim 1670$ | $\sim 266$ |
| 0.8 | 0.372 | $\sim 1440$ | $\sim 229$ |
| 1.0 | 0.333 | $\sim 1290$ | $\sim 205$ |

---

## 3. Detectability

At LISA sensitivity ($\delta\Phi \sim 0.1$ rad), phase lags of $\sim 200$--$400$ cycles are detectable with SNR $> 10^3$, making this the most constraining UQFF prediction for space-based GW detectors.

---

## 4. Source Data

- **File:** smbh_binary_mergers.py
- **Session:** 213
- **CP4 Class:** MergerPhaseLagCalc (#529)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Berti, E. et al. (2006) -- PRD, 73, 064030
3. LISA Consortium (2023) -- arXiv:2402.07571
