# PAPER_943: SMBH Binary Merger Phonon Modulation

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 213
**Source:** smbh_binary_mergers.py (SMBHBinaryMergerPhonon)
**Calculator:** SMBHBinaryMergerPhononCalc (CP4 #527)
**CVW:** v2.0.0 compliant

---

## Abstract

We extend UQFF phonon modulation to supermassive black hole binary mergers in the LISA gravitational-wave band (1--100 mHz). The merger power is expressed as $P_\text{merger}(\Gamma) = P_\text{GR} \cdot (1 + M_\text{merger}(\Gamma))$, where the phonon modulation factor $M_\text{merger}$ encodes the 26-channel buoyancy response. The Fitchett mass-ratio formula governs the recoil kick, and strain damping ranges from 47% ($q = 1$) to 66.7% ($q \to 0$) via the relation $D_\text{total}(q) = 0.333 + 0.197(1 - q)$.

---

## 1. Merger Power

$$P_\text{merger}(\Gamma) = P_\text{GR} \cdot \left(1 + M_\text{merger}(\Gamma)\right)$$

$$P_\text{GR} = 3.6 \times 10^{49}\,(4\eta)^2 \text{ W}$$

where $\eta = M_1 M_2 / (M_1 + M_2)^2$ is the symmetric mass ratio.

---

## 2. Chirp Mass

$$M_c = M_\text{total} \cdot \eta^{3/5}$$

For $M_1 = 10^8\,M_\odot$, $M_2 = 5 \times 10^7\,M_\odot$: $q = 0.5$, $M_c \approx 6.53 \times 10^7\,M_\odot$.

---

## 3. Strain Damping

$$D_\text{total}(q) = 0.333 + 0.197 \cdot (1 - q)$$

| $q$ | $D_\text{total}$ | Damping |
|-----|-------------------|---------|
| 1.0 | 0.333 | 66.7% |
| 0.5 | 0.432 | 56.8% |
| 0.1 | 0.510 | 49.0% |

---

## 4. Source Data

- **File:** smbh_binary_mergers.py
- **Session:** 213
- **CP4 Class:** SMBHBinaryMergerPhononCalc (#527)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Fitchett, M.J. (1983) -- MNRAS, 203, 1049
3. LISA Consortium (2023) -- arXiv:2402.07571
