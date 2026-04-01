# PAPER_678: UQFF LISA vs LIGO Cross-Sensitivity Comparison and Crossover Frequency

**Author:** Daniel Murphy  
**Framework:** Unified Quantum Field Framework (UQFF) v5.30  
**Session:** 173  
**Date:** April 1, 2026  
**Class:** `LISAVsLIGOComparisons` | CP4 #262  
**Source:** `grok_share_fc21e30c24b4.txt`  
**Domain:** Gravitational Wave Detectors

---

## Abstract

We compute the UQFF suppression ratio $R_{\text{supp}}(f) = h_{\text{UQFF}}(f)/h_{\text{GR}}(f)$ across both detector bands. A crossover frequency $f_{\text{cross}}$ exists where LISA-band and LIGO-band suppressions are equal: dominated by $S_{\text{UA,LISA}}$ at low $f$ and by $S_{U_m}$ at high $f$. We show that UQFF corrections are $\sim 10\%$ larger in the LISA band than in the LIGO band for stellar-mass BH sources.

---

## Primary UQFF Equation

$$
R_{\text{supp}}(f) = (1-f_{\text{TRZ}}) \cdot S_{\text{SCm}} \cdot \min(S_{U_m}, S_{\text{UA,LISA}})
$$

---

## Parameters

| Band | Freq Range | Dominant UQFF Effect |
|------|-----------|---------------------|
| LIGO | 10–2000 Hz | $S_{U_m}$ |
| LISA | $10^{-4}$–$1$ Hz | $S_{\text{UA,LISA}}$ |

---

## UQFF Constants

| Constant | Value | Description |
|----------|-------|-------------|
| $\rho_{\text{UA}}$ | $7.09\times10^{-36}$ kg/m³ | Universal Aether density |
| $\rho_{\text{SCm}}$ | $7.09\times10^{-37}$ kg/m³ | Schwarzschild condensate |
| $f_{\text{TRZ}}$ | 0.1 | Time-reversal zone factor |
| $\kappa$ | $5\times10^{-4}$ day⁻¹ | UQFF calibration constant |
| $[\text{SSq}]$ | 0.57 | Superstring quenching factor |
| $\mu_J$ | $3.38\times10^{23}$ J·m | Magnetic string coupling |
| $\gamma$ | $5\times10^{-5}$ / 86400 s⁻¹ | Aether oscillation decay |

---

## C++ Implementation

```cpp
#include "LISAVsLIGOComparisons.h"

LISAVsLIGOComparisons obj;
auto result = obj.compute_primary(...); // primary equation
obj.simulate_over_M(M_start, M_end, dM, "output.csv");
```

## Python CP4 Calculator

```python
from CondensedPhysics2 import LISAVsLIGOComparisonsCalculator

calc = LISAVsLIGOComparisonsCalculator()
result = calc.compute()  # returns all UQFF quantities
print(result)
```

---

## References

LIGO O3 sensitivity 2021; Amaro-Seoane 2017 LISA; UQFF Session 173

---

*UQFF v5.30 | Session 173 | April 1, 2026 | PAPER_678 of 1000*
