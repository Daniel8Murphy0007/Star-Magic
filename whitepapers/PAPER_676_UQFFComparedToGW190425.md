# PAPER_676: UQFF Analysis of GW190425: Heaviest Known NS Binary and Ejecta Constraints

**Author:** Daniel Murphy  
**Framework:** Unified Quantum Field Framework (UQFF) v5.30  
**Session:** 173  
**Date:** April 1, 2026  
**Class:** `UQFFComparedToGW190425` | CP4 #260  
**Source:** `grok_share_fc21e30c24b4.txt`  
**Domain:** Gravitational Waves / Dense Matter

---

## Abstract

GW190425 represents the heaviest binary neutron star system observed (total mass $3.4\\,M_\\odot$). The UQFF framework predicts a post-merger phase shift $\\phi_{\\text{UQFF}} = \\kappa f_{\\text{TRZ}} t_{\\text{merger}}$ and constrains the ejecta mass through the condensate suppression: $M_{\\text{ej,UQFF}} < M_{\\text{ej,GR}} \\cdot (\\rho_{\\text{SCm}}/\\rho_{\\text{UA}}) (1-f_{\\text{TRZ}})$, resulting in an ejecta upper limit nearly two orders of magnitude below GR.

---

## Primary UQFF Equation

$$
M_{\\text{ej,UQFF}} = 0.05 M_{\\text{tot}} \\cdot \\frac{\\rho_{\\text{SCm}}}{\\rho_{\\text{UA}}} (1-f_{\\text{TRZ}})
$$

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| $m_1$ | $1.9\\,M_\\odot$ | NS1 mass |
| $m_2$ | $1.5\\,M_\\odot$ | NS2 mass |
| $d$ | $159\\,\\text{Mpc}$ | Distance |

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
#include "UQFFComparedToGW190425.h"

UQFFComparedToGW190425 obj;
auto result = obj.compute_primary(...); // primary equation
obj.simulate_over_M(M_start, M_end, dM, "output.csv");
```

## Python CP4 Calculator

```python
from CondensedPhysics2 import UQFFComparedToGW190425Calculator

calc = UQFFComparedToGW190425Calculator()
result = calc.compute()  # returns all UQFF quantities
print(result)
```

---

## References

Abbott et al. 2020 (GW190425); UQFF Session 173

---

*UQFF v5.30 | Session 173 | April 1, 2026 | PAPER_676 of 1000*
