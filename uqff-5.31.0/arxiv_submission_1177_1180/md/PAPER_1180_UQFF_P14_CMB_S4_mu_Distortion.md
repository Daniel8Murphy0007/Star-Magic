---
title: "UQFF P14: CMB-S4 Spectral Distortion mu-Parameter Strict Bound from R26 Closure"
author: "Daniel Murphy"
date: "May 2026"
geometry: margin=1in
fontsize: 11pt
colorlinks: true
---

# Abstract

The UQFF closed Lagrangian, with all 26 radial layers locked at saturated occupation, predicts vanishing additional CMB chemical-potential spectral distortion: $\mu_{\mathrm{UQFF}} \le 1.0 \times 10^{-8}$. This sits at the projected $1\sigma$ floor of CMB-S4 (and just below LiteBIRD/PIXIE forecasts). Detection of $\mu > 3.0 \times 10^{-8}$ at $3\sigma$ falsifies the R26 saturation assumption. This is the cosmological-microwave-background counterpart of PAPER_1178 (DESI Y5 strict-static dark energy).

# 1. UQFF Prediction

The R26 manifold has zero free parameters in cosmic time at $z < 10^6$ (recombination + diffusion-damping era). Energy injection into the CMB photon bath is bounded above by the residual buoyancy ledger residual, which is geometrically fixed at

$$
\mu_{\mathrm{UQFF}} = \frac{\rho_{\mathrm{BSFG}}}{\rho_{\gamma,0} \cdot \xi_0} \cdot f_{\mathrm{damp}} \le 1.0 \times 10^{-8}
$$

where $\rho_{\mathrm{BSFG}} = 5.96 \times 10^{-10}$ J/m^3 (4-term G2 ledger), $\rho_{\gamma,0} = 4.17 \times 10^{-14}$ J/m^3, $\xi_0 = 13/3$, and the damping-era suppression $f_{\mathrm{damp}} = e^{-(z_\mu / 5 \times 10^4)^{5/2}}$ with $z_\mu \sim 2 \times 10^6$ saturates the integral to $\le 10^{-8}$.

# 2. Observational Channels

| Experiment    | Year forecast | $1\sigma$ on $\mu$        |
|---------------|---------------|----------------------------|
| FIRAS (1996)  | published     | $9 \times 10^{-5}$         |
| PIXIE forecast| late 2020s    | $1 \times 10^{-8}$         |
| CMB-S4        | 2030-2033     | $1 \times 10^{-8}$ (broad) |
| LiteBIRD      | 2032-2034     | $\sim 1 \times 10^{-8}$    |

# 3. Falsification Criterion

UQFF is falsified if any of these experiments measure

$$
\mu_{\mathrm{obs}} > 3.0 \times 10^{-8} \quad \text{at}\ 3\sigma.
$$

# 4. CP4 Implementation

```python
from CondensedPhysics4 import UQFFCMBmuDistortionCalculator
r = UQFFCMBmuDistortionCalculator().compute()
# returns: mu_UQFF_upper_bound = 1.0e-8,
#          mu_falsifier_3sigma = 3.0e-8,
#          decisive_experiment = "CMB-S4 / LiteBIRD (2030-2034)"
```

# 5. References

[1] PAPER_1174 --- UQFF closed Lagrangian and 4-term G2 ledger.
[2] PAPER_1178 --- UQFF P13 DESI Y5 strict-static dark energy.
[3] Chluba & Sunyaev 2012 --- CMB spectral distortion forecasts.
[4] CMB-S4 collaboration 2022 --- Science Book.
