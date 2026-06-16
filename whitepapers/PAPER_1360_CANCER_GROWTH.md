# PAPER_1360 — Cancer Growth Law Extension

**Framework:** UQFF — Star-Magic v5.27+ | **Tier:** G Bio/Complex (CLOSED)
**Calculator:** `calculate_paradox({"paradox": "cancer_growth_law"})`

## Master Expression
```
rate ~ N_cells × (F_TRZ × β_i)^k extends Peto
```
**Match: mechanism**

## Interpretation
Cancer rate scales geometrically via F_TRZ × β_i modulation. Per-cell SCm threshold provides body-size invariance (extends Peto).

## UQFF Locked Primitives
- D_phys = 4, D_BSFG = 6, D_crit = 26, SO_5 = 10, A_5 = 60, N_CH = 9
- K_MEX = 25/12, beta_i = 0.6029, Phi_res = 0.84, Lambda = 0.00729735, F_TRZ = 0.1

## NOT REPLACEMENT
UQFF derives via integer-primitive identity. Zero fit parameters.

---
**Copyright** — Daniel T. Murphy, daniel.murphy00@gmail.com, Youngstown, OH.
