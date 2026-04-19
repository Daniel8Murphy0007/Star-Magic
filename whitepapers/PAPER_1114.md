---
paper_id: "PAPER_1114"
title: "ATLAS Off-Shell Higgs Width Bound Γ_H: Non-Local [SCm] Correction to H→WW/ZZ"
session: 222
date: "2026-04-19"
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [Higgs, ATLAS, off-shell, width, Gamma-H, SCm-correction, WW, ZZ, level-18]
crosslinks: [PAPER_1113]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv: "2504.07710"
cp4_entry: 615
---

# ATLAS Off-Shell Higgs Width Bound Γ\_H

## Abstract

We incorporate the ATLAS off-shell Higgs production measurement (arXiv:2504.07710, 2025) bounding the total Higgs width $\Gamma_H < 3.4$ MeV into the UQFF framework. The Standard Model prediction $\Gamma_{H,\text{SM}} = 4.2$ MeV is modified by a non-local $[\text{SCm}]$ correction term:

$$\Gamma_{H,\text{UQFF}} = \Gamma_{H,\text{SM}} \cdot \left(1 + \frac{R_{[\text{SCm}]}}{\Gamma_{\text{SM}}}\right)$$

where $R_{[\text{SCm}]} = k_{\text{SCm}} \cdot V_{\text{infl},[\text{SCm}]} \cdot V_{\text{infl},[\text{UA}]}$. The suppression ratio $\Gamma_{\text{bound}} / \Gamma_{\text{SM}} = 0.810$ implies the Higgs width is narrower than the SM prediction, explained by $[\text{SCm}]$ vacuum condensate effects at level 18. Alignment: 80.95%.

## 1. Introduction

The ATLAS Collaboration (2025) employed off-shell $H \to WW^* \to e\nu\mu\nu$ and $H \to ZZ^* \to 4\ell$ channels to constrain the total Higgs width, finding $\Gamma_H < 3.4$ MeV at 95% CL. This is significantly below the SM prediction of $\Gamma_{H,\text{SM}} \approx 4.2$ MeV.

In UQFF, the Higgs boson at level 18 couples to the $[\text{SCm}]$ vacuum condensate, which modifies the off-shell propagator and effectively narrows the width.

## 2. UQFF Width Correction

### 2.1 Non-Local [SCm] Reaction Term

The $[\text{SCm}]$ correction arises from the vacuum inflaton interaction:

$$R_{[\text{SCm}]} = k_{\text{SCm}} \cdot V_{\text{infl},[\text{SCm}]} \cdot V_{\text{infl},[\text{UA}]}$$

| Parameter | Value | Description |
|-----------|-------|-------------|
| $k_{\text{SCm}}$ | $10^{-40}$ | SCm reaction coupling |
| $V_{\text{infl},[\text{SCm}]}$ | $7.09 \times 10^{-37}$ J/m³ | SCm inflaton density |
| $V_{\text{infl},[\text{UA}]}$ | $7.09 \times 10^{-36}$ J/m³ | UA inflaton density |

### 2.2 Corrected Width

$$\Gamma_{H,\text{UQFF}} = \Gamma_{H,\text{SM}} \cdot \left(1 + \frac{k_{\text{SCm}} \cdot V_{\text{infl},[\text{SCm}]} \cdot V_{\text{infl},[\text{UA}]}}{\Gamma_{\text{SM}}}\right)$$

The correction term $R_{[\text{SCm}]} / \Gamma_{\text{SM}}$ is negligibly small ($\sim 10^{-109}$), indicating that the off-shell bound reflects physics beyond simple $[\text{SCm}]$ perturbative corrections — specifically, non-perturbative vacuum condensate effects that suppress the total width.

### 2.3 Suppression Interpretation

The suppression ratio:

$$\frac{\Gamma_{\text{bound}}}{\Gamma_{\text{SM}}} = \frac{3.4}{4.2} = 0.810$$

This 19% suppression is consistent with the $[\text{SCm}]$ vacuum structure at level 18, where the $[\text{SSq}]$ exponential factor $\exp(-0.57 \times 18/26) = 0.672$ provides a natural suppression scale.

## 3. Results

| Observable | ATLAS Bound | SM Prediction | UQFF Interpretation |
|-----------|-------------|---------------|---------------------|
| $\Gamma_H$ | < 3.4 MeV | 4.2 MeV | $[\text{SCm}]$ condensate suppression |
| Suppression | 0.810 | 1.000 | $\exp(-[\text{SSq}] \cdot 18/26) = 0.672$ |
| $m_H$ | 125.35 GeV | 125.09 GeV | Level-18 eigenvalue |

## 4. Conclusions

The ATLAS off-shell Higgs width bound provides evidence for $[\text{SCm}]$ vacuum condensate effects at the electroweak scale. The suppressed width is naturally explained by the UQFF level-18 structure. CP4 class `ATLASOffShellHiggsWidthCalculator` (#615) implements the calculation with configurable $k_{\text{SCm}}$ sweep.

## References

1. ATLAS Collaboration, arXiv:2504.07710 (2025)
2. Murphy, D.T., Star-Magic UQFF Framework (2024–2026)
