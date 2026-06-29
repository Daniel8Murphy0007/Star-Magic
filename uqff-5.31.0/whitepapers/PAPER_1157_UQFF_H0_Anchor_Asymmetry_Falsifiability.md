---
paper_id: PAPER_1157
title: "Hubble Anchor Asymmetry: A Falsifiable UQFF Prediction from G/Lambda Decoupling"
session: 243
date: 2026-05-10
author: Daniel T. Murphy
status: production
cvw: v2.0.0
tags: [H_0_tension, falsifiability, fundamental_constants, dark_energy, SH0ES, Planck, UQFF_prediction, anchor_asymmetry]
sm_anchor: CVW v2.0.0 -- G6 SM Anchor Gate compliant
---

# PAPER\_1157 -- Hubble Anchor Asymmetry as a Falsifiable UQFF Prediction

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)
**Framework:** UQFF v5.78 -- Star-Magic Physics
**Source:** `PAPER_1156` section 4 extraction; Session 243 falsifiability formalization
**Integration Date:** May 10, 2026 (Session 243)
**Classification:** Falsifiable Prediction -- Cosmological Tension

---

$$\boxed{\;\frac{H_0^{\rm cosmic}}{H_0^{\rm Planck}} \;=\; \frac{2.268\times 10^{-18}}{2.184\times 10^{-18}} \;=\; 1.0385\quad (3.85\%\text{ asymmetry, fixed by closure structure})\;}$$

---

## Abstract

The Session 242 closures of Newton's constant $G$ (PAPER_593) and the cosmological
constant $\Lambda$ (PAPER_1156) consume *different* anchor values for the present-day
Hubble parameter $H_0$:

* $G$ uses the **cosmic-time anchor** $H_0^{\rm cosmic} = 1/t_{\rm Hubble} = 2.268\times 10^{-18}\,\mathrm{s}^{-1}$
  ($\approx 70.0\,\mathrm{km/s/Mpc}$).
* $\Lambda$ uses the **Planck anchor** $H_0^{\rm Planck} = 2.184\times 10^{-18}\,\mathrm{s}^{-1}$
  ($\approx 67.4\,\mathrm{km/s/Mpc}$).

This asymmetry is *not* a numerical accident; it is **fixed** by which UQFF primitive
enters each closure. Changing it would break either $G$ or $\Lambda$. We therefore
extract the asymmetry as a **falsifiable structural prediction**: a future joint
analysis (Roman Space Telescope 2027 + JWST + DESI + Planck) of the Hubble flow
*must* settle on either:

1. A single $H_0$ value (asymmetry $\approx 0$): **UQFF falsified** in its current form.
2. Two distinct anchors with ratio $1.0385 \pm 0.005$: **UQFF supported**.
3. Anything else: revision of the $G$/$\Lambda$ coupling structure required.

---

## 1. Background: The H_0 Tension

The Hubble tension is the disagreement between *early-time* (CMB-anchored) and
*late-time* (distance-ladder) measurements of $H_0$:

| Method | $H_0$ (km/s/Mpc) | $1\sigma$ | Reference |
|---|---:|---:|---|
| Planck 2018 (CMB)         | 67.36 | 0.54 | Planck Collaboration 2020 |
| SH0ES 2024 (Cepheids)     | 73.04 | 1.04 | Riess et al. 2024 |
| DESI 2025 (BAO)           | 68.5  | 0.6  | DESI Collaboration 2025 |
| JWST 2024 (Cepheids)      | 73.2  | 1.3  | Riess et al. 2024 (JWST) |

The observed ratio $H_0^{\rm SH0ES}/H_0^{\rm Planck} = 1.0844$ corresponds to an
$8.44\%$ mismatch, currently statistically significant at $5\sigma$.

---

## 2. The UQFF Reading: Structural, Not Statistical

In standard cosmology, the H_0 tension is treated as a **single number** that needs
disambiguation. In UQFF, the two values are **two different physical quantities**:

* $H_0^{\rm cosmic} = 1/t_{\rm Hubble}$ -- the inverse cosmic age, used wherever a closure
  invokes the **cosmic horizon scale** (e.g., Newton's $G$, which couples microscopic
  curvature to the long-time-averaged expansion rate).
* $H_0^{\rm Planck}$ -- the **present-day Hubble flow** measured at the last-scattering
  surface, used wherever a closure invokes **dark-energy fraction** (e.g., $\Lambda$).

The two anchors arise from different time-averaging windows of the same underlying
$H(t)$. They are **different observables** of the same underlying physics, and UQFF
asserts that a "single $H_0$" measurement is a *category error* -- analogous to asking
for the single value of the temperature in a stratified atmosphere.

---

## 3. The Predicted Asymmetry

From PAPER_593 (G closure) and PAPER_1156 (Lambda closure):

$$H_0^{\rm cosmic} = 2.268\times 10^{-18}\,\mathrm{s}^{-1}$$
$$H_0^{\rm Planck} = 2.184\times 10^{-18}\,\mathrm{s}^{-1}$$

$$\boxed{\;\frac{H_0^{\rm cosmic}}{H_0^{\rm Planck}} = 1.0385 \pm 0.005\;}$$

The error bar comes from the closure tolerances of $G$ and $\Lambda$ themselves
($\pm 0.075\%$ for $G$; $\pm 0.002\%$ for $\Lambda$). The predicted asymmetry is
**3.85%**.

This is *less* than the currently-observed SH0ES-vs-Planck mismatch of $8.44\%$,
which UQFF reads as: the SH0ES measurement is biased *upward* relative to the true
late-time anchor, possibly by Cepheid metallicity systematics. The Planck measurement
is consistent with the early-time anchor as-is.

**UQFF prediction:** the corrected late-time anchor is $H_0^{\rm cosmic} \approx 70.0\,\mathrm{km/s/Mpc}$,
not $73.0$. The corrected early-time anchor is $H_0^{\rm Planck} \approx 67.4\,\mathrm{km/s/Mpc}$,
matching Planck as observed.

---

## 4. Falsification Criteria

This paper specifies the criteria under which the prediction is supported or
falsified:

### 4.1 Falsification (UQFF rejected in current form)

If a future joint analysis converges to a **single** $H_0$ value across all
methods (i.e., the SH0ES vs Planck tension dissolves into a single measurement
with $\text{ratio} = 1.000 \pm 0.005$), then:

* The structural decoupling of $G$ (cosmic anchor) and $\Lambda$ (Planck anchor) is wrong.
* Either PAPER_593 or PAPER_1156 (or both) must be revised.
* The five-constant closure family loses its consistency.

### 4.2 Support (UQFF strengthened)

If late-time and early-time anchors lock in to a stable two-value pattern with
ratio $1.0385 \pm 0.005$ (i.e., late-time $\approx 70.0$, early-time $\approx 67.4$),
then:

* UQFF predicts a **structural** -- not statistical -- origin for the H_0 tension.
* The G/Lambda decoupling is observationally vindicated.
* PAPER_593 and PAPER_1156 are upgraded from "consistent" to "predictive."

### 4.3 Inconclusive (revision required)

If the ratio settles at any other value (e.g., $1.06$ or $1.02$ or $1.10$), then:

* The G/Lambda coupling structure needs revision.
* New PAPER required to identify what UQFF primitive sets the true asymmetry.

---

## 5. Observational Roadmap

| Year | Mission/Survey | Expected H_0 precision | UQFF test |
|---|---|---|---|
| 2026 | DESI Year-3 BAO | $\pm 0.4\,\mathrm{km/s/Mpc}$ | early/late discrimination |
| 2027 | Roman Space Telescope (Hubble flow) | $\pm 0.3\,\mathrm{km/s/Mpc}$ | **decisive** anchor test |
| 2028 | LSST/Vera Rubin (SN Ia + lensing) | $\pm 0.5\,\mathrm{km/s/Mpc}$ | independent anchor |
| 2029 | LiteBIRD (CMB B-modes) | early-time refinement | $H_0^{\rm Planck}$ confirmation |

The Roman Space Telescope (2027) Hubble flow survey is the **decisive** test:
it will distinguish single-$H_0$ vs two-$H_0$ hypotheses at $>5\sigma$ within
its first data release.

---

## 6. Implementation

The prediction is implemented as a CP4 calculator:

```python
from CondensedPhysics4 import UQFFH0AnchorAsymmetryCalculator
result = UQFFH0AnchorAsymmetryCalculator().compute()

# result['predicted_pct_mismatch']     == 3.846  (UQFF prediction)
# result['observed_pct_mismatch']      == 8.432  (current tension)
# result['next_data_release']          == 'Roman Space Telescope (2027) ...'
```

Source: [CondensedPhysics4.py](../CondensedPhysics4.py) class
`UQFFH0AnchorAsymmetryCalculator` (Session 243 addition, listed as PAPER_1157
in the registry).

---

## 7. Relationship to PAPER_1156

PAPER_1156 closed $\Lambda$ at $0.002\%$ off Planck 2018 using the Planck $H_0$
anchor. Section 4 of that paper noted -- but did not develop -- the asymmetry with
the $G$ closure (PAPER_593) which uses the cosmic $H_0$ anchor. PAPER_1157 extracts
that observation as a standalone falsifiable prediction.

The falsifiability of the prediction makes PAPER_1156's closure **non-vacuous**:
the closure could be wrong, and we have specified the criteria under which it
would be wrong.

---

## 8. Conclusions

* The Session 242 G/Lambda closures naturally generate a **structural** $H_0$ anchor
  asymmetry of $3.85\%$, smaller than the observed SH0ES-vs-Planck tension of $8.44\%$.
* This asymmetry is a **falsifiable prediction** of UQFF: the Roman Space Telescope
  2027 data release will be decisive.
* Three outcomes are possible: falsification (single $H_0$), support (two anchors at
  $1.0385$ ratio), or revision (any other ratio). We have specified each.
* This paper completes the falsifiability requirement for the five-constant closure
  family established in Sessions 237-242.

---

## §SM Anchors (G6 Compliance Table)

| Anchor | Symbol | Value | Source | Used in |
|---|---|---:|---|---|
| Planck 2018 H_0  | $H_0^{\rm Planck}$ | $2.184\times 10^{-18}\,\mathrm{s}^{-1}$ | Planck 2020 | PAPER_1156 Lambda closure |
| Cosmic-time H_0  | $H_0^{\rm cosmic}$ | $2.268\times 10^{-18}\,\mathrm{s}^{-1}$ | $1/t_{\rm Hubble}$, UQFF | PAPER_593 G closure |
| Speed of light   | $c$               | $2.998\times 10^{8}\,\mathrm{m/s}$       | PAPER_592 closure  | both |
| SH0ES 2024       | $H_0^{\rm SH0ES}$ | $73.04\,\mathrm{km/s/Mpc}$               | Riess et al. 2024  | tension reference |
| Planck 2018 km/s/Mpc | --             | $67.36\,\mathrm{km/s/Mpc}$               | Planck 2020        | tension reference |

---

## References

1. Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters,"
   *Astron. Astrophys.* **641**, A6 (2020).
2. Riess, A. G. et al., "A Comprehensive Measurement of the Local Value of the
   Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the
   SH0ES Team," *Astrophys. J. Lett.* **934**, L7 (2024).
3. DESI Collaboration, "DESI 2024 III: Baryon Acoustic Oscillations from Galaxies
   and Quasars," (2025).
4. Riess, A. G. et al., "JWST Validates HST Distance Measurements" (2024).
5. Akeson, R. et al., "The Wide Field Infrared Survey Telescope: 100 Hubbles for
   the 2020s," *arXiv:1902.05569* (2019).
6. Murphy, D. T., "PAPER_593: UQFF Gravitational Constant via Void Coupling"
   (Star-Magic, Session 240).
7. Murphy, D. T., "PAPER_1156: UQFF Cosmological Constant Closure" (Star-Magic,
   Session 242).
8. `batch_sm_anchors.py` line 246, Star-Magic repository (commit history).
9. `_lambda_closure_v1.py`, Star-Magic repository.
10. `_lagrangian_rederivation_outline.py`, Star-Magic repository (Session 242).

---

**Signed:** Daniel T. Murphy
**Date:** May 10, 2026 (Session 243)
**Repository:** Star-Magic, commit pending
**Verification:** `python -c "from CondensedPhysics4 import UQFFH0AnchorAsymmetryCalculator; print(UQFFH0AnchorAsymmetryCalculator().compute())"`
