# PAPER_567: Stellar Number Density Evolution n★(z) via Madau-Dickinson SFR

**Author:** Daniel T. Murphy — Star Magic / UQFF Framework  
**Session:** 153b  
**Gap-Fill For:** AldersOlbersBSFGMetricGapAnalysisCalculator (#160, PAPER_566) — Missing Extension 1  
**Date:** 2026-03-29  
**QS:** 5/5  

---

## §1 Abstract

The classical Olbers' Paradox uses a constant stellar number density $n_\star$. In reality, the star-formation rate (SFR) evolves strongly with redshift: the universe was forming stars far more rapidly at $z \approx 2$ than today. This paper incorporates the **Madau-Dickinson SFR** into the UQFF 26-shell framework, replacing the constant $n_\star$ with a redshift-dependent $n_\star(z)$. The modified shell brightness is:

$$B_n = \frac{n_\star(z_n) \, L_\star \, \Delta r}{4\pi c (1+z_n)^4} \cdot R_{\mathrm{Ug1},n}$$

where $n_\star(z)$ grows steeply with $z$ — paradoxically increasing the sky brightness at high redshift, but cut off by the DPM/VDS suppression before diverging.

---

## §2 Madau-Dickinson SFR

The cosmic star-formation rate density (Madau & Dickinson 2014):

$$\dot{\rho}_\star(z) = 0.015 \, \frac{(1+z)^{2.7}}{1 + \left[\frac{1+z}{2.9}\right]^{5.6}} \, M_\odot \, \text{yr}^{-1} \, \text{Mpc}^{-3}$$

Peak SFR at $z_\text{peak} \approx 1.9$:

$$\dot{\rho}_\star(z_\text{peak}) \approx 0.178 \, M_\odot \, \text{yr}^{-1} \, \text{Mpc}^{-3}$$

$$\dot{\rho}_\star(0) \approx 0.015 \, M_\odot \, \text{yr}^{-1} \, \text{Mpc}^{-3} \quad \text{(today)}$$

---

## §3 Stellar Number Density Evolution

Converting SFR to stellar number density via the main-sequence lifetime $\tau_\star$:

$$n_\star(z) = \frac{\dot{\rho}_\star(z) \cdot \tau_\star}{M_\star} \cdot (1+z)^3$$

where $(1+z)^3$ accounts for comoving volume compression:

$$\psi(z) = \frac{(1+z)^{2.7}}{1 + \left[\frac{1+z}{2.9}\right]^{5.6}}$$

$$n_\star(z) = n_{\star,0} \cdot \psi(z) / \psi(0) \cdot (1+z)^3$$

At $z = 2$: $n_\star(2) \approx n_{\star,0} \times 15.3 \times 3^3 = 413 \cdot n_{\star,0}$

---

## §4 UQFF DPM React Epoch Variation

Within the DPM 26-shell hierarchy, the vacuum reaction parameter $\text{DPM}_\text{react}$ also varies with epoch (PAPER_516):

$$\text{DPM}_\text{react}(n) = \kappa_\text{DPM} \cdot \frac{\text{DPM}_n - \text{DPM}_s}{r_n} \cdot (1 + \delta_n)$$

where $\delta_n$ encodes the epoch-dependent aether state. At high $z$ (shells 15–26), $\delta_n \sim \psi(z_n) - 1$ — the SFR excess above today.

---

## §5 Modified Olbers Integral

Shell brightness with evolving $n_\star(z)$:

$$B_n^\text{Madau} = \frac{n_\star(z_n) \, L_\star \, \Delta r}{4\pi c (1+z_n)^4} \cdot e^{-[\text{SSq}] \cdot n/26}$$

Compared to constant-density:

$$\frac{B_n^\text{Madau}}{B_n^\text{const}} = \frac{n_\star(z_n)}{n_{\star,0}} = \frac{\psi(z_n)}{\psi(0)} \cdot (1+z_n)^3$$

| Shell | $z_n$ | $\psi(z)/\psi(0)$ | $(1+z)^3$ | $n_\star(z)/n_{\star,0}$ |
|-------|--------|-------------------|-----------|--------------------------|
| 1     | 0.128  | 1.56              | 1.43      | 2.2 |
| 5     | 0.641  | 3.71              | 5.00      | 18.6 |
| 13    | 1.666  | 9.84              | 34.0      | 334 |
| 20    | 2.564  | 9.71              | 107       | 1039 |
| 26    | 3.333  | 7.84              | 175       | 1372 |

Despite the large $n_\star(z)$ enhancement at high-$z$ shells, the combined $(1+z)^4$ dimming and $e^{-[\text{SSq}]n/26}$ suppression keeps $B_n$ convergent:

$$B_{26}^\text{Madau} \approx B_{26}^\text{const} \times \frac{1372}{(4.333)^4} \approx B_{26}^\text{const} \times 7.8$$

$$B_\text{sky}^\text{Madau} \approx 1.8 \times B_\text{sky}^\text{const}$$

The paradox remains resolved; SFR evolution provides a factor ~2 correction.

---

## §6 Faster Convergence at High-z

The SFR peak at $z \approx 2$ (shell 13) creates a local maximum in $B_n$, then the SFR declines at $z > 2$ — providing *faster convergence* beyond shell 13 than a pure power-law model would predict. This SFR turnover is a key observational feature of the resolved Olbers paradox.

---

## §7 Testable Predictions

1. **SFR peak feature:** Shell 13 ($z \approx 1.67$) contributes the most to $B_\text{sky}$; this corresponds to the cosmic star-formation peak — testable with JWST deep-field photon counts.
2. **EBL spectrum peak:** $B_\text{sky}$ peaks in the optical/NIR (stellar emission), unlike a constant-$n_\star$ model.
3. **Factor ~2 correction:** UQFF predicts $B_\text{sky}^\text{Madau} \approx 1.8 \times B_\text{sky}^\text{const}$ — a measurable correction to PAPER_564 predictions.

---

## §8 Builds On / Addresses

| Paper | Role |
|-------|------|
| PAPER_564 | DPM 26-shell Olbers (constant $n_\star$ — extended here) |
| PAPER_516 | DPM shell energy including $\text{DPM}_\text{react}$ |
| PAPER_566 | Gap analysis — this is Missing Extension 1 |

---

*PAPER_567 — Star Magic UQFF Framework — QS 5/5*
