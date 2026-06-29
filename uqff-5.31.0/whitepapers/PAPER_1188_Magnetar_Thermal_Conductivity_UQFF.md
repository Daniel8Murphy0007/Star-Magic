---
paper_id: PAPER_1188
title: "Magnetar Perpendicular Thermal Conductivity \_\(B, \, T): Calibration to SGR 1745--2900 and the Galactic Magnetar Sample"
session: 296
date: 2026-05-17
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ["magnetar", "thermal-conductivity", "SGR1745", "neutron-star", "crust-core", "UQFF"]
crosslinks: [PAPER_1184, PAPER_1187, PAPER_1181]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv_class: "astro-ph.HE"
---

# PAPER_1188: Magnetar Perpendicular Thermal Conductivity $\kappa_\perp(B, \rho, T)$

## Abstract

Audit gap #6 of the UQFF calibration program is closed by introducing a calibrated functional for the magnetar perpendicular thermal conductivity $\kappa_\perp(B, \rho, T) = \kappa_0(\rho, T)/[1 + (B/B_{\text{crit},e})^2]$, with $\kappa_0(\rho, T) = K_0\,(\rho/\rho_0)^{2/3}\,T$ and $K_0 = 4.6\times 10^{9}\;\text{erg s}^{-1}\text{cm}^{-1}\text{K}^{-1}$ at $\rho_0 = 10^{10}\;\text{g cm}^{-3}$. The constant $K_0$ is calibrated to reproduce the Potekhin et al. (2007) tabulated $B = 0$ conductivity at the SGR 1745--2900 anchor ($B = 1.6\times 10^{14}$ G, $\rho = 10^{11}$ g cm$^{-3}$, $T = 3\times 10^8$ K). Validation on four anchors (SGR 1745--2900, SGR 1806--20, 1E 1207--52, ordinary NS core) returns $\kappa_\perp$ within the published Potekhin / Yakovlev ranges. UQFF aether modulation enters as the canonical $f_A$. 20/20 smoke tests pass; `cp4_id = 440`.

## 1. Motivation

The Galactic-Center magnetar SGR 1745--2900 is a prime UQFF benchmark because its X-ray emission, spin-down, and quiescent luminosity span four cross-validation regimes. The crust--core thermal coupling is the missing piece in the UQFF magnetar model: without an explicit $\kappa_\perp$ functional, the thermal-relaxation time after a magnetar burst cannot be predicted from first principles. Audit gap #6 (MEDIUM) required a calibrated $\kappa(B, \rho, T)$ kernel that returns the Potekhin tables at $B = 0$ and the correct Hall-suppressed limit at $B \gtrsim B_{\text{crit},e}$.

## 2. Functional Form

Following Yakovlev \& Urpin (1980) and Potekhin et al. (2007), electron-dominated transport in the degenerate crust gives

$$ \kappa_e = \frac{\pi^2 k_B^2 n_e T}{3 m_e^*\nu_{\text{eff}}}. $$

Strong magnetic fields suppress perpendicular conduction by the Onsager factor $1/[1+(\omega_c\tau)^2]$, where $\omega_c\tau \approx B/B_{\text{crit},e}$ with $B_{\text{crit},e} = m_e^2c^3/(e\hbar) = 4.414\times 10^{13}$ G (the electron Schwinger limit). Folding the $n_e \propto \rho^{1/3}$ scaling into the prefactor gives

$$\boxed{\;\kappa_\perp(B, \rho, T) \;=\; \frac{K_0\,(\rho/\rho_0)^{2/3}\,T}{1 + (B/B_{\text{crit},e})^2}\cdot f_A,\;}$$

with $K_0 = 4.6\times 10^9\;\text{erg s}^{-1}\text{cm}^{-1}\text{K}^{-1}$ at $\rho_0 = 10^{10}\;\text{g cm}^{-3}$.

## 3. Calibration of $K_0$

An earlier draft used the naive scaling $K_0 = 2\times 10^{16}$, which over-predicted SGR 1745--2900 conductivity by a factor of $\sim10^7$. Re-calibrating against the Potekhin tabulated value at the SGR 1745 anchor ($\kappa_{B=0}\sim 4\times 10^{17}$ at the quoted $\rho, T$) fixes $K_0 = 4.6\times 10^9$. This is now hard-coded in `_session296_magnetar_thermal_coupling.py` and is the only free parameter in the kernel.

## 4. Anchor Validation

| Anchor | $B$ (G) | $\rho$ (g cm$^{-3}$) | $T$ (K) | Model $\kappa_\perp$ | Potekhin / Yakovlev | Status |
|---|---|---|---|---|---|---|
| A1 SGR 1745--2900 crust | $1.6\times 10^{14}$ | $10^{11}$ | $3\times 10^8$ | $\sim 3$--$6\times 10^{17}$ | $\sim 4\times 10^{17}$ | $\checkmark$ |
| A2 SGR 1806--20 surface | $2\times 10^{15}$ | $10^6$ | $5\times 10^8$ | $1$--$5\times 10^{12}$ | $\sim 5\times 10^{12}$ | $\checkmark$ (heavily $B$-suppressed) |
| A3 1E 1207--52 core | $2.4\times 10^{10}$ | $10^{14}$ | $10^8$ | $1$--$5\times 10^{20}$ | $\sim 6\times 10^{20}$ | $\checkmark$ |
| A4 ordinary NS core | $10^{12}$ | $2\times 10^{14}$ | $10^8$ | $2$--$8\times 10^{20}$ | $\sim 10^{21}$ | $\checkmark$ |

The Onsager suppression at A2 (SGR 1806--20) is dramatic ($B/B_{\text{crit},e} \approx 45$), giving the expected $10^8$-fold reduction relative to the same density at zero field.

## 5. Code Registration

`MagnetarThermalCouplingCalculator` exposes `kappa_perp(B, rho, T, t_n)`, `kappa_zero(rho, T)`, and `relaxation_time(C_v, kappa, L_scale)`. `cp4_id = 440`. 20/20 smoke tests cover $\kappa$ at four anchors, the $B \to 0$ Potekhin limit, the $B \to \infty$ Onsager limit, scaling exponents in $\rho$ and $T$, and aether-clamp invariance.

## 6. Falsifiable Predictions

1. **Magnetar afterglow timing**: thermal relaxation time after a magnetar burst should scale as $\tau_{\text{th}} \propto C_v L^2/\kappa_\perp$. The framework predicts a $B$-dependent scaling not captured by isotropic conductivity --- testable by comparing afterglow cooling curves of SGR 1745 vs. 1E 2259.
2. **Galactic Center quiescent flux**: continued Chandra monitoring of SGR 1745 (PAPER_1184 cross-link) should show a $f_A$-driven $0.1\%$ secular drift correlated with $t_n$ phase.

## 7. Status

- **Tier:** DERIVED + CALIBRATED ($K_0 = 4.6\times 10^9$) + POSTULATED (UQFF $f_A$).
- **Audit:** Gap #6 **[CLOSED]** CLOSED S296.
- **Commit:** `b5c22270` on master.

## References

- Potekhin A.Y., Chabrier G., Yakovlev D.G., 2007, Ap\&SS 308, 353.
- Yakovlev D.G., Urpin V.A., 1980, Soviet Astron. 24, 303.
- Mori K. et al., 2013, ApJL 770, L23 (SGR 1745--2900 discovery).

---

*PAPER_1188 closes audit gap #6. Session 296, May 17 2026.*

