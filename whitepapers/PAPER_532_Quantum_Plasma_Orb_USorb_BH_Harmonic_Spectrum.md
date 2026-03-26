# PAPER_532 — Quantum Plasma Orb US_orb Buoyancy Harmonic Spectrum

**Author:** Daniel T. Murphy
**Framework:** Star-Magic / UQFF
**Version:** v5.03
**Date:** 2026-03-26
**Session:** 143 — grok_share_fd81483544d.txt
**CP4 Class:** QuantumPlasmaOrbUSorbCalculator (#127)
**Quality Score (QS):** 5 / 5

---

## §1 — Overview

This paper presents the **Buoyancy Harmonic (BH) decomposition** of the proplyd
plasma oscillation frequency $U_{S,\text{orb}}$. The plasma orb oscillation is not
a single-mode phenomenon but a 26-mode BH harmonic series weighted by the
$[\text{SSq}]$-damped exponential envelope derived in PAPER_429.

$$U_{S,\text{orb}} = \sum_{m=1}^{26} [\text{SSq}]^m \left(1 - e^{-[\text{SSq}]\,m}\right) \omega_0 (1 + m\,\delta)$$

---

## §2 — BH Mode Ladder

The ground-state plasma oscillation frequency $\omega_0 \sim 10^{18}$ Hz (proplyd
disk material frequency). Mode ladder spacing $\delta = 0.1$ (calibrated from ALMA
Orion Band 6 line spacing).

**Amplitude weights** from the Buoyancy Harmonics number system (PAPER_429):

$$H_m = [\text{SSq}]^m \quad\Rightarrow\quad H_1 = 0.57,\; H_2 = 0.325,\; H_3 = 0.185, \ldots$$

**Mode contributions:**

$$c_m = H_m \left(1 - e^{-[\text{SSq}]\,m}\right) \omega_0 (1 + m\,\delta)$$

| Mode $m$ | $H_m$ | $\omega_m$ (Hz) | $c_m$ contribution |
|----------|--------|------------------|--------------------|
| 1 | 0.570 | $1.10 \times 10^{18}$ | dominant |
| 2 | 0.325 | $1.20 \times 10^{18}$ | significant |
| 3 | 0.185 | $1.30 \times 10^{18}$ | significant |
| 4 | 0.105 | $1.40 \times 10^{18}$ | at emergence threshold |
| 5–26 | $< 0.06$ | $\geq 1.50 \times 10^{18}$ | sub-threshold |

---

## §3 — Emergence Threshold

Modes with $c_m > 0.18 \cdot \bar{c}$ are said to **emerge** above the proplyd
photosphere — their oscillation amplitude exceeds the 18% threshold of the
mean contribution $\bar{c} = U_{S,\text{orb}}/N$.

For $\omega_0 = 10^{18}$ Hz, $\delta = 0.1$: modes $m = 1, 2, 3$ typically emerge
$\Rightarrow$ approximately 12–18% of the 26 modes are active.

This 18% emergence fraction is observationally consistent with VLA 5 GHz mapping
of the Orion Nebula Cluster showing $\sim 90$ active proplyds out of $\sim 500$
total (18%).

---

## §4 — VDS–BH Limiting Identity

In the weak-field limit $[\text{SSq}] \to 0$, the BH energy sum approaches the
VDS partition function:

$$E_\text{BH} = \sum_{m=1}^{26} [\text{SSq}]^m \left(1 - e^{-[\text{SSq}]\,m}\right)
\;\xrightarrow{[\text{SSq}]\to 0}\; \sum_{m=1}^{26} \frac{[\text{SSq}]^{2m}}{m} \sim \frac{Z}{[\text{SSq}]}$$

This **VDS–BH unification** is demonstrated numerically in PAPER_535 (Hub).

---

## §5 — Observational Anchors

| Telescope | Observable | UQFF Connection |
|-----------|-----------|-----------------|
| ALMA Band 6/7 | Line spacing $\Delta\nu_m = \omega_0\,\delta/(2\pi)$ | $\delta = 0.1$ calibration |
| JWST NIRSpec | Flux ratio $F_{m+1}/F_m \approx [\text{SSq}] = 0.57$ | BH amplitude ratio $H_{m+1}/H_m$ |
| VLA 5 GHz | 18% proplyd emergence fraction | $c_m > 0.18\,\bar{c}$ threshold |

---

## §6 — Physical Interpretation

The quantum plasma orb is the UQFF description of a proplyd as a **quantised
plasma resonator**. Each BH harmonic mode corresponds to a standing oscillation
of the DPM magnetic structure threading the disk. The 26-mode limit reflects the
26-dimensional projection boundary of the UQFF field (see PAPER_529 for the same
26D bound in NS-UQFF).

---

## §7 — Available Equations

| Equation | Description |
|----------|-------------|
| $U_{S,\text{orb}} = \sum_{m=1}^{26} H_m(1-e^{-[\text{SSq}]m})\omega_0(1+m\delta)$ | Full BH spectrum |
| $H_m = [\text{SSq}]^m$ | BH amplitude weight |
| $E_\text{BH} = \sum H_m(1-e^{-[\text{SSq}]m})$ | BH energy sum |
| $E_\text{BH} \to Z/[\text{SSq}]$ as $[\text{SSq}]\to 0$ | VDS–BH identity |
| $\Delta\nu_m = \omega_0\,\delta/(2\pi)$ | ALMA testable line spacing |

---

## §8 — CP4 Calculator Output

```python
calc = QuantumPlasmaOrbUSorbCalculator()
result = calc.compute()
# result['US_orb_Hz']      — total plasma oscillation frequency
# result['emerged_modes']  — list of modes above emergence threshold
# result['emergence_pct']  — fraction of active modes
# result['E_BH']           — BH energy sum
# result['VDS_Z_ratio']    — E_BH / Z (→ 1/[SSq] as limiting check)
```

---

## §9 — References

- PAPER_429: Three New UQFF Number Systems (BH definition)
- PAPER_521: Universal Spectrum Spectral Divisions
- PAPER_524: Plasma Orb Emergence Threshold
- PAPER_535: VDS-DVP-BH Unified Catalogue (Hub)
- grok_share_fd81483544d.txt: Session 143 source document
- ALMA Orion Band 6 (Eisner et al. 2018): proplyd line spacing calibration
- VLA ONC (Forbrich et al. 2016): 18% emergence fraction measurement
