# PAPER_571: t_neg Photon Arrival Timing via Negative Time Delay in DPM Framework

**Author:** Daniel T. Murphy — Star Magic / UQFF Framework  
**Session:** 153b  
**Gap-Fill For:** AldersOlbersBSFGMetricGapAnalysisCalculator (#160, PAPER_566) — Completed Extension 5  
**Date:** 2026-03-29  
**QS:** 5/5  

---

## §1 Abstract

Standard light-travel-time calculations assume photons arrive at precisely $t_\text{obs} = r/c$. UQFF introduces a **negative-time delay** $t_\text{neg}$ (PAPER_519): photons from distant shells experience a DPM-modified arrival time due to vacuum buoyancy lag, giving an adjusted observation time:

$$t_\text{adj} = \frac{t_\text{obs}}{1 + \Delta_\text{dil}} + t_\text{neg}$$

This modifies the Olbers shell brightness: photons that arrive "late" (with $t_\text{neg} < 0$) contribute to a different effective shell, spreading the radiance across the shell hierarchy and further damping $B_\text{sky}$.

---

## §2 Negative Time Delay — PAPER_519

From the ShellRadiancePrototypeEquationCalculator (PAPER_519), the $t_\text{neg}$ timing correction encodes DPM vacuum lag:

$$t_\text{adj} = \frac{t_\text{obs}}{1 + \Delta_\text{dil}} + t_\text{neg}$$

where $\Delta_\text{dil} = [\text{SSq}] \cdot (n/26)^2$ is the DPM dilation factor for shell $n$.

Per-shell $t_\text{neg}$ distribution:

$$\Delta t_{\text{neg},n} = t_\text{neg} \cdot \frac{n}{26}$$

So inner shells (small $n$) have smaller negative delays; outer shells (large $n$) have the full $t_\text{neg}$.

---

## §3 DPM-Modified Light Travel

The radial null geodesic in the DPM vacuum is modified:

$$\frac{dr}{dt}\bigg|_\text{DPM} = c \left(1 - \frac{\kappa_\text{DPM} \, [\text{SSq}]}{r^{1/26}}\right)$$

Integrating over shell $n$:

$$t_n^\text{DPM} = \frac{r_n}{c} + t_\text{neg} \cdot \frac{n}{N} + \int_0^{r_n} \frac{\kappa_\text{DPM} [\text{SSq}]}{c \, r^{1/26}} \, dr$$

The last term gives a logarithmic correction to the classical travel time.

---

## §4 Effect on Shell Brightness

A shell-$n$ photon that arrives at $t_\text{adj}$ instead of $t_\text{obs}$ contributes to an effective redshift:

$$z_n^\text{eff} = z_n + \delta z_n, \qquad \delta z_n = -H_0 \cdot |t_\text{neg}| \cdot \frac{n}{N}$$

Modified shell brightness:

$$B_n^{t_\text{neg}} = \frac{n_\star L_\star \Delta r}{4\pi c (1 + z_n^\text{eff})^4} \cdot R_{\mathrm{Ug1},n}$$

For $t_\text{neg} = -1$ s (a small but non-zero delay), $\delta z_n \approx -2.4 \times 10^{-18} \cdot n$ — negligible individually but cumulative over 26 shells.

---

## §5 $t_\text{neg}$ Gradient Effect

The gradient of $B_n$ with respect to $t_\text{neg}$:

$$\frac{\partial B_n}{\partial t_\text{neg}} = -4 B_n \cdot \frac{H_0}{(1+z_n)} \cdot \frac{n}{N}$$

Summing the gradient correction over all shells:

$$\delta B_\text{sky} = \sum_{n=1}^{26} B_n \cdot \Delta t_{\text{neg},n} \cdot \frac{\partial \ln B_n}{\partial t_\text{neg}}$$

$$= -4 H_0 t_\text{neg} \sum_{n=1}^{26} B_n \cdot \frac{n^2}{N(1+z_n)}$$

This provides a systematic blue/red-shift correction to the total sky brightness.

---

## §6 DPM ProtoH Full Formula

The ProtoH formula from PAPER_519:

$$B_\text{total} = B_\text{sky}^\text{UQFF} + \text{DPM}_\text{react} \cdot P_\text{order} \cdot |t_\text{neg}|$$

In its full shell-explicit form:

$$B_\text{total} = \sum_{n=1}^{26} B_n \left(1 + \frac{\partial B_n}{\partial t_\text{neg}} \cdot \frac{|t_\text{neg}|}{B_n}\right)$$

$$= \sum_{n=1}^{26} B_n \left(1 - \frac{4 H_0 |t_\text{neg}| n^2}{N(1+z_n)}\right)$$

For $|t_\text{neg}| = 1$ s, the correction is of order $10^{-17}$ per shell — negligible at cosmological scales but grows with $|t_\text{neg}| \to t_\text{Hubble}$.

---

## §7 Physical Interpretation

The $t_\text{neg}$ timing effect represents **vacuum buoyancy lag**: photons from distant shells are slightly "retarded" by the DPM vacuum field, arriving later than the classical prediction. This means the universe effectively appears *younger* when observed from a UQFF perspective — reducing the effective sky brightness by delaying photon arrival from high-$z$ shells.

The effect is coupled to the BSFG horizon blinking (PAPER_566): the $\cos(\pi t_n)$ phase in the aether metric creates a periodic $t_\text{neg}$ whose average over many cycles is zero, but whose variance creates an effective line broadening in the Olbers integral.

---

## §8 Testable Predictions

1. **Pulsar timing:** The DPM-modified geodesic predicts nanosecond deviations in pulsar arrival times as a function of $n_\text{shell}$ — testable with PPTA/NANOGrav.
2. **FRB dispersion:** Fast Radio Burst dispersion measures should show a small $t_\text{neg}$ excess at $z > 1$ — encoded in the DPM-modified $\text{DM} \propto \int n_e \, dt^\text{DPM}$.
3. **Integral effect:** The total correction $\delta B_\text{sky} / B_\text{sky} \approx -4 H_0 |t_\text{neg}| \langle n^2/(N(1+z))\rangle$ — effectively $\sim 10^{-17}$ for $|t_\text{neg}| = 1$ s, but larger for $|t_\text{neg}| \sim t_\text{Hubble}$.

---

## §9 Builds On / Addresses

| Paper | Role |
|-------|------|
| PAPER_519 | ShellRadiancePrototypeEquationCalculator — original $t_\text{neg}$ definition |
| PAPER_516 | DPM layered shell energy — $\kappa_\text{DPM}$ coupling |
| PAPER_564 | DPM 26-shell Olbers (extended here with $t_\text{neg}$) |
| PAPER_566 | Gap analysis — this is Missing Extension 5 |

---

*PAPER_571 — Star Magic UQFF Framework — QS 5/5*
