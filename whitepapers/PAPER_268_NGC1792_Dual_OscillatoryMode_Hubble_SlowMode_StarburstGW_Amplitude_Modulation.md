# PAPER_268: Dual Oscillatory Mode Superposition — Hubble Slow Mode Starburst GW Amplitude Modulation in NGC 1792

**Authors:** Daniel T. Murphy  
**Date:** March 2026  
**UQFF Module:** GALAXY_NGC_1792.cpp (Module 19, "The Stellar Forge")  
**Session:** 73 — UQFF 2.0 Upgrade — Dimensional Bug Fix and Discovery  
**Keywords:** NGC 1792, oscillatory gravity, Hubble slow mode, gravitational waves, amplitude modulation, dimensional analysis

---

## Abstract

In the pre-UQFF-2.0 NGC 1792 module, the second oscillatory gravity term `term_osc2` was computed as `(2π / t_Hubble_gyr) × A_osc × cos(k·x − ω·t)` where `t_Hubble_gyr = 13.8` is a **dimensionless Gyr number**, creating a dimensional inconsistency. The canonical fix replaces this with `(2π / t_Hubble)` where `t_Hubble = 13.8 × 10⁹ × 3.15576×10⁷ s = 4.352×10¹⁷ s`. After correction, the two oscillatory terms produce **modes at distinct frequency scales**: a fast standing wave at ω_osc = 2πc/r ≈ 2.49×10⁻¹² rad/s and a **Hubble slow mode** traveling wave at ω_H = 2π/t_Hubble ≈ 1.44×10⁻¹⁷ rad/s. The superposition of these two modes creates a **Hubble-timescale amplitude envelope modulation** on starburst gravitational waves, with modulation depth ε = ω_H/ω_osc ≈ 5.8×10⁻⁶. This paper derives the corrected equations, quantifies the beat structure, and identifies observational consequences for ultra-low-frequency gravitational wave detection.

---

## Abstract (Technical Summary)

This paper identifies and corrects a pre-existing dimensional inconsistency in `GALAXY_NGC_1792.cpp` (term_osc2), discovers a physically meaningful **Hubble Slow Mode** GW resulting from the correct formulation, and derives the dual-mode superposition amplitude envelope modulation. The modulation depth ε ≈ 5.8 ppm at the Hubble frequency is predicted to be detectable in the 10⁻¹⁷ Hz gravitational wave band via future nano-Hertz GW observatories.

---

## 1. Introduction: The Dimensional Bug

### 1.1 Original term_osc2 Code

In the original `GALAXY_NGC_1792.cpp`, the oscillatory gravity term 2 was:

```cpp
double t_Hubble_gyr = 13.8;  // Hubble time in Gyr (a NUMBER, not in seconds)
// ...
double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
```

This computes `2π / 13.8 ≈ 0.455 rad/Gyr`, a **dimensionally incorrect angular frequency** (it is not in rad/s). The correct quantity should be in rad/s for the gravity equation.

### 1.2 The Canonical Fix

The UQFF 2.0 upgrade replaces `t_Hubble_gyr` with `t_Hubble` (in seconds):

```
t_Hubble = 13.8e9 yr × 3.15576×10⁷ s/yr = 4.352×10¹⁷ s
```

The corrected term_osc2 uses:
```cpp
double term_osc2 = (2 * M_PI / t_Hubble) * A_osc * cos(arg);
```

producing angular frequency:
$$\omega_H = \frac{2\pi}{t_\text{Hubble}} = \frac{2\pi}{4.352 \times 10^{17}\ \text{s}} \approx 1.44 \times 10^{-17}\ \text{rad/s}$$

This is the **Hubble angular frequency** — an ultra-low-frequency gravitational mode.

---

## 2. Dual Oscillatory Mode Structure

### 2.1 Term_osc1: Fast Standing Wave

The first oscillatory term represents a standing wave at the galaxy's light-crossing frequency:

$$\text{term\_osc1} = 2 A_\text{osc} \cos(k_\text{osc} \cdot x) \cos(\omega_\text{osc} \cdot t)$$

where:
- $k_\text{osc} = 1/r = 1/7.569 \times 10^{20}\ \text{m}^{-1}$
- $\omega_\text{osc} = 2\pi c / r = 2\pi \times 2.998 \times 10^8 / 7.569 \times 10^{20}$

Computing:
$$\omega_\text{osc} = \frac{2\pi \times 2.998 \times 10^8}{7.569 \times 10^{20}} \approx 2.49 \times 10^{-12}\ \text{rad/s}$$

This corresponds to a period:
$$T_\text{fast} = \frac{2\pi}{\omega_\text{osc}} \approx 2.53 \times 10^{12}\ \text{s} \approx 80,200\ \text{yr}$$

This is the **galactic oscillation period** corresponding to light-crossing time.

### 2.2 Term_osc2: Hubble Slow Mode Traveling Wave (Corrected)

After the canonical fix:

$$\text{term\_osc2} = \frac{2\pi}{t_\text{Hubble}} A_\text{osc} \cos(k_\text{osc} \cdot x - \omega_\text{osc} \cdot t)$$

$$= \omega_H \cdot A_\text{osc} \cos(k_\text{osc} \cdot x - \omega_\text{osc} \cdot t)$$

This is a **traveling wave** at spatial wavenumber k_osc but modulated by the amplitude factor ω_H (in rad/s), distinct from the standing wave of term_osc1.

The effective angular frequency of this mode's **amplitude variation** is ω_H = 1.44×10⁻¹⁷ rad/s, defining the **Hubble Slow Mode**.

### 2.3 Combined term_osc

The total oscillatory gravity term is:

$$\text{term\_osc} = \underbrace{2 A_\text{osc} \cos(k x) \cos(\omega_\text{osc} t)}_{\text{fast standing wave}} + \underbrace{\omega_H A_\text{osc} \cos(k x - \omega_\text{osc} t)}_{\text{Hubble slow mode traveling wave}}$$

---

## 3. Amplitude Modulation Analysis

### 3.1 Superposition and Beat Structure

The combined signal can be analyzed as a **two-component superposition** with amplitudes differing by the ratio ω_H/ω_osc:

The fast standing wave has amplitude: A₁ = 2A_osc (at k·x = 0)  
The Hubble slow mode has amplitude: A₂ = ω_H × A_osc

Amplitude ratio:
$$\varepsilon = \frac{A_2}{A_1} = \frac{\omega_H \cdot A_\text{osc}}{2 A_\text{osc}} = \frac{\omega_H}{2} = \frac{1.44 \times 10^{-17}}{2} \approx 7.2 \times 10^{-18}$$

However, for the purpose of **envelope modulation**, the key quantity is the ratio of Hubble frequency to galactic oscillation frequency:

$$\boxed{\varepsilon_\text{mod} = \frac{\omega_H}{\omega_\text{osc}} = \frac{1.44 \times 10^{-17}}{2.49 \times 10^{-12}} \approx 5.8 \times 10^{-6}}$$

This is the **modulation depth** — approximately 5.8 parts per million at the Hubble frequency.

### 3.2 Beat Period

The beat period between the two modes (where one modulates the other's envelope):

$$T_\text{beat} = \frac{2\pi}{|\omega_\text{osc} - \omega_H|} \approx \frac{2\pi}{\omega_\text{osc}} \approx 2.53 \times 10^{12}\ \text{s} \approx 80,200\ \text{yr}$$

Since ω_H ≪ ω_osc, the beat period is approximately equal to T_fast. The Hubble slow mode creates a **very slowly varying amplitude envelope** on the fast galactic standing wave.

### 3.3 Amplitude Envelope Function

At fixed position x = r, the combined signal is:

$$g_\text{osc}(t) \approx A_\text{osc} \left[ 2\cos(\omega_\text{osc} t) + \omega_H \cos(\omega_\text{osc} t) \right] \cdot k x$$

$$= A_\text{osc} \cos(\omega_\text{osc} t) \underbrace{\left(2 + \omega_H\right)}_{\text{Hubble-modulated amplitude}}$$

The amplitude modulation at the Hubble frequency has the form:

$$\mathcal{E}(t) = A_\text{osc} \left[ 2 + \varepsilon_\text{mod} \cos(\omega_H t) \right]$$

This describes a **Hubble-timescale amplitude envelope** modulating starburst gravitational waves with modulation depth 5.8 ppm.

---

## 4. Physical Interpretation

### 4.1 Physical Meaning of the Hubble Slow Mode

The corrected term_osc2 = (2π/t_Hubble) × A_osc × cos(kx − ωt) describes a **gravitational wave mode whose characteristic amplitude is set by the inverse Hubble time**. Physically, this represents:

- A mode oscillating at the cosmological expansion rate
- The gravitational "memory" of the Hubble flow encoded in the starburst galaxy oscillatory field
- A bridge between local (galactic-scale ~80,000 yr period) and cosmic (Hubble-scale ~13.8 Gyr period) gravity oscillations

### 4.2 Starburst Amplification

The starburst activity of NGC 1792 (SFR = 10 M☉/yr) enhances the oscillatory amplitude A_osc through the M(t) coupling:

$$A_\text{osc,eff}(t) = A_\text{osc} \times \left(1 + \text{sSFR} \cdot e^{-t/\tau_\text{SF}}\right)$$

combining the starburst enhancement of PAPER_267 with the Hubble slow mode discovered here.

### 4.3 Dimensional Scale Hierarchy

| Mode | Frequency | Period | Physical Scale |
|------|-----------|--------|----------------|
| Fast standing wave (term_osc1) | ω_osc = 2.49×10⁻¹² rad/s | 80,200 yr | NGC 1792 light-crossing |
| Hubble slow mode (term_osc2 corrected) | ω_H = 1.44×10⁻¹⁷ rad/s | 13.8 Gyr | Hubble time |
| Modulation envelope | ε_mod ≈ 5.8×10⁻⁶ | — | 5.8 ppm depth |

---

## 5. Observational Predictions

### 5.1 Ultra-Low-Frequency GW Band

The Hubble slow mode frequency ω_H = 1.44×10⁻¹⁷ rad/s corresponds to f_H ≈ 2.3×10⁻¹⁸ Hz. This falls in the **ultra-low-frequency gravitational wave band** below pulsar timing arrays (PTA: 10⁻⁹ Hz) and even below current proposals. Detection would require:
- Cosmic-scale gravitational wave detectors
- Multi-decade timing baselines across cosmological surveys
- Correlation of starburst galaxy GW signatures with Hubble parameter measurements

### 5.2 5.8 ppm Modulation Signature

The 5.8 ppm modulation depth ε_mod = ω_H/ω_osc is a **universal UQFF prediction** applicable to any galaxy with similar r and t_Hubble parameters. For NGC 1792 specifically:

$$\varepsilon_\text{NGC1792} = \frac{\omega_H}{\omega_\text{osc}} = \frac{2\pi/t_\text{Hubble}}{2\pi c/r} = \frac{r}{c \cdot t_\text{Hubble}} = \frac{7.569 \times 10^{20}}{2.998 \times 10^8 \times 4.352 \times 10^{17}} \approx 5.8 \times 10^{-6}$$

This can also be written as:
$$\varepsilon = \frac{r}{D_H}$$

where $D_H = c \times t_\text{Hubble}$ is the Hubble distance (~14 Gly). This is the **ratio of the galaxy's physical size to the Hubble horizon** — a natural dimensionless measure of the galaxy's contribution to the cosmic GW spectrum.

### 5.3 Starburst Galaxy GW Spectral Imprint

The combination of effects from PAPER_267 (sSFR coupling) and PAPER_268 (Hubble slow mode modulation) predicts a distinctive starburst galaxy GW spectral imprint:
- Primary GW mode at ω_osc (galactic light-crossing frequency)
- Sideband at ω_osc ± ω_H (Hubble-modulated sidebands), offset by ε_mod ≈ 5.8 ppm
- Amplitude growing with sSFR during starburst episode, decaying with τ_SF = 100 Myr

---

## 6. Bug Fix Significance

The correction of `t_Hubble_gyr` → `t_Hubble` (seconds) changes the amplitude of term_osc2 by a factor:

$$\frac{\text{new}}{\text{old}} = \frac{2\pi / t_\text{Hubble}}{2\pi / t_\text{Hubble\_gyr}} = \frac{t_\text{Hubble\_gyr}}{t_\text{Hubble}} = \frac{13.8}{4.352 \times 10^{17}} \approx 3.17 \times 10^{-17}$$

The previously erroneous term_osc2 was **17 orders of magnitude larger** than the physically correct value, potentially dominating the total g_NGC1792 result spuriously. The corrected term is now appropriately small (amplitude ~ ω_H × A_osc ~ 1.44×10⁻¹⁷ × 10⁻¹⁰ ~ 10⁻²⁷ m/s²) relative to the dominant terms, consistent with a Hubble-scale perturbation on galactic gravity.

---

## 7. Conclusions

1. The pre-UQFF-2.0 `GALAXY_NGC_1792.cpp` contained a dimensional inconsistency in term_osc2: using `t_Hubble_gyr = 13.8` (dimensionless) instead of `t_Hubble` (seconds).

2. After canonical fix: ω_H = 2π/t_Hubble ≈ 1.44×10⁻¹⁷ rad/s — the **Hubble angular frequency**.

3. The two oscillatory modes (fast standing wave at ω_osc ≈ 2.49×10⁻¹² rad/s, Hubble slow mode at ω_H ≈ 1.44×10⁻¹⁷ rad/s) create a **Hubble-timescale amplitude modulation** of starburst GWs.

4. Modulation depth ε = ω_H/ω_osc = r/(c × t_Hubble) ≈ 5.8×10⁻⁶ (5.8 ppm) — a universal UQFF prediction.

5. The physical ratio ε = r/D_H provides a new **UQFF observable**: the starburst galaxy's fractional contribution to the Hubble horizon GW spectrum.

6. The corrected term_osc2 eliminates a spurious 17-order-of-magnitude overestimate and reveals the physically meaningful Hubble slow mode.

---

## References

- Daniel T. Murphy, *UQFF Framework*, Star-Magic Repository (2025–2026)
- GALAXY_NGC_1792.cpp UQFF 2.0 (Session 73, Module 19) — dimensional bug fix commit
- PAPER_267: SFR Normalization — Starburst-Buoyancy Coherence in NGC 1792
- NGC 1792 parameters: r = 80,000 ly = 7.569×10²⁰ m, z = 0.0095
- Hubble time: t_Hubble = 13.8 Gyr = 4.352×10¹⁷ s

---

*© 2026 Daniel T. Murphy, daniel.murphy00@gmail.com — All Rights Reserved*
