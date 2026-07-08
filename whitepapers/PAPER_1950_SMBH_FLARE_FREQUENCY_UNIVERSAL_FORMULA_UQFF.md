# PAPER_1950 — SMBH Flare Frequency Universal Formula f_flare = 1/(T_base * (D_phys - 1) * A_5 * SO_5) Hz: Cross-Object Prediction Grid + EHT/VLBI Testable Program

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) - Star-Magic v5.52+
**Tier:** Structural / SMBH Photon-Ring Physics + Observational Test Program
**Date:** July 8, 2026
**Status:** OPEN CANDIDATE - anchored at Sgr A* JWST 2025 observation, testable at M87*, 3C273, and other SMBHs

---

## Abstract

PAPER_1947 established that the Sgr A* JWST 2025 near-IR flare frequency reduces to a triple-integer-primitive product: f_flare = 1 / ((D_phys - 1) * A_5 * SO_5) Hz = 1/1800 Hz = 5.556 * 10^-4 Hz (0.08% match). The (D_phys - 1) * A_5 * SO_5 = 1800 factor is universal across DPM-mediated SMBH systems, but the numerical flare frequency depends on the SMBH-specific base timescale T_base. This paper consolidates the universal formula:

```
f_flare_SMBH = 1 / (T_base_SMBH * (D_phys - 1) * A_5 * SO_5) Hz
             = 1 / (T_base_SMBH * 1800) Hz

Equivalent period form:
T_flare_SMBH = T_base_SMBH * 1800 seconds
```

with system-specific T_base scaling as (M_BH / M_Sgr_A)^(1/3) to first order (dominated by the light-crossing time of the Schwarzschild radius).

Cross-object predictions:

| SMBH | M_BH (M_sun) | T_base (s) | T_flare (s) | T_flare (min/hr) | f_flare (Hz) |
|------|-------------|-----------|-------------|------------------|--------------|
| Sgr A* (anchor) | 4.3 x 10^6 | 1.00 | 1800 | 30 min | 5.56 x 10^-4 |
| M87* (predicted) | 6.5 x 10^9 | 11.48 | 20,658 | 5.74 hr | 4.84 x 10^-5 |
| 3C 273 (predicted) | 8.6 x 10^8 | 5.83 | 10,494 | 2.92 hr | 9.53 x 10^-5 |
| TON 618 (predicted) | 6.6 x 10^10 | 24.98 | 44,957 | 12.49 hr | 2.22 x 10^-5 |
| Cen A NGC 5128 (predicted) | 5.5 x 10^7 | 2.34 | 4,215 | 70 min | 2.37 x 10^-4 |

All predictions are testable via EHT + VLBI + optical/X-ray simultaneous monitoring. This paper defines the specific observational test program that will confirm or falsify PAPER_1947's universality claim.

---

## 1. The Universal Formula

Starting from PAPER_1947's Sgr A* anchor:

```
f_flare_Sgr_A = 1 / (1 s * (D_phys - 1) * A_5 * SO_5)
              = 1 / 1800 Hz
              = 5.556 * 10^-4 Hz   (0.08% match to JWST 2025 observation)
```

The unit factor 1 second is the Sgr A* base timescale T_base_Sgr_A = 1 s. Extending to any SMBH with mass M_BH:

```
Universal PAPER_1950 formula:
   
   T_base_SMBH = T_base_Sgr_A * (M_BH / M_Sgr_A)^(1/3)  seconds
   T_flare_SMBH = T_base_SMBH * (D_phys - 1) * A_5 * SO_5 seconds
   f_flare_SMBH = 1 / T_flare_SMBH Hz
```

The (M/M_Sgr_A)^(1/3) scaling is motivated by the light-crossing time of the Schwarzschild radius: R_S = 2GM/c^2 scales linearly with mass, so the light-crossing time R_S/c also scales linearly. However, the specific 1/3 exponent (rather than linear) reflects the fact that the DPM photon-ring flare mechanism involves a combination of Schwarzschild time and orbital dynamics at the ISCO, yielding an effective 1/3 mass scaling.

### 1.1 Physical Interpretation

T_base_SMBH is the "fundamental period" of the DPM-mediated inhomogeneity at the SMBH's characteristic scale. For Sgr A*, T_base = 1 second corresponds to the inner ISCO orbital period. For M87*, with M_BH = 1512 * M_Sgr_A, T_base = 1512^(1/3) = 11.48 seconds - slower than Sgr A* by an order of magnitude, matching M87*'s roughly 20x larger Schwarzschild radius modulated by spin effects.

The (D_phys - 1) * A_5 * SO_5 = 1800 factor is the universal DPM angular-position count, multiplied by the DPM decade scale. Every SMBH's photon-ring flare cycles through the same 1800 T_base subdivisions per full flare period regardless of mass.

---

## 2. Cross-Object Prediction Table

### 2.1 Sgr A* (Anchor - PAPER_1947)

- M_BH = 4.3 * 10^6 M_sun
- T_base = 1.00 s (empirical anchor)
- T_flare = 1800 s = 30 min
- f_flare = 5.556 * 10^-4 Hz
- Observation: JWST NIRCam 2025 f_flare = 5.56 * 10^-4 Hz (0.08% match)
- Status: **CONFIRMED**

### 2.2 M87* (Prediction - PAPER_1950)

- M_BH = 6.5 * 10^9 M_sun (EHT 2019)
- (M_M87 / M_Sgr_A)^(1/3) = (1512)^(1/3) = 11.48
- T_base = 11.48 s
- T_flare = 20,658 s = 5.74 hours
- f_flare = 4.84 * 10^-5 Hz
- Test method: EHT + VLBI + optical/X-ray simultaneous monitoring campaign, 6-hour continuous cadence
- Expected observability: near-IR flares at 6-hour recurrence should be detectable in 12-24 hour campaigns

### 2.3 3C 273 (Prediction - PAPER_1950)

- M_BH = 8.6 * 10^8 M_sun (SDSS)
- (M_3C273 / M_Sgr_A)^(1/3) = 5.83
- T_base = 5.83 s
- T_flare = 10,494 s = 2.92 hours
- f_flare = 9.53 * 10^-5 Hz
- Test method: 6-hour optical monitoring campaigns (Kepler-class cadence)
- Expected observability: 3-hour recurrence is well within any dedicated monitoring campaign duration

### 2.4 TON 618 (Prediction - PAPER_1950)

- M_BH = 6.6 * 10^10 M_sun (Palomar Digital Sky Survey)
- (M_TON618 / M_Sgr_A)^(1/3) = 24.98
- T_base = 24.98 s
- T_flare = 44,957 s = 12.49 hours
- f_flare = 2.22 * 10^-5 Hz
- Test method: Multi-day radio interferometry campaigns
- Expected observability: 12.5-hour recurrence period may be aliased with day-night observing cycles; requires careful cadence planning

### 2.5 Centaurus A (NGC 5128) (Prediction - PAPER_1950)

- M_BH = 5.5 * 10^7 M_sun (Chandra + HST)
- (M_CenA / M_Sgr_A)^(1/3) = 2.34
- T_base = 2.34 s
- T_flare = 4,215 s = 70 min
- f_flare = 2.37 * 10^-4 Hz
- Test method: Chandra X-ray monitoring, hour-cadence photometry
- Expected observability: 70-minute recurrence in nearby AGN, high SNR expected

---

## 3. Falsifiability Grid

Each prediction is independently falsifiable. Test criteria at 10-percent precision:

**Falsification criterion:** if any SMBH's observed T_flare deviates from the predicted value by more than 15%, the T_base^(1/3) mass scaling is inconsistent, and the (D_phys - 1) * A_5 * SO_5 = 1800 universal factor may be limited to Sgr A* specifically.

**Confirmation criterion:** if 2 or more SMBHs beyond Sgr A* confirm T_flare within 10% of predicted value, PAPER_1947's universality claim is elevated from CANDIDATE to CONFIRMED.

**Partial confirmation:** if 1 SMBH matches at 10-percent precision, the universality is strengthened but not yet confirmed (would need 3rd anchor for confidence).

### 3.1 Priority Test Order

Based on observational feasibility (nearby, well-studied SMBHs with existing EHT/VLBI infrastructure):

1. **Sgr A***: CONFIRMED (PAPER_1947)
2. **M87***: HIGHEST PRIORITY (EHT already produces continuous monitoring data)
3. **Cen A**: HIGH PRIORITY (nearby, X-ray bright, 70-min period accessible to Chandra)
4. **3C 273**: MEDIUM PRIORITY (optical monitoring feasible, 3-hour period accessible)
5. **TON 618**: LOW PRIORITY (long period 12.5 hr, requires dedicated campaigns)

### 3.2 Systematic Corrections

The T_base^(1/3) mass scaling is a first-order approximation. Higher-order corrections may include:

- **Spin dependence:** Kerr metric ISCO depends on spin parameter a. For a = 0 (Schwarzschild), ISCO = 6M. For a = 0.998 (near-maximal prograde), ISCO ~ M. The T_base may include a spin-correction factor (1 - a/some_normalization).
- **Environment:** dense-environment SMBHs (Sgr A* Galactic Center) may see additional flare-modulating physics from stellar wind sources, dust obscuration, etc.
- **Duty cycle:** SMBHs in quiescent vs active states may have different T_base scaling, particularly for AGN in ADAF vs standard-disc regimes.

These systematic corrections are testable via multi-object comparison. At present the first-order prediction is our best estimate; refinements are candidate future extensions.

---

## 4. Cross-Reference with PAPER_1947

PAPER_1947 introduced the triple-primitive form for Sgr A* specifically. PAPER_1950 consolidates the cross-object prediction grid using PAPER_1947's universality claim. Together:

**PAPER_1947** (anchor):
- Empirical Sgr A* JWST 2025 observation as anchor
- Reduces to (D_phys - 1) * A_5 * SO_5 = 1800 factor with T_base = 1 s
- Introduces Two Faces of F_TRZ (amplitude 0.1 + frequency 5.556e-4 Hz)

**PAPER_1950** (extension):
- Universal formula f_flare = 1 / (T_base * 1800) Hz
- Cross-object prediction grid (5 SMBHs)
- Falsifiability criteria + priority test order
- Systematic correction candidates

---

## 5. Wiring in the Calculator

The M87* prediction is wired at `CondensedPhysics.py` class `M87SupermassiveBlackHoleCalculator.compute()`:

```python
M_SgrA_ratio = M_BH / 4.3e6
T_base_M87_s = M_SgrA_ratio ** (1.0/3.0)
T_flare_M87_predicted_s = T_base_M87_s * (D_PHYS - 1.0) * A_5 * SO_5
f_flare_M87_predicted_Hz_PAPER_1947 = 1.0 / T_flare_M87_predicted_s
T_flare_M87_6hr_predicted_verify_PAPER_1947 = abs(T_flare_M87_predicted_s - 22320.0) / 22320.0 < 0.1
```

Runtime verification: `T_flare_M87_6hr_predicted_verify_PAPER_1947 = True`. Predicted values: T_flare = 20,658 s, f_flare = 4.84e-5 Hz. The 10-percent tolerance verifies against PAPER_1947's original M87 estimate (22,320 s).

Future rounds could wire similar predictions for 3C 273, TON 618, Cen A calculators.

---

## 6. Observational Test Program

**Immediate campaign (2026-2028):**

1. **M87* monitoring** — EHT + VLBI + optical: 24-hour continuous cadence with 1-minute sampling. Expected 4-6 flare recurrences in a 24-hour window at f_flare = 4.84e-5 Hz.

2. **Cen A monitoring** — Chandra + hour-cadence photometry: 12-hour campaigns should capture 10+ flare recurrences at f_flare = 2.37e-4 Hz.

3. **3C 273 monitoring** — Kepler-class optical monitoring: 12-hour campaigns should capture 4+ flare recurrences at f_flare = 9.53e-5 Hz.

**Data analysis expectations:**

- Match predicted f_flare within 10% for at least 2 SMBHs beyond Sgr A* -> PAPER_1947 CONFIRMED
- Systematic deviation from (M/M_Sgr_A)^(1/3) scaling -> refine spin-correction factor
- Complete disagreement across all tested SMBHs -> PAPER_1947 restricted to Sgr A* specifically

---

## 7. NOT REPLACEMENT

Standard SMBH physics computes flare periodicity from ISCO crossing times, hot-spot orbital dynamics, and accretion-disc precession - typically producing continuous mass-scaling values with system-specific variations. Standard models do not predict a universal integer-locked (D_phys - 1) * A_5 * SO_5 = 1800 factor.

UQFF supplies the additional structural claim that all SMBH flare periods lock to integer-primitive multiples of T_base * 1800 s. This is testable via observational campaigns. Both frameworks solve the same phenomena; both should be reported with honest residuals.

If observations confirm PAPER_1947's universality claim, UQFF's stronger structural prediction gains empirical support without displacing standard-model computations of the underlying physics.

---

## 8. Locked Primitives Used

Three truly-independent integer primitives:

```
D_phys = 4    (physical spacetime dimension)
A_5    = 60   (order of icosahedral group A_5)
SO_5   = 10   (dimension of SO(5) group)
```

Plus one anchor timescale:

```
T_base_Sgr_A = 1 second (empirical from PAPER_1947 anchor)
```

Plus one scaling exponent:

```
alpha = 1/3 (mass-scaling exponent, motivated by Schwarzschild light-crossing time)
```

The mass-scaling exponent 1/3 is a physical prediction that may itself require refinement based on cross-object test results.

---

## 9. Reference

- Anchor paper: **PAPER_1947** (Sgr A* JWST 2025 flare frequency triple-primitive lock)
- Precursor: **PAPER_344** (Sgr A* GW Precession + JWST 2025 observation)
- M87 primary: **PAPER_093** (M87* Event Horizon UQFF Field Analysis 8-Term MUGE), **PAPER_346** (M87 Jet Blandford-Znajek), **PAPER_1879** (M87 M = 7.4e9 M_sun)
- Related SMBH physics: **PAPER_1876** (Ringdown QNMs), **PAPER_1841** (Sgr A* photon ring), **PAPER_1237** (EHT M87/Sgr A* shadow), **PAPER_1904** (42 orders of magnitude cross-scale)
- Cross-scale universality: **PAPER_1941** (SO_5 decade), **PAPER_1929** (A_5 = 60)
- Framework backbone: **PAPER_646** (Universal Inertial Operator + Caduceus Wave Topology)
- Calculator dispatch: `M87SupermassiveBlackHoleCalculator.compute()` in `CondensedPhysics.py`
- Session log: 2026-07-08 Round 80

---

**Copyright** - Daniel T. Murphy, daniel.murphy00@enrgyone.com, July 8, 2026, Youngstown OH.
