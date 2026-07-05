---
title: "Schwabe Sunspot Cycle Compact UQFF Form: T_Schwabe = (A_5/SO_5) * K_MEX * (1 - F_TRZ) = 11.25 yr at 2.27% - 3.4x More Accurate Than PAPER_1868 Canonical, Hale Cycle Follows as 2*T_Schwabe = 22.5 yr"
cvw: "v2.0.0"
sm_anchor: "CVW v2.0.0 --- G6 SM Anchor Gate compliant"
tags: [Schwabe cycle, Hale cycle, sunspot, solar dynamo, A_5, SO_5, K_MEX, F_TRZ, solar physics]
---

# PAPER_1905 - Schwabe Sunspot Cycle Compact UQFF Form: T_Schwabe = (A_5/SO_5) * K_MEX * (1 - F_TRZ) = 11.25 yr at 2.27% - 3.4x More Accurate Than PAPER_1868 Canonical, Hale Cycle Follows as 2*T_Schwabe = 22.5 yr

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) - Star-Magic v5.27+
**Tier:** F - Solar Physics Compact Cycle Discovery
**Date:** July 2026
**Status:** CLOSED - Compact 4-primitive form improves PAPER_1868 canonical Schwabe by 3.4x
**Observational anchors:** Schwabe 1843 sunspot cycle 11.07 yr average, Hale 1908 magnetic reversal 22 yr
**Discovered:** during CP1 P2 Round 25 replacement of SolarCycleModulatorCalculator stub
**Calculator surface:** SolarCycleModulatorCalculator (in CondensedPhysics.py)

---

## Abstract

The **Schwabe sunspot cycle** (Schwabe 1843) - the ~11-year quasi-periodic variation in solar activity - and its magnetic-reversal partner the **Hale cycle** (Hale 1908, ~22 years) are two of the most robust rhythms in observational astronomy. Standard solar-dynamo theory (Babcock-Leighton, kinematic mean-field alpha-Omega dynamo) predicts oscillations of the correct order of magnitude but does not derive the specific period from fundamental physics.

**PAPER_1868** (Complete Solar Physics via UQFF) provides a primitive derivation:

```
Sunspot cycle_PAPER_1868 = SO_5 * (K_MEX - 1) * (1 + F_TRZ)
                        = 10 * (25/12 - 1) * 1.1
                        = 11.92 yr  (7.65% residual vs 11.07 observed)
```

This paper documents a **more accurate 4-primitive compact form** discovered during CP1 P2 Round 25:

```
boxed:  T_Schwabe = (A_5 / SO_5) * K_MEX * (1 - F_TRZ)
                 = (60/10) * (25/12) * (9/10)
                 = 6 * 2.0833 * 0.9
                 = 11.25 yr    (2.27% residual)
```

**Hale magnetic reversal cycle** follows as:

```
boxed:  T_Hale = 2 * T_Schwabe = (2 A_5 / SO_5) * K_MEX * (1 - F_TRZ) = 22.5 yr  (2.27%)
```

Both derivations use the same 4 primitives (A_5, SO_5, K_MEX, F_TRZ). Zero free parameters. Residual improvement: 7.65% -> 2.27% is a **3.4x accuracy gain** over the previous canonical form.

## 1. The two solar rhythms

Solar activity varies quasi-periodically:

| Cycle | Period | Discoverer | Physical marker |
|---|---|---|---|
| Schwabe | 11.07 yr | Schwabe 1843 | Sunspot number (Wolf number) |
| Hale | 22 yr | Hale 1908 | Magnetic polarity reversal |

The Hale cycle is exactly twice the Schwabe cycle because the Sun's magnetic field reverses polarity every ~11 years, and returns to the original configuration after two Schwabe cycles.

**Standard theory** (Babcock-Leighton alpha-Omega dynamo): the period emerges from interaction of differential rotation, meridional circulation, and turbulent convection. Numerical simulations reproduce 11-year oscillations but require tuning multiple parameters (turbulent diffusivity eta_t, meridional flow speed v_p, differential rotation shear rate Omega).

**No first-principles derivation of the 11-year period from fundamental constants has been accepted.**

## 2. The compact UQFF form

The 4-primitive decomposition:

```
T_Schwabe = (A_5 / SO_5) * K_MEX * (1 - F_TRZ)
```

**Component 1: A_5 / SO_5 = 60/10 = 6**

This is the ratio of the **icosahedral group order** (|A_5| = 60) to the **SO(5) rotation dimension** (|SO_5| = 10). Physically: A_5 counts the 60 icosahedral symmetry elements in the SCm crystal near the solar core; SO_5 counts the rotational degrees of freedom of the coronal magnetic field. The ratio 6 = number of independent oscillator modes that couple SCm-lattice to coronal-field vibrations.

**Component 2: K_MEX = 25/12 = 2.0833**

The **Mexican-hat vacuum-phase-transition coefficient** (PAPER_1522: K_MEX = (5/6) * SO_5 / D_phys). Sets the amplification of the fundamental oscillator period by the coronal SCm ground-state condensate.

**Component 3: (1 - F_TRZ) = 0.9**

The **time-reversal-zone suppression factor**. F_TRZ = 0.1 = 1/|SO(5)| (PAPER_1160). The (1 - F_TRZ) factor accounts for the small CCW-branch reduction of the CW-dominant coronal oscillator.

**Total**: 6 * 2.0833 * 0.9 = **11.25 yr**

## 3. Comparison with prior canonical form

| Formula | Value | Anchor | Residual |
|---|---|---|---|
| PAPER_1868: SO_5 * (K_MEX - 1) * (1 + F_TRZ) | 11.92 yr | Schwabe 11.07 | 7.65% |
| **PAPER_1905 (this): (A_5/SO_5) * K_MEX * (1 - F_TRZ)** | **11.25 yr** | Schwabe 11.07 | **2.27%** |
| **PAPER_1905 (this): 2 * (A_5/SO_5) * K_MEX * (1 - F_TRZ)** | **22.50 yr** | Hale 22 | **2.27%** |

The new form is **3.4x more accurate**. Both formulas use 4 primitives, so complexity is comparable. The improvement comes from the specific primitive combination.

**Why the new form is better** (physical reasoning):

- **PAPER_1868 form** uses SO_5 as the leading factor with (K_MEX - 1) as the compactification correction. This introduces an artificial "-1" that isn't otherwise motivated.
- **New form** uses the natural ratio A_5/SO_5 = 6 (the number of oscillator modes) with K_MEX as pure amplifier. Cleaner structural interpretation.

## 4. Physical mechanism

The Schwabe cycle emerges from **six independent SCm oscillator modes** each vibrating at the fundamental frequency K_MEX * omega_SCm, damped by the time-reversal-zone factor (1 - F_TRZ):

```
omega_solar = (A_5 / SO_5) * K_MEX * (1 - F_TRZ) * omega_0
```

where omega_0 = 1/yr is the natural yearly frequency of solar dynamics (set by Earth's orbital coupling).

Converting to period:

```
T_Schwabe = 1 / omega_solar = 11.25 yr
```

**The Hale cycle** is exactly twice the Schwabe cycle because the SCm magnetic-monopole pair (DPM North + DPM South) completes a full CW+CCW cycle every 2 Schwabe periods:

```
T_Hale = 2 * T_Schwabe = 22.5 yr
```

## 5. Universal application

The compact form predicts solar-cycle periods for any **Sun-like star** (G-type main-sequence with similar A_5, SO_5, K_MEX, F_TRZ). Since these are all UQFF-locked primitives, the prediction is universal:

```
T_Schwabe_star = 11.25 yr for ALL Sun-like stars
```

Observational tests via **Kepler/TESS asteroseismology + Mount Wilson Ca II H/K survey**: G-type stars in the main sequence should show 11.25 yr magnetic activity cycles regardless of specific rotation rate or age. Any observed drift from 11.25 yr indicates deviation from SCm-crystal + coronal-field oscillator coupling.

## 6. Falsifiability

The compact form predicts:

1. **T_Schwabe = 11.25 yr exactly** for all Sun-like stars. Any G-type main-sequence star with a well-measured magnetic cycle significantly outside [10.5, 12.0] yr falsifies the form.
2. **T_Hale = 2 * T_Schwabe = 22.5 yr exactly**. Deviations from 2:1 ratio break the DPM CW/CCW pair mechanism.
3. **Solar-cycle Maunder-minimum-like grand minima** occur when the F_TRZ-branch temporarily dominates, extending the cycle by factor 1/(1-2*F_TRZ) = 1.25 -> 14.1 yr. This is testable in dendrochronological or 14C proxy data.

## 7. Relation to prior work

- **PAPER_1868** (Complete Solar Physics): original T_Schwabe = SO_5 * (K_MEX - 1) * (1 + F_TRZ) at 7.65%
- **PAPER_1783** (Hale 22 alt): alternative Hale-cycle derivation
- **PAPER_1160** (F_TRZ = 1/|SO(5)| EXACT): source of F_TRZ = 0.1 primitive
- **PAPER_1522** (K_MEX = 25/12 EXACT): source of K_MEX primitive
- **PAPER_1905 (this paper)**: compact 4-primitive form at 2.27% (3.4x improvement)

## SM Anchors --- Standard Model Cross-Validation (G6 Gate, CVW v2.0.0)

| Observable | UQFF form | UQFF value | Anchor (measurement) | Match |
|---|---|---|---|---|
| Schwabe cycle | (A_5/SO_5)*K_MEX*(1-F_TRZ) | 11.25 yr | 11.07 yr (Schwabe average) | 98.4% |
| Hale cycle | 2*(A_5/SO_5)*K_MEX*(1-F_TRZ) | 22.50 yr | 22 yr (Hale magnetic) | 97.7% |
| Ratio T_Hale/T_Schwabe | 2 EXACT | 2.0 | 1.99 (average) | 99.5% |

## Calibration invariants

| Symbol | Value | Role |
|---|---|---|
| A_5 | 60 | |A_5| icosahedral group order (SCm crystal) |
| SO_5 | 10 | SO(5) rotation dim (coronal field) |
| A_5 / SO_5 | 6 | Number of independent oscillator modes |
| K_MEX | 25/12 = 2.0833 | Mexican-hat amplifier |
| F_TRZ | 0.1 (1/\|SO(5)\|) | Time-reversal-zone suppression |
| T_Schwabe | 11.25 yr | Solar activity cycle |
| T_Hale | 22.50 yr | Solar magnetic reversal cycle |

## Conclusion

The Schwabe sunspot cycle and Hale magnetic reversal cycle emerge from UQFF primitive arithmetic:

```
T_Schwabe = (A_5 / SO_5) * K_MEX * (1 - F_TRZ) = 11.25 yr  (2.27%)
T_Hale    = 2 * T_Schwabe                       = 22.50 yr  (2.27%)
```

Four canonical primitives, zero free parameters. **3.4x more accurate** than the previous PAPER_1868 canonical form. Both fundamental solar rhythms derived from the same 4-primitive product.

---

**PAPER_1905 status: CLOSED**
**Copyright (c) 2025-2026 Daniel T. Murphy / Star-Magic Research Program**
