# PAPER_319: Compact HII SFR Gravitational Binding Phase Transition
## t_cross = 67,730 yr | sSFR = 5×10⁻⁴ yr⁻¹ (50× Lagoon) | m_factor(t_age) = 151
### FIRST UQFF Compact HII Region SFR Runaway Gravitational Binding Phase Transition

**Session:** 91  
**Module:** ORION_UQFF_MODULE.cpp (33rd C++ module)  
**System:** Orion Nebula M42/NGC 1976 — compact HII region with active star formation  
**Watermark:** Copyright — Daniel T. Murphy, Session 91, March 2026  

---


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
## Abstract

This paper derives the SFR-driven gravitational binding phase transition for the Orion Nebula. With SFR = 1 M_sun/yr acting on an initial mass M = 2000 M_sun, the specific star formation rate **sSFR = 5×10⁻⁴ yr⁻¹ is 50× that of the Lagoon Nebula** (1×10⁻⁵ yr⁻¹, PAPER_305). The system is born wind-dominated (unbound), but as SFR continuously adds mass, the effective gravitational acceleration amplified by the SFR mass factor m_factor(t) = 1 + SFR×t_yr/M_sun_count crosses the growing wind ram pressure at **t_cross = 67,730 yr**, transitioning the nebula from unbound to gravitationally bound. By t_age = 300 kyr, binding_ratio = g_SFR/a_wind = 2.654. This is the FIRST UQFF compact HII SFR runaway gravitational binding phase transition.

---

## System Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| M | 2000 M_sun = 3.978×10³³ kg | Initial mass |
| SFR | 1 M_sun/yr = 6.303×10²² kg/s | Star formation rate |
| sSFR | 5×10⁻⁴ yr⁻¹ | Specific SFR = SFR/M |
| t_age | 3×10⁵ yr | Nebula age |
| v_wind | 8×10³ m/s | HII ionization front |
| g_base | 1.907×10⁻¹¹ m/s² | Base gravity (PAPER_317) |
| a_wind(t=0) | 5.424×10⁻¹⁰ m/s² | Initial wind acceleration |

---

## Key Equations

### SFR Mass Factor
$$M_{\rm sf}(t) = \frac{\rm SFR_{yr} \times t_{\rm yr}}{M_{\rm sun,count}}, \qquad m_{\rm factor}(t) = 1 + M_{\rm sf}(t)$$

### SFR-Amplified Gravitational Acceleration
$$g_{\rm SFR}(t) = g_{\rm base} \times m_{\rm factor}(t) = g_{\rm base}\left(1 + \frac{{\rm SFR_{yr}} \cdot t_{\rm yr}}{M_{\rm sun,count}}\right)$$

### Wind Acceleration (time-evolving)
$$a_{\rm wind}(t) = \frac{v_{\rm wind}^2}{r}\left(1 + \frac{t}{t_{\rm age}}\right)$$

### Crossover Time t_cross (PAPER_319 Key Result)
Setting g_SFR(t) = a_wind(t):
$$g_{\rm base}\left(1 + {\rm sSFR} \cdot t\right) = a_{\rm wind,0}\left(1 + \frac{t}{t_{\rm age,yr}}\right)$$

$$t_{\rm cross} = \frac{a_{\rm wind,0} - g_{\rm base}}{g_{\rm base} \cdot {\rm sSFR} - a_{\rm wind,0}/t_{\rm age,yr}} = \boxed{67{,}730\ \text{yr}}$$

### Specific SFR
$${\rm sSFR} = \frac{\rm SFR_{yr}}{M_{\rm sun,count}} = \frac{1}{2000} = 5\times10^{-4}\ \text{yr}^{-1} = 50\times {\rm Lagoon}$$

### Gas Depletion Timescale
$$t_{\rm consume} = \frac{M_{\rm sun,count}}{\rm SFR_{yr}} = \frac{2000}{1} = 2000\ \text{yr} \quad \text{(without OMC-1 replenishment)}$$

---

## Results Summary

| Quantity | Value | Significance |
|----------|-------|--------------|
| sSFR | **5×10⁻⁴ yr⁻¹** | 50× Lagoon Nebula |
| t_cross | **67,730 yr** | Unbound → bound transition |
| m_factor(t_age=300 kyr) | **151** | SFR amplification |
| g_SFR(t_age) | 2.878×10⁻⁹ m/s² | 151× g_base |
| a_wind(t_age) | 1.085×10⁻⁹ m/s² | |
| binding_ratio(t_age) | **2.654** | Gravitationally bound |
| m_factor(1 Myr) | **501** | |
| binding_ratio(1 Myr) | **4.069** | Increasingly bound |
| t_consume | **2000 yr** | Gas depletion (no replenishment) |

---

## Phase Transition Diagram

```
t=0        t_cross=67.7 kyr    t_age=300 kyr     t=1 Myr
|               |                    |               |
UNBOUND         TRANSITION           BOUND 2.65×     BOUND 4.07×
a_wind>g_SFR    a_wind = g_SFR       g_SFR>a_wind    g_SFR>>a_wind
η=28.47 wind    ─── crossover ───    η_B=2.654       η_B=4.069
```

---

## Comparison: UQFF SFR Class

| Module | SFR (M_sun/yr) | M (M_sun) | sSFR (yr⁻¹) | PAPER |
|--------|---------------|-----------|-------------|-------|
| Lagoon | 0.1 | 10,000 | 1×10⁻⁵ | PAPER_305 |
| **Orion** | **1.0** | **2000** | **5×10⁻⁴** | **PAPER_319** |
| M16 (Eagle) | ~0.01× | 1200 | — | PAPER_284 |

Orion sSFR is 50× Lagoon — UQFF identifies this as the "ultra-compact HII" class with rapid gravitational binding crossover. Lagoon has t_consume = 100 kyr; Orion has t_consume = **2000 yr**, the shortest gas depletion time in the UQFF series, sustained only by continuous OMC-1 giant molecular cloud inflow.

---

## UQFF Physical Interpretation

The phase transition at t_cross = 67,730 yr is a structural boundary in Orion's evolution:

- **Before t_cross:** The system is wind-dominated (η_wind > 1). UV photoionization and ram pressure drive a champagne flow outward; newly formed stars are subject to wind erosion.
- **After t_cross:** SFR has added sufficient mass that g_SFR > a_wind. The system becomes self-gravitating with respect to its own stellar formation. The cluster proceeds to form stars under its own gravity — consistent with the Orion OB1 association framework where the cluster is now 300 kyr old and gravitationally bound.

Within ORION_UQFF_MODULE.cpp, the SFR mass factor m_factor(t) enters the base gravity term via `G*M*(1+M_sf(t))/r^2`, registered as `WOLFRAM_TERM_ORION_BASE` and `WOLFRAM_TERM_ORION_SFR_BINDING`.

---

## WOLFRAM_TERM Registration

```cpp
#define WOLFRAM_TERM_ORION_BASE(val)        (val)
// g_base = (G*M*(1+M_sf(t))/r^2)*(1+Hz*t)*(1-B/B_crit)*(1+f_TRZ)
// M_sf(t) = SFR_yr*t_yr / M_sun_count  [PAPER_319 SFR amplification]

#define WOLFRAM_TERM_ORION_SFR_BINDING(val) (val)
// wind = v_wind^2/r * (1+t/t_age)  [crosses SFR gravity at t_cross=67730 yr]
```

*Series first: FIRST UQFF compact HII SFR runaway gravitational binding phase transition. Establishes sSFR as a new UQFF classification axis for HII regions: compact (Orion, sSFR=5×10⁻⁴ yr⁻¹) vs extended (Lagoon, sSFR=1×10⁻⁵ yr⁻¹).*
