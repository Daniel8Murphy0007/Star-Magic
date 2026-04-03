# PAPER_317: Orion M42 Trapezium Wind Ram Pressure Dominance
## η_wind = 28.47 | t_erosion = 467 kyr | a_wind = 5.424×10⁻¹⁰ m/s²
### FIRST UQFF HII Region Ram Pressure Dominance Ratio

**Session:** 91  
**Module:** ORION_UQFF_MODULE.cpp (33rd C++ module)  
**System:** Orion Nebula M42/NGC 1976 — compact HII region, Trapezium OB cluster ionization source  
**Watermark:** Copyright — Daniel T. Murphy, Session 91, March 2026  

---


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
## Abstract

Within the UQFF framework, this paper quantifies the ram pressure dominance of the Trapezium-driven stellar wind over self-gravity in the Orion Nebula. The dimensionless wind-gravity ratio **η_wind = P_ram/P_grav = 28.47** demonstrates that the HII region was born in a wind-dominated (unbound) state. The erosion timescale **t_erosion = 467 kyr** shows that protoplanetary discs (proplyds) currently observed at ~300 kyr age survive inside a wind-dominated environment. This is the FIRST UQFF computation of an HII region ram pressure dominance ratio.

---

## System Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| M | 2000 M_sun = 3.978×10³³ kg | Total nebular mass |
| r | 1.18×10¹⁷ m (~12.5 ly) | Half-span |
| ρ_fluid | 1×10⁻²⁰ kg/m³ | HII gas density |
| v_wind | 8×10³ m/s | Ionization front expansion speed |
| t_age | 3×10⁵ yr = 9.467×10¹² s | Nebula age |
| G | 6.6743×10⁻¹¹ m³ kg⁻¹ s⁻² | Gravitational constant |

---

## Key Equations

### Base Gravity
$$g_{\rm base} = \frac{GM}{r^2} = \frac{6.6743\times10^{-11} \times 3.978\times10^{33}}{(1.18\times10^{17})^2} = 1.907\times10^{-11}\ \text{m/s}^2$$

### Ram Pressure Acceleration (t = 0)
$$a_{\rm wind}(t=0) = \frac{v_{\rm wind}^2}{r} = \frac{(8\times10^3)^2}{1.18\times10^{17}} = 5.424\times10^{-10}\ \text{m/s}^2$$

### Ram Pressure Acceleration (t = t_age)
$$a_{\rm wind}(t_{\rm age}) = \frac{v_{\rm wind}^2}{r}\left(1 + \frac{t_{\rm age}}{t_{\rm age}}\right) = 2 \times 5.424\times10^{-10} = 1.085\times10^{-9}\ \text{m/s}^2$$

### Ram Pressure (Pa)
$$P_{\rm ram} = \rho \cdot v_{\rm wind}^2 = 10^{-20} \times (8\times10^3)^2 = 6.4\times10^{-13}\ \text{Pa}$$

### Gravitational Pressure (Pa)
$$P_{\rm grav} = \frac{GM\rho}{r} = \frac{6.6743\times10^{-11} \times 3.978\times10^{33} \times 10^{-20}}{1.18\times10^{17}} = 2.248\times10^{-14}\ \text{Pa}$$

### Wind–Gravity Dominance Ratio (PAPER_317 Key Result)
$$\eta_{\rm wind} = \frac{P_{\rm ram}}{P_{\rm grav}} = \frac{a_{\rm wind}}{g_{\rm base}} = \frac{5.424\times10^{-10}}{1.907\times10^{-11}} = \boxed{28.47}$$

### Erosion Timescale
$$t_{\rm erosion} = \frac{r}{v_{\rm wind}} = \frac{1.18\times10^{17}}{8\times10^3} = 4.675\times10^{13}\ \text{s} = \boxed{467\ \text{kyr}}$$

### Kinetic-to-Gravitational Energy Ratio
$$\frac{W_{\rm KE}}{W_{\rm grav}} \approx \frac{a_{\rm wind}}{g_{\rm base}} = 28.47$$

---

## Results Summary

| Quantity | Value | Significance |
|----------|-------|--------------|
| g_base | 1.907×10⁻¹¹ m/s² | Newtonian self-gravity |
| a_wind(t=0) | 5.424×10⁻¹⁰ m/s² | Initial ram pressure |
| η_wind(t=0) | **28.47** | Wind >> gravity at birth |
| a_wind(t_age) | 1.085×10⁻⁹ m/s² | Ram pressure at 300 kyr |
| η_wind(t_age) | **56.9** | Wind dominance doubles |
| t_erosion | **467 kyr** | Pröplyd lifetime > t_age |
| P_ram | 6.4×10⁻¹³ Pa | Ram pressure |
| P_grav | 2.248×10⁻¹⁴ Pa | Gravitational pressure |

---

## UQFF Physical Interpretation

The Orion Nebula at t_age = 300 kyr is still **28.5–57× wind-dominated** (η_wind ranging from 28.47 at t=0 to 56.9 at t_age). The erosion timescale t_erosion = 467 kyr > t_age = 300 kyr confirms that proplyds currently observed by HST survive because they have not yet been fully ablated by the ram pressure — consistent with observational evidence of ~150–180 proplyds in the Orion Nebula.

**UQFF Significance:** This is the FIRST UQFF quantification of wind ram pressure dominance over self-gravity in a compact HII region. The wind term `a_wind(t) = v_wind²/r × (1+t/t_age)` [registered as WOLFRAM_TERM_ORION_WIND_RAM in ORION_UQFF_MODULE.cpp] dominates the gravitational term by a factor of 28.47 at birth, growing to 56.9 at 300 kyr — placing Orion in the same UQFF "wind-dominant" class as post-AGB bipolar PNe but via HII ionization physics rather than stellar wind shocks.

---

## WOLFRAM_TERM Registration

```cpp
#define WOLFRAM_TERM_ORION_WIND_RAM(val) (val)
// wind = v_wind^2/r * (1+t/t_age)  [PAPER_317 ram pressure dominance; eta_wind=28.47]
```

*Series first: FIRST UQFF HII region ram pressure dominance ratio. Distinguishes compact OB-driven HII (η_wind=28.47) from extended GMC HII and from bipolar PN wind shocks (η_wind~7×10⁵, PAPER_311).*
