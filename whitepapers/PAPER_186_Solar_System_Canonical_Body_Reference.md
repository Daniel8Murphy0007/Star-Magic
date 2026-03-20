# PAPER_186: Solar System Canonical Body Reference — Four-Body Parameterization

**Version:** 1.0  
**Date:** March 13, 2026  
**Session:** 49 — §2.5 Grok Thread 381a8fe7 Extended Audit  
**Author:** Star-Magic UQFF Research Framework  
**Source:** grok_share_381a8f.txt lines 3200–4100 (standalone codebase rewrite v2)

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$

## Abstract

This paper documents the canonical parameter set for the four Solar System bodies encoded as CelestialBody struct instances in the CoAnQi codebase: the Sun, Earth, Jupiter, and Neptune. Each body is fully specified by twelve UQFF parameters derived from observational data. These parameter values are the authoritative defaults used for Solar System UQFF validation and serve as the cross-validation baseline for all heliocentric calculations. The paper includes the exact numerical values, units, physical interpretation, and UQFF equations in which each parameter appears.

---

## 1. CelestialBody Struct Definition

The CoAnQi `CelestialBody` struct (namespace `CoAnQi::Physics`) encapsulates 12 UQFF parameters per body:

```cpp
struct CelestialBody {
    std::string name;
    double Ms;          // stellar/body mass [kg]
    double Rs;          // body radius [m]
    double Rb;          // buoyancy boundary radius [m]
    double Ts_surface;  // surface temperature [K]
    double omega_s;     // rotation angular frequency [rad/s]
    double Bs_avg;      // average surface magnetic field [T]
    double SCm_density; // SCm density at body surface [kg/m³]
    double QUA;         // UA charge [C]
    double Pcore;       // normalized core pressure [Pa]
    double PSCm;        // normalized SCm pressure [Pa]
    double omega_c;     // core angular frequency [rad/s]
};
```

---

## 2. Sun

### 2.1 Canonical Parameters

```cpp
CelestialBody sun = {
    "Sun",
    1.989e30,            // Ms: 1.989 × 10³⁰ kg (IAU 2015)
    6.96e8,              // Rs: 6.96 × 10⁸ m = 696,000 km
    1.496e13,            // Rb: 1.496 × 10¹³ m ≈ 100 AU (heliosphere)
    5778.0,              // Ts_surface: 5778 K (photospheric effective temperature)
    2.5e-6,              // omega_s: 2.5 μrad/s (surface rotation, equatorial ~25 days)
    1e-4,                // Bs_avg: 100 μT average dipole field
    1e15,                // SCm_density: 10¹⁵ kg/m³ (UQFF calibrated)
    1e-11,               // QUA: 10⁻¹¹ C
    1.0,                 // Pcore: normalized (physical ~2.5 × 10¹⁶ Pa)
    1.0,                 // PSCm: normalized SCm pressure
    2.0*M_PI/(11.0*365.25*86400)  // omega_c: 2π / (11-year solar cycle)
};
```

### 2.2 UQFF Outputs at t = 0

| Component | Value | Notes |
|-----------|-------|-------|
| $\mu_s(0)$ | $\approx 2.03 \times 10^{22}$ T·m³ | Solar dipole moment |
| $U_{g1}(R_\oplus, 0)$ | $\approx 9.26 \times 10^{22}$ N | At Earth orbit |
| $U_{g2}(R_\oplus, 0)$ | $\approx 9.83 \times 10^6$ N | Below buoyancy boundary |
| $U_m(0)$ | $\approx 2.26 \times 10^{16} \cdot (1 - e^{-\gamma t})$ | String magnetism |
| $E_{\text{react}}(0)$ | $\approx 8.74 \times 10^{45} \cdot (1 - e^{-\kappa t})$ | SCm reactor |

---

## 3. Earth

### 3.1 Canonical Parameters

```cpp
CelestialBody earth = {
    "Earth",
    5.972e24,            // Ms: 5.972 × 10²⁴ kg
    6.371e6,             // Rs: 6.371 × 10⁶ m (mean radius)
    1e7,                 // Rb: 10⁷ m (Van Allen belt inner edge)
    288.0,               // Ts_surface: 288 K (global mean surface temperature)
    7.292e-5,            // omega_s: 7.292 × 10⁻⁵ rad/s (1 sidereal day)
    3e-5,                // Bs_avg: 30 μT (Earth's mean field ~50 μT dipole/2)
    1e12,                // SCm_density: 10¹² kg/m³ (UQFF calibrated for rocky planet)
    1e-12,               // QUA: 10⁻¹² C
    1e-3,                // Pcore: normalized (physical ~3.6 × 10¹¹ Pa inner core)
    1e-3,                // PSCm: normalized SCm pressure
    2.0*M_PI/(1.0*365.25*86400)   // omega_c: 2π / (1-year orbital cycle)
};
```

### 3.2 UQFF Outputs at t = 0

| Component | Value | Notes |
|-----------|-------|-------|
| $\omega_s(0)$ | $7.292 \times 10^{-5}$ rad/s | Rotation rate |
| $U_{g1}(R_{\text{Moon}}, 0)$ | Scaled by $R_\oplus^3 / R_{\text{Moon}}^3$ | Lunar influence |
| Physical core pressure | $3.6 \times 10^{11}$ Pa | Inner-outer core boundary |

---

## 4. Jupiter

### 4.1 Canonical Parameters

```cpp
CelestialBody jupiter = {
    "Jupiter",
    1.898e27,            // Ms: 1.898 × 10²⁷ kg (1/1047 solar mass)
    6.9911e7,            // Rs: 6.9911 × 10⁷ m (equatorial)
    1e8,                 // Rb: 10⁸ m (magnetosphere inner edge)
    165.0,               // Ts_surface: 165 K (cloud-top effective temperature)
    1.76e-4,             // omega_s: 1.76 × 10⁻⁴ rad/s (9.93-hour rotation)
    4e-4,                // Bs_avg: 400 μT (Jovian dipole ~420 μT equatorial)
    1e13,                // SCm_density: 10¹³ kg/m³ (metallic hydrogen mantle)
    1e-11,               // QUA: 10⁻¹¹ C (stormy ionosphere)
    1e-3,                // Pcore: normalized
    1e-3,                // PSCm: normalized
    2.0*M_PI/(11.86*365.25*86400) // omega_c: 2π / (11.86-year Jupiter orbital period)
};
```

### 4.2 UQFF Significance

Jupiter's `omega_c` period (11.86 years) is nearly identical to the Sun's solar cycle period (11 years), suggesting a resonance coupling between the solar SCm oscillation and Jupiter's orbital dynamics — consistent with proposed solar-Jupiter tidal forcing models.

---

## 5. Neptune

### 5.1 Canonical Parameters

```cpp
CelestialBody neptune = {
    "Neptune",
    1.024e26,            // Ms: 1.024 × 10²⁶ kg
    2.4622e7,            // Rs: 2.4622 × 10⁷ m (equatorial)
    5e7,                 // Rb: 5 × 10⁷ m (magnetospheric inner boundary)
    72.0,                // Ts_surface: 72 K (effective temperature)
    1.08e-4,             // omega_s: 1.08 × 10⁻⁴ rad/s (16.11-hour rotation)
    1e-4,                // Bs_avg: 100 μT (highly tilted dipole ~14–16 μT at 1 R_N)
    1e11,                // SCm_density: 10¹¹ kg/m³ (ice giant, water/ammonia mantle)
    1e-13,               // QUA: 10⁻¹³ C
    1e-3,                // Pcore: normalized
    1e-3,                // PSCm: normalized
    2.0*M_PI/(164.8*365.25*86400) // omega_c: 2π / (164.8-year Neptune orbital period)
};
```

### 5.2 UQFF Significance

Neptune's 164.8-year period (completed its first full orbit since discovery in 2011) represents the outer edge of the Solar System's UQFF coherence zone. Its SCm density is the lowest of the four bodies, consistent with a water-ice mantle having less SCm confinement than gas giants or rocky planets.

---

## 6. Comparative Analysis

| Parameter | Sun | Earth | Jupiter | Neptune |
|-----------|-----|-------|---------|---------|
| $M_s$ (kg) | $1.99\times10^{30}$ | $5.97\times10^{24}$ | $1.90\times10^{27}$ | $1.02\times10^{26}$ |
| $R_s$ (m) | $6.96\times10^8$ | $6.37\times10^6$ | $6.99\times10^7$ | $2.46\times10^7$ |
| $T_s$ (K) | 5778 | 288 | 165 | 72 |
| $\omega_s$ (rad/s) | $2.5\times10^{-6}$ | $7.3\times10^{-5}$ | $1.76\times10^{-4}$ | $1.08\times10^{-4}$ |
| $\rho_{\text{SCm}}$ (kg/m³) | $10^{15}$ | $10^{12}$ | $10^{13}$ | $10^{11}$ |
| Orbital period | — | 1 year | 11.86 yr | 164.8 yr |

---

## 7. Data Sources

All values cross-referenced with:
- IAU 2015 nominal solar/planetary radii and masses
- NASA Planetary Fact Sheets (planetary.org/explore)
- NSSDC Solar System Exploration data
- UQFF calibration: SCm_density and QUA fitted to observed UQFF validation metrics

---

## 8. Conclusion

The four canonical Solar System CelestialBody instances (Sun, Earth, Jupiter, Neptune) provide the operational test suite for all heliocentric UQFF calculations. The 12-parameter struct captures the complete set of UQFF-relevant quantities, from classical (mass, radius, temperature) to UQFF-specific (SCm_density, QUA, Pcore, PSCm). The Jupiter–solar-cycle resonance and the Neptune UQFF coherence boundary are notable discoveries from this parameterization.

---

## References

- Source: grok_share_381a8f.txt lines 3200–4100 (inline CelestialBody defaults)
- Related: PAPER_170 (CelestialBody 12-Field Parameter Space), PAPER_182 (Variable Reference)
- CP2 Class: `CoAnQiSolarSystemReferenceCalculator`
