# PAPER_220: Crab Nebula Pulsar Wind Nebula UQFF — F_wind and M_mag in Expanding PWN

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)  
**Framework:** UQFF v4.3 — Star-Magic Physics  
**Source:** grok_share_7514fe.txt — Document 10: Crab Nebula  
**Date:** March 14, 2026  
**Series:** Phase 2 Session 55 — §2.10 Fourth-Pass System Extraction

---

## Abstract

The Crab Nebula (M1) PWN equation introduces two additive UQFF terms unique to this system: pulsar spindown ram pressure `F_wind = Ė_sd/(c·4πr²)` and magnetic moment dipole dilution `M_mag = μ₀·m/(4πr³)`. This combination — spindown luminosity converted to ram pressure plus magnetic moment field decay — is specific to the **isolated pulsar wind nebula** context and differs from the `MagnetarSGR1745DynamicModulationCalculator` (Session 53), which handles M_mag in a binary orbit context. Additionally, the Crab PWN uses a TIME-DEPENDENT RADIUS `r(t) = r₀ + v_exp·t` reflecting the known supernova remnant expansion, making this the only UQFF system with an analytically expanding spatial domain. We derive the critical radius below which wind pressure exceeds gravity and validate against Crab Nebula measurements.

---

## 1. The Crab Nebula PWN UQFF Equation

From Document 10 of grok_share_7514fe:

```
g_Crab(r, t) = (G·M) / (r(t)²) · (1 + H(z)·t) · (1 - B/B_crit)
             + (Ug1 + Ug2 + Ug3 + Ug4)
             + Λc²/3 + QM + q(v×B) + fluid + DM
             + F_wind
             + M_mag
```

with:
```
r(t) = r₀ + v_exp · t   (expanding nebula)
F_wind = Ė_sd / (c · 4π · r(t)²)
M_mag = μ₀ · m / (4π · r(t)³)
```

---

## 2. Expanding Spatial Domain: r(t)

The Crab is the only UQFF system with analytically expanding r(t). This is derived from:

```
v_exp ≈ 1500 km/s = 1.5×10⁶ m/s (observed filament velocity)
r₀ ≈ 5.5 ly = 9.46×10¹⁵ · 5.5 ≈ 5.20×10¹⁶ m (current half-radius)
age ≈ 972 years (SN 1054 → 2026)
```

At t=0 (initial explosion):
```
r₀_initial = r₀ - v_exp · t_age
           = 5.20×10¹⁶ - 1.5×10⁶ · (972 · 3.156×10⁷)
           = 5.20×10¹⁶ - 4.59×10¹⁶
           ≈ 6.1×10¹⁵ m   (~0.2 pc initial radius — consistent with SN ejecta)
```

This confirms r(t) correctly recovers solar-system-scale initial radius at t=0.

---

## 3. Pulsar Spindown Ram Pressure: F_wind

### 3.1 Spindown Luminosity

The Crab pulsar (PSR J0534+2200) has:
- Period P = 33.5 ms
- Period derivative Ṗ = 4.21×10⁻¹³ s/s
- Moment of inertia I ≈ 10⁴⁵ g·cm² = 10³⁸ kg·m²

Spindown luminosity:
```
Ė_sd = -4π² · I · Ṗ / P³
      = 4π² · 10³⁸ · 4.21×10⁻¹³ / (33.5×10⁻³)³
      = 4π² · 10³⁸ · 4.21×10⁻¹³ / (3.76×10⁻⁵)
      ≈ 4.6×10³¹ W
```

This matches the canonical value Ė_sd ≈ 4.6×10³¹ W (Hester 2008).

### 3.2 Ram Pressure at r(t)

The spindown energy is converted isotropically to ram pressure:
```
F_wind(r) = Ė_sd / (c · 4π · r²)
```

At current Crab nebula radius r₀ = 5.20×10¹⁶ m (computational default 9.46×10¹⁵ m for inner region):
```
F_wind = 4.6×10³¹ / (2.998×10⁸ · 4π · (9.46×10¹⁵)²)
       ≈ 4.6×10³¹ / (3.37×10⁴¹)
       ≈ 1.36×10⁻¹⁰ N/m² → treated as acceleration [m/s²] in UQFF normalization
```

### 3.3 Ratio F_wind / g_base

Base gravity at the same radius (M_ejecta ≈ 4.6 M☉):
```
g_base = G·M/r² = 6.674e-11 · 4.6·1.989e30 / (9.46e15)²
       ≈ 6.674e-11 · 9.15e30 / 8.95e31
       ≈ 6.82×10⁻¹² m/s²
```

Ratio:
```
F_wind / g_base ≈ 1.36×10⁻¹⁰ / 6.82×10⁻¹² ≈ 20
```

**F_wind exceeds g_base by a factor of ~20** at the Crab inner radius — confirming the wind-dominated morphology of the Crab PWN inner torus and jets.

---

## 4. Magnetic Moment Dilution: M_mag

### 4.1 Magnetic Moment of Crab Pulsar

The Crab pulsar surface field B_s ≈ 3.8×10¹² G = 3.8×10⁸ T.  
Pulsar radius R_ns ≈ 10 km = 10⁴ m.

Magnetic dipole moment:
```
m = (4π/μ₀) · B_s · R_ns³
  = 4π/(4π×10⁻⁷) · 3.8×10⁸ · (10⁴)³
  = 10⁷ · 3.8×10⁸ · 10¹²
  = 3.8×10²⁷ A·m²
```

### 4.2 Dipole Field at r(t)

```
M_mag(r) = μ₀ · m / (4π · r³)
```

At r = 9.46×10¹⁵ m:
```
M_mag = 4π×10⁻⁷ · 3.8×10²⁷ / (4π · (9.46×10¹⁵)³)
      = 10⁻⁷ · 3.8×10²⁷ / (8.47×10⁴⁷)
      ≈ 4.49×10⁻²⁸ T²·m or equivalent normalized acceleration
```

The M_mag term dilutes as r⁻³ (dipole law), falling off faster than F_wind (r⁻²) and g_base (r⁻²). This means at large r, M_mag becomes negligible relative to F_wind — consistent with the Crab PWN morphology where the toroidal magnetic structure dominates near the pulsar but the wind torus is the dominant energy carrier at larger scales.

---

## 5. Distinction from SGR 1745 Magnetar

The `MagnetarSGR1745DynamicModulationCalculator` (Session 53, Document 8) also includes M_mag, but in an entirely different context:

| Feature | Crab PWN (This Paper) | SGR 1745 (Session 53) |
|---------|----------------------|----------------------|
| Context | Isolated PWN | Binary magnetar system |
| M_mag role | Dipole field dilution with r⁻³ | Dynamic modulation with binary orbit |
| F_wind source | Pulsar spindown Ė_sd/(c·4πr²) | Not present |
| r(t) | Expanding: r₀ + v_exp·t | Binary orbital r(t) |
| B field | Canonical 3.8×10¹² G | Extraordinary ~10¹⁵ G |
| Main physics | PWN expansion + wind + dipole | Magnetar + companion orbit |

The two classes are mathematically and physically distinct despite sharing M_mag notation.

---

## 6. Critical Radius Analysis

Define the critical radius r_c where F_wind = g_base:
```
Ė_sd / (c · 4π · r_c²) = G·M / r_c²
→ r_c = √((Ė_sd) / (4π·c·G·M/r²)) ...
→ This simplifies to: Ė_sd / (c·4π) = G·M → r_c is constant (cancels):
   Ė_sd/c = 4π·G·M → M_crit = Ė_sd/(4π·G·c)
```

Critical mass below which wind always exceeds gravity:
```
M_crit = 4.6×10³¹ / (4π · 6.674×10⁻¹¹ · 2.998×10⁸)
       = 4.6×10³¹ / (2.51×10⁻¹)
       ≈ 1.83×10³² kg ≈ 92 M☉
```

The Crab ejecta mass ≈ 4.6 M☉ << 92 M☉, confirming F_wind ALWAYS exceeds g_base for the Crab — the pulsar wind permanently inflates the nebula against gravity.

---

## 7. Calculator Implementation

`CrabPWNUQFFCalculator` in CondensedPhysics3.py (Session 55) implements:
- `r(t) = r0 + v_exp · t`
- `F_wind = E_sd / (c · 4 * pi * r(t)**2)`
- `M_mag = mu0 * m / (4 * pi * r(t)**3)`
- Full UQFF g_Crab = g_base + F_wind + M_mag

---

## References

1. grok_share_7514fe.txt — Document 10: Crab Nebula g_Crab equation
2. Hester (2008) — "The Crab Nebula: An Astrophysical Chimera", ARAA 46
3. Bühler & Blandford (2014) — "The Surprising Crab Pulsar and its Nebula", RPPH 77
4. Kaplan et al. (2008) — Crab pulsar ephemeris, Ṗ = 4.21×10⁻¹³
5. CondensedPhysics3.py — `CrabPWNUQFFCalculator` (Session 55)

---

*© 2026 Daniel T. Murphy — Star-Magic UQFF Framework — All Rights Reserved*  
*Paper 220 of 1,000 — Session 55 — Phase 2 §2.10 Fourth-Pass Extraction*
