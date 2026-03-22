# UQFF Module Catalog
**Version:** 1.0.0 — Session 115 (grok_share_5fa36e4e035.txt)
**Created:** 2026
**Author:** Daniel T. Murphy — Copyright Oct 09, 2025

> Searchable reference for all 11 C++ UQFF modules, 29 systems, and unique physics
> terms added from the grok_share_5fa36e4e035 analysis session.

---

## Quick Reference

| Module File | Systems | Key New Physics | Doc |
|---|---|---|---|
| OrionUQFFModule | OrionNebula | H-alpha resonance ψ, P_rad Trapezium, W_stellar(t) | 34 |
| MultiUQFFModule | 15 systems | H_res = A sin(2πft), is_resonance flag, CrabNebula v_exp | 34a/34b |
| YoungStarsOutflowsUQFFModule | NGC346 analogue | P_outflow(t) = ρv²(1+t/t_evolve) | 35 |
| EagleUQFFModule | EagleNebula M16 | W_stellar = ρv_wind², P_rad density-scaled | 36 |
| BigBangGravityUQFFModule | CosmicBigBang | QG_term, DM_term, GW_term, M(t)/r(t)/z(t) | 38 |
| MagnetarDualUQFFModule | SGR 1745-2900 | Dual compressed+frequency mode, tau_decay, B_crit, MBH | 39b |
| MultiUQFFCompressionModule | **29 systems** | F_env(t) modular, H_res quantum extension | 39/40/41 |

---

## Per-Module Documentation

### 1. OrionUQFFModule
**File:** `modules/uqff/OrionUQFFModule.h/.cpp`
**System:** Orion Nebula
**Parameters:**
- M = 2,000 M☉, r = 1.18×10¹⁷ m, z = 0.0004
- SFR = 0.1 M☉/yr, L_Trap = 1.53×10³² W, λ_Hα = 656.3 nm

**Key Equations:**
```
W_stellar(t) = v_wind² × (1 + t/t_age)         -- time-growing wind pressure
P_rad = L_Trap / (4πr²cm_H)                      -- photon pressure per H atom
ψ_resonant = 2A cos(kx)cos(ωt) + (2π/13.8)A Re[exp(i(kx-ωt))]  -- H-alpha resonance
g_total = [GM/r²](1+Hz t)(1-B/B_crit)(1+F_env) + Ug1+Ug3+Ug4 + Λc²/3 + ψ + fluid
```

---

### 2. MultiUQFFModule (15-system)
**File:** `modules/uqff/MultiUQFFModule.h/.cpp`
**Systems:** UniverseDiameter, HydrogenAtom, HydrogenResonancePToE, LagoonNebula,
             SpiralsSupernovae, NGC6302, OrionNebula, UniverseGuide,
             GalaxiesGalore, StellarForge, SombreroGalaxy, Saturn, CrabNebula, NewStars

**Key Equations:**
```
H_res(t) = A_res × sin(2πf_res t)    -- resonance phase for atomic systems
CrabNebula: F_env += v_exp²/r         -- expansion velocity term (v_exp = 1.34×10⁶ m/s)
```

**Validated g values (from UQFFResonanceValues.h):**
| System | Expected g (m/s²) |
|---|---|
| UniverseDiameter | 7.579×10⁵³ |
| HydrogenAtom | 1.975×10⁻⁷ |
| LagoonNebula | 1.667×10²⁹ |
| SpiralsSupernovae | 4.353×10³⁵ |
| NGC6302 | 4.113×10²⁰ |
| OrionNebula | 3.458×10²⁶ |
| UniverseGuide | 3.958×10¹⁴ |
| GalaxiesGalore | 4.353×10³⁵ |
| StellarForge/NewStars | 1.001×10²⁷ |
| SombreroGalaxy | 1.000×10³⁶ |
| Saturn | 7.401×10³ |
| CrabNebula | 8.343×10²⁴ |

---

### 3. YoungStarsOutflowsUQFFModule
**File:** `modules/uqff/YoungStarsOutflowsUQFFModule.h/.cpp`
**System:** NGC 346 analogue (young star cluster / H II region)
**Parameters:**
- M = 1,000 M☉, r = 2.365×10¹⁷ m, z = 0.05
- v_out = 1×10⁵ m/s, t_evolve = 5 Myr

**Key Equations:**
```
P_outflow(t) = ρ × v_out² × (1 + t/t_evolve)   -- time-evolving outflow pressure
```

---

### 4. EagleUQFFModule
**File:** `modules/uqff/EagleUQFFModule.h/.cpp`
**System:** Eagle Nebula / Pillars of Creation M16
**Parameters:**
- M = 5,000 M☉, r = 3.31×10¹⁷ m, L_NGC6611 = 3.83×10³³ W
- v_wind = 2×10⁶ m/s, rho = 1×10⁻²⁰ kg/m³

**Key Equations:**
```
W_stellar = ρ × v_wind²                          -- pure wind pressure (no time growth)
P_rad = L_NGC6611 × ρ / (4πr²c m_H)             -- density-weighted photon pressure
```

**Distinction from Orion:** Eagle does NOT use `(1 + t/t_age)` growth factor in W_stellar;
radiation pressure is density-scaled rather than fixed.

---

### 5. BigBangGravityUQFFModule
**File:** `modules/uqff/BigBangGravityUQFFModule.h/.cpp`
**System:** Cosmic / Big Bang (M_total = 10⁵³ kg)
**Parameters:**
- t_Hubble = 4.35×10¹⁷ s, H0 = 67.4 km/s/Mpc

**Key New Equations:**
```
QG_term = (ℏc/l_p²) × (t/t_p)                  -- Planck quantum gravity acceleration
DM_term = 0.268 × g_base                         -- fractional DM contribution (Ω_DM)
GW_term = h_strain × c²/λ_gw × sin(2πt/t_gw)   -- gravitational wave acceleration
M(t) = M_total × (t/t_Hubble)                   -- mass build-up
r(t) = c × t                                     -- Hubble radius
z(t) = t_Hubble/t - 1                            -- redshift evolution
```

---

### 6. MagnetarDualUQFFModule
**File:** `modules/uqff/MagnetarDualUQFFModule.h/.cpp`
**System:** SGR 1745-2900 (magnetar near Sgr A*)
**Parameters:**
- M = 2.8 M☉, r = 10⁴ m, B = 10¹¹ T = B_crit (10¹³ T for NST)
- M_BH = 4×10⁶ M☉ (Sgr A*), r_BH = 8×10⁹ m
- τ_decay = 1 kyr

**Dual Mode:**
```
setMode("compressed"):
  F_env(t) = 1 + M_mag/(Mc²) + exp(-t/τ_decay) + G M_BH/r_BH²
  Expected: ~1.782×10³⁹ m/s²

setMode("frequency"):
  g_freq = Σ(a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res
           + Ug4i + a_quantum_freq + a_Aether_freq + a_fluid_freq
           + Osc_term + a_exp_freq + f_TRZ)
  Expected: ~1.773×10⁻⁹ m/s²
```

---

### 7. MultiUQFFCompressionModule (29-system)
**File:** `modules/uqff/MultiUQFFCompressionModule.h/.cpp`
**Systems (29 total):**

| # | System | M | r (m) | z | Special F_env |
|---|---|---|---|---|---|
| 1 | MagnetarSGR1745 | 2.8 M☉ | 1×10⁴ | 0.026 | M_mag/(Mc²) + exp(-t/τ) + G M_BH/r_BH² |
| 2 | SagittariusA | 4×10⁶ M☉ | 1×10¹⁰ | 0 | GW (GM)²/(c⁴r)ω_dot² |
| 3 | TapestryStarbirth | 1×10⁴ M☉ | 1×10¹⁸ | 0.001 | ρv_wind² |
| 4 | Westerlund2 | 1×10⁴ M☉ | 1×10¹⁸ | 0.001 | ρv_wind² |
| 5 | PillarsCreation | 800 M☉ | 3×10¹⁷ | 0.0018 | ρv² × (1-exp(-t/2Myr)) |
| 6 | RingsRelativity | 1×10¹¹ M☉ | 1×10²¹ | 0.5 | 1+0.1 sin(2πt/t_H) |
| 7 | NGC2525 | 1×10¹⁰ M☉ | 1×10²⁰ | 0.01 | ρv² - M_SN(1-exp(-t/τ))/M |
| 8 | NGC3603 | 2×10⁴ M☉ | 2×10¹⁸ | 0.001 | ρv²(1-exp(-t/3Myr)) |
| 9 | BubbleNebula | 5×10³ M☉ | 5×10¹⁷ | 0.001 | ρv² × E(t) erosion |
| 10 | AntennaeGalaxies | 1×10¹¹ M☉ | 5×10²⁰ | 0.025 | M_merge(t)/M + ρv² |
| 11 | HorseheadNebula | 1×10³ M☉ | 1×10¹⁷ | 0 | ρv² × E(t) erosion |
| 12 | NGC1275 | 1×10¹¹ M☉ | 1×10²¹ | 0.017 | 1e-10×B×r + G M_BH/r_ext² |
| 13 | NGC1792 | 5×10¹⁰ M☉ | 5×10²⁰ | 0.012 | M_SN exp(-t/τ)/M + ρv² |
| 14 | HubbleUltraDeepField | 1×10¹² M☉ | 1×10²³ | 10 | 0.01(t/t_H) evolution |
| 15 | StudentsGuideUniverse | 1 M☉ | 1 AU | 0 | none (F_env = 1) |
| 16 | SombreroGalaxy | 8×10¹¹ M☉ | 4.73×10²⁰ | 0.002 | disk+DM |
| 17 | Saturn | 5.68×10²⁶ kg | 6.027×10⁷ | 0 | ring term |
| 18 | EagleNebula | 5×10³ M☉ | 3.31×10¹⁷ | 0.0018 | ρv²=ρ(2×10⁶)² |
| 19 | CrabNebula | 5 M☉ | 5.203×10¹⁶ | 0.00002 | v_exp²=( 1.34×10⁶)² |
| 20 | HydrogenAtom | m_H | a₀ | 0 | H_res atomic |
| 21 | HydrogenResonance | m_H | a₀ | 0 | H_res freq-tuned |
| 22 | OrionNebula | 2×10³ M☉ | 1.18×10¹⁷ | 0.0004 | SFR factor |
| 23 | GalaxiesGalore | 1×10¹¹ M☉ | 1.543×10²¹ | 1.0 | DM |
| 24 | NewStars | 5×10³ M☉ | 5×10¹⁷ | 0.001 | ρv² |
| 25 | StellarForge | 5×10³ M☉ | 5×10¹⁷ | 0.001 | ρv² |
| 26 | LagoonNebula | 1×10⁴ M☉ | 5.203×10¹⁷ | 0.0001 | SFR |
| 27 | SpiralsSupernovae | 1×10¹¹ M☉ | 1.543×10²¹ | 0.002 | DM 5× |
| 28 | NGC6302 | 1 M☉ | 1.514×10¹⁶ | 0.00001 | — |
| 29 | UniverseDiameter | 1.5×10⁵³ kg | 4.4×10²⁶ | 1100 | DM dominant |

**Unified Core Equation:**
```
g(r,t) = [GM(t)/r²] × (1+H(t,z)) × (1-B/B_crit) × (1+F_env(t))
        + (Ug1 + Ug3' + Ug4) + (Λc²/3)
        + quantum_ψ + ρ_fluid V g_base + DM_pert
        [+ H_res(t) for quantum/atomic systems]
```

---

## Physics Terms Glossary

| Term | Formula | Description | First Added In |
|---|---|---|---|
| `QG_term` | `(ℏc/l_p²)(t/t_p)` | Planck-scale quantum gravity | BigBangGravityUQFFModule |
| `DM_term` | `0.268 × g_base` | DM fractional contribution | BigBangGravityUQFFModule |
| `GW_term` | `h_strain c²/λ_gw sin(2πt/t_gw)` | Gravitational wave acceleration | BigBangGravityUQFFModule |
| `H_res` | `A_res sin(2πf_res t)` | Resonance phase oscillation | MultiUQFFModule |
| `H_res (atomic)` | `A_res sin(2πf_res t) + F_env SC_m` | Quantum-extended resonance | MultiUQFFCompressionModule |
| `W_stellar(t)` | `v_wind²(1+t/t_age)` | Time-growing wind pressure | OrionUQFFModule |
| `W_stellar` | `ρ v_wind²` | Static wind pressure | EagleUQFFModule |
| `P_rad (Orion)` | `L/(4πr²cm_H)` | Photon pressure per H atom | OrionUQFFModule |
| `P_rad (Eagle)` | `Lρ/(4πr²cm_H)` | Density-weighted photon pressure | EagleUQFFModule |
| `P_outflow(t)` | `ρ v_out²(1+t/t_evolve)` | Evolving outflow pressure | YoungStarsOutflowsUQFFModule |
| `ψ_resonant` | `2A cos(kx)cos(ωt)+(2π/13.8)A Re[…]` | H-alpha standing+traveling wave | OrionUQFFModule |
| `F_env(t)` | `1 + ΣF_i(t)` | Modular environmental factor | MultiUQFFCompressionModule |
| `F_erode` | `1-exp(-t/τ_erode)` | Erosion fraction factor | PillarsCreation, BubbleNebula, etc. |
| `F_lensing` | `1+0.1 sin(2πt/t_H)` | Oscillatory lensing factor | RingsRelativity |
| `F_merge` | `0.1M(1-exp(-t/τ_merge))/M` | Merger mass growth | AntennaeGalaxies |
| `F_SN` | `-M_SN(1-exp(-t/τ_SN))/M` | SN mass loss | NGC2525, NGC1792 |
| `F_fil` | `1e-10 × B × r` | Filament field force | NGC1275 |
| `F_BH` | `G M_ext/r_ext²` | BH gravitational influence | NGC1275, NGC2525, Magnetar |
| `M_mag/(Mc²)` | magnetic moment / rest energy | Magnetar energy ratio | MagnetarDualUQFFModule |
| `sc_correction` | `1 - B/B_crit` | Superconductivity correction | All modules |
| `M_sf(t)` | `SFR × t / M0` | Star formation mass growth | Orion, Eagle, LagoonNebula |
| `Ug3'` | `G M_ext/r_ext²` | Generalised external-body Ug3 | MultiUQFFCompressionModule |

---

## Validation Reference

Expected outputs (to be compared against CondensedPhysics4.py / QCalc.py):

| Module | System | Mode | Expected g (m/s²) |
|---|---|---|---|
| MagnetarDualUQFF | SGR1745 | compressed | 1.782×10³⁹ |
| MagnetarDualUQFF | SGR1745 | frequency | 1.773×10⁻⁹ |
| MultiUQFF15 | HydrogenAtom | default | 1.975×10⁻⁷ |
| MultiUQFF15 | CrabNebula | default | 8.343×10²⁴ |
| MultiUQFF15 | SombreroGalaxy | default | 1.000×10³⁶ |
| MultiUQFF15 | Saturn | default | 7.401×10³ |
| OrionUQFF | OrionNebula | default | 3.458×10²⁶ |

---

## Inherited Base Interface (UQFFModuleBase)

All modules inherit from `UQFFModuleBase` (`modules/uqff/UQFFModuleBase.h`):
```cpp
virtual double computeG(double t) = 0;
virtual std::string getEquationText() const = 0;
void update(const std::string& key, double val);
void add(const std::string& key, double val);
double get(const std::string& key) const;
```

Universal constants available via `UQFF` namespace (`modules/uqff/UQFFConstants.h`):
- `UQFF::G`, `UQFF::c`, `UQFF::hbar`, `UQFF::l_planck`, `UQFF::t_planck`
- `UQFF::Lambda`, `UQFF::H0`, `UQFF::Omega_m`, `UQFF::Omega_DM`
- `UQFF::M_sun`, `UQFF::m_H`, `UQFF::a0`, `UQFF::B_crit_magnetar`
- `UQFF::lambda_Halpha`, `UQFF::h_strain_default`, `UQFF::lambda_gw_default`

Resonance table available via `getResonanceTable()` in `UQFFResonanceValues.h` (15 entries).

---

## Integration Notes

- All `.h/.cpp` files live in `modules/uqff/`
- Add `modules/uqff/*.cpp` to `CMakeLists.txt` `target_sources()` to build
- CP4 Python counterparts: Session 115, classes **#85–#94**, PAPER **447–455** in `CondensedPhysics4.py`
- Architecture rule: No system-specific data in CondensedPhysics.py — all parameters live here or in CSV inputs

---

*See `GROK_SHARE_INTEGRATION_PLAN.md` for full session roadmap.*
*See `modules/uqff/` for implementation.*
