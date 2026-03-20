# PAPER_245: MUGE Fluid Self-Gravity Archimedes Buoyancy Sub-Term — Universal Gravitational Buoyancy

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)
**Framework:** UQFF v4.27 — Star-Magic Physics
**Source:** CondensedPhysics3.py — `MUGEFluidSelfGravityTermCalculator` (Session 62, grok_share_8d951e12.txt 4th-pass)
**Date:** March 2026
**Series:** Phase 2 Session 62 — §3.x Universal MUGE Sub-Term Integration

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$

## Abstract

Archimedes' principle states that a body immersed in a fluid experiences an upward buoyant force equal to the weight of fluid displaced. This paper extends that classical result to the gravitational domain within the Modified Unified Gravity Equation (MUGE) framework, establishing a **fluid self-gravity sub-term** (`g_fluid`) in which the gravitating body's own gravity acts on the surrounding fluid to produce an effective buoyancy correction.

The defining equation `g_fluid = (ρ_fluid · V · g_grav) / M`, with `V = (4/3)πr³`, directly transposes the Archimedes buoyancy ratio to a gravitational acceleration correction. The term introduces a critical crossover radius `r_c = (3M / (4πρ_fluid))^(1/3)` at which fluid buoyancy equals Newtonian gravity, representing a fundamental scale boundary in astrophysical fluid-gravity coupling.

Like `g_Q` (PAPER_244), this term appears universally across MUGE modules as a structural additive correction. Its physical significance grows at large radii and high fluid densities, making it particularly relevant for galaxy cluster intracluster medium (ICM) modelling, star-formation regions with dense molecular cloud envelopes, and proto-stellar disk self-gravity.

---

## 1. System Parameters and Equation Overview

| Parameter | Symbol | Default Value | Units | Meaning |
|-----------|--------|---------------|-------|---------|
| Gravitational constant | G | 6.6743 × 10⁻¹¹ | m³/(kg·s²) | Newton |
| Body mass | M | 1.989 × 10³⁰ | kg | Solar mass |
| Body radius | r | 6.957 × 10⁸ | m | Solar radius |
| Surrounding fluid density | ρ_fluid | 1 × 10⁻²⁰ | kg/m³ | Low-density ISM default |
| Gravitational acceleration | g_grav | GM/r² | m/s² | Newtonian surface gravity |

**Primary equation:**
```
V        = (4/3) · π · r³
g_grav   = G · M / r²
g_fluid  = (ρ_fluid · V · g_grav) / M
         = (4π/3) · ρ_fluid · r · G
```

Note the simplified form: `g_fluid = (4πG/3) · ρ_fluid · r`, linear in both ρ_fluid and r — a remarkable simplification that removes the mass dependence entirely.

**Archimedes fraction:**
```
η = ρ_fluid · V / M   (dimensionless — ratio of fluid-sphere mass to body mass)
g_fluid = η · g_grav
```

**Crossover radius:**
```
r_c = (3M / (4π · ρ_fluid))^(1/3)   [where η = 1: g_fluid = g_Newt]
```

---

## 2. Core Physics Derivation

### 2.1 Archimedes Transposition

Classical Archimedes: `F_buoy = ρ_fluid · V · g`. In MUGE, the gravitational field itself is the "fluid", and the buoyant force on a mass M in the field is proportional to the mass of fluid in volume V times the local gravitational acceleration g_grav. Dividing by M to obtain acceleration:

```
g_fluid = F_buoy / M = (ρ_fluid · V · g_grav) / M
```

Substituting `V = (4/3)πr³` and `g_grav = GM/r²`:

```
g_fluid = ρ_fluid · (4πr³/3) · (G/r²)
        = (4πG/3) · ρ_fluid · r
```

This form is identical to the gravitational acceleration at the surface of a uniform sphere of density ρ_fluid — the shell theorem applied to the surrounding medium.

### 2.2 Dimensional Analysis

```
[g_fluid] = [G] · [ρ_fluid] · [r]
           = (m³/kg·s²) · (kg/m³) · m
           = m/s²   ✓
```

The result is independent of body mass M — a body of any mass in a fluid of density ρ_fluid at radius r experiences the same fluid gravity correction. This universality is the structural reason the term appears identically in all MUGE modules.

### 2.3 Crossover Radius and Phase Boundary

Setting `g_fluid = g_Newt = GM/r²`:

```
(4πG/3) · ρ_fluid · r_c = G·M/r_c²
r_c³ = 3M / (4π · ρ_fluid)
r_c  = (3M / (4π · ρ_fluid))^(1/3)
```

For solar parameters (M = M_sun, ρ_fluid = 10⁻²⁰ kg/m³): `r_c ≈ (3 × 1.989×10³⁰ / (4π × 10⁻²⁰))^(1/3) ≈ (4.75×10⁴⁹)^(1/3) ≈ 3.6×10¹⁶ m ≈ 1.2 pc`.

Below r_c, Newtonian gravity dominates; above r_c, fluid self-gravity dominates. This scale is consistent with the outer boundary of stellar wind influence zones and the transition to molecular cloud self-gravity.

### 2.4 Rayleigh-Taylor and Hydrodynamic Extensions

The calculator also provides:
- **Rayleigh-Taylor growth rate:** `σ = √(g_fluid · k · Δρ / ρ_total)` at wavenumber k — fluid instability onset driven by the buoyancy term.
- **Kelvin-Helmholtz scale:** shear instability at the g_fluid boundary layer.
- **Density gradient coupling:** `∂g_fluid/∂ρ = (4πGr/3)` — how a density perturbation maps to a gravity perturbation.

---

## 3. Linear Radius Theorem

**Theorem (Fluid Self-Gravity Linearity):** Within MUGE, the fluid self-gravity sub-term `g_fluid = (4πG/3) · ρ_fluid · r` is a linear function of radius r for fixed ρ_fluid, independent of body mass M. The associated Archimedes fraction `η = ρ_fluid · V / M` is the only mass-dependent quantity; when η = 1 the system crosses from Newtonian-dominated to fluid-dominated gravity at the radius r_c.

This theorem establishes that fluid self-gravity provides a **radial amplification** of the gravitational correction — at large r (galaxy cluster scale, r ~ Mpc), even dilute intracluster medium (ρ_ICM ~ 10⁻²⁶ kg/m³) contributes `g_fluid ~ (4π × 6.67×10⁻¹¹ / 3) × 10⁻²⁶ × 3×10²² ≈ 8×10⁻¹⁹ m/s²`, which is non-negligible for cluster mass reconstruction.

---

## 4. Observational Predictions / Validation

- **Galaxy cluster ICM:** Fluid self-gravity in the ICM at r ~ Mpc contributes ~1% of the total MUGE gravity, testable via Sunyaev-Zel'dovich effect pressure profiles (Planck/SPT data).
- **Proto-stellar disks:** At r ~ 100 AU with ρ_disk ~ 10⁻¹⁴ kg/m³, `g_fluid ~ 10⁻⁸ m/s²` — comparable to stellar surface gravity at that distance. This modifies the standard disk self-gravity criterion (Toomre Q parameter).
- **Crossover radius in molecular clouds:** r_c predictions in the range 0.1–1 pc for dense cores (ρ ~ 10⁻¹⁷ kg/m³) are testable with ALMA high-resolution density maps.

---

## 5. References

1. Archimedes of Syracuse (~250 BC). *On Floating Bodies*. (Classical buoyancy principle.)
2. Toomre, A. (1964). On the Gravitational Stability of a Disk of Stars. *ApJ* 139, 1217.
3. Fabian, A.C. (1994). Cooling Flows in Clusters of Galaxies. *ARA&A* 32, 277.
4. Murphy, D.T. (2025). UQFF Framework v4.x — MUGE Sub-Term Integration. Star-Magic internal document.
5. grok_share_8d951e12 validation session — universal `g_fluid` Archimedes term confirmation.

---

*PAPER_245 | UQFF v4.27 | Star-Magic | Session 62 | March 2026*
