# PAPER_284: M16 Eagle Nebula UQFF — Dual Mass Co-Action Product (Φ_dm)
## SFR Growth × Photoevaporation Erosion Multiplicative Coupling

**Classification:** UQFF 2.0 Gravitational Physics — Nebular Mass Dynamics  
**System:** M16 Eagle Nebula (IC 4703), Eagle Nebula Star-Forming Region  
**Session:** 80 | **Module:** M16_UQFF_MODULE.cpp (22nd C++ UQFF module)  
**Author:** Daniel T. Murphy | **Date:** March 2026

---

## 1. Abstract

This paper introduces the **Dual Mass Co-Action Product** (Φ_dm) — a UQFF gravity modulation factor that couples star-formation-driven mass accumulation and radiation-driven photoevaporation erosion through a **multiplicative** product rather than the previously used additive form. For M16 (Eagle Nebula), with star formation rate SFR = 1 M☉/yr over initial gas mass M₀ = 1200 M☉, and maximum photoevaporation fraction E₀ = 0.3 (30%), the multiplicative form produces a 24.3% reduction in Φ_dm relative to the additive approximation at t = 5 Myr. This is the **first UQFF module** to simultaneously apply an additive-gain and saturation-subtractive product on the same gravity term.

---

## 2. Physical Motivation

In active star-forming nebulae, two competing processes drive mass evolution:

1. **Star Formation Accretion** — molecular gas accretes onto protostars, increasing the effective gravitational mass fraction by SFR_rate × t.
2. **Photoevaporation Erosion** — UV radiation from newly formed massive stars erodes the surrounding gas, progressively reducing the effective mass by a saturating fraction E_rad(t).

The **additive approximation** (used in prior UQFF modules):
$$\Phi_{dm}^{add} = g_{base} \times (1 + M_{sf}) - E_{rad}$$

linearly superposes the two effects, implicitly treating them as independent processes acting on separate mass reservoirs.

The **multiplicative form** (this paper):
$$\Phi_{dm}(t) = (1 + \text{SFR\_rate} \times t) \times (1 - E_{rad}(t))$$

correctly encodes that **the eroded mass is drawn from the same growing reservoir** — the fraction lost to photoevaporation scales with the mass being accreted, not the original quiescent mass. This is physically accurate for pillar-geometry star formation (e.g., M16's "Pillars of Creation").

---

## 3. Mathematical Formulation

### 3.1 Parameters (M16 Eagle Nebula)

| Parameter | Value | Description |
|-----------|-------|-------------|
| M₀ | 2.387 × 10³³ kg (1200 M☉) | Initial nebula gas mass |
| SFR | 1 M☉/yr | Star formation rate |
| SFR_rate | 2.639 × 10⁻¹¹ s⁻¹ | = SFR / (M₀/M☉) / (3.156×10⁷ s/yr) |
| τ_erode | 9.468 × 10¹³ s (3 Myr) | Photoevaporation e-folding time |
| E₀ | 0.3 | Maximum photoevaporation fraction |
| g_base | 1.454 × 10⁻¹² m/s² | G × M / r² at r = 3.31 × 10¹⁷ m |

### 3.2 Dual Co-Action Product

$$M_{sf}(t) = \text{SFR\_rate} \times t$$

$$E_{rad}(t) = E_0 \left(1 - e^{-t/\tau}\right)$$

$$\boxed{\Phi_{dm}(t) = (1 + M_{sf}) \times (1 - E_{rad})}$$

### 3.3 Dynamic Gravity Term

$$g_{dyn}(t) = g_{base} \times \Phi_{dm}(t)$$

### 3.4 Multiplicative–Additive Gap

The gap relative to the additive form:
$$\Delta_{gap} = \Phi_{dm}^{mult} - \Phi_{dm}^{add} = -(M_{sf} \times E_{rad})$$

This cross-term is always **negative** — the multiplicative form predicts lower gravity than the additive approximation whenever both SFR accumulation and erosion are simultaneously active.

---

## 4. Numerical Results at t = 5 Myr

t = 5 Myr = 1.578 × 10¹⁴ s

| Quantity | Value |
|----------|-------|
| M_sf_frac | SFR_rate × t = 2.639×10⁻¹¹ × 1.578×10¹⁴ = **4164.8** |
| E_rad | E₀ × (1 − exp(−5/3)) = 0.3 × 0.8110 = **0.2433** |
| Φ_dm (multiplicative) | (1 + 4164.8) × (1 − 0.2433) = 4165.8 × 0.7567 = **3151.9** |
| Φ_dm (additive) | (1 + 4164.8) − 0.2433 = **4165.6** |
| **gap_mult_add** | −(4164.8 × 0.2433) = **−1013.3** (24.3% less) |
| g_dyn(5 Myr) | 1.454×10⁻¹² × 3151.9 = **4.583 × 10⁻⁹ m/s²** |

The multiplicative gap of −1013.3 confirms that treating erosion as acting on the growing mass reservoir (not the static initial mass) produces a **measurable 24.3% reduction** compared to the additive approximation.

---

## 5. Connection to UQFF 2.0 g_total

In the full M16 UQFF 2.0 equation:

$$g_{total}(r, t) = \left[g_{dyn}(t) + U_{g,sum}(26) + \Lambda + Q + L + F + g_{exp}\right] \times \text{corr}_{SC}$$

The Φ_dm product modulates only the dynamic base gravity term g_dyn. The 26-layer Triadic (U_g,sum), cosmological Λ, quantum, Lorentz, fluid, and Friedmann expansion terms are all independent of Φ_dm — the modulation is cleanly scoped to the time-evolving mass component.

---

## 6. UQFF Historical Distinction

| Module | SFR Term | Erosion Term | Form |
|--------|----------|-------------|------|
| Session 55 CP3 M16EagleNebulaRadiationSFR | g_base×(1+M_sf) | −E_rad | **Additive** |
| **This paper (PAPER_284)** | (1+M_sf) | ×(1−E_rad) | **Multiplicative** |

This is the **first UQFF module** to use the multiplicative dual co-action form, correctly encoding the coupled feedback between star formation accretion and pillar photoevaporation for M16-class nebulae.

---

## 7. Wolfram KB Term

```
M16UQFF:Phi_dm=(1+SFR_rate*t)*(1-E_0*(1-Exp[-t/tau])); SFR_rate=2.639e-11/s; M_sf_frac=SFR_rate*t [PAPER_284]
```

---

## 8. Cross-References

- **PAPER_285:** Erosion Saturation Half-Time (t_half, ΔgMax)
- **PAPER_286:** Nebular Friedmann Redshift (κ_neb, z=0.0015)
- **M16_UQFF_MODULE.cpp:** Full UQFF 2.0 C++ implementation (22nd module)
- **CondensedPhysics3.py:** `M16DualMassCoActionProductCalculator`

---

*Copyright — Daniel T. Murphy, Session 80, March 2026. UQFF 2.0.*
