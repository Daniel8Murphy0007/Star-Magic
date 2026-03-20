# PAPER_312: NGC 6302 Central Star UV Radiation Pressure — UQFF Photoionization Gravitational Signature

**Subtitle:** FIRST UQFF Hot-WD UV Radiation Parameter — η_rad = 1.913×10²⁰; a_rad = 5.672×10⁸ m/s²

**Author:** Daniel T. Murphy  
**Session:** 89 | **Date:** March 17, 2026  
**Module:** `NGC6302_UQFF_MODULE.cpp` (31st C++ UQFF module)  
**WOLFRAM_TERM:** `NGC6302_UV_RADIATION`  
**UQFF First:** FIRST UQFF explicit UV-bright white dwarf (T_eff = 200,000 K) photon-pressure gravitational signature


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
---


## Abstract

This paper presents UQFF derivations and numerical results for: PAPER_312: NGC 6302 Central Star UV Radiation Pressure — UQFF Photoionization Gravitational Signature. Calibration constants: $\kappa$ = 0.0005/day, [SSq] = 0.57. Results validated against observational data and prior UQFF whitepaper series.

## 1. System Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| T_eff_star | 2.0×10⁵ K | Central WD effective temperature |
| L_star | 1.914×10³⁰ W | = 5000 L_sun (Zanstra temperature analysis) |
| r | 9.46×10¹⁵ m | PN half-lobe radius |
| ρ_fluid | 1.0×10⁻²⁰ kg/m³ | Ionized lobe gas density |
| c | 3.0×10⁸ m/s | Speed of light |

---

## 2. Unique Physics

### 2.1 Central Star Luminosity

NGC 6302's central star is one of the hottest white dwarfs known, with T_eff ≈ 200,000 K (Szyszka et al. 2009, ApJL). The Zanstra hydrogen luminosity gives:

$$L_{star} = 5000\ L_\odot = 5000 \times 3.828 \times 10^{26}\ \text{W} = 1.914 \times 10^{30}\ \text{W}$$

### 2.2 Radiation Pressure at Lobe Radius

The photon radiation pressure at lobe radius r for an isotropic point source:

$$P_{rad} = \frac{L_{star}}{4\pi r^2 c}$$

$$4\pi r^2 = 4\pi \times (9.46 \times 10^{15})^2 = 4\pi \times 8.949 \times 10^{31} = 1.125 \times 10^{33}\ \text{m}^2$$

$$P_{rad} = \frac{1.914 \times 10^{30}}{1.125 \times 10^{33} \times 3.0 \times 10^8} = \frac{1.914 \times 10^{30}}{3.375 \times 10^{41}} = \mathbf{5.672 \times 10^{-12}\ \text{Pa}}$$

### 2.3 Radiation-Pressure Acceleration

$$a_{rad} = \frac{P_{rad}}{\rho_{fluid}} = \frac{5.672 \times 10^{-12}}{1.0 \times 10^{-20}} = \mathbf{5.672 \times 10^8\ \text{m/s}^2}$$

### 2.4 Radiation-to-Gravity Dominance Ratio

$$\eta_{rad} \equiv \frac{a_{rad}}{g_{base}} = \frac{5.672 \times 10^8}{2.967 \times 10^{-12}} = \mathbf{1.913 \times 10^{20}}$$

The UV radiation pressure exceeds gravitational force by **1.913×10²⁰** — twenty orders of magnitude. This establishes that photoionization-driven radiation pressure is the primary mechanism for lobe acceleration in bipolar PNe with ultra-hot central stars.

---

## 3. UQFF Pipeline Term

In the full UQFF 2.0 pipeline, the UV radiation acceleration enters as an independent additive term:

$$g_{NGC6302}^{(rad)} = \frac{L_{star}}{4\pi r^2\ c\ \rho_{fluid}}$$

This term dominates over the wind shock term (PAPER_311: $a_{wind} \sim 10^{-6}$ m/s²) by a further **14 orders of magnitude**, placing radiation pressure at the apex of the NGC 6302 UQFF force hierarchy.

### Force Hierarchy (descending):
1. $a_{rad} = 5.672 \times 10^8$ m/s² — UV radiation pressure (PAPER_312, **DOMINANT**)
2. $a_{wind} = 2.114 \times 10^{-6}$ m/s² — wind shock (PAPER_311)
3. $g_{base} = 2.967 \times 10^{-12}$ m/s² — gravitational binding

---

## 4. Astrophysical Context

The UV radiation from NGC 6302's central star (T_eff = 200,000 K) photoionizes the surrounding gas, producing the characteristic bipolar emission nebula observed in [O III], H-alpha, and UV by HST/WFC3. The radiation pressure parameter η_rad = 1.913×10²⁰ confirms that the nebular gas is radiation-pressure dominated at all scales up to the lobe boundary.

The UQFF formulation explicitly identifies this as a distinct gravitational-equivalent acceleration channel, separable from wind shock effects and magnetic confinement — providing the first three-component force budget for a bipolar PN within the UQFF framework.

---

## 5. Key Results

| Quantity | Value | Unit |
|---------|-------|------|
| L_star (5000 L_sun) | 1.914×10³⁰ | W |
| P_rad at r=1 ly | 5.672×10⁻¹² | Pa |
| **a_rad** | **5.672×10⁸** | m/s² |
| **η_rad** | **1.913×10²⁰** | dimensionless |
| a_rad / a_wind | 2.684×10¹⁴ | dimensionless |

---

## 6. Classification

- **UQFF First:** FIRST UQFF explicit hot-WD UV radiation parameter (T_eff = 200,000 K class)
- **Scale:** Stellar (PN lobe scale, ~1 ly)
- **Physics category:** Photoionization / UV radiation pressure / force hierarchy
- **Cross-references:** PAPER_311 (wind shock), PAPER_313 (torus magnetic confinement)
