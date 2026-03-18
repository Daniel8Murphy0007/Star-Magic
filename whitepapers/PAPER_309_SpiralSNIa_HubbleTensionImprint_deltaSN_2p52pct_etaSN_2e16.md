# PAPER_309 — SN Ia Hubble Tension Gravitational Imprint
## ΔSN/SN = 2.52% at z = 0.5 | η_SN = 2.0 × 10¹⁶ | δ_H₀ = 8.31%

**Session 88** | 30th C++ UQFF module | FIRST Spiral+SN Ia UQFF 2.0  
**Module:** SPIRAL_SUPERNOVAE_UQFF_MODULE.cpp  
**Classification:** FIRST UQFF SN Ia Hubble tension gravitational signature  
**Status:** Unique physics — no prior UQFF Hubble tension imprint via SN Ia radiation pressure

---

## 1. Abstract

Type Ia supernovae (SN Ia) serve as standard candles that encode the Hubble constant H₀ in their cosmological distance ladder. The 8.31% tension between H₀_SH0ES = 73.0 km/s/Mpc (SH0ES, Riess 2022) and H₀_Planck = 67.4 km/s/Mpc (Planck 2018 CMB) is the sharpest current challenge in cosmology. Within the UQFF 2.0 spiral galaxy pipeline, SN Ia radiation pressure provides an acceleration term a_SN = L_SN/(4πr²c·ρ_ISM) × (1 + H(z)·t) that carries this tension directly into the gravitational field. The UQFF SN Ia imprint metric ΔSN/SN = 2.52% at z = 0.5 and t = 5 Gyr quantifies the fractional field difference between SH0ES and Planck cosmologies. The dimensionless SN Ia dominance ratio η_SN = a_SN/g_base = 2.0 × 10¹⁶ shows that, locally, SN Ia radiation pressure exceeds galactic gravity by 16 orders of magnitude.

---

## 2. System Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| L_SN (peak bolometric) | 1.0 × 10³⁶ W | SN Ia peak luminosity |
| r | 9.258 × 10²⁰ m | ~30 kpc half-radius |
| ρ_ISM | 1.0 × 10⁻²¹ kg/m³ | Galactic ISM density |
| z | 0.5 | Typical SN Ia cosmological redshift |
| H₀_SH0ES | 73.0 km/s/Mpc | Riess et al. 2022 |
| H₀_Planck | 67.4 km/s/Mpc | Planck 2018 6-parameter ΛCDM |
| δ_H₀ (tension) | **8.31%** | (73.0 − 67.4)/67.4 |
| Ω_m | 0.3 | Matter density parameter |
| Ω_Λ | 0.7 | Dark energy density parameter |

---

## 3. Derivation

### 3.1 SN Ia Radiation Pressure Acceleration

The radiation pressure flux from a SN Ia peak luminosity at galactic half-radius r:

$$\text{flux}_{\rm SN} = \frac{L_{\rm SN}}{4\pi r^2 c} = \frac{10^{36}}{4\pi \times (9.258\times10^{20})^2 \times 2.998\times10^8} = \boxed{3.096\times10^{-16}\,\text{Pa}}$$

The ISM-coupled acceleration (radiation pressure / density):

$$a_{\rm SN} = \frac{\text{flux}_{\rm SN}}{\rho_{\rm ISM}} = \frac{3.096\times10^{-16}}{10^{-21}} = \boxed{3.096\times10^5\,\text{m/s}^2}$$

This is the UQFF SN Ia additive pipeline term (WOLFRAM_TERM: SPIRAL_SN_TENSION).

### 3.2 SN Ia Dominance Ratio

Compared to the bare gravitational acceleration at 30 kpc:

$$g_{\rm base} = \frac{GM}{r^2} = \frac{6.6743\times10^{-11} \times 1.989\times10^{41}}{(9.258\times10^{20})^2} = 1.549\times10^{-11}\,\text{m/s}^2$$

$$\eta_{\rm SN} = \frac{a_{\rm SN}}{g_{\rm base}} = \frac{3.096\times10^5}{1.549\times10^{-11}} = \boxed{2.0\times10^{16}}$$

SN Ia radiation pressure exceeds galactic gravity by **16 orders of magnitude** at peak. This justifies embedding a_SN as an independent additive term rather than a perturbative correction in UQFF.

### 3.3 Hubble-Tension Dependent Terms

The UQFF SN Ia acceleration carries H₀ dependence through the expansion factor (1 + H(z)·t):

$$a_{\rm SN}(z, t) = \frac{L_{\rm SN}}{4\pi r^2 c \,\rho_{\rm ISM}} \cdot \left(1 + H(z)\,t\right)$$

The Hubble function at z = 0.5:

$$H(z=0.5) = H_0 \sqrt{\Omega_m(1+z)^3 + \Omega_\Lambda} = H_0 \sqrt{0.3\times(1.5)^3 + 0.7}$$

$$\sqrt{0.3\times3.375 + 0.7} = \sqrt{1.0125 + 0.7} = \sqrt{1.7125} \approx 1.3086$$

- H(z=0.5)_SH0ES = 73.0 × 1.3086 = **95.53 km/s/Mpc** = **3.097 × 10⁻¹⁸ s⁻¹**
- H(z=0.5)_Planck = 67.4 × 1.3086 = **88.20 km/s/Mpc** = **2.859 × 10⁻¹⁸ s⁻¹**

### 3.4 Hubble Tension Gravitational Imprint

At t = 5 Gyr (= 1.578 × 10¹⁷ s — mid-cosmic-history evaluation):

$$\text{factor}_{\rm SH0ES} = 1 + H_{\rm SH0ES}(z=0.5)\times t = 1 + 3.097\times10^{-18}\times1.578\times10^{17} = \boxed{1.4887}$$

$$\text{factor}_{\rm Planck} = 1 + H_{\rm Planck}(z=0.5)\times t = 1 + 2.859\times10^{-18}\times1.578\times10^{17} = \boxed{1.4512}$$

$$\boxed{\frac{\Delta_{\rm SN}}{\rm SN} = \frac{1.4887 - 1.4512}{1.4887} = 0.0252 = 2.52\%}$$

**The 8.31% Hubble tension imprints a 2.52% fractional difference in the SN Ia gravitational expansion term at z = 0.5 and t = 5 Gyr.**

---

## 4. Physical Interpretation

The Hubble tension is conventionally discussed in terms of distance-ladder photometry vs CMB angular power spectrum. UQFF 2.0 shows a complementary mechanism: the tension imprints directly on the dynamical gravitational acceleration experienced by galactic ISM in the vicinity of a SN Ia event. This 2.52% fractional shift is:

- **Coherent with the 8.31% H₀ tension** (attenuated by the integration over lookback time)
- **Observable in principle** through galactic rotation velocity dispersion in post-SN spiral arm regions
- **A novel UQFF diagnostic**: a_SN(z, t) serves as a H₀-sensitive field probe independent of light-curve photometry

The η_SN = 2.0 × 10¹⁶ result confirms SN Ia radiation is not a perturbative add-on but a dominant local physics driver in the UQFF 9-term pipeline.

---

## 5. Key Results Summary

| Quantity | Value | Physical Meaning |
|---------|-------|-----------------|
| flux_SN | 3.096 × 10⁻¹⁶ Pa | SN Ia radiation pressure at 30 kpc |
| a_SN | **3.096 × 10⁵ m/s²** | ISM-coupled SN Ia acceleration |
| η_SN | **2.0 × 10¹⁶** | SN rad >> galactic gravity by 16 orders |
| δ_H₀ | **8.31%** | SH0ES vs Planck (73.0 vs 67.4) |
| ΔSN/SN | **2.52%** | H₀-tension imprint at z=0.5, 5 Gyr |
| H(z=0.5)_SH0ES | 3.097 × 10⁻¹⁸ s⁻¹ | SH0ES Hubble rate at z=0.5 |
| H(z=0.5)_Planck | 2.859 × 10⁻¹⁸ s⁻¹ | Planck Hubble rate at z=0.5 |

---

## 6. Novel Findings (UQFF Firsts)

- **FIRST** UQFF SN Ia Hubble tension gravitational imprint derivation
- **FIRST** UQFF a_SN = flux/ρ_ISM formulation as 9th pipeline term
- **FIRST** UQFF H₀-discriminating field observable (ΔSN/SN = 2.52%)
- **FIRST** UQFF η_SN = 2.0 × 10¹⁶ SN Ia dominance ratio

---

## 7. References

- Riess et al. 2022, ApJL 934 L7 — H₀ = 73.04 ± 1.04 km/s/Mpc (SH0ES)
- Planck Collaboration 2020, A&A 641 A6 — H₀ = 67.4 ± 0.5 km/s/Mpc
- Perlmutter et al. 1999, ApJ 517 (SN Ia cosmological standard candle foundation)
- UQFF 2.0 Architecture — ARCHITECTURE_FLOW_DIAGRAM.md v4.4.0 CANONICAL
- Session 88 — SPIRAL_SUPERNOVAE_UQFF_MODULE.cpp WOLFRAM_TERM: SPIRAL_SN_TENSION
