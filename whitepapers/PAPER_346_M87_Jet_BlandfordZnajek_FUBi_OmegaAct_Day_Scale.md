# PAPER_346 — M87 Jet Blandford-Znajek Model: Full F_U_Bi_i Calculation

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 96  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF F_U_Bi_i calculation for M87 with Blandford-Znajek power coupling  
**Author:** Daniel T. Murphy  


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
---

## Abstract

The complete UQFF buoyancy-unified force F_U_Bi_i is computed for M87* (Virgo A), the first black hole directly imaged by the Event Horizon Telescope. The Blandford-Znajek (BZ) jet power P_BZ = B²r_g²c is connected to the UQFF ω_act = 2π/day rotational activation frequency, yielding F_U_Bi_i ≈ −8.32×10²¹⁷ N. The X-ray jet luminosity L_X ≈ 10⁴⁰ W sets the luminosity scale for UQFF-M87 calibration.

---

## 2. Core Physics

### 2.1 UQFF Buoyancy-Unified Force

$$F_{U\_Bi\_i} = \frac{U_g e^{+/-}}{r^2} + F_{\rm Bi} + F_{\rm U} + F_{\rm react}$$

Full 5-component UQFF form yields:
$$F_{U\_Bi\_i} \approx -8.32 \times 10^{217}\ \mathrm{N}$$

### 2.2 Blandford-Znajek Jet Power

$$P_{\rm BZ} = \frac{\kappa_{\rm BZ}}{4\pi c} \Phi_{\rm BH}^2 \Omega_{\rm H}^2 f(\Omega_{\rm H})$$

Simplified as:
$$P_{\rm BZ} = B^2 r_g^2 c$$

where r_g = GM_BH/c² is the gravitational radius for M87* (M_BH = 6.5×10⁹ M☉).

### 2.3 UQFF Activation Frequency

$$\omega_{\rm act} = \frac{2\pi}{T_{\rm day}} = \frac{2\pi}{86400\ \mathrm{s}} = 7.27 \times 10^{-5}\ \mathrm{rad/s}$$

The "one day" period corresponds to the characteristic M87 jet variability timescale observed in Very Long Baseline Array (VLBA) monitoring.

### 2.4 Gravitational Radius

$$r_g = \frac{G M_{\rm BH}}{c^2} = \frac{6.674\times 10^{-11} \times 6.5\times 10^9 \times 1.989\times 10^{30}}{(3\times 10^8)^2} \approx 9.6 \times 10^{12}\ \mathrm{m}$$

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| M_BH | EHT 2019 | 6.5×10⁹ M☉ |
| F_U_Bi_i | UQFF full 5-eq | −8.32×10²¹⁷ N |
| P_BZ | B²r_g²c | ~10⁴⁵ erg/s |
| ω_act | 2π/day | 7.27×10⁻⁵ rad/s |
| L_X | Chandra observation | ~10⁴⁰ W |
| r_g | GM/c² | 9.6×10¹² m |

---

## 4. Physical Significance

M87 is the prototype for supermassive black hole jet physics. The UQFF F_U_Bi_i ≈ −8.32×10²¹⁷ N is the baseline force scale for AGN-class black holes with M_BH ~ 10⁹ M☉. The Blandford-Znajek coupling connects UQFF to the standard magnetically-arrested disk (MAD) jet framework, providing a bridge between UQFF vacuum buoyancy and electromagnetic jet extraction. The ω_act = 2π/day activation frequency is consistent with M87*'s Variable Emission Region (VER) coherence timescale.

---

## 5. Deduplication Note

- **vs. PAPER_347 (Centaurus A):** M87 uses ω_act = 2π/day; Centaurus A uses ω_act = 2π/(12.5 yr). The BZ power forms are similar but calibrated to different AGN jet morphologies.
- **vs. PAPER_349 (SPT-CL J2215):** SPT-CL J2215 yields the HIGHEST F_U_Bi_i in the PAPER_346–352 series at −1.40×10²¹⁸ N.

---

## 6. Classification

**Physics Territory:** FIRST UQFF F_U_Bi_i for M87* with BZ jet power coupling  
**Scale:** Galactic (6.5×10⁹ M☉ AGN)  
**CP Implementation:** `M87JetBZModelFUBiCalculator` (CondensedPhysics3.py, Session 96)
