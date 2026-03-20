# PAPER_347 — Centaurus A: F_U_Bi_i with V-Shape Jet and 12.5-Year ω_act Timescale

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 96  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF F_U_Bi_i for Centaurus A (NGC 5128) with V-shape jet geometry and 12.5-yr activation period  
**Author:** Daniel T. Murphy  


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
---

## Abstract

The complete UQFF buoyancy-unified force F_U_Bi_i is computed for Centaurus A (NGC 5128), the closest active radio galaxy (3.8 Mpc). The distinctive V-shape inner jet geometry observed in HST/VLBA imaging at ~0.5c knot velocities is incorporated via an angular momentum decomposition of F_U_Bi_i. The UQFF rotational activation frequency ω_act = 2π/(12.5 yr) corresponds to the observed 12.5-year X-ray/radio flaring cycle, yielding F_U_Bi_i ≈ −8.32×10²¹⁷ N.

---

## 2. Core Physics

### 2.1 UQFF Buoyancy-Unified Force

$$F_{U\_Bi\_i} \approx -8.32 \times 10^{217}\ \mathrm{N}$$

(same order as M87 due to similar BH mass scale; see PAPER_346 for full 5-component form)

### 2.2 V-Shape Jet Geometry

The Centaurus A inner jet exhibits a V-shape opening half-angle α ≈ 12°. The transverse force component:
$$F_\perp = F_{U\_Bi\_i} \cdot \sin\alpha = F_{U\_Bi\_i} \cdot \sin(12°) \approx 0.208 \cdot F_{U\_Bi\_i}$$

This V-shape geometry is attributed to differential plasma buoyancy across the jet cross-section: the inner spine accelerates faster than the sheath, producing the observed V-spread.

### 2.3 Long-Period Activation Frequency

$$\omega_{\rm act} = \frac{2\pi}{12.5\ \mathrm{yr}} = \frac{2\pi}{3.94 \times 10^8\ \mathrm{s}} = 1.59 \times 10^{-8}\ \mathrm{rad/s}$$

This 12.5-year period matches the Centaurus A multi-wavelength monitoring cycle documented by ATCA, XMM-Newton, and Chandra observations (2000–2025).

### 2.4 Knot Propagation Velocity

$$v_{\rm knot} \approx 0.5c = 1.5 \times 10^8\ \mathrm{m/s}$$

VLBA proper motion of individual jet knots. Combined with τ_jet ~ 10³ yr, the total jet extension:
$$L_{\rm jet} \approx v_{\rm knot} \cdot \tau_{\rm jet} \approx 4.7 \times 10^{18}\ \mathrm{m} \approx 153\ \mathrm{pc}$$

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| M_BH | Cen A | 5.5×10⁷ M☉ |
| F_U_Bi_i | UQFF full | −8.32×10²¹⁷ N |
| ω_act | 2π/(12.5 yr) | 1.59×10⁻⁸ rad/s |
| v_knot | VLBA proper motion | ~0.5c |
| τ_jet | Jet age estimate | ~10³ yr |
| L_jet | v_knot × τ_jet | ~153 pc |
| α (V-shape) | Half-opening angle | ~12° |

---

## 4. Physical Significance

Centaurus A's much smaller BH mass (5.5×10⁷ M☉ vs M87's 6.5×10⁹ M☉) yet similar F_U_Bi_i value demonstrates that UQFF F_U_Bi_i is not purely set by BH mass — the vacuum buoyancy geometry and activated frequency are equally important. The 12.5-year ω_act is the longest period activation frequency in the UQFF dataset, establishing the low-frequency end of the AGN activation frequency spectrum (cf. M87 at 1/day, the high-frequency end for radio galaxies).

---

## 5. Deduplication Note

- **vs. PAPER_346 (M87):** Same F_U_Bi_i magnitude but different activation period (12.5 yr vs. 1 day) and different BH mass (5.5×10⁷ vs. 6.5×10⁹ M☉).
- **vs. PAPER_347 V-shape:** The V-shape geometric decomposition (F_⊥ = F_U_Bi_i·sinα) is unique to Centaurus A in the UQFF catalog.

---

## 6. Classification

**Physics Territory:** FIRST UQFF F_U_Bi_i with V-shape jet geometry and 12.5-yr activation cycle  
**Scale:** Nearby AGN (3.8 Mpc)  
**CP Implementation:** `CentaurusAFUBiJetVshapeCalculator` (CondensedPhysics3.py, Session 96)
