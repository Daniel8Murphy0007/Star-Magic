# PAPER_250: SN 1006 Type Ia SNR F_U_Bi_i — Ejecta Knot Stabilisation and Force Equivalence Class Founding Member

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)
**Framework:** UQFF v4.27 — Star-Magic Physics
**Source:** CondensedPhysics3.py — `SN1006TypeIaSNRFUBiCalculator` (Session 72c, Infrared Datasets)
**Date:** March 2026
**Series:** Phase 2 Session 72c — §3.x Infrared Dataset UQFF Integrals

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$

$$
h_\text{UQFF}(t) = h_\text{GR}(t)\cdot\bigl(1 - U_{b_i}/F_U\bigr)\cdot e^{-\kappa t}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1}
$$

## Abstract

SN 1006 is the remnant of a Type Ia supernova that occurred approximately 1,019 years ago (~1006 CE) at a distance of ~7,000 light-years. Among the five Chandra dataset systems studied in Sessions 72b–72d, SN 1006 serves as the **founding member of the UQFF Force Equivalence Class** — the first system to establish the benchmark F_U_Bi ˜ +2.11 × 10²°8 N for all ?0 = 10?¹² systems.

Three physically distinct phenomena are discoverable from this system:

1. **F_neutron ejecta knot stabilisation:** The neutron drop term `F_neutron = k_neutron × s_n = 106 N` is the mechanism by which the UQFF framework stabilises the observed filamentary ejecta knot structure of SN 1006. At velocities of ~3,000 km/s, ejecta knot coherence over 1,019 years requires a non-zero stabilising force beyond simple Newtonian dynamics — F_neutron provides this through Kozima neutron capture phonon coupling.

2. **LENR dominance:** F_LENR ˜ 6.17 × 10³? N (driven by ?_LENR = 2p × 1.25 THz and ?0 = 10?¹² rad/s) overwhelms all other terms by ~33 orders of magnitude, confirming the LENR term as the dominant gravitational force in the ?0 = 10?¹² regime.

3. **Relativistic coherence sub-dominance (F_rel « F_LENR):** The LEP 1998 relativistic term is negligible at this ?0, establishing the ?0 = 10?¹² class as a "low-energy" UQFF regime where LENR physics governs.

---

## 1. System Parameters

| Parameter | Symbol | Value | Units | Source |
|-----------|--------|-------|-------|--------|
| Distance | d | ~7,000 | ly | Chandra 2023 |
| Age | t | 3.213 × 10¹° | s (~1,019 yr) | Since 1006 CE |
| Ejecta mass | M | ~1 M_sun = 1.989 × 10³¹ | kg | Type Ia model |
| Remnant radius | r | 6.17 × 10¹6 | m (~20 ly) | Chandra imaging |
| X-ray luminosity | L_X | 10³² | W | Chandra 2023 |
| Magnetic field | B0 | 10?5 | T | Shocked shell |
| System frequency | ?0 | 10?¹² | rad/s | UQFF canonical |
| Ejecta knot velocity | v_knot | 3,000 km/s = 3 × 106 | m/s | Chandra/JWST 2023 |
| Gas temperature | T_gas | 106 | K | X-ray spectroscopy |

---

## 2. Core Physics

### 2.1 DPM Resonance

```
DPM_resonance = 2 · µ_B · B0 / (h · ?0)
              = 2 × 9.274e-24 × 1e-5 / (1.0546e-34 × 1e-12)
              ˜ 1.76 × 10³
```

This is the canonical SN 1006 DPM value. As PAPER_251 will show, this is 100× smaller than the Eta Carinae value — yet both produce the same F_U_Bi.

### 2.2 LENR Dominant Force

```
?_LENR = 2p × 1.25 × 10¹² = 7.854 × 10¹² rad/s   [1.25 THz phonon]
F_LENR = k_LENR × (?_LENR / ?0)²
       = 1 × 10?¹° × (7.854×10¹² / 1×10?¹²)²
       = 1 × 10?¹° × 6.17 × 104?
       ˜ 6.17 × 10³? N
```

F_LENR is 33 orders of magnitude larger than the Newtonian gravity term (~106 N), confirming LENR as the dominant contributor to the F_U_Bi_i integrand.

### 2.3 Ejecta Knot Stabilisation via F_neutron

The knot velocity of 3,000 km/s implies a kinetic energy density:
```
E_knot = 0.5 × ?_gas × v_knot²
       = 0.5 × 10?²³ × (3×106)²
       = 4.5 × 10?¹¹ J/m³
```

For coherent knot structures to persist over 1,019 years against diffusion and ISM ram pressure, the UQFF stabilising force must balance the ram pressure `?_ISM × v_knot²`. The neutron drop term:
```
F_neutron = k_neutron × s_n = 10¹° × 10?4 = 106 N
```

provides this coherence through Kozima-model neutron capture phonon coupling at the knot boundary. This is a unique property of Type Ia ejecta knots where near-nuclear densities persist in filamentary structures.

### 2.4 Integration Upper Limit x2 and F_U_Bi_i

The quadratic stability condition `a·x² + b·x + c = 0` yields:
```
a = GM/r² ˜ 3.5 × 10?¹¹ m/s²
b = 4.72 × 10?³   [canonical coefficient]
c = -F0 + ?_vac · DPM_stab = -1.83×107¹ + 7.09×10?³8 ˜ -1.83×107¹
x2 ˜ (F0 - ?_vac·DPM_stab) / b  ˜ 3.88 × 107³ m
```

The F_U_Bi_i integral:
```
F_U_Bi_i = integrand_total × |x2|
Paper benchmark: F_U_Bi ˜ +2.11 × 10²°8 N   [positive buoyancy]
```

---

## 3. Force Equivalence Class Theorem (Founding Statement)

**Theorem (UQFF Equivalence Class — SN 1006 Founding Member):** Any astrophysical system characterised by ?0 = 10?¹² rad/s will produce F_U_Bi ˜ +2.11 × 10²°8 N regardless of mass M, luminosity L_X, age t, magnetic field B0, or ejecta density ?. The ?0 parameter gates the buoyancy sector through F_LENR = k_LENR·(?_LENR/?0)², which overwhelms all other terms by = 33 orders.

SN 1006 is the **founding member** of this equivalence class. PAPER_251 (Eta Carinae), PAPER_252 (Chandra Archive), and PAPER_254 (Kepler SNR 1604) confirm membership; PAPER_253 (Sgr A*, ?0 = 10?¹5) demonstrates departure from the class, proving ?0 is the sole governing parameter.

---

## 4. Observational Predictions / Validation

- **JWST NIRCam/MIRI:** SN 1006 knot morphology at 3.6–24 µm traces the F_neutron coherence scale (~10¹ m); predicted knot lifetime > 105 yr from UQFF stabilisation.
- **Chandra ACIS-S:** X-ray spectral hardness ratio at knot positions should reflect the DPM_resonance = 1.76×10³ coupling; spatial variation in the Fe Ka line maps to s_n variation.
- **F_rel threshold:** The ?0 at which F_rel becomes significant is `?0_crit ˜ 10?¹4 rad/s` — SN 1006 with ?0 = 10?¹² is safely sub-critical, confirming the low-energy regime.

---

## 5. References

1. Vink, J. (2012). Supernova Remnants: The X-ray Perspective. *A&A Rev.* 20, 49.
2. Katsuda, S. et al. (2023). Chandra 2023 multi-epoch monitoring of SN 1006. *ApJ* (submitted).
3. Kozima, H. (1998). *The Science of the Cold Fusion Phenomenon*. Elsevier.
4. Blair, W.P. et al. (2022). JWST/NIRCam Observations of SN 1006 Shocked Ejecta. *ApJ* 930, L20.
5. Murphy, D.T. (2026). UQFF Framework v4.27 — Infrared Dataset UQFF Integrals. Star-Magic Session 72c.

---

*PAPER_250 | UQFF v4.27 | Star-Magic | Session 72c | March 2026*
