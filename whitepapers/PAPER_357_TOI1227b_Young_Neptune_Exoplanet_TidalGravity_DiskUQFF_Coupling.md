# PAPER_357 — TOI-1227b: Young Neptune Exoplanet with Tidal Gravity and Disk-UQFF Coupling

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 97  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF calculation for a young exoplanet (T_age = 11 Myr) with tidal + disk force coupling  
**Author:** Daniel T. Murphy  


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
---

## Abstract

TOI-1227b is a rare young sub-Neptune (T_age = 11 Myr, P_orb = 11 days) still embedded in the debris disk of its M-dwarf host. UQFF computes the tidal gravitational acceleration g_tide = GM_star/a_orb² at the orbital radius, a disk-UQFF force F_disk = ρ_disk·v_disk²·(1+H₀t)·SC_m incorporating Hubble expansion and superconductive modifier, and the full F_U_Bi_i. TOI-1227b provides a benchmark for UQFF at typical planetary mass and orbital scales during the planet formation epoch.

---

## 2. Core Physics

### 2.1 Tidal Gravitational Acceleration

$$g_{\rm tide} = \frac{G M_\star}{a_{\rm orb}^2}$$

At a_orb from P_orb = 11 days via Kepler's third law:
$$a_{\rm orb} = \left(\frac{G M_\star P_{\rm orb}^2}{4\pi^2}\right)^{1/3}$$

For M_star ≈ 0.17 M☉ (M-dwarf), P_orb = 11 d:
$$a_{\rm orb} \approx 0.05\ \mathrm{AU} = 7.5 \times 10^9\ \mathrm{m}$$
$$g_{\rm tide} \approx \frac{6.674\times 10^{-11} \times 0.17 \times 1.989\times 10^{30}}{(7.5\times 10^9)^2} \approx 0.40\ \mathrm{m/s}^2$$

### 2.2 Disk-UQFF Force Coupling

$$F_{\rm disk} = \rho_{\rm disk} \cdot v_{\rm disk}^2 \cdot (1 + H_0 t) \cdot {\rm SC}_m$$

where:
- ρ_disk = protoplanetary disk density at a_orb (~10⁻⁹ kg/m³ at 11 Myr)
- v_disk = disk gas velocity (~1 km/s at 0.05 AU)
- SC_m = superconductive modifier
- (1 + H₀t) = Hubble expansion factor over 11 Myr

For t = 11 Myr = 3.47×10¹⁴ s:
$$(1 + H_0 t) \approx 1 + 2.27\times 10^{-18} \times 3.47\times 10^{14} \approx 1.0008$$

### 2.3 UQFF F_U_Bi_i (Planetary Scale)

$$F_{U\_Bi\_i}^{\rm planet} = F_{U\_Bi\_i}(M_{\rm planet}, r_{\rm orbit}, \omega_{\rm act})$$

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| T_age | TESS age estimate | 11 Myr |
| P_orb | TESS + RV | 11 days |
| a_orb | Kepler 3rd law | ~0.05 AU |
| g_tide | GM_star/a² | ~0.40 m/s² |
| F_disk | ρ_disk·v²·SC_m | disk-epoch force |
| (1+H₀t) | Hubble correction | 1.0008 |

---

## 4. Physical Significance

TOI-1227b is exceptional because it is young enough that the UQFF disk-coupling term F_disk is non-negligible: the protoplanetary disk provides a dense medium that couples the SC_m modifier to the local vacuum field. This is impossible for older planetary systems where the disk has dissipated. The UQFF prediction is that disk-embedded planets during the 1–100 Myr epoch receive a systematic F_disk force that contributes at the ~0.1% level to their orbital energy budget, producing a measurable shift in transit timing variations (TTVs) detectable by PLATO/Ariel.

---

## 5. Deduplication Note

- **vs. PAPER_357 vs. Solar System UQFF papers:** All prior UQFF planet papers used mature (>1 Gyr) systems; TOI-1227b at 11 Myr is the youngest planet in the UQFF dataset.
- **Unique:** Disk coupling F_disk = ρ_disk·v²·(1+H₀t)·SC_m is unique to young-planet UQFF.

---

## 6. Classification

**Physics Territory:** FIRST UQFF young exoplanet (11 Myr) with disk-UQFF SC_m coupling  
**Scale:** Planetary (sub-Neptune, 0.05 AU orbit)  
**CP Implementation:** `TOI1227bYoungNeptuneExoplanetFUBiCalculator` (CondensedPhysics4.py, Session 97)
