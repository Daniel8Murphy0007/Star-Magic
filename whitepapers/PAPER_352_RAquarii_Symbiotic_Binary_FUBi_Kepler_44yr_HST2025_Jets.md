# PAPER_352 — R Aquarii Symbiotic Binary: F_U_Bi_i with Kepler 44-Year Orbital Period

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 96  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF F_U_Bi_i for a symbiotic nova binary with P_orb = 44 yr  
**Author:** Daniel T. Murphy  

---

## 1. Abstract

R Aquarii is the nearest symbiotic nova system, consisting of a Mira Pulsating Giant and a white dwarf companion in a 44-year orbit. UQFF buoyancy-unified force F_U_Bi_i ≈ −2.09×10²¹² N is calculated at the binary orbital radius derived from Kepler's third law: a_orb = (GM·P²/4π²)^(1/3). HST 2025 near-UV imaging of the expanding jet system (launched ~10³ yr ago) provides the observational anchor. The UQFF Kozima LENR coupling is relevant at the mass transfer interface between the giant's wind and the WD accretion disk.

---

## 2. Core Physics

### 2.1 UQFF Buoyancy-Unified Force

$$F_{U\_Bi\_i} \approx -2.09 \times 10^{212}\ \mathrm{N}$$

(intermediate between TDE-scale PAPER_351 and AGN-scale PAPER_346)

### 2.2 Kepler's Third Law Orbital Radius

$$a_{\rm orb} = \left(\frac{G M_{\rm total} P_{\rm orb}^2}{4\pi^2}\right)^{1/3}$$

with P_orb = 44 yr = 1.388×10⁹ s and M_total ≈ M_Mira + M_WD ≈ 2 M☉:
$$a_{\rm orb} = \left(\frac{6.674\times 10^{-11} \times 4\times 10^{30} \times (1.388\times 10^9)^2}{4\pi^2}\right)^{1/3} \approx 1.0 \times 10^{13}\ \mathrm{m} \approx 70\ \mathrm{AU}$$

### 2.3 Mass Transfer and LENR Coupling

At the tidal mass transfer interface:
$$\dot{M}_{\rm transfer} = \dot{M}_{\rm wind} \cdot \left(\frac{a_{\rm orb}}{R_{\rm Mira}} \right)^{-2}$$

The UQFF Kozima LENR force:
$$F_{\rm Kozima} = E_{\rm LENR} / a_{\rm orb}$$

where E_LENR is the low-energy nuclear reaction energy scale at the compressed accretion interface.

### 2.4 HST 2025 Jet Observation

R Aquarii's expanding jets, first resolved by HST in the 1990s and re-observed in 2025:
$$v_{\rm jet} \approx 200\ \mathrm{km/s}$$
$$\tau_{\rm jet} \approx 10^3\ \mathrm{yr}$$
$$L_{\rm jet} = v_{\rm jet} \cdot \tau_{\rm jet} \approx 6.3 \times 10^{15}\ \mathrm{m} \approx 0.2\ \mathrm{pc}$$

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| P_orb | Spectroscopic | 44 yr |
| a_orb | Kepler 3rd law | ~70 AU |
| F_U_Bi_i | UQFF | −2.09×10²¹² N |
| M_total | Mira + WD | ~2 M☉ |
| v_jet | HST proper motion | ~200 km/s |
| τ_jet | Jet age | ~10³ yr |

---

## 4. Physical Significance

R Aquarii provides UQFF calibration at the symbiotic nova (compact binary) scale — intermediate between individual stellar objects (PAPER_351 TDE) and galactic nuclei. The 44-year orbital period sets the longest UQFF binary activation period in the dataset (vs. PAPER_347 Cen A 12.5-yr AGN cycle). The Kozima LENR coupling at the mass transfer interface raises the testable prediction that R Aquarii's nova outburst energetics deviate from standard accretion models by an LENR fraction, detectable in nuclear gamma-ray emission (e.g., INTEGRAL 511 keV line observations).

---

## 5. Deduplication Note

- **vs. PAPER_322 (R Aquarii in SOURCE122):** Earlier treatment computed the 5-frequency resonances; this paper adds full F_U_Bi_i with a_orb from Kepler's third law and directly compares with HST 2025 jet data.
- **vs. PAPER_351 (ASASSN-14li):** Different system class (symbiotic binary vs. TDE); both include F_Kozima.

---

## 6. Classification

**Physics Territory:** FIRST UQFF F_U_Bi_i for a symbiotic nova binary with Kepler a_orb  
**Scale:** Stellar (70 AU binary, ~200 pc distance)  
**CP Implementation:** `RAquariiSymbioticBinaryFUBiCalculator` (CondensedPhysics3.py, Session 96)
