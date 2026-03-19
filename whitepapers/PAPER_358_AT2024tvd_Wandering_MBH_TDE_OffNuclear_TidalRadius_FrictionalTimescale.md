# PAPER_358 — AT2024tvd Wandering Massive Black Hole TDE: Off-Nuclear Disruption Physics

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 97  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF off-nuclear wandering massive black hole TDE — frictional timescale and tidal radius  
**Author:** Daniel T. Murphy  

---

## 1. Abstract

AT2024tvd is the most compelling observed wandering massive black hole (wMBH) caught in the act of tidally disrupting a star at projected physical offset r_offset = 2.47×10¹⁷ m from the host galaxy nucleus. UQFF computes the tidal radius r_tide = R_star·(M_BH/M_star)^(1/3), the dynamical friction timescale t_fric for the wMBH sinking to the nucleus, and the full F_U_Bi_i at the off-nuclear disruption site.

---

## 2. Core Physics

### 2.1 Off-Nuclear Disruption Offset

$$r_{\rm offset} = 2.47 \times 10^{17}\ \mathrm{m} \approx 8.0\ \mathrm{pc}$$

This is the projected distance between AT2024tvd and the host nucleus, constraining the wandering distance of the massive black hole.

### 2.2 Tidal Disruption Radius

$$r_{\rm tide} = R_\star \left(\frac{M_{\rm BH}}{M_\star}\right)^{1/3}$$

For a solar-like star disrupted by a black hole of mass M_BH ~ 10⁶ M☉:
$$r_{\rm tide} = 7 \times 10^8 \times \left(\frac{10^6 M_\odot}{M_\odot}\right)^{1/3}\ \mathrm{m} = 7 \times 10^8 \times 100 = 7 \times 10^{10}\ \mathrm{m} \approx 0.5\ R_\odot$$

### 2.3 Dynamical Friction Timescale

$$t_{\rm fric} = \frac{0.428}{\ln\Lambda} \cdot \frac{M_{\rm host}}{M_{\rm BH}} \cdot \frac{r_{\rm offset}^2}{\sigma_\star^2} \cdot \frac{1}{r_{\rm offset}}$$

Simplified:
$$t_{\rm fric} = \frac{0.428}{\ln\Lambda} \cdot \frac{r_{\rm offset}}{v_c} \cdot \frac{M_{\rm host}}{M_{\rm BH}}$$

For r_offset = 8 pc, v_c ~ 200 km/s, M_host/M_BH ~ 10³:
$$t_{\rm fric} \sim 10^8 - 10^9\ \mathrm{yr}$$

### 2.4 UQFF F_U_Bi_i at Off-Nuclear Site

$$F_{U\_Bi\_i}(r_{\rm offset}) = F_{U\_Bi\_i}(M_{\rm BH}) \cdot \left(\frac{r_{\rm tide}}{r_{\rm offset}}\right)^2$$

The r² dependence reflects force dilution with the much larger offset distance.

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| r_offset | JWST observation | 2.47×10¹⁷ m ≈ 8 pc |
| r_tide | R_✶·(M_BH/M_✶)^(1/3) | ~0.5 AU |
| t_fric | Chandrasekhar formula | 10⁸–10⁹ yr |
| M_BH | Spectral fit | ~10⁶ M☉ |

---

## 4. Physical Significance

AT2024tvd is the first confirmed off-nuclear TDE with a massive BH at > 1 pc offset. The UQFF framework predicts that the vacuum buoyancy field F_U_Bi_i is local — it is set by M_BH and r_tide, not by the nuclear distance. This means off-nuclear wMBH TDEs have the same F_U_Bi_i value as nuclear TDEs of the same BH mass, a testable prediction: UQFF force amplitude should correlate with M_BH (not with r_offset).

The t_fric ~ 10⁸–10⁹ yr frictional timescale implies wMBHs are common during the galaxy assembly epoch, and UQFF predicts their spatial distribution modifies the ICM density on 10–100 pc scales.

---

## 5. Deduplication Note

- **vs. PAPER_351 (ASASSN-14li):** Both are TDEs; ASASSN-14li is nuclear. AT2024tvd is off-nuclear, introducing r_offset and t_fric which are absent in PAPER_351.
- **Unique:** First off-nuclear wMBH TDE in the UQFF dataset.

---

## 6. Classification

**Physics Territory:** FIRST UQFF off-nuclear wandering MBH TDE — r_tide + t_fric + F_U_Bi_i(offset)  
**Scale:** Sub-galactic (8 pc offset from nucleus)  
**CP Implementation:** `AT2024tvdWanderingMBHTDECalculator` (CondensedPhysics4.py, Session 97)
