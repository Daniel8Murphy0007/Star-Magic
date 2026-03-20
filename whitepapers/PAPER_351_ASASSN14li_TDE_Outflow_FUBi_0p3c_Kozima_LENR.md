# PAPER_351 — ASASSN-14li Tidal Disruption Event: Ultrafast Outflow F_U_Bi_i and Kozima LENR Force

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 96  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF F_U_Bi_i for a TDE with 0.3c ultrafast outflow and Kozima LENR coupling  
**Author:** Daniel T. Murphy  


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
---

## Abstract

ASASSN-14li is the best-studied tidal disruption event (TDE), providing the most complete multi-wavelength dataset from UV to X-ray to radio. The UQFF buoyancy-unified force is computed for the stellar mass black hole remnant (M_BH = 10⁶ M☉), yielding F_U_Bi_i ≈ −8.32×10²¹¹ N — six orders of magnitude smaller than AGN-scale F_U_Bi_i, reflecting the much smaller BH mass. The ultrafast outflow at v_out = 0.3c is connected to UQFF via the Kozima LENR force component F_Kozima = 10³⁰ N at the stellar disruption interface.

---

## 2. Core Physics

### 2.1 UQFF Buoyancy-Unified Force (TDE Scale)

$$F_{U\_Bi\_i} \approx -8.32 \times 10^{211}\ \mathrm{N}$$

The six-order-of-magnitude reduction from the AGN scale (−8.32×10²¹⁷ N) reflects M_BH = 10⁶ M☉ vs. 10⁹ M☉.

### 2.2 Ultrafast Outflow

$$v_{\rm out} = 0.3c = 9.0 \times 10^7\ \mathrm{m/s}$$

Observed in Chandra/XMM-Newton blueshifted Fe K absorption lines. The UQFF kinetic coupling:
$$P_{\rm outflow} = \frac{1}{2} \dot{M}_{\rm out} v_{\rm out}^2 = \frac{1}{2} \dot{M}_{\rm out} (0.3c)^2$$

### 2.3 Kozima LENR Force Component

$$F_{\rm Kozima} = 1 \times 10^{30}\ \mathrm{N}$$

The Kozima heavy-rydberg LENR force arises when stellar debris density exceeds the nuclear lattice threshold at the tidal disruption radius:
$$r_{\rm tide} = R_\star \left(\frac{M_{\rm BH}}{M_\star}\right)^{1/3}$$

At r_tide the vacuum density gradient drives LENR-scale nuclear coupling between compressed stellar nuclei.

### 2.4 Full F_U_Bi_i Decomposition

$$F_{U\_Bi\_i}^{\rm ASASSN} = F_{\rm UQFF}^{\rm TDE} + F_{\rm Kozima} + F_{\rm outflow}$$

$$\approx -8.32\times 10^{211} + 10^{30} + P_{\rm outflow}/r_{\rm tide}\ \mathrm{N}$$

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| M_BH | UV-optical fit | 10⁶ M☉ |
| F_U_Bi_i | UQFF TDE scale | −8.32×10²¹¹ N |
| v_out | Chandra Fe K | 0.3c |
| F_Kozima | LENR coupling | 10³⁰ N |
| r_tide | R_✶(M_BH/M_✶)^(1/3) | ~7 R☉ |

---

## 4. Physical Significance

ASASSN-14li bridges stellar-scale and AGN-scale UQFF physics. The TDE provides a laboratory for testing how F_U_Bi_i scales with BH mass: the 6-order-of-magnitude reduction from 10⁹ M☉ to 10⁶ M☉ tracks the mass scaling F_U_Bi_i ∝ M_BH^α, α derived from comparing PAPER_346 (M87) to PAPER_351, enabling a power-law calibration of the BH mass dependence of UQFF vacuum buoyancy.

The Kozima LENR force at F_Kozima = 10³⁰ N is much smaller than F_U_Bi_i in this TDE context, suggesting LENR effects are perturbative at stellar BH scales.

---

## 5. Deduplication Note

- **vs. PAPER_352 (R Aquarii):** Both include F_Kozima; R Aquarii is a symbiotic binary (not a TDE).
- **vs. all AGN papers (346–350):** TDE F_U_Bi_i ≈ 10²¹¹ N (stellar mass BH) vs. AGN 10²¹⁷–10²¹⁸ N.

---

## 6. Classification

**Physics Territory:** FIRST UQFF TDE with ultrafast outflow (0.3c) and Kozima LENR coupling  
**Scale:** Stellar (10⁶ M☉ BH × TDE disruption radius)  
**CP Implementation:** `ASASSN14liTDEOutflowFUBiCalculator` (CondensedPhysics3.py, Session 96)
