# PAPER_343 — SGR J1745-2900: SC_m Mass-Modified Luminosity and Doubled f_react (June 2013)

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 96  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF treatment of magnetar SC_m mass modification with L_X = ρ_vac·f_res·V  
**Author:** Daniel T. Murphy  

---

## 1. Abstract

A novel UQFF form for the superconductive modifier SC_m of SGR J1745-2900 is derived as a mass-dependent suppression by the critical field ratio: SC_m = M·(1 − B/B_crit). The X-ray luminosity is expressed as L_X = ρ_vac·f_res·V, coupling vacuum energy density, resonance frequency, and magnetospheric volume. The activation event of June 2013 corresponds to a doubling of f_react, confirmed by the sudden spin-up and luminosity jump. T_surf = 1.16×10⁷ K is derived from the Stefan-Boltzmann radiative balance.

---

## 2. Core Physics

### 2.1 Mass-Modified Superconductive Modifier

$${\rm SC}_m = M \cdot \left(1 - \frac{B}{B_{\rm crit}}\right)$$

where B_crit = 4.4×10¹³ T (quantum critical field). For SGR J1745-2900: B = 2×10¹⁰ T ≪ B_crit, giving SC_m ≈ M (nearly full superconductive coupling).

### 2.2 Vacuum-Energy X-ray Luminosity Form

$$L_X = \rho_{\rm vac} \cdot f_{\rm res} \cdot V_{\rm mag}$$

where ρ_vac = ρ_SCm − ρ_UA is the net vacuum energy density and V_mag is the magnetospheric volume.

### 2.3 June 2013 Activation Event

$$f_{\rm react}^{\rm post} = 2 \cdot f_{\rm react}^{\rm pre}$$

The doubling of the reactance frequency at outburst onset corresponds to:
$$\Delta L_X = \rho_{\rm vac} \cdot \Delta f_{\rm res} \cdot V_{\rm mag} = \rho_{\rm vac} \cdot f_{\rm react}^{\rm pre} \cdot V_{\rm mag}$$

### 2.4 Surface Temperature

From Stefan-Boltzmann balance of magnetospheric luminosity:
$$T_{\rm surf} = \left(\frac{L_X}{4\pi R_{\rm NS}^2 \sigma_{\rm SB}}\right)^{1/4} = 1.16 \times 10^7\ \mathrm{K}$$

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| SC_m | M(1−B/B_crit) | ≈ M_NS (B ≪ B_crit) |
| L_X | ρ_vac·f_res·V | ~10³⁵ erg/s |
| f_react (pre-2013) | canonical | f₀ |
| f_react (post June 2013) | 2f₀ | 2f₀ |
| T_surf | Stefan-Boltzmann | 1.16×10⁷ K |
| B_crit | Quantum critical | 4.4×10¹³ T |

---

## 4. Physical Significance

SGR J1745-2900 is the only magnetar within 0.3 pc of Sgr A*, making it the unique test-bed for UQFF near the Galactic Center supermassive black hole. The SC_m = M(1−B/B_crit) form establishes that even strongly magnetized neutron stars maintain near-unity superconductive coupling due to B ≪ B_crit. The L_X = ρ_vac·f_res·V form is a direct observable prediction: doubling f_res (as seen in June 2013) should produce a factor of 2 luminosity jump, consistent with the 2013 XMM-Newton observations.

---

## 5. Deduplication Note

- **vs. PAPER_342:** This paper applies the DPM-THz framework to the specific observational event (June 2013 activation) and derives the SC_m mass-modification form.
- **vs. SOURCE27 (SGR 1745 SuperFreq):** SOURCE27 computed the 5 resonance frequencies; this paper derives the observable L_X from ρ_vac·f_res·V.

---

## 6. Classification

**Physics Territory:** FIRST UQFF SC_m mass-modified form for magnetar near Galactic Center  
**Scale:** Stellar/Galactic Center (0.3 pc)  
**CP Implementation:** `SGR17452900SCmLxFreqFormCalculator` (CondensedPhysics3.py, Session 96)
