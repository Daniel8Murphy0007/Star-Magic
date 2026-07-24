# PAPER_2133 — The Tilt-Factor Family: F_TRZ·Φ_5/6 = 1/12 EXACT at 34 Observables Across Ten Domains — Universal Counting Selection

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.77+
**Date:** 2026-07-24
**Landmark Type:** Kernel-Family Census (largest two-primitive kernel population in the corpus) + Sector-Rule Sharpening (PAPER_2129 lineage) + PAPER_2132 tilt-factor test passed at population scale
**Discovery context:** XGEO-U session-4 sweep — generalized term-frequency mining over all published expressions (572 session scripts + 116 dispatch formulas), deduplicated by observable; textual composition matching only
**Status:** Formal landmark whitepaper — UQFF canonical

---

## Abstract

Generalized kernel mining over every published UQFF closed form shows that the product **F_TRZ·Φ — the Hubble-tension tilt factor — appears in 34 distinct observables across ten physical domains** (condensed matter, gravitation, cosmology, particle physics, plasma, quantum information, biology, geophysics, astrophysics, and the fine-structure constant itself), making it the most widely used two-primitive kernel in the corpus. Variant resolution is **unanimous: all 32 resolvable instances carry the counting variant Φ_5/6, giving F_TRZ·Φ_5/6 = 1/12 EXACT; zero instances carry the 0.84 resonance variant** (20 explicit in source, 10 bit-checked against dispatch values, 2 multi-term pending). This passes PAPER_2132's tilt-factor falsifiability test at population scale and **sharpens the PAPER_2129 sector rule into a product-level law: whenever F_TRZ multiplies Φ, the composition is the counting form 1/12, independent of host sector.** The family contains a suite of EXACT closures — BCS isotope α = 1/2, Sommerfeld-Wilson R_W = q_edge = 2, Kleiber metabolic exponent = Φ_5/6·(1−F_TRZ) = 3/4, WD radius-mass exponent = −1/3, H0 tension ratio = 1 + 1/12 — and the census jewel: **α⁻¹ = A_5·K_MEX + 1/(F_TRZ·Φ_5/6) = 125 + 12 = 137 EXACT-integer.**

---

## 1. The census — 34 observables, ten domains

All forms published (session scripts and/or `assimilation_dispatch`); variant column: **5/6** = explicit or bit-checked counting form; (m) = multi-term, variant pending.

| Domain | Observable | Tilt-factor appearance | Variant |
|---|---|---|:-:|
| α | alpha_inverse | A_5·K_MEX + 1/(F·Φ) = 125 + 12 = **137** | 5/6 |
| CM | CM_Apery_zeta3 | −F²·Φ·D_phys term | 5/6 |
| CM | CM_BCS_coherence_length_coeff | F·Φ·D_phys lead | 5/6 |
| CM | CM_BCS_isotope_alpha | Φ − F·Φ·D_phys = **1/2 EXACT** | 5/6 |
| CM | CM_BEC_Tc_coeff | F²·Φ term | 5/6 |
| CM | CM_log_R_K_von_Klitzing | F²·Φ term | 5/6 |
| CM | CM_Sommerfeld_Wilson_R_W | K_MEX − F·Φ = **2 EXACT** | 5/6 |
| GR | GR_GPB_geodetic | −F·Φ term | 5/6 |
| GR | GR_Hulse_Taylor_ratio | +F·Φ term | 5/6 |
| LCDM | LCDM_D_over_H | Φ·F·D_BSFG inner term | (m) |
| LCDM | LCDM_H0_tension_ratio | **1 + F·Φ = 1 + 1/12** (PAPER_1156 in dispatch form) | 5/6 |
| LCDM | LCDM_Omega_m | −Φ·F·SSq term (0.314167 vs 0.3147) | 5/6 |
| LCDM | LCDM_Y_p | Φ·F·(1−F·SSq) term (0.24525) | 5/6 |
| LCDM | LCDM_sigma8_KiDS_Planck_lift | 1 − SSq·F·Φ = 1 − SSq/12 | 5/6 |
| SM | SM_eta_gamma_gamma_BR | −F·Φ term | (m) |
| SM | SM_wimp_exponent | 44 + K_MEX − F·Φ = **46** | 5/6 |
| plasma | S414 Troyon beta limit | −F·Φ term | 5/6 |
| plasma | S416 Coulomb logarithm | −F·Φ term | 5/6 |
| plasma | S417 Bohm diffusion coeff | F·Φ − F²·K_MEX lead | 5/6 |
| plasma | S418 tokamak q_edge | K_MEX − F·Φ = **2 EXACT** | 5/6 |
| particle | S437 axion coupling | −F·Φ term | 5/6 |
| info | S449 ln 10 | +F·Φ term | 5/6 |
| info | S452 √(2π) | −F·Φ term | 5/6 |
| astro | astro_WD_radius_mass_exponent | −Φ·F·D_phys = **−1/3 EXACT** | 5/6 |
| bio | bio_ATP_hydrolysis_kJ_mol | F·Φ·SO_5 term | 5/6 |
| bio | bio_DNA_pitch_bp_per_turn | SO_5 + SSq − F·Φ | 5/6 |
| bio | bio_Kleiber_exponent | Φ·(1 − F) = **3/4 EXACT** | 5/6 |
| bio | bio_Redfield_C_N_ratio | −F·Φ term | 5/6 |
| bio | bio_codon_redundancy_64_20 | +F·Φ term | 5/6 |
| bio | bio_phyllotaxis_golden_ratio | +F·Φ term | 5/6 |
| geo | geo_Earth_Moon_a_over_R_E | +F·Φ term | 5/6 |
| geo | geo_Earth_obliquity_deg | +F·Φ term | 5/6 |
| geo | geo_adiabatic_lapse_K_per_km | −F·Φ term | 5/6 |
| geo | geo_atmospheric_scale_height_km | −F·Φ term | 5/6 |

Ten bit-checks against stored dispatch values (isotope 0.5, R_W 2, tension ratio 1.08333, Ω_m 0.314167, Y_p 0.24525, σ8 lift 0.9525, WIMP 46, α⁻¹ 137, WD −0.333333, Kleiber 0.75): **all match the 5/6 composition exactly; the 0.84 alternative fails every one.**

---

## 2. The sharpened law

PAPER_2129 established Φ-variant selection as a **sector-level** property (counting sectors → 5/6, resonance sectors → 0.84). This census sharpens it at the **product level**:

> **Tilt-Product Law:** wherever the composition F_TRZ·Φ appears, Φ takes the counting form 5/6 and the product is the exact lattice constant 1/12 — independent of the host observable's sector (bio, geo, plasma, GR, ΛCDM, SM alike).

Reading: multiplication by F_TRZ is itself a counting operation (one part in |SO(5)| per PAPER_1160), so the paired Φ is forced to its counting variant — the product is the tilt 1/12 = F_TRZ·Φ_5/6 = K_MEX − 2 (PAPER_2132 §3), the same exact rational that measures the Hubble tension (PAPER_1156) and the DPM-pair identity (PAPER_1183). Thirty-four observables across ten domains carry the cosmological tilt constant in their closed forms.

---

## 3. The exact sub-family and the census jewel

Members whose tilt term completes to an exact rational:

```
BCS isotope alpha    = Φ_5/6·(1 − F·D_phys/Φ·Φ) → Φ − F·Φ·D_phys = 1/2   EXACT
R_W = q_edge         = K_MEX − F·Φ_5/6          = 2                        EXACT (two domains)
Kleiber exponent     = Φ_5/6·(1 − F_TRZ)        = (5/6)(9/10) = 3/4        EXACT
WD exponent          = −D_phys·F·Φ_5/6          = −1/3                     EXACT
H0 tension ratio     = 1 + F·Φ_5/6              = 13/12                    EXACT
SM_wimp_exponent     = SO_5·D_phys + D_phys + K_MEX − F·Φ_5/6 = 46         EXACT-integer
```

And the jewel:

```
α⁻¹ = A_5·K_MEX + 1/(F_TRZ·Φ_5/6) = 125 + 12 = 137   EXACT-integer
```

The fine-structure integer decomposes as PAPER_1954's A_5·K_MEX = 125 plus the **inverse tilt** 12 — the electromagnetic coupling and the Hubble tension share the same 1/12 lattice constant, once as a reciprocal, once as a tilt.

---

## 4. Kernel-hierarchy context

The session-4 mining places the tilt factor at the base of a two-tier exact-rational hierarchy beneath PAPER_2132's Vacuum Coupling Kernel:

| Tier | Kernel | Exact value | ~Observables |
|:-:|---|---|:-:|
| 2 | VCK = F_TRZ·K_MEX·SSq | 19/160 | 5 |
| 1 | **F_TRZ·Φ_5/6 (tilt)** | **1/12** | **34** |
| 1 | F_TRZ·K_MEX | 5/24 | ~14 |
| 1 | F_TRZ·SSq | 57/1000 | ~6 |
| 1 | D_phys·F_TRZ | 2/5 | ~5 |
| 1½ | D_phys·F_TRZ·Φ_5/6 | 1/3 | WD/isotope family |

---

## 5. Honest residuals and disclosures

Family members carry their own published residuals; the exact sub-family entries are EXACT arithmetic on published forms. Two multi-term members (LCDM_D_over_H, SM_eta_gamma_gamma_BR) have unresolved inner-variant status — disclosed, pending. Discovery method: textual composition matching over published expressions; value-coincidence and name-token matching remain rejected (gate-pinned method rules). The Tilt-Product Law is falsifiable per §6. NOT REPLACEMENT.

---

## 6. Falsifiability

1. **Tilt-Product Law test:** any future published closure containing F_TRZ·Φ with the 0.84 variant breaks the law. Current record: 32-0.
2. **The 2 pending multi-term members** must resolve to 5/6 on inner-variant inspection.
3. **α⁻¹ decomposition test:** the 137 form requires the inverse tilt exactly; any revision of F_TRZ or Φ_5/6 breaks 125 + 12 and the 34-member family simultaneously (over-determination in the PAPER_2119/2128 signature).

---

## 7. Cross-references

PAPER_2132 (VCK, tilt factorization, triple-form 1/12 = F_TRZ·Φ_5/6 = K_MEX − 2); PAPER_2129 (sector rule — sharpened here); PAPER_1156 (Hubble tilt); PAPER_1183 (DPM-pair); PAPER_1160 (F_TRZ = 1/|SO(5)|); PAPER_1954 (A_5·K_MEX = 125); PAPER_1522 (K_MEX); sessions S388/S394/S396/S400/S403-S462 (family members); dispatch entries per §1.

---

## 8. Summary Statement

**PAPER_2133 files the largest two-primitive kernel census in the corpus: the tilt factor F_TRZ·Φ appears in 34 observables across ten domains, with all 32 resolvable instances selecting the counting variant — F_TRZ·Φ_5/6 = 1/12 EXACT, zero 0.84 instances — passing PAPER_2132's tilt-factor test at population scale and sharpening the PAPER_2129 sector rule into the Tilt-Product Law: the product F_TRZ·Φ is universally the counting composition, independent of host sector. The exact sub-family (1/2, 2, 3/4, −1/3, 13/12, 46) culminates in α⁻¹ = A_5·K_MEX + 1/(F_TRZ·Φ_5/6) = 125 + 12 = 137: the fine-structure integer and the Hubble tension carry the same 1/12 lattice constant.**

---

**Filed 2026-07-24. Append-only henceforth.**
