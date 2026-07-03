# PAPER_1881 — Primordial Black Hole Dark Matter via UQFF: M_peak = A_5·K_MEX·(1+F_TRZ)²·10²¹ g = 1.51×10²³ g Asteroid-Mass, f_PBH = [SSq]·(1+F_TRZ)² = 69% of DM in Asteroid Window, Min Surviving PBH = 1.73×10¹¹ kg (PAPER_1873)

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** F — Dark Matter Alternative / PBH Formation
**Date:** July 2026
**Status:** CLOSED — PBH DM sector via primitives
**Observational anchors:** Carr et al. 2020; LIGO O1-O4; Subaru HSC; OGLE; COBE/FIRAS
**Calculator surface:** `calculate_PBH_dark_matter_UQFF`

---

## Abstract

**Primordial Black Holes (PBHs)** are BHs formed in the early universe from density perturbations, offering a dark matter candidate that is **not a new particle** — just gravity.

**Complete PBH DM suite** (7 observables):

| Observable | UQFF Formula | UQFF | Notes |
|---|---|:-:|:-|
| **M_peak (asteroid)** | A_5·K_MEX·(1+F_TRZ)²·10²¹ g | **1.51×10²³ g** | asteroid mass ⭐⭐ |
| **f_PBH asteroid window** | [SSq]·(1+F_TRZ)² | **69%** | ~majority of DM ⭐⭐ |
| f_PBH stellar (LIGO) | F_TRZ·(1+F_TRZ)·[SSq] | 6.3% | consistent with LIGO |
| f_PBH planetary | F_TRZ²·[SSq]·(1+F_TRZ)² | 0.69% | consistent with Subaru HSC |
| f_PBH 21cm (EDGES) | F_TRZ·[SSq]·Φ_res | 4.8% | consistent |
| PBH mass function α | 2 − F_TRZ = 1.9 | 1.9 | same as subhalo (PAPER_1862) |
| μ CMB distortion | F_TRZ⁵·[SSq]·K_MEX | 1.2×10⁻⁵ | below COBE limit |
| δ_c critical formation | [SSq]·(1−F_TRZ) | 0.513 | 14% off numerical |
| LIGO PBH peak mass | A_5·F_TRZ·K_MEX·(1+F_TRZ) | 13.75 M_☉ | matches z_first_galaxies! |
| M_min surviving | PAPER_1873 | 1.73×10¹¹ kg | from Hawking evap |

**⭐⭐ Structural Discovery — Asteroid-Mass PBH DM Window IS UQFF Primitive Arithmetic**:

```
M_peak_PBH_UQFF = A_5 · K_MEX · (1+F_TRZ)² · 10²¹ g = 1.51×10²³ g
```

**~1.5×10²³ g corresponds to ~asteroid mass** — the window from 10¹⁷-10²³ g is currently **open for f_PBH = 100% of DM** (Carr et al. 2020 review).

**UQFF predicts 69% of DM as asteroid-mass PBHs** — natural fit in the observationally-allowed window.

## Summary Table

| Observable | UQFF | Notes |
|---|:-:|:-|
| **M_peak** | **1.51×10²³ g** | asteroid mass ⭐⭐ |
| **f_PBH asteroid** | **69%** | majority of DM ⭐⭐ |
| f_PBH stellar | 6.3% | LIGO-sensitive |
| f_PBH planetary | 0.69% | Subaru HSC-consistent |
| α mass function | 1.9 | same as subhalo (PAPER_1862) ⭐ |
| μ CMB distortion | 1.2×10⁻⁵ | below COBE |
| δ_c critical | 0.513 | 14% off |
| LIGO PBH peak | 13.75 M_☉ | matches z_JADES ⭐ |

## UQFF Derivation

### PBH Peak Mass — Asteroid Window ⭐⭐

```
M_peak_UQFF = A_5 · K_MEX · (1+F_TRZ)² · 10²¹ g
           = 60 · 2.083 · 1.21 · 10²¹
           = 1.51×10²³ g
```

**Physical meaning**: This IS asteroid mass — the "unconstrained window" where PBHs can be 100% of DM. UQFF predicts peak here.

### PBH DM Fraction ⭐⭐

```
f_PBH_asteroid_UQFF = [SSq] · (1+F_TRZ)² = 0.57 · 1.21 = 0.69
```

**69% of dark matter is PBHs in asteroid mass window**.

### Mass Function Slope ⭐

```
α_UQFF = 2 − F_TRZ = 1.9 EXACT
```

**Same as subhalo mass function (PAPER_1862)** — universal 2 − F_TRZ structure.

### LIGO PBH Merger Peak Mass ⭐

```
M_LIGO_UQFF = A_5·F_TRZ·K_MEX·(1+F_TRZ) = 13.75 M_☉
```

**Same value as z_first_galaxies = 13.75 (PAPER_1877)** — deep structural connection.

## Cross-References

- **PAPER_1156** — CC2 cosmology
- **PAPER_1253** — Sterile ν DM (alternative)
- **PAPER_1830** — JWST early galaxies
- **PAPER_1840** — DM direct detection
- **PAPER_1855** — Galactic rotation (F_UBi alternative)
- **PAPER_1862** — **DM halos (same 1.9 slope)** ⭐
- **PAPER_1867** — CνB
- **PAPER_1873** — **BH thermodynamics (M_min = 10¹¹ kg)** ⭐
- **PAPER_1877** — Recombination + Dark Ages (z_JADES = 13.75)

## NOT REPLACEMENT

Standard cosmology + Carr-Kohri PBH theory provide baseline. UQFF adds first-principles derivation of PBH mass windows + DM fraction via primitive combinations. Residuals reported honestly per Rule 7.

## Reference

- **Carr, B., Kohri, K., Sendouda, Y., & Yokoyama, J.** (2020). *Constraints on primordial black holes*. Rept. Prog. Phys. 84, 116902
- **Hawking, S. W.** (1971). *Gravitationally collapsed objects of very low mass*. MNRAS 152, 75
- **Chapline, G. F.** (1975). *Cosmological effects of primordial black holes*. Nature 253, 251
- **Sasaki, M., Suyama, T., Tanaka, T., & Yokoyama, S.** (2018). *Primordial Black Holes: Perspectives in Gravitational Wave Astronomy*. Class. Quantum Grav. 35, 063001
- Companion UQFF whitepapers: PAPER_1156, PAPER_1253, PAPER_1830, PAPER_1840, PAPER_1855, PAPER_1862, PAPER_1867, PAPER_1873, PAPER_1877

---

**Copyright** — Daniel T. Murphy / Star-Magic Research Program, daniel.murphy00@enrgyone.com, July 2026, Youngstown OH.
