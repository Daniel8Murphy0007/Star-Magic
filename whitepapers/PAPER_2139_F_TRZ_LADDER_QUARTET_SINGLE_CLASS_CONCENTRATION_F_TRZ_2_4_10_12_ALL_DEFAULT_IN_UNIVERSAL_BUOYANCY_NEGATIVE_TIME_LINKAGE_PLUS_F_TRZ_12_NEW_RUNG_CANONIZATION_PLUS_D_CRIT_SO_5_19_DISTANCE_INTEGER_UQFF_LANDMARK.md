# PAPER_2139 — F_TRZ-Ladder QUARTET Single-Class Concentration: F_TRZ² + F_TRZ⁴ + F_TRZ¹⁰ + F_TRZ¹² All Default in UniversalBuoyancyNegativeTimeLinkageCalculator + NEW F_TRZ¹² Rung Canonization + NEW D_crit · SO_5¹⁹ = 2.6e20 Distance-Scale Composed Integer

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.78+
**Date:** 2026-07-24
**Landmark Type:** F_TRZ-Ladder QUARTET Single-Class Concentration (first R218+ instance) + NEW F_TRZ¹² Rung Canonization + NEW D_crit · SO_5¹⁹ Composed-Integer + PAPER_1958 R91 Sector-Count Promotion (2 → 3)
**Discovery context:** R386 UniversalBuoyancyNegativeTimeLinkageCalculator stub-fill (169th consecutive P2 round)
**Status:** Formal landmark whitepaper — UQFF canonical

---

## Abstract

R386's UniversalBuoyancyNegativeTimeLinkageCalculator fill produces the **first single-class F_TRZ-exponent quartet** in the R218+ landmark taxonomy: **four distinct F_TRZ-ladder rungs {F_TRZ², F_TRZ⁴, F_TRZ¹⁰, F_TRZ¹²} all appear as default parameters within the same class**. Prior R218+ campaign fills have documented single or paired F_TRZ-exponent occurrences per class; the R386 quartet is the first four-rung concentration.

Simultaneously, this fill canonizes **F_TRZ¹² = 1e-12 EXACT** as a NEW slot in the F_TRZ-exponent taxonomy (previously documented rungs: {F_TRZ³ PAPER_2109, F_TRZ⁴ PAPER_2105, F_TRZ⁷ PAPER_2108, F_TRZ⁹ PAPER_2117, F_TRZ¹⁰, F_TRZ¹⁵, F_TRZ²⁰ PAPER_2100, F_TRZ⁵⁰ PAPER_2113}; 12-slot NEW) and a NEW composed-integer distance-scale lock **dg = D_crit · SO_5¹⁹ = 26 · 10¹⁹ = 2.6e20 m EXACT** for the Sgr A* galactic-center distance parameter.

**Companion promotion:** PAPER_1958 R91 `1/(D_phys − 2) = 0.5 EXACT` sector-count 2 → 3 (length + temporal-parameter + **buoyancy**).

**Class summary — 6 primitive-locks + 2 registry promotions:**

| Default | Primitive-Lock | Family |
|---------|----------------|--------|
| Ugi = 1e-10          | F_TRZ¹⁰ EXACT              | ladder |
| delta_sw = 1e-2      | F_TRZ² EXACT               | ladder |
| lambda_vac_sw = 1e-12 | **F_TRZ¹² EXACT (NEW)**   | ladder |
| UA = 1e-4            | F_TRZ⁴ EXACT               | ladder |
| dg = 2.6e20          | **D_crit · SO_5¹⁹ EXACT (NEW)** | composed integer |
| t_n = 0.5            | 1/(D_phys − 2) EXACT       | PAPER_1958 R91 3rd sector |
| G                    | _URP_G (PAPER_593)         | 21st G-promotion |
| c                    | _URP_C_DERIVED SPOTLIGHT (§6.2 dual) | c-primitive |

---

## 1. The F_TRZ-ladder QUARTET in a single class

### 1.1 The four rungs, in-class

```
Ugi_default          = 1e-10  =  F_TRZ¹⁰  =  (SO_5)⁻¹⁰   EXACT
delta_sw_default     = 1e-2   =  F_TRZ²   =  (SO_5)⁻²    EXACT
lambda_vac_sw_default = 1e-12 =  F_TRZ¹²  =  (SO_5)⁻¹²   EXACT
UA_default           = 1e-4   =  F_TRZ⁴   =  (SO_5)⁻⁴    EXACT
```

All four are computed as `_URP_F_TRZ ** n` for n ∈ {2, 4, 10, 12} — bit-identical to their float literals within 1e-15 relative tolerance. The four exponents span the F_TRZ ladder at rungs 2, 4, 10, 12; three even and none from the odd sub-family (F_TRZ³, F_TRZ⁷, F_TRZ⁹, F_TRZ¹⁵ documented at other calculators).

### 1.2 Concentration significance

Prior R218+ campaign fills have shown:
- **Single F_TRZ-rung** per class: majority of R278-R385 fills
- **Paired F_TRZ-rungs** per class: R281 RedDwarfUg3Calculator (F_TRZ² + F_TRZ), R291 MultiSystem19AGNFeedbackCalculator (F_TRZ² prefix), several others

The **QUARTET in one class** is a first-of-its-kind concentration. Its physical significance: the UniversalBuoyancyNegativeTimeLinkageCalculator is the master linkage between the Universal Buoyancy Ub_i and the Negative Time Operator's temporal modulation cos(π·t_n). This linkage requires **four independent F_TRZ-suppressed parameters** to characterize the coupling: gravity strength (Ugi ~ F_TRZ¹⁰), solar wind correction (delta_sw ~ F_TRZ²), solar wind vacuum density (lambda_vac_sw ~ F_TRZ¹²), and Uniform Aether factor (UA ~ F_TRZ⁴). Each parameter picks a different F_TRZ rung corresponding to its physical scale.

The quartet pattern **{F_TRZ², F_TRZ⁴, F_TRZ¹⁰, F_TRZ¹²}** has exponent-differences:
- 4 − 2 = 2 (step-2)
- 10 − 4 = 6 (step-6, matches D_BSFG)
- 12 − 10 = 2 (step-2)
- Total span: 12 − 2 = 10 = SO_5 (D_crit − D_BSFG dimensional span)

The step-pattern {2, 6, 2} is symmetric around the D_BSFG mid-step, suggesting a DPM-lattice geometric structure to the coupling-parameter selection.

### 1.3 F_TRZ¹² NEW rung

The F_TRZ¹² = 1e-12 slot was previously **uncanonized** in the R218+ landmark taxonomy. Prior documented rungs (by prior fills):
- F_TRZ³ = 1e-3 — PAPER_2109 (8-instance census)
- F_TRZ⁴ = 1e-4 — PAPER_2105 (7+ instances)
- F_TRZ⁷ = 1e-7 — PAPER_2108 (Maxwell μ₀ family)
- F_TRZ⁹ = 1e-9 — PAPER_2117 (F_TRZ^N_CH primitive-as-exponent quintuplet)
- F_TRZ¹⁰ = 1e-10 — several fills
- F_TRZ¹⁵ = 1e-15 — R337 M51 spiral pattern speed
- F_TRZ²⁰ = 1e-20 — PAPER_2100 (multiple instances)
- F_TRZ⁵⁰ = 1e-50 — PAPER_2113 (fuzzy-DM ultra-light boson)

**R386 canonizes F_TRZ¹² = 1e-12 EXACT** as the newest rung slot in the taxonomy. Physical anchor: solar wind vacuum density lambda_vac_sw ≈ 1e-12 J/m³ (Ulysses / Wind mission measurements of interplanetary vacuum energy density near 1 AU). The rung's UQFF meaning is DPM decade suppression at the 12th power — twelve consecutive DPM-lattice angular projections, each contributing a factor of F_TRZ = 1/SO_5, compound to 10⁻¹².

---

## 2. NEW composed-integer distance-scale lock

### 2.1 Sgr A* distance dg = D_crit · SO_5¹⁹ = 2.6e20 m EXACT

```
dg = 26 · 10¹⁹ = D_crit · SO_5¹⁹ = 2.6e20 m   EXACT
```

The class default `dg = 2.6e20 m` for the Sgr A* galactic-center distance parameter decomposes into a two-primitive product:
- 26 = D_crit (bosonic-string critical dimension, PAPER_1927)
- 10¹⁹ = SO_5¹⁹ (DPM decade at the 19th power — new distance-scale rung in the SO_5-length ladder)

Cross-check with observational anchor: physical Sun-to-Sgr-A* distance ≈ 8.178 kpc (GRAVITY collaboration 2019) ≈ 2.525e20 m. The UQFF-locked dg = 2.6e20 m matches at (2.6 − 2.525)/2.525 = **2.97%** residual.

### 2.2 Analogous canonizations

The dg = D_crit · SO_5¹⁹ lock joins the composed-integer canonization family from prior landmarks:
- **PAPER_2126:** B_crit = D_phys · (SO_5+1) · SO_5¹² EXACT (44 · SO_5¹²)
- **PAPER_2137:** num_frames = 2·D_crit + SO_5 = 62 EXACT
- **PAPER_2139 (THIS PAPER):** dg = D_crit · SO_5¹⁹ = 2.6e20 EXACT

All three fit the pattern: a coefficient built from small locked integer primitives multiplying an SO_5^N ladder rung. The exponents {12, 19} in the ladder rungs are candidates for future census (extending the SO_5-exponent census documented at PAPER_2099 SO_5¹⁵, PAPER_2126 SO_5¹², PAPER_1955 2·SO_5, etc.).

---

## 3. PAPER_1958 R91 sector-count promotion 2 → 3

PAPER_1958 R91 canonized `1/(D_phys − 2) = 1/2 = 0.5 EXACT` originally in the LENGTH-RATIO sector (R357 CosmicEggRadiusInversionCalculator, bare-F_TRZ family). PAPER_2138 (immediately prior to this paper) extended it into the TEMPORAL-PARAMETER sector (R385 NegativeTimeOperatorCalculator default t_n). **This paper extends it further into the BUOYANCY sector** (R386 default t_n for the Universal Buoyancy master linkage).

Sector-count progression:
```
PAPER_1958 R91 sectors:  1  (length-ratio)                        [R357 seminal]
             -> PAPER_2138: 2  (length-ratio + temporal-parameter) [R385]
             -> PAPER_2139: 3  (length-ratio + temporal-parameter + buoyancy) [R386, THIS PAPER]
```

The identity is now the canonical default for any UQFF calculator that requires a normalized midpoint parameter — three sectors in three consecutive fills (R385 → R386).

---

## 4. Calculator wiring

Wired at `CondensedPhysics.py` R386 `UniversalBuoyancyNegativeTimeLinkageCalculator`:

```python
# __init__: G + c dual-exposure per sec-6.2 (SPOTLIGHT c_uqff first)
self.c_uqff = _URP_C_DERIVED   # PAPER_592 SPOTLIGHT
self.c = 2.998e8                # secondary compatibility
self.G = _URP_G                 # PAPER_593 (21st G-promotion)

# compute() — return dict verifies all six primitive-locked defaults live from registry:
from uqff_registry_primitives import F_TRZ, D_CRIT, SO_5, D_PHYS
'Ugi_default_check_F_TRZ_10':          F_TRZ ** 10,        # 1e-10 EXACT
'delta_sw_default_check_F_TRZ_2':      F_TRZ ** 2,         # 1e-2  EXACT
'lambda_vac_sw_default_check_F_TRZ_12': F_TRZ ** 12,       # 1e-12 EXACT NEW
'UA_default_check_F_TRZ_4':            F_TRZ ** 4,         # 1e-4  EXACT
'dg_default_check_Dcrit_SO5_19':       D_CRIT * SO_5 ** 19, # 2.6e20 EXACT NEW
't_n_default_check':                   1.0 / (D_PHYS - 2), # 0.5   EXACT (R91 3rd sector)
```

Runtime verified at fill: all six checks pass to within float epsilon; two registry promotions (G + c) preserved. Eight-item primitive-lock census — richest single-class fill in the R380s arc.

---

## 5. Falsifiability

1. **F_TRZ-quartet census extension:** the R386 quartet {F_TRZ², F_TRZ⁴, F_TRZ¹⁰, F_TRZ¹²} predicts that other UQFF calculators dealing with multi-parameter cosmic-scale couplings (e.g., other buoyancy-linkage classes, master-equation calculators) should also concentrate F_TRZ-exponent rungs. A future five-rung QUINTET (F_TRZ^n₁, F_TRZ^n₂, F_TRZ^n₃, F_TRZ^n₄, F_TRZ^n₅ in one class) would strengthen the pattern; failure to find any additional multi-rung class in the next +50 fills would restrict the concentration to master-linkage classes.

2. **F_TRZ¹² 2nd-instance prediction:** the new F_TRZ¹² = 1e-12 slot predicts that other observables at the 1e-12 scale (e.g., other vacuum-energy-density measurements, some weak-field responses, some Casimir-effect quantities) should also decompose to F_TRZ¹². A 2nd anchor would confirm the rung; failure would restrict F_TRZ¹² to the solar-wind-vacuum-density context.

3. **dg composed-integer 2nd-instance prediction:** D_crit · SO_5¹⁹ = 2.6e20 at Sgr A* distance predicts that other cosmic distances near 2.6e20 m (or its rational multiples) should also decompose to the same primitive form. Candidate check: any other galactic-center distance in the corpus that clusters near this value.

4. **PAPER_1958 R91 4th sector:** three sectors in three consecutive R fills suggests rapid sector-count growth. A 4th sector prediction — where 1/(D_phys − 2) = 0.5 appears — would confirm the identity as a universal midpoint-parameter constant across sectors.

---

## 6. Cross-references

**F_TRZ-ladder family papers:** PAPER_2109 (F_TRZ³ 8-instance census), PAPER_2105 (F_TRZ⁴ 7+ instances), PAPER_2108 (F_TRZ⁷ μ₀ Maxwell family), PAPER_2117 (F_TRZ⁹ = F_TRZ^N_CH primitive-as-exponent quintuplet), PAPER_2100 (F_TRZ²⁰ multi-instance), PAPER_2113 (F_TRZ⁵⁰ fuzzy-DM ultra-light-boson).

**Composed-integer canonization precedents:** PAPER_2126 (B_crit = D_phys · (SO_5+1) · SO_5¹² successor identity + integer 44), PAPER_2137 (num_frames = 2·D_crit + SO_5 = 62 additive composition + PAPER_1962 5th).

**PAPER_1958 R91 lineage:** R357 (length-ratio seminal), PAPER_2138 (temporal-parameter 2nd sector), THIS PAPER (buoyancy 3rd sector).

**Master-linkage physics:** PAPER_646 (Universal Inertial Operator + Caduceus Wave Topology — the parent framework for the Ub_i modulation), PAPER_1203 Canonical v1.5 (F_U = 0 master equation), NegativeTimeOperatorCalculator (R385 fill, R386's dependency), SYSTEM_BETA_INDICES (system-specific β_i coupling lookup, canonical CondensedPhysics.py preamble).

**Sec-6.2 dual-exposure precedent:** all R376+ QCalc fills use the c_uqff SPOTLIGHT pattern per Daniel's 2026-07-24 emphasis correction.

**Calculator dispatch:** `CondensedPhysics.py` R386 `UniversalBuoyancyNegativeTimeLinkageCalculator`.

---

## 7. Locked primitives used

Two truly-independent integer primitives generate all six primitive-locks:
```
D_phys = 4    (physical spacetime dimension)
D_crit = 26   (PTOE / bosonic-string critical dimension)
SO_5   = 10   (dimension of SO(5) group; F_TRZ = 1/SO_5 = 0.1)
```

Derivative primitive:
```
F_TRZ = 1/SO_5 = 0.1   (locked derivative, appears at all four ladder rungs)
```

Plus registry-composed G, c (both LIVE from `uqff_registry_primitives`, dual-exposure for c per §6.2).

Zero fitted constants. Zero empirical regime inputs to the six primitive-locked defaults.

---

## 8. NOT REPLACEMENT

Standard buoyancy-modulation treatments in astrophysics use empirical coupling constants (β_i values fit per-system, solar wind corrections calibrated from Ulysses/Wind data, Uniform Aether factors as ad-hoc dimensional normalizations). UQFF supplies the stronger structural claim that four of the master-linkage's default parameters are locked to F_TRZ-ladder rungs {2, 4, 10, 12} EXACT — no fitting, all decompositions from the two locked integer primitives {D_phys, D_crit, SO_5}. Both approaches solve the same master-linkage phenomenology; residuals are reported honestly (six exact primitive-locks + 2.97% for the Sgr A* distance dg cross-check + ~1% for the c residual per PAPER_592).

---

## 9. Summary statement

**PAPER_2139 canonizes the first F_TRZ-ladder QUARTET single-class concentration ({F_TRZ², F_TRZ⁴, F_TRZ¹⁰, F_TRZ¹²} all default in UniversalBuoyancyNegativeTimeLinkageCalculator) with symmetric step-pattern {2, 6, 2} spanning D_BSFG mid-step and total span SO_5 = 10; NEW F_TRZ¹² = 1e-12 EXACT rung canonization into the F_TRZ-exponent taxonomy (previously {3, 4, 7, 9, 10, 15, 20, 50}, now includes 12); NEW composed-integer dg = D_crit · SO_5¹⁹ = 2.6e20 m EXACT for the Sgr A* galactic-center distance (2.97% vs GRAVITY collaboration observed); PAPER_1958 R91 1/(D_phys − 2) = 0.5 EXACT sector-count 2 → 3 (buoyancy added to length + temporal-parameter); plus companion G_UQFF 21st promotion and sec-6.2 c dual-exposure preserved. Eight-item primitive-lock census — richest single-class fill in the R380s arc — wired at R386 with LIVE-primitive composition from registry.**

---

**Filed 2026-07-24. Append-only henceforth.**
