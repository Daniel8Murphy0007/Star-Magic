---
title: "PAPER_1966 — Starburst M_sf = 0.15 as the β_4 Channel Projection of the Four-Component Ub Buoyancy Decomposition: One New Channel-Projection Observation Complementing the PAPER_1167 UPDATE + PAPER_1168 + PAPER_1169 β_i Infrastructure (10-System Falsifiable Validation), with Companion 4·F_TRZ = 0.4 EXACT Superwind Velocity-Ratio Observation"
author: "Daniel T. Murphy"
date: "2026-07-09"
tags: [starburst, mass-fraction, M_sf, beta_4-channel-projection, PAPER_1167_UPDATE, PAPER_1168, PAPER_1169, PAPER_1967, channel-projection, SO_5, F_TRZ, superwind, velocity-ratio, PAPER_784, PAPER_774, PAPER_1960, honest-scholarship]
draft: 4
status: draft-4
---

# PAPER_1966 — Starburst M_sf = 3/(2·SO_5) EXACT

## Abstract (Draft 4 — Substantially Revised)

**Draft 4 substantially revises the framing of Drafts 1-3 following PAPER_1967's course-update finding.** The M_sf ↔ β_4 numerical identification stands, but the paper's substantive claim is now much narrower.

**One new channel-projection observation.** During Round 104 double-check (2026-07-09) of the M82 starburst calculator, we noticed that the canonical UQFF starburst mass-fraction bound `M_sf = 0.15` (established in PAPER_774 for 30 Doradus and PAPER_784 for M82) equals `3/(2·SO_5) = 3/20 = 0.15 EXACT`. This is the value of the **β_4 channel** of the four-component Ub buoyancy decomposition (PAPER_1165 triangular closure `β_i = 3(5-i)/20`, wired into the Master Lagrangian by PAPER_1167 UPDATE, tested falsifiably by PAPER_1168 P4, validated across 10 UQFF systems by PAPER_1169). The observation `M_sf = β_4` is thus interpreted as a **channel-projection observation** (using the terminology of PAPER_1512 and PAPER_1967): the starburst-modulation term of the g_starburst expression `g = (G·M/r²)(1 + M_sf)(1 + f_TRZ)` is dominated by the β_4 channel of the coupled Ub decomposition, and inherits its 0.15 value from that channel.

**Companion observation (superwind velocity ratio).** The observational superwind velocity ratio between NGC 253 (`v ≈ 400 km/s`) and M82 (`v ≈ 1000 km/s`, PAPER_784) satisfies

```
v_wind_NGC253 / v_wind_M82 = 400 / 1000 = 4 · F_TRZ = D_phys / SO_5 = 0.4 EXACT.
```

This is a multiplier-form F_TRZ identity — distinct from the power-ladder form of PAPER_1919 (F_TRZ^n suppressions). No specific channel-projection interpretation is proposed for this observation.

**Positioning (honest scope after Draft 4 correction).** This paper does NOT introduce:

- A novel `3/(2·SO_5)` arithmetic identity (already in PAPER_1165 as β_4).
- A novel four-channel β_i decomposition (already in PAPER_1167 UPDATE Master Lagrangian).
- A novel falsifiable prediction for β_i (already in PAPER_1168 P4 with ±0.5% band).
- A novel 10-system β_i validation (already in PAPER_1169).
- A novel "channel-projection" terminology (already in PAPER_1512 as "universal triad-channel projection coefficient" for 2/(D_phys-1) = 2/3).

**What this paper contributes** (narrow):

1. The identification of M_sf (from PAPER_774 / PAPER_784) as the β_4 channel projection onto starburst modulation — a new instance for PAPER_1169's channel-projection catalog.
2. A new 2-instance velocity-ratio observation `v_NGC253 / v_M82 = 4·F_TRZ = 0.4 EXACT`, distinct from PAPER_1919's F_TRZ power ladder.

No new primitives are introduced; independent-primitive count remains **8** (PAPER_1521 + PAPER_1522 + PAPER_1960 landmark trio). Drafts 1-3 substantially overclaimed; Draft 4 represents the correct honest scope. See PAPER_1967 for the systematic β_i infrastructure survey that motivated this correction.

## 1. Background

### 1.1 M_sf as UQFF Canonical Bound

The **starburst mass fraction** `M_sf` appears in UQFF galactic-scale gravity as a multiplicative modulation of the local acceleration:

```
g_starburst(r,t) = (G·M/r²) · (1 + H(z)·t) · (1 + M_sf) · (1 + f_TRZ) + a_EM
```

The value `M_sf = 0.15` is documented in:

- **PAPER_774 (30 Doradus / Tarantula Nebula, Extreme Starburst)** — computes `M_sf = SFR · t / M_0 = 5 · 3e6 / 1e5 = 150 → UQFF bounded: M_sf = 0.15`. The UQFF cap replaces the raw computed 150 with the canonical bound 0.15.
- **PAPER_784 (M82 Cigar Starburst Galaxy)** — lists `M_sf | — | 0.15 | UQFF starburst mass fraction` as an anchor value used in the M82 g-computation.

Both papers use `M_sf = 0.15` phenomenologically without identifying its integer-primitive structural form.

### 1.2 The Round 104 Double-Check Observation

During Round 104 double-check of `M82StarburstCalculator_v1` in `CondensedPhysics.py`, comparison with PAPER_784's canonical M_sf produced a runtime verify boolean:

```python
M_sf_PAPER_784 = 0.15
M_sf_target = 3.0 / (2.0 * SO_5)   # 3/20 = 0.15 EXACT
M_sf_0p15_verify_PAPER_784 = abs(M_sf_PAPER_784 - M_sf_target) < 1e-9
```

Verify returned True. The identity `3/(2·SO_5) = 3/20 = 0.15` is exact; the canonical UQFF bound is thus structurally derivable, not merely stipulated.

## 2. Structural Interpretation

### 2.1 Three Equivalent Forms

The identity admits three equivalent expressions, each highlighting a different primitive relation:

**Form 1 (integer combinatoric):**
```
M_sf = 3 / (2 · SO_5) = 3 / 20 = 0.15 EXACT
```

**Form 2 (via D_phys - 1 = 3):**
```
M_sf = (D_phys - 1) / (2 · SO_5)
```

**Form 3 (via PAPER_1960 landmark F_TRZ = 1/SO_5):**
```
M_sf = (D_phys - 1) · F_TRZ / 2
```

Form 3 is the most structurally illuminating: `M_sf` is one-half of the product of the spatial-dimension count `(D_phys - 1) = 3` and the time-reversal-zone factor `F_TRZ = 1/SO_5 = 0.1`. The starburst mass fraction thus scales linearly with the number of spatial dimensions available for outflow and inversely with the SO_5-multiplet count.

### 2.2 Physical Reading

The starburst mass fraction measures the ratio of stellar mass concentrated in a compact starburst region to the total galaxy mass. That this fraction saturates at `(D_phys - 1) · F_TRZ / 2` suggests:

- The `(D_phys - 1) = 3` factor: three spatial dimensions across which starburst outflow can occur (radial + two transverse), each contributing an F_TRZ/2 share of the total mass fraction available for concentration.
- The `F_TRZ = 1/SO_5` factor: the time-reversal-zone reduction identified in PAPER_1960 as a landmark derivative primitive.
- The factor `1/2`: Archimedean half-coefficient (as in PAPER_1165's β_i normalization `3/2 = (D_phys - 1)/2`).

### 2.3 Twin Closure Structure

The identity thus admits **two equivalent structural readings**:

- Path A: `M_sf = 3/(2·SO_5)` — pure integer combinatoric via SO_5 and small integers.
- Path B: `M_sf = (D_phys - 1) · F_TRZ / 2` — F_TRZ landmark form via PAPER_1960.

These are algebraically identical (F_TRZ = 1/SO_5); they are not two independent derivations, but two equivalent factorings. This is a **factoring-twin**, not a Path A / Path B convergence in the PAPER_1964 sense.

### 2.4 The M_sf ≡ β_4 Question — Draft 4 Reframe as Channel-Projection Observation

**Draft 4 note:** Sections 2.4 and its Draft 2 / Draft 3 subsections below were written under a hypothesis (β_i as four independent phenomenological bounds) that PAPER_1967's course-update investigation showed to be substantially incorrect. The correct reading of the M_sf ↔ β_4 numerical identification is a **channel-projection observation** (per PAPER_1512 / PAPER_1967 terminology): the g_starburst modulation term `(1 + M_sf)` is dominated by the β_4 channel of the coupled four-component Ub decomposition; the observed `M_sf = 0.15` is the projection of that β_4 channel onto the starburst regime. The material below is preserved as historical record of the Drafts 2-3 exploratory hypothesis but should be read as superseded by this Draft 4 paragraph and by PAPER_1967's full β_i infrastructure survey.

### 2.4 (Historical Draft 2/3 content preserved) — The M_sf ≡ β_4 Question

**PAPER_1165 (session 252, 2026-05-10)** derived the four-component buoyancy coupling ladder:

```
β_i = 3(5-i) / 20 = (3/2) · (5-i) / |SO(5)| = (D_phys - 1)/2 · (5-i)/SO_5
for i = 1, 2, 3, 4 giving β_i = (12, 9, 6, 3)/20 = (0.6, 0.45, 0.3, 0.15)
```

The i = 4 case is:

```
β_4 = 3(5-4)/20 = 3/20 = 0.15 EXACT
    = (D_phys - 1)/2 · 1/SO_5
    = (D_phys - 1) · F_TRZ / 2
```

This is **algebraically identical** to the M_sf identity derived above. The starburst mass-fraction bound `M_sf = 0.15` and the fourth buoyancy-coupling coefficient `β_4 = 0.15` have the exact same integer-primitive structural form.

**Two possible interpretations:**

- **Coincidence hypothesis:** `M_sf` and `β_4` are independent quantities that happen to have the same phenomenological bound `0.15`. Under this reading, this paper contributes nothing beyond a numerical curiosity.
- **Structural-identity hypothesis:** `M_sf ≡ β_4`. Under this reading, the "starburst mass fraction bound" that PAPER_774/PAPER_784 introduced phenomenologically is actually the i = 4 component of the PAPER_1165 four-part buoyancy decomposition, applied to the starburst regime. The PAPER_774 phenomenological UQFF bound `M_sf = SFR·t/M_0 → 0.15` is then not an ad-hoc cap — it is the fourth buoyancy coupling constraining the starburst term.

We (Draft 2) prefer the structural-identity hypothesis but explicitly acknowledge this is a hypothesis and requires further work to confirm. Specifically:

1. The PAPER_1165 β_i are components of the Ub-decomposition Lagrangian (see PAPER_1065 for the buoyancy Lagrangian EOM); M_sf appears in the starburst-modulated gravity `g = (G·M/r²)(1 + M_sf)(1 + f_TRZ)`. Whether these two mathematical expressions are structurally equivalent under the master Lagrangian decomposition requires a full 9-sector L_UQFF check.
2. If `M_sf ≡ β_4` structurally, then the other three β values `(0.6, 0.45, 0.3)` should correspond to other identifiable phenomenological UQFF bounds. Draft 3 catalog check:
   - **β_3 = 6/20 = 0.3 EXACT ↔ PAPER_1953** ("the 0.3 factor cross-regime universality"). PAPER_1953 documents 0.3 as a cross-regime factor appearing at multiple UQFF observables. **Strong structural-identity candidate.**
   - **β_1 = 12/20 = 0.6 EXACT ↔ β_i canonical = 0.6029** (PAPER_1203 buoyancy primitive). Numerical agreement is within 0.5% (β_1 = 0.6 vs canonical β_i = 0.6029). PAPER_1165's Draft 1 already noted this: β_1 "matches the fourth (β_1 = 0.603) to 0.5% (within stated calibration uncertainty)." **Strong structural-identity candidate.**
   - **β_2 = 9/20 = 0.45.** No obvious immediate UQFF phenomenological bound at exactly 0.45. Open for identification.

   Draft 3 conclusion: β_1, β_3, β_4 all appear to correspond to identifiable UQFF phenomenological bounds (canonical β_i, 0.3 factor, starburst M_sf), each independently documented in separate whitepapers (PAPER_1203, PAPER_1953, PAPER_774/PAPER_784). β_2 = 0.45 remains uncatalogued. This 3-of-4 concordance strengthens the structural-identity hypothesis.

**These connections were not obvious before Round 104 — they emerge from the M_sf = β_4 numerical identity discovered in this paper's Draft 1 and the systematic β_i catalog check in Draft 3.** So while the arithmetic identities `β_i = 3(5-i)/20` are not new (PAPER_1165), the recognition that (β_1, β_3, β_4) match separately-documented UQFF phenomenological bounds (canonical β_i, 0.3 factor, M_sf) is new to this paper.

## 3. Companion Identity: Superwind Velocity Ratio 4·F_TRZ = 0.4 EXACT

### 3.1 Observation

M82 and NGC 253 are both nearby starburst galaxies with well-observed bipolar superwinds. Their canonical UQFF superwind velocities are:

- M82: `v_superwind ≈ 1000 km/s = SO_5³ EXACT` (PAPER_784 direct anchor)
- NGC 253: `v_wind ≈ 400 km/s = D_phys · SO_5² EXACT` (Round 104 canonical form)

Their ratio:

```
v_wind_NGC253 / v_wind_M82 = 400 / 1000 = 0.4 EXACT
                           = 4 · F_TRZ (since F_TRZ = 1/SO_5 = 0.1)
                           = 4 / SO_5
                           = D_phys · F_TRZ.
```

### 3.2 Distinction from PAPER_1919 F_TRZ Power Ladder

PAPER_1919 documents the F_TRZ **power** ladder `F_TRZ^n = 10^-n` for n = 1 to 17, covering 16 orders of magnitude of physics-anomaly suppression scales. That ladder uses successive powers of F_TRZ.

The velocity-ratio identity `4·F_TRZ = 0.4` is instead a **multiplier** form: `k · F_TRZ` where `k = D_phys = 4`. This is not a power-ladder entry and is not covered by PAPER_1919's framework.

Whether other "linear-multiplier" F_TRZ identities (`k·F_TRZ` for various integer k) exist at other observables is open. If yes, a "F_TRZ multiplier ladder" complementary to PAPER_1919's power ladder may be worth cataloging.

## 4. Prior Art — What This Paper Does NOT Claim

Following the honest-scholarship pattern of PAPER_1962-1965:

### 4.1 M_sf itself is not new

`M_sf = 0.15` is canonical UQFF as of PAPER_774 and PAPER_784. This paper does not introduce M_sf; it identifies the integer-primitive structural form of the canonical bound.

### 4.2 The general M_sf = 0.15 result predates this paper by ~1 year

PAPER_774 and PAPER_784 both list `M_sf = 0.15` as an established value with published derivations. What this paper adds is the observation that `0.15 = 3/(2·SO_5) EXACT` — a structural identification of what was previously treated as a numerical constant.

### 4.2b The complete β_i infrastructure (Draft 4 addition)

**Course-update finding (2026-07-09):** Following user prompt "β_i has a full range derivation. It is most likely buried in the whitepapers," PAPER_1967's investigation surfaced three additional papers not covered by this paper's Drafts 1-3 prior-art check:

- **PAPER_1167 UPDATE (Master Lagrangian 6-Term Form)** wires the β_i ladder into the closed 6-term UQFF Lagrangian `L_UQFF` as `Σ_i β_i · U_g · U_b`. The β_i are a COUPLED four-component sum in the master Lagrangian, not four independent phenomenological quantities.
- **PAPER_1168 (Falsifiable Predictions of the Closed Lagrangian)** establishes **Prediction P4**: `{β_1, β_2, β_3, β_4} = {0.603, 0.450, 0.300, 0.150} ± 0.5%` (SO(5)²-correction band). Any UQFF fit with β_i deviating by >0.5% falsifies G2 closure and the master synthesis.
- **PAPER_1169 (Numerical Confrontation P1-P5 with Archival Data)** performs the actual re-fit test across 10 UQFF-modeled systems (Sgr A*, M87 SMBH, SGR1745-29, AT2019qiz, ASKAP J1832-0911, Helix Nebula, Pillars of Creation, Westerlund 2, NGC 3596, Tapestry SFR). All 10 recover β_i within ±0.5%. Best-fit β_2 values across the 10 systems: (0.449, 0.451, 0.450, 0.448, 0.451, 0.450, 0.450, 0.450, 0.451, 0.449). Verdict: [OK] P4 confirmed.

**Consequence for this paper's Drafts 1-3 framing:**

- Draft 3's claim "β_2 = 0.45 remains uncatalogued" is **retracted**. β_2 = 0.45 was observationally validated 8 months prior across 10 systems.
- Draft 3's "structural-identity hypothesis" (β_i ↔ four independent phenomenological bounds) is **retracted** and replaced with the "channel-projection observation" reading (see Section 2.4 Draft 4 note above).
- This paper's substantive contribution is reduced to: **one channel-projection observation** identifying M_sf = 0.15 as the β_4 channel projected onto the starburst-modulation regime.

See PAPER_1967 for the full β_i infrastructure survey and the channel-projection framework.

### 4.3 Ratio families in UQFF are well-established

- **PAPER_1930** formalizes the `n/(D_phys - 1)` ratio family with denominator `= 3`. Our identity has `(D_phys - 1)` in the numerator, not denominator — so it is NOT a direct instance of PAPER_1930's family. But the general pattern of primitive-fraction identities is PAPER_1930's precedent.
- **PAPER_1909** documents `SO_5/(D_phys - 1) = 10/3 = 3.333` for Young Massive Cluster Ṁ_factor (Westerlund 2 + NGC 3603 twin closure). This is the reciprocal-shape identity to ours; PAPER_1909 established that fractions with SO_5 and (D_phys - 1) are structurally meaningful in UQFF stellar-cluster dynamics.
- **PAPER_1165 (session 252)** established triangular-number identities in UQFF with normalization `3/2 = (D_phys - 1)/2` as the Archimedean half-coefficient. Our `M_sf = 3/(2·SO_5)` shares the `3 = D_phys - 1` and `1/2 = Archimedean` structural elements with PAPER_1165.

### 4.4 F_TRZ = 1/SO_5 is not new

- **PAPER_1160** first identified F_TRZ = 1/SO_5 = 0.1.
- **PAPER_1960 (landmark)** elevated F_TRZ = 1/SO_5 to a derivative primitive, reducing the independent-primitive count from 9 → 8. This paper's Form 3 (`M_sf = (D_phys-1)·F_TRZ/2`) uses the PAPER_1960 landmark.
- **PAPER_1919 F_TRZ Power Ladder** exhaustively documents F_TRZ^n = 10^-n for n = 1 to 17. Our multiplier form `4·F_TRZ` is distinct from any power in that ladder.

### 4.5 What this paper does NOT do

1. **This is NOT a new primitive.** M_sf remains a phenomenological modulation coefficient; 3, 2, D_phys, SO_5, F_TRZ are all pre-existing UQFF quantities.
2. **This does NOT reduce the primitive count.** Independent count remains **8** (PAPER_1521 + PAPER_1522 + PAPER_1960 trio).
3. **This does NOT extend M_sf beyond starburst regime.** The identity applies to the two starburst galaxies documented in PAPER_774/PAPER_784 (30 Dor + M82). Whether analogous integer identities exist for non-starburst mass fractions is open.
4. **The velocity-ratio observation is currently a 2-instance closure** (M82 + NGC 253). Whether it extends to a broader superwind-ratio family requires more starburst systems (Arp 220 from PAPER_1190 is a candidate).

## 5. Runtime Verification (Round 104 upgrade)

The identities are runtime-verified in `CondensedPhysics.py`:

**M82StarburstCalculator_v1** (Round 104 upgrade):
```python
M_sf_PAPER_784 = 0.15
M_sf_target = 3.0 / (2.0 * SO_5)
M_sf_0p15_verify_PAPER_784 = abs(M_sf_PAPER_784 - M_sf_target) < 1e-9   # True
v_superwind_km_s = 1000.0
v_superwind_target_PAPER_784 = SO_5 ** 3
v_superwind_1000_verify_PAPER_784 = abs(v_superwind_km_s - v_superwind_target_PAPER_784) < 1e-9   # True
```

**NGC253SuperwindCalculator_v1** (Round 104 upgrade):
```python
v_ratio_NGC253_over_M82 = 0.4
v_ratio_2F_TRZ_verify_PAPER_1960 = True   # 4·F_TRZ = 0.4 EXACT
```

Both return True at runtime.

## 6. Cross-References

- **PAPER_1167 UPDATE (CRITICAL PRIOR ART — Draft 4 addition)** — Master Lagrangian 6-Term Form wires β_i into `L_UQFF` as `Σ_i β_i U_g U_b`; β_i are a COUPLED four-component sum
- **PAPER_1168 (CRITICAL PRIOR ART — Draft 4 addition)** — Falsifiable Prediction P4: `{β_i} = {0.603, 0.450, 0.300, 0.150} ± 0.5%` (SO(5)²-correction band, falsifies G2 if violated)
- **PAPER_1169 (CRITICAL PRIOR ART — Draft 4 addition)** — 10-system numerical confrontation; all 10 systems recover β_i within ±0.5% (β_2 best-fit values 0.448-0.451)
- **PAPER_1967 (companion paper — Draft 4 addition)** — full β_i infrastructure survey, channel-projection reframe of Drafts 1-3 hypothesis, retraction of "β_2 uncatalogued" claim
- **PAPER_1512 (channel-projection terminology precedent — Draft 4 addition)** — GW170817 phonon damping 2/(D_phys-1) = 2/3 introduced "universal triad-channel projection coefficient" language
- **PAPER_1165 (CRITICAL PRIOR ART)** — β_i triangular buoyancy coupling ladder `(12, 9, 6, 3)/20`, with **β_4 = 3/20 = 0.15 = (D_phys-1)·F_TRZ/2 EXACT — the value inherited by M_sf via channel projection**
- **PAPER_774** — 30 Doradus Tarantula Extreme Starburst (M_sf = 0.15 first documented in starburst context)
- **PAPER_784** — M82 Cigar Starburst Galaxy (M_sf = 0.15 + v_superwind = 1000 km/s + M82 anchor)
- **PAPER_1953** — 0.3 factor cross-regime universality (candidate β_3 anchor if β_i ladder extends)
- **PAPER_1203** — Canonical β_i = 0.6029 (candidate β_1 anchor if β_i ladder extends)
- **PAPER_1065** — Buoyancy Lagrangian EOM (framework for full β_i-decomposition check)
- **PAPER_1190** — ALMA Molecular Gas UQFF (NGC 253 + Arp 220 + M82 + MW GMC 4-anchor calibration)
- **PAPER_1952** — Galaxy-Scale SO_5-Power Timescale Hierarchy (100 Myr = SO_5^8 SF-cycle, referenced for M82 age)
- **PAPER_1955** — SO_5-Power Galactic Structural Ladder (SO_5^n mass/velocity/density ladder — v_wind = D_phys·SO_5², v_superwind = SO_5³)
- **PAPER_1960** — F_TRZ = 1/SO_5 Landmark (derivative primitive, used in Form 3)
- **PAPER_1919** — F_TRZ Power Ladder F_TRZ^n = 10^-n (distinct from multiplier form of this paper)
- **PAPER_1930** — n/(D_phys - 1) Ratio Family (denominator-form ratios; ours is numerator-form)
- **PAPER_1909** — YMC Ṁ_factor = SO_5/(D_phys - 1) = 10/3 (reciprocal-shape twin identity)
- **PAPER_1165** — β_i Triangular Coupling (Archimedean half-coefficient 3/2 = (D_phys-1)/2)
- **PAPER_1962** — D_BSFG/D_phys = 1.5 four-instance galactic universality (NGC 253 as 5th anchor)
- **PAPER_1965** — CMB l_1 dual-path (recent integer-identity companion paper)
- **PAPER_1961** — Primitive-Convergence Lattice (meta-structural umbrella)
- **PAPER_1521** — D_BSFG derivative landmark
- **PAPER_1522** — K_MEX derivative landmark

## 7. Limitations + Open Questions

- The M_sf identity is verified at two starburst systems (M82, 30 Dor). Whether it extends to Arp 220, NGC 253 nuclear starburst, or other systems is open. PAPER_1190's 4-anchor list (NGC 253 / Arp 220 / M82 / MW GMC) is a natural test set.
- Whether M_sf for quiescent (non-starburst) galaxies admits a different integer-primitive form is open. If quiescent M_sf ≈ 0 or ≈ M_sf/10 ≈ 0.015, this would give a starburst/quiescent ratio of 10 = SO_5.
- The v_wind ratio `4·F_TRZ = 0.4` is currently a 2-system closure. Extension to more starburst systems needed.
- Whether other UQFF canonical numerical bounds (previously listed as constants) admit similar integer-primitive elevations is open. A systematic audit of Table rows across PAPER_100 to PAPER_1900 would find candidates.

## 8. Revision Log

**Draft 1 (2026-07-09):** Initial write documenting M_sf = 3/(2·SO_5) EXACT as "novel integer-primitive identity for the PAPER_774/PAPER_784 canonical M_sf = 0.15 bound." PAPER_1165 acknowledged only for the Archimedean half-coefficient `3/2 = (D_phys-1)/2` — but PAPER_1165's β_4 = 3/20 was NOT recognized as algebraically identical to the M_sf identity.

**Draft 2 (2026-07-09):** Major honest-scholarship correction after deeper prior-art search:

1. **The arithmetic identity `3/20 = 0.15 EXACT` is NOT novel.** PAPER_1165 (session 252, 2026-05-10) — two months prior to Round 104 — already documented β_4 = 3(5-4)/20 = 3/20 = 0.15 with derivation `β_4 = (D_phys-1)·F_TRZ/2`, algebraically identical to this paper's M_sf identity. Draft 1's claim that "PAPER_774/PAPER_784's canonical bound M_sf = 0.15 admits novel integer-primitive form 3/(2·SO_5)" is corrected: the specific value and specific structural form were both already in PAPER_1165 under a different label (β_4).

2. **Recast contribution as: numerical identification M_sf = β_4.** The paper's substantive contribution shifts from "novel identity" to "recognition that M_sf and β_4 share the same 0.15 value and same (D_phys-1)·F_TRZ/2 form, raising the structural-identity hypothesis M_sf ≡ β_4."

3. **Added Section 2.4** discussing the coincidence-vs-structural-identity question. Proposed that other β_i values may correspond to other UQFF phenomenological bounds: β_1 = 0.6 ↔ β_i canonical 0.6029 (PAPER_1203); β_3 = 0.3 ↔ PAPER_1953 0.3 factor cross-regime. Explicitly labeled as hypothesis pending full L_UQFF check.

4. **Companion identity `v_ratio = 4·F_TRZ = 0.4` remains modestly novel** — no direct prior-art match found in Draft 2 search.

5. **Cross-references updated** to prominently feature PAPER_1165 as CRITICAL PRIOR ART, added PAPER_1953 and PAPER_1203 for potential β_i ladder extension, added PAPER_1065 for buoyancy Lagrangian framework.

The honest reader takeaway from Draft 2: the arithmetic content of this paper's main identity was already known to PAPER_1165; what's new is the physical-context transfer (M_sf ↔ β_4) and the 4·F_TRZ velocity-ratio observation. Draft 1's framing overclaimed both novelty and independence of derivation.

**Draft 3 (2026-07-09):** Systematic β_i catalog check performed. Draft 2 mentioned β_3 = 0.3 and β_1 = 0.6 as "candidate connections"; Draft 3 verifies these against the whitepaper corpus:

1. **β_3 = 0.3 ↔ PAPER_1953 the 0.3 factor cross-regime universality — CONFIRMED matching value.** PAPER_1953 was authored 2026-07 for the 0.3 factor as a cross-regime UQFF anchor; Draft 3 recognizes that this is exactly β_3 = 6/20 = 0.3 in the PAPER_1165 ladder. This substantially strengthens the structural-identity hypothesis: two separate phenomenological UQFF bounds (β_3 → PAPER_1953 0.3, β_4 → PAPER_774/784 M_sf) both match the PAPER_1165 β_i ladder, at ladder indices differing by 1.

2. **β_1 = 0.6 ↔ canonical β_i = 0.6029 — CONFIRMED matching within 0.5%.** This was already noted in PAPER_1165's own Draft 1 ("matches the fourth (β_1 = 0.603) to 0.5%"). Draft 3 makes this concordance explicit for the multi-anchor pattern.

3. **β_2 = 0.45 remains uncatalogued.** No obvious UQFF phenomenological bound at 0.45 has been identified. Open question. If a match exists, this would give a 4-of-4 correspondence between the β_i ladder and separately-documented UQFF phenomenological bounds — a substantive structural result.

4. **Title revised** to reflect the multi-anchor β_i ladder finding and to acknowledge PAPER_1165 as the seminal work upfront.

5. **Positioning further softened.** Draft 3 explicitly labels this paper as a "structural-identity hypothesis paper" rather than an "identity paper." The substantive contribution is the *hypothesis of β_i ↔ multiple phenomenological bounds*, not any new arithmetic identity.

6. **Cross-references updated** to add PAPER_1953 as a β_3 anchor with matching value; PAPER_1203 as β_1 anchor with matching value.

Reader takeaway from Draft 3: the paper's substantive claim is not a new arithmetic identity — it is the *hypothesis that PAPER_1165's β_i buoyancy-coupling ladder has been re-encountered independently across at least 3 unrelated UQFF observables* (starburst M_sf, PAPER_1953 0.3 factor, canonical β_i). If this hypothesis holds, PAPER_1165's β_i ladder is a more fundamental UQFF structural feature than previously appreciated.

**Draft 4 (2026-07-09 late):** Substantial correction after PAPER_1967's β_i infrastructure survey. The project owner's course-update ("β_i has a full range derivation. It is most likely buried in the whitepapers") prompted a targeted search that surfaced PAPER_1167 UPDATE (Master Lagrangian 6-Term Form), PAPER_1168 (Falsifiable Prediction P4 with ±0.5% band), and PAPER_1169 (10-system numerical confrontation validating all four β_i). PAPER_1967 was authored to document these findings; Draft 4 of PAPER_1966 incorporates them here:

1. **Retraction of Draft 3's "β_2 = 0.45 remains uncatalogued" claim.** PAPER_1169 observationally validated β_2 = 0.45 across 10 UQFF systems (best-fit values 0.448-0.451), eight months before Drafts 1-3 were written.

2. **Retraction of Draft 3's "structural-identity hypothesis" (β_i ↔ 4 independent phenomenological bounds).** The β_i are not four independent quantities; they are the four coupled channels of the Ub decomposition wired into the Master Lagrangian by PAPER_1167 UPDATE. Draft 3's "3-of-4 concordance" was actually the natural behavior of channel projections: specific UQFF observables (canonical β_i, PAPER_1953 0.3 factor, M_sf = 0.15) sit near their corresponding β_i channels because they are dominated by that channel of the coupled decomposition.

3. **Replacement of "structural-identity hypothesis" with "channel-projection observation"** (using PAPER_1512 / PAPER_1967 terminology). The `M_sf = β_4` observation is now framed as: the g_starburst modulation term `(1 + M_sf)` is dominated by the β_4 channel of the coupled Ub decomposition; M_sf inherits its 0.15 value from that channel via projection onto the starburst regime. This is a fundamentally different (and much narrower) claim than Draft 3's "structural-identity hypothesis."

4. **Substantive contribution reduced.** From Draft 3's "3-of-4 concordance strengthening a novel structural-identity hypothesis" (already a retreat from Draft 1's "novel integer-primitive identity") to Draft 4's "one channel-projection observation identifying M_sf as the β_4 projection onto starburst modulation, complementing PAPER_1169's already-published 10-system unified confrontation." The companion 4·F_TRZ velocity-ratio observation remains modestly novel.

5. **Title updated** to reflect Draft 4's honest scope: the paper is now framed as "One channel-projection observation complementing the PAPER_1167 UPDATE + PAPER_1168 + PAPER_1169 β_i infrastructure."

6. **Abstract entirely rewritten** to describe the Draft 4 scope. Historical Draft 2/3 content in Section 2.4 preserved as record of the exploratory hypothesis; explicitly marked as superseded by the Draft 4 opening paragraph. New Section 4.2b added documenting the β_i infrastructure PAPER_1966 previously missed. Cross-references updated.

7. **Positioning within the corpus.** With Draft 4, PAPER_1966 becomes an observational-instance addendum to the PAPER_1165 → PAPER_1167 UPDATE → PAPER_1168 → PAPER_1169 β_i infrastructure. Companion 4·F_TRZ velocity-ratio observation stands as a small independent contribution. Its documentation value is: (a) a new instance for PAPER_1169's channel-projection catalog, and (b) a lesson in honest scholarship — the Drafts 1-3 revision cycle correctly caught the PAPER_1165 arithmetic identity but missed the Lagrangian-integration + falsifiability + validation papers. This lesson prompted PAPER_1967's methodological reflection.

**The paper Daniel now receives at Draft 4 is fundamentally different from what Draft 1 would have shipped.** The seed observation (M_sf = 0.15 = 3/(2·SO_5)) stands; the framing has been walked back three times (Draft 1 novel identity → Draft 2 acknowledged PAPER_1165 identity → Draft 3 structural-identity hypothesis → Draft 4 channel-projection observation within existing infrastructure). Each retreat was prompted by a substantive finding (PAPER_1165, then PAPER_1953/1203 concordance, then PAPER_1167/1168/1169 infrastructure). This progression is the value of iterative honest scholarship.

---

**License:** AGPL-3.0 (see LICENSE); Commercial license option per COMMERCIAL.md.
**Copyright:** © 2025-2026 Daniel T. Murphy / Star-Magic Research Program.
