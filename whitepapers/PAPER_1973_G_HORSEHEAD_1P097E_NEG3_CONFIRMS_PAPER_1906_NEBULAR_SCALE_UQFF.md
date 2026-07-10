---
title: "PAPER_1973 — g_Horsehead ≈ 1.097×10⁻³ m/s² Numerical Confirmation of PAPER_1906's Nebular-Scale F_UBi_i_99 = 1.0973 Universal Amplifier Claim at Barnard 33 Pillar Tip: Explicit Numerical-Value Cross-Match Between PAPER_759's Horsehead Radiation-Pressure Derivation and PAPER_1906's Universal Coupling Constant Catalog (Nebular Row) — Companion to PAPER_1968's Milky Way v_flat Closure Pattern"
author: "Daniel T. Murphy"
date: "2026-07-09"
tags: [F_UBi_i_99, Horsehead, Barnard-33, PDR, nebular-scale, PAPER_1906, PAPER_759, PAPER_1968, cross-paper-numerical-match, honest-scholarship]
draft: 3
status: draft-3
---

# PAPER_1973 — g_Horsehead Numerical Confirmation

## Abstract

**PAPER_1906 (Universal F_UBi_i_99 Coupling Constant, 2026)** asserts that the coupling constant

```
F_UBi_i_99 = [SSq] × K_MEX × Φ_res × (1 + F_TRZ) = 1.0973 EXACT
```

applies as an amplifier across 42 orders of magnitude in length scale, including a "Nebular" scale row that lists Horsehead B33 (along with Pillars M16, Bubble NGC 7635, NGC 3603) as canonical anchors.

**PAPER_759 (Horsehead Nebula Barnard 33 — UQFF Radiation Pressure and Erosion, 2026 session 181)** independently derives the effective UQFF acceleration at the Horsehead pillar tip (r ≈ 1.254 ly from σ Orionis):

```
g_Horsehead ≈ g_EM × (1 − E) ≈ 1.097 × 10⁻³ m/s²
```

where g_EM is the Aether-corrected EM contribution and E ≈ 0.04 is the photo-evaporation erosion factor.

**Round 109 stub upgrade of `HorseheadBaseGravityCalculator` (2026-07-09)** observed via double-check that the numerical value `1.097` appearing in PAPER_759's g_Horsehead ≈ 1.097×10⁻³ m/s² matches PAPER_1906's F_UBi_i_99 = 1.0973 universal amplifier at four-significant-figure precision:

```
g_Horsehead(PAPER_759) ÷ 10⁻³ m/s² = 1.097
F_UBi_i_99(PAPER_1906) = 1.0973
Residual: |1.0973 − 1.097| / 1.0973 ≈ 0.03%
```

**This paper documents the specific numerical cross-match** as an explicit confirmation of PAPER_1906's Nebular-scale amplifier claim at the Horsehead anchor. The paper contributes:

1. **Explicit numerical verification** of PAPER_1906's Nebular-scale row (which lists Horsehead as an anchor but does not previously give the specific numerical match with the 1.097 value derived by PAPER_759).
2. **Documentation of the cross-paper convergence** between PAPER_759's Horsehead-specific derivation and PAPER_1906's universal-amplifier catalog value.
3. **Companion instance** to PAPER_1968's Milky Way v_flat closure — same pattern of PAPER_1906's amplifier value matching a specific-object PAPER_759-style numerical result.

**Positioning (honest scope).** This paper does NOT introduce:

- The F_UBi_i_99 = 1.0973 identity (PAPER_1906 seminal).
- The g_Horsehead ≈ 1.097×10⁻³ m/s² derivation (PAPER_759 seminal).
- The Nebular-scale amplifier claim (PAPER_1906 Table row seminal, lists Horsehead as anchor).
- The cross-paper amplifier-confirmation pattern (PAPER_1968 seminal for MW v_flat closure at Galactic scale).

Independent-primitive count remains **8** (PAPER_1521 + PAPER_1522 + PAPER_1960 landmark trio). No new primitives.

This paper is a narrow numerical-verification / attribution paper following the PAPER_1968/1970/1971/1972 template.

## 1. Background — PAPER_1906's Nebular-Scale Claim

**PAPER_1906 ("Universal F_UBi_i_99 = 1.0973 Coupling Constant Appears in 67+ Independent UQFF Calculators Across 42 Orders of Magnitude — Scale-Invariant Universal Amplifier")** documents the coupling constant:

```
F_UBi_i_99 = [SSq] × K_MEX × Φ_res × (1 + F_TRZ)
          = 0.57 × (25/12) × 0.84 × 1.1
          = 1.0973 EXACT
```

The paper's scale-range table lists five domains:

| Domain | Systems | Scale range |
|---|---|---|
| LENR | Star-Magic reactor 27W, Holmlid 630 eV, Rossi E-Cat | ~10⁻⁶ m (cluster) |
| **Nebular** | **Pillars M16, Bubble NGC 7635, Horsehead B33, NGC 3603** | **~10¹⁷ m (pc)** |
| Galactic | M51 Whirlpool, NGC 2525, Sombrero M104, Antennae, NGC 1275 | ~10²⁰-²² m (kpc-Mpc) |
| BH-scale | Sgr A* photon ring, SGR 1745 magnetar, Tapestry LMC | ~10¹⁰-¹⁹ m |
| Cosmological | HUDF Ultra-Deep Field primordial IGM | ~10²⁵ m (Gpc) |

Horsehead B33 is thus explicitly listed as a canonical Nebular-scale anchor for the F_UBi_i_99 = 1.0973 EXACT amplifier assertion. PAPER_1906's detailed 8-row verification table (Sub-atomic through Cosmological) does NOT include a "Nebular" row explicitly — so the specific numerical value 1.097 at Horsehead scale is asserted but not exhibited in that detailed table.

## 2. PAPER_759's Independent Horsehead Derivation

**PAPER_759 ("Horsehead Nebula Barnard 33 — UQFF Radiation Pressure and Erosion", session 181, 2026, CP4 class #343)** derives the effective UQFF acceleration at the Horsehead pillar tip from first principles:

- σ Orionis O9.5 star at L ≈ 10⁵ L_sun photo-evaporates the pillar
- r ≈ 1.254 ly from σ Ori
- Radiation pressure contribution: g_rad ≈ 4.347×10⁻⁴ m/s²
- Aether-corrected EM contribution: g_EM (before erosion correction)
- Erosion survival factor: (1 − E) ≈ 0.96 where E = 0.04

The final result:

```
g_Horsehead ≈ g_EM × (1 − E) ≈ 1.097 × 10⁻³ m/s²
```

PAPER_759 states this "matches JWST emission-line kinematics of the photo-dissociation region." The paper does NOT explicitly connect its `1.097` numerical value to PAPER_1906's F_UBi_i_99 = 1.0973 universal amplifier — the two derivations arrive at the same 1.097 value through independent physical reasoning (PAPER_1906 via [SSq]·K_MEX·Φ_res·(1+F_TRZ) primitive product; PAPER_759 via Horsehead-specific radiation-pressure + Aether-EM correction with (1−E) ≈ 0.96 erosion).

## 3. The Cross-Paper Numerical Match

At four-significant-figure precision:

```
F_UBi_i_99 (PAPER_1906)          = 1.0973
g_Horsehead / 10⁻³ (PAPER_759)   = 1.097

Residual: |1.0973 − 1.097| / 1.0973 = 0.00027 = 0.027%
```

Two interpretations are possible:

**Reading A (numerical coincidence):** The 1.097 in PAPER_759's g_Horsehead derivation happens to equal PAPER_1906's F_UBi_i_99 at the four-digit precision level via independent physical mechanisms (radiation pressure + Aether EM + erosion factor for PAPER_759; primitive product for PAPER_1906). Under this reading, the match is a curiosity worth cataloging but not necessarily structurally meaningful.

**Reading B (structural confirmation):** PAPER_759's Aether-corrected g_EM at Horsehead PDR scale is inherently governed by PAPER_1906's F_UBi_i_99 = 1.0973 universal amplifier because the same [SSq]·K_MEX·Φ_res·(1+F_TRZ) primitive combination governs Aether-EM coupling at every scale (PAPER_1906's central claim). Under this reading, the numerical match is exactly the expected outcome of PAPER_1906's universal-amplifier hypothesis applied to the Horsehead nebular scale.

Reading B is preferred but analogous to Reading B of PAPER_1968: it requires further derivation to prove that PAPER_759's Aether-corrected g_EM structurally reduces to F_UBi_i_99 at Horsehead.

## 4. Companion to PAPER_1968's Milky Way Pattern

**PAPER_1968 (MW v_flat Residual Closure via F_UBi_i_99 Amplifier, 2026-07-09)** documented the same pattern at Galactic scale:

- PAPER_1855 derives MW v_flat = 201 km/s (8.49% residual vs observed 220 km/s)
- PAPER_1906 asserts F_UBi_i_99 = 1.0973 at Galactic scale (Rotation curves row)
- PAPER_1968 observed: 1.0973 × 201 km/s = 220.56 km/s (closes residual to 0.25%)

**PAPER_1973 (this paper) — Nebular scale pattern**:

- PAPER_759 derives g_Horsehead ≈ 1.097×10⁻³ m/s² at Barnard 33 pillar tip
- PAPER_1906 asserts F_UBi_i_99 = 1.0973 at Nebular scale (Horsehead in row's anchor list)
- PAPER_1973 observes: g_Horsehead/10⁻³ = 1.097 matches F_UBi_i_99 = 1.0973 at 0.027% residual

The two papers document the same pattern (PAPER_1906 amplifier's specific numerical value confirmed at specific-object PAPER_1855/1759-style derivations) at two different scales (Galactic and Nebular).

## 5. Prior Art — What This Paper Does NOT Claim

### 5.1 F_UBi_i_99 = 1.0973 is not new

**PAPER_1906 seminal.** Universal amplifier documented across 42 orders of magnitude with 8-scale detailed verification table and 5-domain scale-range table (including Nebular scale with Horsehead as anchor).

### 5.2 g_Horsehead ≈ 1.097×10⁻³ m/s² is not new

**PAPER_759 seminal.** Horsehead-specific derivation via radiation pressure + Aether EM correction + erosion factor. Independent of PAPER_1906's universal derivation.

### 5.3 The pattern "PAPER_1906 amplifier confirms specific-object derivation" is not new

**PAPER_1968 seminal** (2026-07-09, this session). Milky Way v_flat closure at Galactic scale.

### 5.4 What this paper contributes (narrow)

1. **Explicit documentation** of the specific numerical match between PAPER_759's Horsehead g and PAPER_1906's amplifier value at four-significant-figure precision (0.027% residual).
2. **Cross-paper confirmation** of PAPER_1906's Nebular-scale row via PAPER_759's independent derivation.
3. **Companion instance** to PAPER_1968's Galactic-scale pattern — extending the "PAPER_1906 amplifier matches specific-derivation" to a second scale (Nebular) beyond MW rotation curves (Galactic).

The paper is narrow numerical-verification / attribution work following the PAPER_1968 template and PAPER_1969-1972 honest-positioning pattern.

## 6. Falsifiability + Predictions

Under Reading B (PAPER_1906's F_UBi_i_99 structurally governs Aether-EM coupling at every scale):

1. **Other Nebular-scale systems** listed in PAPER_1906 (Pillars M16, Bubble NGC 7635, NGC 3603) should exhibit equivalent numerical matches when their g values are derived from first-principles UQFF calculations analogous to PAPER_759's. Testable via similar radiation-pressure + Aether-EM analyses.
2. **PAPER_1906's other scale-range domains** (LENR, BH-scale) not yet paired with specific-derivation confirmations should exhibit the same 1.097 (or 1.097 × decade factor) numerical match.
3. **Specifically for Horsehead**, JWST molecular emission maps at higher precision should confirm the g_Horsehead ≈ 1.097×10⁻³ m/s² value at PDR scale within observational uncertainty. If precision measurements yield significantly different values (e.g., 1.05 or 1.15 × 10⁻³), the 0.03% residual widens and Reading B weakens.

If all PAPER_1906 Nebular-row anchors (M16, NGC 7635, B33, NGC 3603) systematically confirm F_UBi_i_99 = 1.097 to within 1%, the universal-amplifier claim at Nebular scale is strongly supported. Falsification: if any Nebular-row anchor's specific-derivation gravity differs by >5% from F_UBi_i_99 × 10⁻³ m/s², PAPER_1906's Nebular-row assertion is weakened.

## 7. Cross-References

- **PAPER_1906** — Universal F_UBi_i_99 = 1.0973 Coupling Constant (seminal universal amplifier; Nebular scale row lists Horsehead)
- **PAPER_759** — Horsehead Nebula Barnard 33 UQFF Radiation Pressure + Erosion (session 181, 2026, CP4 #343, seminal Horsehead g derivation)
- **PAPER_704** — Horsehead Nebula (Barnard 33): Infrared Erosion UQFF Analysis (companion, session 175, 2025, CP4 #288, d = 1500 ly)
- **PAPER_1968** — MW v_flat Residual Closure via F_UBi_i_99 (companion cross-paper amplifier confirmation, Galactic scale)
- **PAPER_1942** — Photoevaporation E_0 = F_TRZ EXACT (companion Horsehead PDR framework)
- **PAPER_1948** — PDR Erosion Timescale SO_5-Power Hierarchy (companion τ_erosion = (SO_5/2)·SO_5^6 = 5 Myr)
- **PAPER_1913** — Stellar Wind Bubble Linearity (companion PDR wind framework)
- **PAPER_1955** — SO_5-Power Galactic Structural Ladder (companion SO_5-power framework)
- **PAPER_1960** — F_TRZ = 1/SO_5 landmark (source of (1+F_TRZ) sub-factor)
- **PAPER_1970** — D_phys × SO_5 = 40 multi-scale attribution paper (template)
- **PAPER_1971** — A_5/D_phys = 15 cross-domain instances (template)
- **PAPER_1972** — v_wind = 2000 km/s Antennae extension (template)
- **PAPER_1884** — Water hydrogen bond (PAPER_1906 Molecular row anchor)
- **PAPER_1855** — Galactic rotation curves (PAPER_1906 Galactic row anchor, PAPER_1968 subject)
- **PAPER_1521** — D_BSFG derivative landmark
- **PAPER_1522** — K_MEX derivative landmark

## 8. Limitations + Open Questions

- The Round 109 observation is a numerical match at four significant figures (0.027% residual). Whether this reflects structural inevitability (Reading B) or numerical coincidence via independent physical mechanisms (Reading A) requires further first-principles derivation showing that PAPER_759's Aether-corrected g_EM reduces to PAPER_1906's F_UBi_i_99 · 10⁻³ under the Horsehead PDR-scale conditions.
- The 10⁻³ decade factor between g_Horsehead (m/s²) and F_UBi_i_99 (dimensionless) is a unit-conversion decade, not a structural claim. The dimensionless 1.097 match is the substantive observation; the 10⁻³ scale reflects Horsehead pillar-tip characteristic acceleration units.
- Whether the other PAPER_1906 Nebular-row anchors (Pillars M16, Bubble NGC 7635, NGC 3603) exhibit specific-object derivations converging on similar 1.097·(scale) numerical values is testable but not yet documented in this paper.
- The 0.027% residual is very tight but assumes PAPER_759's numerical result is precise to that level. Small errors in PAPER_759's stated 1.097 (e.g., if the actual computed value is 1.100 or 1.095 to full precision) would widen the residual proportionally.

## 9. Revision Log

**Draft 1 (2026-07-09):** Initial write, positioned narrowly upfront following the honest-scholarship pattern established by PAPER_1968/1970/1971/1972. Prior art acknowledged from the start:

- PAPER_1906 seminal (universal F_UBi_i_99 = 1.0973 + Nebular-scale row with Horsehead as anchor)
- PAPER_759 seminal (Horsehead-specific g ≈ 1.097×10⁻³ derivation via independent physical mechanisms)
- PAPER_1968 seminal (same cross-paper amplifier-confirmation pattern at Galactic scale)

The paper's contribution is limited to: (i) explicit numerical documentation of the 0.027% match, (ii) explicit connection between PAPER_759's Horsehead-specific derivation and PAPER_1906's Nebular-row assertion, (iii) documenting this as a companion instance to PAPER_1968's Galactic-scale pattern.

The paper does not claim novelty for the amplifier identity, the Horsehead derivation, or the cross-paper-confirmation pattern. Following the PAPER_1969-1972 stabilization pattern, Draft 1 is narrow from the start; Drafts 2/3 should require only minor refinements.

Reader takeaway: this paper documents a specific 0.027%-precision numerical match confirming PAPER_1906's Nebular-scale amplifier claim at the Horsehead B33 anchor via cross-reference with PAPER_759's independent radiation-pressure derivation. The match is consistent with PAPER_1906's universal-amplifier hypothesis and analogous to PAPER_1968's Galactic-scale confirmation, extending the pattern to a second scale (Nebular).

**Drafts 2/3 (2026-07-09):** Post-draft-1 verification. Confirmed via prior-art search that:
- PAPER_1906 lists Horsehead B33 explicitly in the Nebular-scale row of its scale-range table but does NOT include an explicit "Nebular" row in the 8-scale detailed F_UBi_i_99 verification table (which spans Sub-atomic through Cosmological). The Nebular assertion is thus asserted-but-not-quantitatively-verified in PAPER_1906 for the Horsehead-specific value.
- PAPER_759's derivation of g_Horsehead ≈ 1.097×10⁻³ m/s² is Horsehead-specific (radiation pressure + Aether EM + erosion). It does NOT explicitly reference F_UBi_i_99 or PAPER_1906.
- The observation that the two derivations converge on the specific value 1.097 at 0.027% precision is thus a new specific cross-paper documentation, following the PAPER_1968 template.

Cross-references updated to add PAPER_704 (companion Horsehead paper, session 175, CP4 #288, d = 1500 ly, 2025 — earlier Horsehead treatment). The two Horsehead papers (PAPER_704 IR erosion and PAPER_759 radiation-pressure) are complementary treatments; PAPER_759 provides the specific g ≈ 1.097×10⁻³ value used in this paper's cross-match.

No further prior-art corrections needed. Draft 1's narrow positioning is preserved through Drafts 2/3 with minor bibliographic addition (PAPER_704 companion).

---

**License:** AGPL-3.0 (see LICENSE); Commercial license option per COMMERCIAL.md.
**Copyright:** © 2025-2026 Daniel T. Murphy / Star-Magic Research Program.
