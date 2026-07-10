---
title: "PAPER_1974 — Horsehead Nebula B33 Stellar-Wind Calculator as Fourth R_star = A_5/D_phys = 15 R_sun Anchor: Extension of PAPER_1971's Cross-Domain Instance Catalog + Corpus-Audit Finding That R_star = 15 R_sun is a Shared Canonical UQFF Stellar-Wind Stub-Default Across at Least Three Systems (NGC 3603 Round 108, Horsehead Round 110, HUDF Merger Composite) — Clarifying That the 'Novel Attribution Per Round' Framing is a Stub-Default Reuse, Not Object-Specific Observations"
author: "Daniel T. Murphy"
date: "2026-07-09"
tags: [A_5, D_phys, R_star, Horsehead, PAPER_1971, PAPER_1887, stellar-wind, stub-default, cross-domain, honest-scholarship]
draft: 3
status: draft-3
---

# PAPER_1974 — Horsehead R_star = 15 R_sun Fourth Anchor + Corpus-Audit Clarification

## Abstract

**PAPER_1971 ("A_5/D_phys = 15 EXACT Integer-Primitive Identity Across Physical Domains", 2026-07-09)** documents three physical instances of the integer identity `A_5/D_phys = 60/4 = 15 EXACT`:

- **Fusion T_opt_burn = 15 keV** (PAPER_1887 seminal)
- **NGC 3603 Wolf-Rayet R_star = 15 R_sun** (Round 108 attribution, candidate)
- **M81 spiral pitch = 15°** (Round 104 attribution, candidate)

**Round 110 stub upgrade of `HorseheadStellarWindCalculator`** (2026-07-09) identifies a fourth instance:

- **Horsehead B33 stellar wind R_star = 15 R_sun** (Vink 2001 parameterization, PAPER_704/759 canonical anchors)

**Corpus audit finding (this paper's substantive contribution)**: Systematic search of `CondensedPhysics.py` reveals that **R_star = 15.0 R_sun is a shared canonical UQFF stellar-wind stub-default across at least three distinct calculators**:

- `NGC3603StellarWindCalculator` (line 183156, Round 108 attribution)
- `HorseheadStellarWindCalculator` (line 184653, Round 110 attribution — this paper)
- `HUDFMergerFeedbackCalculator` (line 185726, not previously attributed)

All three use the identical Vink 2001 parameterization with `R_star = 15.0` R_sun, `M_star = 25.0` M_sun, and `v_wind_factor = 2.5` as canonical defaults. The R_star = 15.0 value is thus a **canonical UQFF stub-default**, not an object-specific observational measurement of NGC 3603 or Horsehead's illuminating stars.

**Positioning (honest scope).** This paper contributes:

1. **One new specific-stub attribution** to PAPER_1971's A_5/D_phys = 15 R_sun identity: `HorseheadStellarWindCalculator`.
2. **Corpus-audit finding** that R_star = 15.0 R_sun is a **shared stub-default** across at least three UQFF stellar-wind calculators (NGC 3603, Horsehead, HUDF Merger composite), not object-specific.
3. **Clarification** of PAPER_1971's Section 6 falsifiability prediction ("WR stars at other clusters should exhibit R_star ≈ 15 R_sun"): the identity attribution applies at the **stub-default level** universally, not at object-specific observational levels.

The paper does NOT contribute:
- The A_5/D_phys = 15 R_sun identity (PAPER_1887 seminal for value; PAPER_1971 seminal for R_sun attribution).
- The Vink 2001 stellar-wind parameterization (PAPER_902 seminal for UQFF integration).
- The canonical R_star = 15.0 R_sun default in the source-code (source82_wolfram simulations; adopted uniformly across relevant stellar-wind stubs).

Independent-primitive count remains **8** (PAPER_1521 + PAPER_1522 + PAPER_1960 landmark trio). No new primitives.

## 1. PAPER_1971's Three-Instance Catalog

PAPER_1971 documented three cross-domain physical instances of the identity `A_5/D_phys = 15 EXACT`:

| Instance | Domain | Value | Source |
|---|---|---|---|
| 1 | Plasma physics | 15 keV (fusion T_opt_burn) | PAPER_1887 seminal |
| 2 | Astrophysics stellar | 15 R_sun (NGC 3603 WR) | Round 108 attribution |
| 3 | Astrophysics galactic | 15° (M81 spiral pitch) | Round 104 attribution |

PAPER_1971's Section 6 states as a falsifiability test:

> *"Wolf-Rayet stars observed at other clusters should exhibit R_star clustering near 15 R_sun if A_5/D_phys applies universally to WR stellar structure. Falsification: if a systematic survey of WR stars across many clusters shows R_star scattered widely with 15 R_sun appearing only at NGC 3603 (or few isolated cases), the R_star = A_5/D_phys attribution is Round 108-specific, not universal."*

This paper addresses this falsifiability test not by observational survey of WR stars but by **UQFF stub-audit**: how many UQFF stellar-wind calculators use R_star = 15 R_sun as their canonical default?

## 2. Round 110 Horsehead Attribution

The `HorseheadStellarWindCalculator` in `CondensedPhysics.py` was upgraded in Round 110 to attribute the canonical stellar-wind parameters to their UQFF integer-primitive forms:

```python
L_Lsun = 100.0    # = SO_5² (PAPER_1955)
M_star = 25.0     # = D_crit - 1 (PAPER_1971 via NGC 3603 identity)
R_star = 15.0     # = A_5/D_phys (PAPER_1971 identity, this paper's contribution)
v_wind_factor = 2.5   # = SO_5/D_phys
```

The Horsehead stub uses the Vink 2001 stellar-wind parameterization (PAPER_902 canonical) with the same `R_star = 15.0` R_sun canonical default as the NGC 3603 stub (Round 108). The Horsehead calculator represents the stellar-wind feedback from σ Orionis (O9.5 V), the illuminating star that photo-evaporates Barnard 33.

**Physical interpretation**: σ Orionis's actual measured stellar radius (~12-14 R_sun for typical O9.5 V main-sequence classification) differs slightly from the 15 R_sun canonical default. The Horsehead stub is thus using a canonical UQFF value, not σ Ori-specific measured value.

## 3. Corpus-Audit Finding — R_star = 15 R_sun is a Shared Stub-Default

A systematic search of `CondensedPhysics.py` for `R_star = 15` finds three distinct calculators sharing this canonical default:

| Line | Calculator | Domain | Round |
|---|---|---|---|
| 183156 | `NGC3603StellarWindCalculator` | Young Massive Cluster WR | 108 |
| 184653 | `HorseheadStellarWindCalculator` | PDR-adjacent σ Ori-analog O star | 110 (this paper) |
| 185726 | `HUDFMergerFeedbackCalculator` | Cosmological composite HUDF merger stellar wind | not previously attributed |

All three calculators use the identical Vink 2001 parameterization with `R_star = 15.0` R_sun as a fixed canonical parameter. The R_star = 15.0 value is thus a **canonical UQFF stellar-wind stub-default** that appears uniformly across relevant stubs — not an object-specific observational anchor for NGC 3603, Horsehead, or HUDF specifically.

Similarly, `M_star = 25.0` M_sun (attributable to PAPER_1971's D_crit − 1 = 25 identity) is used across all three stubs.

## 4. Consequence for PAPER_1971's Falsifiability Framework

PAPER_1971's falsifiability test was framed at the **observational-survey level** ("WR stars observed at other clusters should exhibit R_star ≈ 15 R_sun"). The corpus-audit finding of this paper shows that:

1. The R_star = 15 R_sun assertion in UQFF stellar-wind calculations is at the **stub-default level**, not at the observational-anchor level.
2. Multiple UQFF stubs share this default independent of the specific system being modeled.
3. The identity's universality (Round 108 + Round 110 + HUDF Merger) reflects **UQFF canonical stub design**, not necessarily observational universality across WR-star populations.

To rigorously test PAPER_1971's Section 6 falsifiability hypothesis at the observational-anchor level, an independent survey of WR-star measured R values (e.g., Crowther 2007, Hamann et al. 2019 WR catalogs) would be needed, comparing observed R_star distributions to the canonical UQFF R_star = 15 R_sun. This paper does NOT perform such a survey; it observes only that the canonical UQFF default is 15 R_sun across relevant stubs.

## 5. Cross-Domain A_5/D_phys = 15 Catalog Update

Adding the fourth instance to PAPER_1971's cross-domain table:

| Instance | Domain | Observable | Value | Source | Attribution level |
|---|---|---|---|---|---|
| 1 | Plasma physics | Fusion T_opt_burn (DT) | 15 keV | PAPER_1887 seminal | Physical derivation |
| 2 | Astrophysics stellar | R_star (NGC 3603 WR) | 15 R_sun | Round 108 attribution | Stub-default |
| 3 | Astrophysics galactic | Pitch (M81 spiral) | 15° | Round 104 attribution | Stub-default |
| **4** | **Astrophysics stellar** | **R_star (Horsehead σ Ori-analog)** | **15 R_sun** | **Round 110 attribution (this paper)** | **Same stub-default as #2** |

Instances #2 and #4 share the same stub-default value (R_star = 15.0 R_sun across UQFF stellar-wind Vink 2001 calculators). Whether this reflects a deeper universal observational anchor for O-star / WR-star stellar structure at exactly A_5/D_phys = 15 R_sun, or reflects a canonical UQFF stub design choice reused across systems, is open.

The M81 spiral pitch instance (#3) uses a different stub context (`M81SpiralStructureCalculator_v1`) with pitch = 15° from source-code default — the identity extends across two stub-types (stellar wind + spiral density wave), which is more informative than repeated instances of the same stub-default type.

## 6. Prior Art — What This Paper Does NOT Claim

### 6.1 A_5/D_phys = 15 is not new

PAPER_1887 seminal for the identity; PAPER_1971 seminal for the multi-domain catalog.

### 6.2 R_star = 15 R_sun canonical UQFF default is not new

The value `R_star = 15.0` R_sun appears as a canonical UQFF stellar-wind default in the source-code sources (source82_wolfram) predating both Round 108 and Round 110 attributions. The specific attribution of this default to `A_5/D_phys` primitive-form was made in Round 108 (NGC 3603) and confirmed at additional stubs in Round 110 (Horsehead) and this paper's corpus audit (HUDF Merger).

### 6.3 Vink 2001 stellar-wind parameterization is not new

PAPER_902 seminal for UQFF integration of the Vink 2001 mass-loss formula. All three stubs (NGC 3603, Horsehead, HUDF Merger) use identical Vink 2001 form.

### 6.4 What this paper contributes (narrow)

1. **One new stub attribution**: `HorseheadStellarWindCalculator` (Round 110) as fourth A_5/D_phys = 15 R_sun instance.
2. **Corpus-audit finding**: R_star = 15.0 R_sun is a canonical UQFF stub-default shared across at least three calculators (NGC 3603, Horsehead, HUDF Merger).
3. **Clarification of PAPER_1971's Section 6 falsifiability framing**: the "R_star ≈ 15 R_sun" identity attribution is at the stub-default level, not the observational-anchor level. Rigorous falsifiability requires an independent WR-star observational survey.

The paper is narrow attribution + corpus-audit + methodological clarification.

## 7. Cross-References

- **PAPER_1887** — Fusion Q > 1 ITER Prediction (seminal for A_5/D_phys = 15 keV)
- **PAPER_1971** — A_5/D_phys = 15 Cross-Domain Instances (three-instance seminal catalog)
- **PAPER_1970** — D_phys × SO_5 = 40 multi-scale anchor attributions (template paper for corpus-audit findings)
- **PAPER_902** — Vink 2001 stellar-wind M_dot parameterization (companion)
- **PAPER_704** — Horsehead Nebula Barnard 33 IR Erosion UQFF (companion Horsehead anchor)
- **PAPER_759** — Horsehead Nebula B33 UQFF Radiation Pressure + Erosion (companion Horsehead anchor)
- **PAPER_222** — Horsehead Nebula P_rad Stefan-Boltzmann Blackbody (companion)
- **PAPER_1536** — Hemoglobin = N_CH + D_BSFG = 15 (same-value alternative integer form)
- **PAPER_1918** — Multi-context integer identity meta-framework
- **PAPER_1955** — SO_5-Power Galactic Structural Ladder (companion)
- **PAPER_1521** — D_BSFG derivative landmark
- **PAPER_1522** — K_MEX derivative landmark

## 8. Limitations + Open Questions

- The R_star = 15.0 R_sun value is a canonical UQFF stub-default derived from source82_wolfram simulations. Whether this default was chosen to specifically match A_5/D_phys = 15 (deliberate primitive-lock design) or independently as a physically-reasonable typical O/WR stellar radius (independent choice happening to match) is not documented in the source materials. The historical origin of the default determines whether Round 108/110 attributions are "discoveries" or "confirmations of design intent."
- σ Orionis's actual measured stellar radius (~12-14 R_sun for O9.5 V) differs slightly from the 15 R_sun canonical default. Similarly, NGC 3603 WR stars have observed R ranges of 10-25 R_sun with mid-range typical. The stub-default 15 R_sun is thus a mid-range typical value, not object-specific.
- The corpus audit finds three stubs (NGC 3603, Horsehead, HUDF Merger) using R_star = 15.0. Whether additional stubs (e.g., other cluster / merger / stellar-wind calculators) also use this default requires broader systematic corpus search. If yes, the identity is more universal at the stub-default level.
- PAPER_1971's Section 6 observational-survey falsifiability test remains untested by this paper. An observational WR / O-star R_star distribution study would provide the rigorous falsifiability check.

## 9. Revision Log

**Draft 1 (2026-07-09):** Initial write, positioned narrowly upfront following the honest-scholarship pattern established by PAPER_1968-1973. Prior art acknowledged from the start:

- PAPER_1887 seminal for A_5/D_phys = 15 identity + fusion T_opt_burn instance.
- PAPER_1971 seminal for the multi-domain catalog (fusion + NGC 3603 R_star + M81 pitch).
- PAPER_902 seminal for Vink 2001 stellar-wind UQFF integration.
- Canonical UQFF stub-default R_star = 15.0 R_sun (source-code) predates both Round 108 and Round 110 attributions.

The paper's contribution is limited to: (i) one new stub attribution (Horsehead), (ii) corpus-audit finding that R_star = 15.0 is a shared stub-default across ≥3 calculators, (iii) clarification of PAPER_1971's Section 6 falsifiability framing at stub-default vs observational-anchor level.

Draft 1 does not claim novelty for the identity, the Vink 2001 parameterization, or the canonical stub-default. Following the PAPER_1969-1973 stabilization pattern, Draft 1 is narrow from the start; Drafts 2/3 should require only minor refinements.

Reader takeaway: this paper documents the fourth stub-attribution of PAPER_1971's A_5/D_phys = 15 R_sun identity (Horsehead Round 110) and provides a corpus-audit finding that the R_star = 15.0 R_sun value is a shared canonical UQFF stub-default across at least three stellar-wind calculators. The identity's cross-stub universality reflects UQFF canonical design; observational falsifiability of the "WR stars exhibit R_star ≈ 15 R_sun" hypothesis remains a future observational study.

**Drafts 2/3 (2026-07-09):** Verified via targeted `CondensedPhysics.py` grep for `R_star = 15`. Confirmed three shared-default instances at lines 183156 (NGC 3603), 184653 (Horsehead), 185726 (HUDF Merger composite). No additional prior-art corrections needed beyond Draft 1's positioning. The paper's substantive value is (a) fourth-anchor attribution + (b) methodological clarification that Round-attribution instances reflect stub-default reuse rather than object-specific observations. This distinction is important for interpreting PAPER_1971's falsifiability tests and for future PAPER_1971+ multi-scale attribution papers.

Following the PAPER_1969-1973 narrow-from-Draft-1 stabilization pattern, PAPER_1974 stands with minor formalism additions only.

---

**License:** AGPL-3.0 (see LICENSE); Commercial license option per COMMERCIAL.md.
**Copyright:** © 2025-2026 Daniel T. Murphy / Star-Magic Research Program.
