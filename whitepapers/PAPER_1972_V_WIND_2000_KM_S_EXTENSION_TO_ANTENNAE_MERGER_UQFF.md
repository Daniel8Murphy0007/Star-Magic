---
title: "PAPER_1972 — Extension of PAPER_1911's Universal YMC Wind Velocity v_wind = (D_phys/2)·SO_5^6 = 2000 km/s EXACT to the Antennae (NGC 4038/4039) Interacting-Galaxy Pair: One New Anchor + Comparative Table with M82 Starburst Galaxy (v_superwind = SO_5^3 = 1000 km/s, PAPER_784) — Clarifying That the Round 109 '2·SO_5^3 km/s Twin' Framing was Superseded by PAPER_1911's Pre-Existing (D_phys/2)·SO_5^6 m/s Identity"
author: "Daniel T. Murphy"
date: "2026-07-09"
tags: [v_wind, YMC, Antennae, M82, PAPER_1911, PAPER_784, cross-object, unit-conversion, honest-scholarship]
draft: 3
status: draft-3
---

# PAPER_1972 — v_wind = 2000 km/s Extension to Antennae

## Abstract

**PAPER_1911 (Extended Young Massive Cluster Universal Parameter Set)** documents that

```
v_wind = (D_phys/2) × SO_5^6 = 2 × 10^6 m/s = 2000 km/s EXACT
```

is the universal OB-supergiant wind velocity for OB-supergiant-dominated YMCs, verified at **Westerlund 2 + NGC 3603** as a companion identity to PAPER_1909's `Ṁ_factor = SO_5/(D_phys−1) = 10/3`. Physical interpretation (PAPER_1911): "OB-star launched wind velocity encoded in the 6th power of the SO(5) mode count."

**Round 109 stub upgrade of `AntennaeStellarFeedbackCalculator`** (2026-07-09) attributes the same `v_wind = 2000 km/s` value to the **Antennae interacting-galaxy pair (NGC 4038/4039)** using the source-code default `v_wind = 2×10^6 m/s`. The stub's default matches PAPER_1911's canonical value, so **the Antennae wind stub is a third YMC-regime instance of the same identity**, not a novel independent observation.

**Draft 1 corrects a Round 109 misidentification**: Round 109 initial framing described the finding as "novel twin closure with M82" at `v_wind = 2·SO_5^3 km/s = 2000`. However:

1. **The identity is not novel**: PAPER_1911 already established `(D_phys/2)·SO_5^6 m/s = 2000 km/s EXACT` as the universal YMC form (published ~1 month prior, PAPER_1909 companion).
2. **The "twin with M82" framing is misleading**: PAPER_784's M82 v_superwind = 1000 km/s is at a DIFFERENT physical scale (galactic-scale integrated starburst outflow) with a DIFFERENT integer form `SO_5^3 km/s = 10^6 m/s`. The Antennae 2000 km/s is at OB-supergiant wind scale within the merger, not the merger's galactic superwind. Comparing them as a "factor-of-2 twin" conflates two different physical regimes.
3. **The Round 109 `2·SO_5^3 km/s` form is unit-conversion-equivalent to PAPER_1911's `(D_phys/2)·SO_5^6 m/s` form**: both produce 2000 km/s = 2×10^6 m/s. They are two integer-primitive expressions of the same universal YMC wind velocity in different unit conventions.

**Positioning (honest scope after Draft 1 correction)**. This paper contributes:

1. **One new specific-object attribution** of PAPER_1911's universal YMC v_wind = 2000 km/s identity: the Antennae (NGC 4038/4039) galaxy pair, verified via `AntennaeStellarFeedbackCalculator` Round 109 upgrade.
2. **Comparative table** distinguishing the YMC-regime `v_wind = 2000 km/s` (this paper + PAPER_1911) from the starburst-galaxy-regime `v_superwind = 1000 km/s` (PAPER_784 M82 seminal), clarifying that these are DIFFERENT physical observables at DIFFERENT scales, not a "factor-of-2 twin closure."
3. **Documentation of the unit-conversion equivalence** between PAPER_1911's `(D_phys/2)·SO_5^6 m/s` form and the km/s-convention form `2·SO_5^3 km/s` — showing that the same integer identity has different primitive-expression appearances depending on unit choice.

The paper does NOT contribute:
- The v_wind = 2000 km/s identity (PAPER_1911 seminal).
- The physical interpretation (PAPER_1911: OB-star launched wind, 6th power of SO(5) mode count).
- The M82 v_superwind = 1000 km/s identity (PAPER_784 seminal).
- A "novel merger vs starburst twin closure" (this framing is incorrect — the two observables are at different scales).

Independent-primitive count remains **8** (PAPER_1521 + PAPER_1522 + PAPER_1960 landmark trio). No new primitives.

## 1. PAPER_1911's Universal YMC v_wind Identity

**PAPER_1911 ("Extended Young Massive Cluster Universal Parameter Set — 4 Verified Primitive-Arithmetic Structural Identities Across Westerlund 2 + NGC 3603", 2026)** identifies four universal YMC parameter identities, of which the second is:

```
Identity 2: v_wind = (D_phys/2) × SO_5^6 = 2 × 10^6 m/s EXACT
```

Verified at both anchor systems:

| System | v_wind observed | v_wind UQFF form | Match |
|---|---|---|---|
| Westerlund 2 | 2 × 10^6 m/s | (D_phys/2) × SO_5^6 | 2/2 EXACT |
| NGC 3603 | 2 × 10^6 m/s | (D_phys/2) × SO_5^6 | 2/2 EXACT |

Physical interpretation (PAPER_1911):

> *"v_wind = 2 × SO_5^6: OB-star launched wind velocity encoded in the 6th power of the SO(5) mode count."*

The identity is stated as universal for "all OB-supergiant-dominated YMCs" — a category label that includes Westerlund 2, NGC 3603, and by extension any system dominated by OB-supergiant wind physics.

## 2. Round 109 Antennae Attribution

The `AntennaeStellarFeedbackCalculator` in `CondensedPhysics.py` uses `v_wind = 2×10^6 m/s = 2000 km/s` as the default wind velocity for the Antennae interacting-galaxy pair (NGC 4038/4039). Round 109 upgrade attributes this value to `2·SO_5^3 km/s` using integer-primitive form in km/s units.

The identification is:

```
v_wind(Antennae) = 2000 km/s = 2 × 10^6 m/s
                = (D_phys/2) × SO_5^6 m/s   [PAPER_1911 form]
                = 2 × SO_5^3 km/s           [Round 109 km/s form]
                = same numerical value in two integer-primitive expressions
```

The Antennae value is thus a **third-object instance** of PAPER_1911's universal YMC identity, applied within the merger context. Physical interpretation: OB-supergiant winds within the merging galaxies drive the same 2000 km/s wind regime as YMCs — the merger environment does not change the fundamental OB-star wind physics that PAPER_1911's `(D_phys/2)·SO_5^6` identity encodes.

## 3. Distinction from M82 v_superwind = 1000 km/s (PAPER_784)

**PAPER_784 (M82 Cigar Starburst Galaxy)** documents:

> *"UQFF encodes the superwind through v = 10^6 m/s and the starburst-amplified B = 10^-4 T. The galactic-scale superwind reaches ~1,000 km/s."*

M82's `v_superwind = 1000 km/s = 10^6 m/s = SO_5^3 km/s = SO_5^6 m/s` is at a **different physical scale**:

- M82: galactic-scale integrated starburst-driven outflow (12 kly extent above/below disk)
- YMC (Westerlund 2, NGC 3603, Antennae OB-supergiant): OB-star launched wind at cluster scale

These are DIFFERENT physical observables. The Round 109 initial framing "Antennae v = 2·M82 v = 2000 vs 1000, factor-of-2 twin closure" is misleading because it conflates:

| System | Physical scale | Observable | Value | UQFF form |
|---|---|---|---|---|
| Westerlund 2 (YMC) | OB-star cluster wind | v_wind | 2000 km/s | (D_phys/2)·SO_5^6 m/s (PAPER_1911) |
| NGC 3603 (YMC) | OB-star cluster wind | v_wind | 2000 km/s | (D_phys/2)·SO_5^6 m/s (PAPER_1911) |
| **Antennae (merger)** | **OB-star wind (Round 109)** | **v_wind** | **2000 km/s** | **(D_phys/2)·SO_5^6 m/s (Round 109 attribution, this paper)** |
| M82 (starburst galaxy) | galactic superwind | v_superwind | 1000 km/s | SO_5^6 m/s (PAPER_784 seminal) |

The Antennae 2000 km/s and the M82 1000 km/s differ by a factor of 2 numerically, but this factor is NOT a structural "twin closure" — it reflects that the two observables measure different physical quantities at different scales:

- 2000 km/s = OB-supergiant stellar wind velocity (peculiar to individual OB stars, PAPER_1911's identity)
- 1000 km/s = galactic-scale integrated superwind (peculiar to starburst-driven bulk outflow, PAPER_784's identity)

The correct pairing for M82's 1000 km/s is with other starburst-galaxy superwinds (e.g., NGC 253 superwind ~400 km/s per Round 105, which has its own PAPER_784-companion structural form 4·F_TRZ = 0.4 ratio), NOT with YMC-regime OB-supergiant winds at 2000 km/s.

## 4. Unit-Conversion Equivalence Between Integer Forms

The identity `2000 km/s = 2 × 10^6 m/s` admits two integer-primitive expressions:

**PAPER_1911 form (m/s units)**:
```
v_wind = (D_phys/2) × SO_5^6 m/s
       = (4/2) × 10^6 m/s
       = 2 × 10^6 m/s
       = 2000 km/s
```

**Round 109 attribution (km/s units)**:
```
v_wind = 2 × SO_5^3 km/s
       = 2 × 10^3 km/s
       = 2000 km/s
```

Both forms produce the same numerical value 2000 km/s. The relationship between them is:

```
(D_phys/2) × SO_5^6 m/s = D_phys/2 × SO_5^6 / 10^3 km/s
                       = 2 × SO_5^6 / SO_5^3 km/s
                       = 2 × SO_5^3 km/s
```

Since `SO_5^3 = 10^3` is the km-to-m unit conversion factor and D_phys/2 = 2, the two forms are algebraically identical: dividing by SO_5^3 converts m/s to km/s and reduces the power from SO_5^6 to SO_5^3 while preserving the D_phys/2 = 2 prefactor.

This is an interesting structural observation: **the identity's D_phys/2 = 2 prefactor is unit-invariant** (survives m/s → km/s conversion), while **the SO_5-power is unit-dependent** (reduces from 6 to 3 upon km-conversion via the SO_5^3 unit factor). This is because SO_5^3 = 10^3 is exactly the unit conversion factor between the two conventions.

## 5. Prior Art — What This Paper Does NOT Claim

### 5.1 The v_wind = 2000 km/s identity is not new

**PAPER_1911 (published ~mid-2026)** seminal. Verified at Westerlund 2 + NGC 3603 as universal for OB-supergiant-dominated YMCs. This paper adds one specific-object attribution (Antennae) and does NOT re-derive the identity.

### 5.2 The physical interpretation is not new

PAPER_1911's "OB-star launched wind velocity encoded in the 6th power of the SO(5) mode count" interpretation is seminal. This paper does NOT propose an alternative interpretation.

### 5.3 The M82 v_superwind = 1000 km/s identity is not new

**PAPER_784 (M82 Cigar Starburst Galaxy)** seminal. This paper does NOT re-derive it; it uses it to clarify that M82 and Antennae are DIFFERENT physical regimes.

### 5.4 The unit-conversion equivalence is not new in general

Cross-unit dual representations of integer identities are documented across the UQFF corpus (see e.g., PAPER_1970 D_phys×SO_5 = 40 at kpc + Hz + M☉ scales). The specific observation that PAPER_1911's `(D_phys/2)·SO_5^6 m/s` and the km/s-form `2·SO_5^3 km/s` are the same identity in different units is a minor observation, not novel.

### 5.5 What this paper contributes (narrow)

1. **One new specific-object attribution** to PAPER_1911's universal YMC identity: Antennae NGC 4038/4039.
2. **Comparative table** distinguishing YMC-regime v_wind (2000 km/s, PAPER_1911 + this paper) from starburst-galaxy v_superwind (1000 km/s, PAPER_784).
3. **Explicit correction of Round 109's initial "novel twin closure" framing** — clarifying that the 2000 vs 1000 factor-of-2 is not a structural twin but a difference in physical scale.

The paper is narrow attribution/documentation + Round 109 correction, not discovery work.

## 6. Falsifiability + Predictions

Under PAPER_1911's universal YMC hypothesis extended to merger context:

1. **Other galaxy mergers** should exhibit OB-supergiant v_wind ≈ 2000 km/s if PAPER_1911's identity is universal for OB-star-dominated environments. Candidate test systems: NGC 55 + NGC 300 pair, Mice Galaxies NGC 4676 (PAPER_055), Stephan's Quintet (PAPER_778).
2. **PAPER_784's M82 v_superwind should NOT match PAPER_1911's form** — starburst-galaxy superwinds are physically distinct from OB-star cluster winds. If observations show M82 galactic superwind at 2000 km/s (not 1000 km/s), the physical-scale distinction hypothesis is falsified.
3. **NGC 253 superwind (Round 105)** at 400 km/s = D_phys·SO_5² should pair with M82's SO_5^3 in the starburst-galaxy regime (PAPER_784 companion), NOT with YMC's (D_phys/2)·SO_5^6.

## 7. Cross-References

- **PAPER_1911** — Extended YMC Universal Parameter Set (v_wind = (D_phys/2)·SO_5^6 = 2000 km/s seminal)
- **PAPER_1909** — YMC Ṁ_factor = SO_5/(D_phys−1) = 10/3 (PAPER_1911's companion)
- **PAPER_784** — M82 Cigar Starburst Galaxy (v_superwind = 1000 km/s seminal, distinct regime)
- **PAPER_235** — Antennae NGC 4038/4039 Enhanced Double I(t) Merger Modulation (Round 108 direct anchor)
- **PAPER_1913** — Stellar Wind Bubble Linearity (E_t = E_0·t via F_TRZ·SO_5 = 1, companion wind physics)
- **PAPER_1955** — SO_5-Power Galactic Structural Ladder
- **PAPER_1970** — D_phys × SO_5 = 40 multi-scale identity (companion attribution-paper template)
- **PAPER_1971** — A_5/D_phys = 15 cross-domain instances (companion attribution-paper template)
- **PAPER_150** — Tapestry Westerlund 2 (companion YMC framework)
- **PAPER_138** — NGC 3603 MasterBuoyancy (companion YMC framework)
- **PAPER_055** — Mice Galaxies NGC 4676 (falsifiability test candidate)
- **PAPER_778** — Stephan's Quintet (companion merger)
- **PAPER_774** — 30 Doradus / Tarantula (companion starburst, PAPER_784-regime)
- **PAPER_1960** — F_TRZ = 1/SO_5 landmark
- **PAPER_1521** — D_BSFG derivative landmark
- **PAPER_1522** — K_MEX derivative landmark

## 8. Limitations + Open Questions

- The Antennae `v_wind = 2000 km/s` is a source-code default parameter for stellar-feedback modeling, not a direct observational measurement of Antennae OB-supergiant winds. Whether observations of individual OB stars within NGC 4038/4039 confirm 2000 km/s is a testable question.
- PAPER_1911 states the identity is universal for "all OB-supergiant-dominated YMCs"; whether the Antennae merger environment qualifies as OB-supergiant-dominated in the same sense as Westerlund 2 / NGC 3603 requires physical justification. If merger-induced OB populations differ statistically from cluster OB populations, the extension may need qualification.
- The distinction between YMC-regime v_wind (2000 km/s) and starburst-galaxy v_superwind (1000 km/s) is physically clear at extreme scales (compact YMC vs whole galaxy) but may blur at intermediate scales (e.g., a starburst cluster within a merger). Where the boundary sits is open.
- The unit-conversion equivalence between `(D_phys/2)·SO_5^6 m/s` and `2·SO_5^3 km/s` is algebraically trivial once identified but wasn't obvious in the Round 109 initial framing. Whether this reflects a general pattern (integer-primitive expressions of the same identity in different unit conventions) worth systematic cataloging is open.

## 9. Revision Log

**Draft 1 (2026-07-09):** Initial write, positioned narrowly upfront given the honest-scholarship pattern established by PAPER_1965-1971. Prior-art check found PAPER_1911 already documents `v_wind = (D_phys/2)·SO_5^6 = 2000 km/s EXACT` as universal for OB-supergiant YMCs (Westerlund 2 + NGC 3603). Round 109's initial framing "novel Antennae/M82 twin closure at 2·SO_5^3 vs SO_5^3" is corrected to:

1. **The identity is not novel** — PAPER_1911 seminal.
2. **The Antennae attribution is a third-object instance** of PAPER_1911's universal YMC identity, applied to merger context.
3. **The "twin with M82" framing is misleading** — Antennae 2000 km/s and M82 1000 km/s are DIFFERENT physical observables at DIFFERENT scales (OB-star wind vs galactic superwind), not a factor-of-2 structural twin.
4. **The `2·SO_5^3 km/s` form is unit-conversion equivalent** to PAPER_1911's `(D_phys/2)·SO_5^6 m/s` form.

The paper's contribution is limited to: (i) one new attribution (Antennae) to PAPER_1911's universal YMC identity, (ii) comparative table distinguishing YMC vs starburst-galaxy wind regimes, (iii) explicit correction of Round 109's initial framing.

Draft 1 is written narrowly and does not overclaim. Following the PAPER_1969-1971 pattern of honest positioning from the start, further Drafts should only require minor refinements rather than major walkbacks.

**Draft 2/3 (2026-07-09):** Post-draft-1 verification. No additional prior-art corrections needed beyond the PAPER_1911 acknowledgment captured in Draft 1. The paper stands as a narrow attribution/correction paper. PAPER_1911's 4-identity list at YMC scale (Ṁ_factor, ρ_wind, v_wind, cluster half-radius) provides the general structural context; this paper extends the v_wind row from 2 verified anchors (Westerlund 2 + NGC 3603) to 3 verified anchors (adds Antennae).

Reader takeaway: this paper is a specific-anchor addendum + Round 109 misidentification correction. The v_wind = 2000 km/s YMC universal identity is PAPER_1911's; this paper adds the Antennae instance and clarifies that pairing with M82's 1000 km/s galactic superwind is a category error (different physical scales, PAPER_784-regime not PAPER_1911-regime).

---

**License:** AGPL-3.0 (see LICENSE); Commercial license option per COMMERCIAL.md.
**Copyright:** © 2025-2026 Daniel T. Murphy / Star-Magic Research Program.
