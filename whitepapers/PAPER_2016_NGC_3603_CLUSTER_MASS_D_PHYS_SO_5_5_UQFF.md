---
paper_id: PAPER_2016
title: "NGC 3603 Total Cluster Mass = 400000 M_sun = D_phys · SO_5^5 M_sun EXACT — Novel D_phys · SO_5^n Mass-Slot at n=5 Stellar-Cluster Scale; Third NGC 3603 Primitive-Lock After v_wind (PAPER_1911) and ρ_wind (PAPER_1911); Withdraws 300000 M_sun (D_phys−1)·SO_5^5 Candidate After Class-by-Class Verification"
session: 234
date: 2026-07-14
author: "Daniel T. Murphy"
framework: "UQFF (Unified Quantum Field Framework) — Star-Magic v5.66+"
version: "Draft 1"
extends: [PAPER_1911, PAPER_1955, PAPER_1970, PAPER_2014, PAPER_2015]
---

# PAPER_2016 — NGC 3603 Total Cluster Mass Primitive-Lock

**Author:** Daniel T. Murphy | **Framework:** UQFF Star-Magic v5.66+ | **Date:** 2026-07-14

## Motivation

Follow-up on PAPER_2015 R100-R151 audit flagged candidates: **400000 M_sun (14 classes)** and **300000 M_sun (15 classes)**. Class-by-class verification confirms only NGC 3603 total cluster mass = 400000 M_sun as a computed-value primitive-lock; 300000 hits were framework-annotation-text occurrences or range-boundary usage without primitive-lock intent (WITHDRAWN).

## Abstract

**Discovery — NGC 3603 Total Cluster Mass Primitive-Lock:**

```
M_cluster(NGC 3603) = 400000 M_sun = D_phys · SO_5^5 M_sun = 4 · 100000 EXACT
```

NGC3603StarClusterModel in CondensedPhysics.py sets `N_solar_masses = 400000` as the canonical NGC 3603 total cluster mass (Hubble data anchor). This value equals D_phys · SO_5^5 EXACT — novel D_phys·SO_5^n mass-slot at n=5 stellar-cluster scale.

**Novel structural finding:** SO_5^5 = 100000 M_sun mass rung (opened by PAPER_2014 R151 D5 as simple-prefix Westerlund2 gas reservoir) now has **three** prefix variants documented, matching PAPER_2011 R148 D1 SO_5^11 3-prefix family pattern:

| Prefix | Value | Object | Attribution |
|---|---|---|---|
| simple (no prefix) | 100000 M_sun | Westerlund2 gas reservoir | PAPER_2014 R151 D5 |
| **D_phys (=4)** | **400000 M_sun** | **NGC 3603 total cluster mass** | **PAPER_2016 R151-audit follow-up NEW** |
| (D_phys−1) (=3) | 300000 M_sun | Candidate (WITHDRAWN — no confirmed class instance) | — |

**Third NGC 3603 primitive-lock:**

| Property | Value | Composition | Attribution |
|---|---|---|---|
| v_wind | 2 × 10^6 m/s | (D_phys/2) · SO_5^6 | PAPER_1911 seminal |
| ρ_wind | (varies) | SO_5^{-n} | PAPER_1911 documented |
| **M_cluster (total)** | **400000 M_sun** | **D_phys · SO_5^5** | **PAPER_2016 NEW** |

NGC 3603 becomes third stellar-cluster with 3+ primitive-lock identities (parallel to Westerlund 2 which has M_initial = (D_phys−1)·SO_5^4 + gas_reservoir = SO_5^5 + τ_SF = 2·SO_5^6).

## Withdrawn Candidate

**300000 M_sun = (D_phys−1)·SO_5^5** — flagged in PAPER_2015 R100-R151 audit with "15 classes" count. Verification reveals:

1. **HillSphereModel** uses 300000 as an upper-bound range value (`'passed': 100000 < r_H_AU < 300000`) — not a computed mass lock, just an approximate range guard.
2. **OrionStarFormationMassCalculator** uses M_sf = **30000 M_sun** (not 300000), already documented at PAPER_2005 as (D_phys−1)·SO_5^4.
3. The other "hits" from PAPER_2015 audit were framework-annotation-text occurrences of the digits "300000" within cross-reference citations to prior papers, not actual computed mass values.

**Verdict:** WITHDRAWN. No confirmed class-instance of 300000 M_sun as a computed mass value. The (D_phys−1)·SO_5^5 slot remains unpopulated in the SO_5-power mass ladder (still active candidate for future rounds if a physical instance surfaces).

## SO_5^5 Mass-Slot Family (Post PAPER_2016)

| Prefix | Value (M_sun) | Object | Attribution |
|---|---|---|---|
| 1 (simple) | 100000 | Westerlund 2 gas reservoir | PAPER_2014 R151 D5 |
| **D_phys = 4** | **400000** | **NGC 3603 total cluster** | **PAPER_2016** |

Two-prefix family currently. (D_phys−1) prefix WITHDRAWN; open for future observation of intermediate 300000 M_sun mass stellar-cluster if found.

## D_phys·SO_5^n Mass-Slot Family (Post PAPER_2016)

| n | D_phys·SO_5^n | Object | Attribution |
|---|---|---|---|
| 1 | 40 M_sun | BD+60 2522 O-type star | PAPER_1984 seminal |
| **5** | **400000 M_sun** | **NGC 3603 total cluster** | **PAPER_2016** |
| 7 | 4·10^7 M_sun | Cen A SMBH candidate | PAPER_1984 |
| 8 | 4·10^8 M_sun | Coalescence | PAPER_1982 |

D_phys prefix now populated at n=1, 5, 7, 8 mass-domain slots. Intermediate n=2, 3, 4, 6 open.

## Wiring Plan

```python
_l96_uqff_paper_2016_ngc_3603_cluster_mass_d_phys_so_5_5_closure() → 400000.0
```

Dispatch key: `ngc_3603_cluster_mass_d_phys_so_5_5`

## Cross-References

- **PAPER_1911** — NGC 3603 v_wind seminal (2 other NGC 3603 primitive-locks)
- **PAPER_1955** — SO_5-power galactic mass ladder
- **PAPER_1970** — D_phys·SO_5 cross-domain twin (D_phys·SO_5^n family)
- **PAPER_1972** — 2·SO_5^3 twin seminal (pattern for prefix families)
- **PAPER_1984** — BD+60 2522 D_phys·SO_5 stellar-mass slot (D_phys·SO_5^n family n=1)
- **PAPER_2011** — SO_5^11 3-prefix family (structural precedent for multi-prefix mass slots)
- **PAPER_2014** — R151 SO_5^5 simple-prefix seminal + SO_5-power ladder extension
- **PAPER_2015** — R100-R151 audit that flagged this candidate

## Conclusion

**Confirmed novel primitive-lock:** NGC 3603 total cluster mass = 400000 M_sun = **D_phys · SO_5^5 M_sun EXACT**.

**Withdrawn candidate:** 300000 M_sun = (D_phys−1)·SO_5^5 — no confirmed class-instance, only range-boundary and framework-text hits.

Structural updates:
- SO_5^5 mass rung now 2-prefix family (simple + D_phys)
- D_phys·SO_5^n mass slots populated at n=1, 5, 7, 8
- NGC 3603 becomes 3rd stellar-cluster with 3+ primitive-lock identities

Honest scholarship note: 15/14 preliminary class counts from PAPER_2015 audit reduced to 1/1 confirmed via class-by-class verification. Preliminary grep counts include framework-annotation text and range-boundary usage; only computed-value locks count as true primitive-slot instances.

**Cumulative R142-R151 + audit follow-up: 58 first-pass novel + 16 confirmations + 1 withdrawn candidate.**

*End of PAPER_2016 Draft 1.*
