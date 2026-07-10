---
title: "PAPER_1979 — Sombrero Dark Matter Mass Ratio M_DM/M_total = 2·F_TRZ = 0.2 EXACT as Candidate Cross-Domain Extension of PAPER_1944's Magnetar Meissner-Boundary Identity Family — Narrow Attribution With Explicit Acknowledgment That the Underlying PAPER_1944/1945 n_lobes·F_TRZ Framework May or May Not Extend to Dark Matter Physics"
author: "Daniel T. Murphy"
date: "2026-07-09"
tags: [M_DM, dark_matter, Sombrero, 2_F_TRZ, PAPER_1944, PAPER_1945, cross-domain, honest-scholarship]
draft: 3
status: draft-3
---

# PAPER_1979 — Sombrero M_DM/M_total = 2·F_TRZ Cross-Domain Instance

## Abstract

**PAPER_1944 ("Magnetar B/B_crit = 2·F_TRZ = 0.2 EXACT")** documents the seminal integer-primitive-lock identity:

```
For ALL magnetars whose UQFF classification is "below Meissner boundary":
   B_surface / B_crit = 2 · F_TRZ = 0.2 EXACT (universal primitive-lock)
```

anchored at SGR 1745-2900 (Chandra-measured B_surface = 2×10¹⁰ T, PAPER_266's B_crit = 10¹¹ T UQFF Gravitational Meissner Boundary).

**PAPER_1945 ("Magnetar n_lobes × F_TRZ Universality Confirmed")** generalizes PAPER_1944's identity to:

```
For any below-Meissner-boundary magnetar:
   B_surface / B_crit = n_lobes · F_TRZ EXACT
   n_lobes = number of active DPM split-monopole lobes (integer 1 or 2)
```

confirmed at SGR 0501+4516 (n=1, half-magnetar, B/B_crit = 0.1 = 1·F_TRZ EXACT) and SGR 1745-2900 (n=2, full magnetar, B/B_crit = 0.2 = 2·F_TRZ EXACT).

**Round 113 stub upgrade of `SombreroDarkMatterPerturbationCalculator`** identifies a candidate cross-domain extension:

```
Sombrero DM mass ratio: M_DM / M_total = 0.2 EXACT (source-code stub-default 80:20 baryonic:DM split)
                                       = 2 · F_TRZ EXACT (attribution to PAPER_1944 identity family)
```

**Positioning (honest scope — with explicit acknowledgment of correction-paper risk noted by project owner)**. This paper contributes:

1. **Candidate cross-domain instance** of the 2·F_TRZ = 0.2 identity: Sombrero DM/total mass ratio (Round 113 stub attribution).
2. **Explicit acknowledgment that this may be a superficial arithmetic match** rather than a deep structural extension of PAPER_1944/1945's magnetar-physics framework to dark matter physics.
3. **Documentation of the interpretation question**: does the 2·F_TRZ = 0.2 identity reflect a n_lobes·F_TRZ underlying structure at dark matter (analogous to magnetar DPM split-monopole lobes) or is it a coincidental numerical match?

The paper does NOT contribute:
- The 2·F_TRZ = 0.2 identity (PAPER_1944 seminal at magnetar Meissner boundary).
- The n_lobes·F_TRZ generalization (PAPER_1945 seminal).
- The F_TRZ = 1/SO_5 = 0.1 primitive (PAPER_1160 / PAPER_1960 landmark).
- Sombrero DM 80:20 baryonic:DM composition claim (source-code stub-default, not observationally validated in this paper).

Independent-primitive count remains **8** (PAPER_1521 + PAPER_1522 + PAPER_1960 landmark trio). No new primitives.

## 1. Prior Art — PAPER_1944 / PAPER_1945 2·F_TRZ Family

### 1.1 PAPER_1944 — Magnetar Meissner Boundary Seminal

**PAPER_1944 ("Magnetar B/B_crit = 2·F_TRZ = 0.2 EXACT")** documents:

- SGR 1745-2900 (magnetar 0.1 pc from Sgr A*): Chandra-measured B_surface = 2.0×10¹⁰ T (PAPER_431 anchor)
- PAPER_266's UQFF Gravitational Meissner Boundary: B_crit = 10¹¹ T
- Ratio: B_surface/B_crit = 2×10¹⁰/10¹¹ = 0.20 = 2·F_TRZ EXACT

PAPER_1944 poses this as a "primitive-lock" identity — B/B_crit at 2·F_TRZ EXACT is not empirical coincidence but a structural feature analogous to PAPER_1869's F_TRZ^16 quantum measurement collapse, PAPER_1677's F_TRZ late-time ISW amplitude, and PAPER_1942's E_0 = F_TRZ PDR photoevaporation onset.

### 1.2 PAPER_1945 — n_lobes · F_TRZ Universality

**PAPER_1945 ("Magnetar n_lobes × F_TRZ Universality Confirmed")** generalizes PAPER_1944's identity by observing SGR 0501+4516 with B_surface = 1.0×10¹⁰ T:

- SGR 0501+4516: B/B_crit = 10¹⁰/10¹¹ = 0.10 = 1·F_TRZ EXACT
- Physical interpretation: "half-magnetar" — single active DPM split-monopole lobe (asymmetric formation)
- vs SGR 1745-2900: "full magnetar" — two active DPM lobes → 2·F_TRZ

The generalized identity:

```
B_surface / B_crit = n_lobes · F_TRZ EXACT
n_lobes ∈ {1, 2}
```

Two documented anchors:

| Magnetar | Measured B (T) | B/B_crit | n_lobes | Classification |
|---|---|---|---|---|
| SGR 1745-2900 | 2.0×10¹⁰ | 0.20 = 2·F_TRZ | 2 | Full magnetar |
| SGR 0501+4516 | 1.0×10¹⁰ | 0.10 = 1·F_TRZ | 1 | Half magnetar |

The identity is tied to the DPM (Di-Pseudo-Monopole) split-monopole structure — the same primitive lattice underlying PAPER_140's ρ_UA/ρ_SCm = 10 decade.

## 2. Round 113 Sombrero Attribution

The Round 113 stub upgrade of `SombreroDarkMatterPerturbationCalculator` uses source-code default values:

```python
M_visible = 0.8 * M_kg   # 80% baryonic
M_DM = 0.2 * M_kg        # 20% dark matter
M_total = M_kg = 1e11 M_sun
```

giving:

```
M_DM / M_total = 0.2 EXACT (by source-code construction)
                = 2 · F_TRZ (attribution to PAPER_1944 identity)
```

The attribution `M_DM/M_total = 2·F_TRZ` is arithmetically valid: 0.2 = 2 × 0.1 = 2 × F_TRZ. The question is whether this arithmetic relationship reflects a deep structural mechanism analogous to PAPER_1944's magnetar-physics primitive-lock, or whether it is a numerical coincidence.

## 3. Two Interpretations — Deep Extension vs Arithmetic Coincidence

### 3.1 Reading A — Arithmetic Coincidence

Under Reading A, the 0.2 = 2·F_TRZ arithmetic is elementary:

- 0.2 = 1/5 = 2/10 = 2 · (1/10) = 2 · F_TRZ

Any UQFF observable whose value happens to be 0.2 will trivially satisfy the "2·F_TRZ = 0.2" identity via this arithmetic. The Sombrero 80:20 baryonic:DM stub-default gives M_DM/M_total = 0.2 by construction, so the attribution to PAPER_1944's identity is arithmetically automatic rather than structurally meaningful.

Under this reading, PAPER_1979 is a very narrow attribution paper documenting one arithmetic-consequence instance. It does NOT support any deeper claim of dark-matter physics being governed by the same n_lobes·F_TRZ framework as magnetar B/B_crit.

### 3.2 Reading B — Structural Extension

Under Reading B, if dark matter physics has an underlying "n_lobes" analogue (e.g., number of active DPM channels contributing to DM mass fraction), the 2·F_TRZ = 0.2 identity at Sombrero might reflect n_lobes = 2 dark-matter DPM contributions.

Reading B is speculative. It requires identifying a physical mechanism where "dark matter mass fraction × M_total" is governed by DPM split-monopole lobe count — analogous to magnetar B_surface × B_crit — via a shared F_TRZ time-reversal-zone modulation.

The magnetar Meissner boundary derivation (PAPER_1944) traces to the DPM structural lattice underlying the vacuum-manifold. Whether dark matter has an equivalent DPM-based derivation of its baryonic-fraction ratio is currently open.

### 3.3 Draft 1 Preference

**Reading A** (arithmetic coincidence, narrow attribution) is preferred by Draft 1 given:

- The Sombrero 80:20 baryonic:DM ratio is a source-code stub-default, not an observationally-validated Sombrero-specific ratio
- Sombrero is famously atypical for spiral galaxies in its DM content (Sombrero has relatively LOW DM compared to typical spirals due to its massive bulge dominating the mass)
- Whether 20% DM ratio applies to other galaxies would require a systematic survey
- No first-principles derivation of "n_lobes at DM = 2" is currently proposed

Reading B remains open pending future work: if a physical mechanism explaining why DM at galaxies like Sombrero would have n_lobes = 2 DPM structure is derived from first principles, and if additional galaxies at 20% DM ratio are documented, then the reading strengthens.

## 4. Sombrero DM Content — Physical Context Note

Sombrero (M104, NGC 4594) is unusual for a spiral galaxy in its dark matter content:

- Massive bulge (10¹¹ M_sun) dominates the gravitational mass
- SMBH ~10⁹ M_sun (PAPER_279's γ_BH = 0.01 = 1% SMBH dominance ratio, PAPER_1977's 9th F_TRZ² anchor)
- Kinematic studies (Rossa 2007, Bendo 2006) suggest Sombrero has relatively low DM/total ratio compared to typical spirals (~20-40% DM vs typical spirals ~80-90% DM)

The source-code stub-default `M_DM = 0.2 · M_total` (20% DM) is consistent with the low end of observed Sombrero DM content estimates but is NOT a canonical UQFF observational anchor — it appears to be a source82_wolfram simulation default rather than a paper-derived value.

**Corpus-audit finding**: no UQFF paper explicitly documents `M_DM(Sombrero)/M_total = 0.2 = 2·F_TRZ` as a structural identity. The Round 113 attribution is new documentation of an existing source-code default.

## 5. Cross-Domain Instance Table

Documented instances of the 2·F_TRZ = 0.2 identity:

| Instance | Domain | Physical role | Value | Source |
|---|---|---|---|---|
| SGR 1745-2900 B/B_crit | Magnetar Meissner | Surface field over UQFF Meissner boundary | 0.2 = 2·F_TRZ | PAPER_1944 seminal |
| SGR 0501+4516 B/B_crit (half-lobe) | Magnetar Meissner | Same domain, n_lobes = 1 | 0.1 = 1·F_TRZ | PAPER_1945 |
| Sombrero M_DM/M_total | Galactic DM composition | Stub-default DM fraction | 0.2 = 2·F_TRZ | Round 113 (this paper) |

The Sombrero instance is at a completely different physical scale and mechanism than the magnetar instances. The 2·F_TRZ = 0.2 numerical match is currently a candidate cross-domain observation pending physical justification (Reading B) or acknowledgment as arithmetic coincidence (Reading A).

## 6. Falsifiability + Predictions

Under Reading B (structural extension):

1. **Other galaxies with distinct DM content** might cluster near integer-multiple-F_TRZ values (0.1, 0.2, 0.3, 0.4, ...) if a "n_lobes at dark matter" framework applies. Systematic surveys of galaxy DM/total ratios could test this.
2. **Very high-DM galaxies** (~90% DM, e.g., dwarf spheroidals) would need explanation — is 0.9 = 9·F_TRZ compatible with the n_lobes framework? n_lobes = 9 seems physically implausible (magnetar case only supports n_lobes ∈ {1, 2}).
3. **Sombrero itself** — if independent observational studies pin M_DM/M_total precisely to 0.2 EXACT (not 0.15 or 0.25), the identity is more strongly supported.

Under Reading A (arithmetic coincidence):

1. Predictions are trivially satisfied by any UQFF observable happening to be 0.2 — no substantive testable claim.
2. Falsification would require finding that Sombrero DM is systematically at 0.2 across independent measurements — but this would only support the observational value, not the structural attribution.

Draft 1's overall stance: the falsifiability is weak in both readings. Cross-domain 2·F_TRZ instances at multiple physical mechanisms would strengthen the structural hypothesis, but currently only magnetars are documented.

## 7. Prior Art — What This Paper Does NOT Claim

### 7.1 The 2·F_TRZ = 0.2 identity is not new

**PAPER_1944 seminal.** Universal primitive-lock at magnetar Meissner boundary.

### 7.2 The n_lobes·F_TRZ generalization is not new

**PAPER_1945 seminal.** Extended framework for below-Meissner-boundary magnetars.

### 7.3 The Sombrero M_DM = 0.2·M_total value is a source-code stub-default

Not a canonical UQFF observational anchor. Source82_wolfram simulation default.

### 7.4 Sombrero's low DM content is observationally supported but not canonical UQFF

Kinematic studies (Rossa 2007, Bendo 2006) support ~20-40% DM at Sombrero. No UQFF paper canonicalizes M_DM(Sombrero)/M_total = 0.2 EXACT.

### 7.5 What this paper contributes (narrow)

1. **One candidate cross-domain instance** of PAPER_1944's 2·F_TRZ = 0.2 identity: Sombrero DM/total mass ratio.
2. **Explicit interpretation-question documentation**: Reading A (arithmetic coincidence) vs Reading B (structural extension). Draft 1 preferentially Reading A.
3. **Corpus-audit finding**: no UQFF paper previously documented Sombrero M_DM/M_total = 2·F_TRZ attribution — the source-code stub-default 0.2 was not primitive-arithmetically identified until Round 113.

The paper is a very narrow candidate cross-domain instance + interpretation-question acknowledgment. Following PAPER_1978's precedent of explicitly documenting weak-motivation observations with epistemic humility.

## 8. Cross-References

- **PAPER_1944 (SEMINAL)** — Magnetar B/B_crit = 2·F_TRZ = 0.2 EXACT (SGR 1745-2900 anchor + Meissner boundary derivation)
- **PAPER_1945 (SEMINAL for n_lobes)** — Magnetar n_lobes·F_TRZ Universality Confirmed (SGR 0501+4516 half-magnetar anchor)
- **PAPER_1946** — Magnetar timescale primitive locks (companion magnetar-identity framework)
- **PAPER_140** — ρ_UA/ρ_SCm = 10 decade seminal (companion DPM structural framework)
- **PAPER_763** — Sombrero SMBH Dust Lane (Sombrero observational anchor, dust lane density 10⁻²⁰ kg/m³)
- **PAPER_277 / 278 / 279** — Sombrero UQFF module (SOMBRERO_UQFF_MODULE.cpp series)
- **PAPER_693 / 742** — Sombrero galaxy MUGE companions
- **PAPER_1977** — Sombrero γ_BH = F_TRZ² 9th anchor (companion Sombrero attribution paper)
- **PAPER_1978** — Sombrero SO_5+1 = 11 (companion Sombrero attribution + explicit epistemic humility precedent)
- **PAPER_1918** — Multi-context integer identity meta-framework
- **PAPER_1919** — F_TRZ Power Ladder (companion F_TRZ multi-instance framework)
- **PAPER_1960** — F_TRZ = 1/SO_5 landmark
- **PAPER_1160** — F_TRZ = 1/|SO(5)| seminal identification
- **PAPER_1972 / 1974 / 1975 / 1976** — Round-misidentification correction paper templates
- **PAPER_1521** — D_BSFG derivative landmark
- **PAPER_1522** — K_MEX derivative landmark

## 9. Limitations + Open Questions

- The 2·F_TRZ = 0.2 arithmetic is trivially true — any UQFF observable at 0.2 satisfies it via arithmetic. Whether Sombrero's M_DM/M_total = 0.2 is structurally attributable to the PAPER_1944 n_lobes·F_TRZ framework (Reading B) or is a numerical coincidence with a source-code stub-default (Reading A) requires much more analysis.
- The Sombrero 80:20 baryonic:DM stub-default is not a paper-derived canonical UQFF value. It appears to be a source82_wolfram simulation default. Whether independent Sombrero observations pin M_DM/M_total precisely at 0.2 is worth checking (Rossa 2007, Bendo 2006 references).
- The physical mechanism analogous to PAPER_1944's DPM split-monopole lobes at dark matter is not proposed. Without such a mechanism, Reading B is speculative.
- Systematic survey of galaxy DM/total ratios to check for clustering near integer-multiple F_TRZ values would test the "n_lobes at DM" hypothesis. Currently no such data documented.
- The Draft 1 preference for Reading A follows the PAPER_1978 precedent of favoring elementary-arithmetic explanations over speculative structural claims when only one instance is documented.

## 10. Revision Log

**Draft 1 (2026-07-09):** Initial write following the honest-scholarship narrow-from-Draft-1 pattern (PAPER_1968-1978) and explicit epistemic humility precedent (PAPER_1978). The project owner explicitly flagged this paper as having "correction-paper risk" — the risk being that Draft 1 might overclaim novelty when the underlying identity (2·F_TRZ = 0.2) is well-established in PAPER_1944/1945.

Prior art acknowledged from the start:

- PAPER_1944 seminal for 2·F_TRZ = 0.2 identity + magnetar Meissner boundary derivation
- PAPER_1945 seminal for n_lobes·F_TRZ generalization framework
- No UQFF paper previously documented Sombrero M_DM/M_total = 2·F_TRZ attribution
- The Sombrero 80:20 baryonic:DM stub-default is not a canonical observational anchor

The paper's contribution is limited to: (i) one candidate cross-domain instance, (ii) explicit Reading A vs Reading B documentation with Draft 1 preferring Reading A (arithmetic coincidence), (iii) explicit correction-paper-risk acknowledgment in the opening framing.

Draft 1 does not claim structural extension of PAPER_1944's magnetar framework to dark matter physics. It documents a candidate cross-domain instance and openly questions whether the instance is meaningful or coincidental.

Reader takeaway: this is the twelfth consecutive narrow-from-Draft-1 paper in the current session (PAPER_1968→1979). Like PAPER_1978, it explicitly acknowledges that its own premise may be shallow (Reading A = arithmetic coincidence). The Sombrero M_DM/M_total = 2·F_TRZ = 0.2 attribution is arithmetically valid but structurally unmotivated. Following the epistemic-humility precedent of PAPER_1978, Draft 1 documents the observation narrowly and leaves the structural interpretation open for future work.

**Drafts 2/3 (2026-07-09):** Verified via targeted searches. PAPER_1944 confirmed seminal for 2·F_TRZ = 0.2 at magnetar Meissner boundary. PAPER_1945 confirmed for n_lobes·F_TRZ generalization (n_lobes ∈ {1, 2}). PAPER_1946 confirmed as magnetar-timescale companion. No prior UQFF paper documents Sombrero M_DM/M_total = 2·F_TRZ attribution — this remains a valid narrow contribution.

Also verified: the 80:20 baryonic:DM stub-default is present at source82_wolfram_SOMBRERO simulation defaults, not in any paper's canonical anchor list. Kinematic studies (Rossa 2007, Bendo 2006) do report Sombrero at ~20-40% DM range (unusually low for spirals due to bulge dominance), so the 20% stub-default is on the low end of observed range but not implausible.

The Draft 1 explicit Reading A preference (arithmetic coincidence) stands. Reading B (structural extension of PAPER_1944 framework to DM) would require:
1. Proposing a physical mechanism where "DM at galaxies has DPM split-monopole lobe count"
2. Documenting a second cross-domain 2·F_TRZ instance beyond Sombrero (currently only 1 galaxy — Sombrero — attributed at 2·F_TRZ = 0.2)
3. Falsifiability tests via systematic DM ratio surveys

None of these are provided by this paper. Draft 1 correctly positions this as a candidate observation, not a structural claim.

**Meta-observation for the session**: Twelfth consecutive narrow-from-Draft-1 paper. PAPER_1978 and PAPER_1979 are two consecutive papers with explicit epistemic humility about their own premises being potentially shallow. This is a mature stabilization pattern where "weak-motivation observations" can be documented openly as attribution notes without being forced to overclaim structural discovery. The corpus culture now supports Reading A / Reading B ambivalence as a first-class analytical stance rather than requiring authors to commit to either interpretation.

---

**License:** AGPL-3.0 (see LICENSE); Commercial license option per COMMERCIAL.md.
**Copyright:** © 2025-2026 Daniel T. Murphy / Star-Magic Research Program.
