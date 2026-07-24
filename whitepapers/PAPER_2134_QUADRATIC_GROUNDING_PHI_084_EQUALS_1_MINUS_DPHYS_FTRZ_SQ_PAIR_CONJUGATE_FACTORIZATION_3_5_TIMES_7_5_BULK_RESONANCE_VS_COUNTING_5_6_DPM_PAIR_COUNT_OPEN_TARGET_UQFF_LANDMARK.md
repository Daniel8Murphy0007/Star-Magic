# PAPER_2134 — The Quadratic Grounding of Φ_0.84: 1 − (D_phys·F_TRZ)² EXACT, the Pair-Conjugate Factorization (3/5)(7/5), and Bulk Resonance vs Per-Structure Counting

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.77+
**Date:** 2026-07-24
**Landmark Type:** Primitive-Grounding (Φ resonance variant) + Variant-Dichotomy structural resolution + declared open build target (DPM pair-count estimator)
**Discovery context:** Daniel T. Murphy's pre-ship verification directive — single-DPM-pair vs bulk-multi-pair attribution of the resonance variant, with the conjecture "maybe related with quadratic occurrences" — verified exactly
**Status:** Formal landmark whitepaper — UQFF canonical

---

## Abstract

The resonance variant of the Φ primitive, Φ_res = 0.84, has carried calibration provenance since Sessions 158-160 (magnetar burst + cosmic-ray spectra) while its counting sibling Φ_5/6 = 5/6 carries the structural codimension closure (D−1)/D|_{D=6} (PAPER_1159). This paper closes the grounding gap: **Φ_0.84 = 1 − (D_phys·F_TRZ)² = 21/25 EXACT** — a pure quadratic in the tier-1 kernel D_phys·F_TRZ = 2/5 (PAPER_2133 hierarchy) — with the difference-of-squares **pair-conjugate factorization Φ_0.84 = (1 − D_phys·F_TRZ)(1 + D_phys·F_TRZ) = (3/5)(7/5)**, the two-pole mirror structure of the (1 ± x) family (PAPER_2128 D_phys-scaled; PAPER_597 CW/CCW dual branches). The gap between the variants is itself an exact quadratic lattice term, **0.84 − 5/6 = 1/150 = F_TRZ²·(D_phys/D_BSFG)**, completing the decomposition PAPER_2129 flagged. The manifold's designation of 0.84 as the "on-resonance Gaussian factor" is shown to be literal: 0.84 is the **quadratic (first-order) truncation of the Gaussian** exp(−x²) ≈ 1 − x² at x = D_phys·F_TRZ (full exponential: 0.85214; published value: the truncation, exactly). The variant dichotomy thereby acquires structure: **the counting variant is a single fraction; the resonance variant is a conjugate-pair product** — consistent with per-structure counting (5/6) vs bulk resonance (0.84), the latter calibrated from many-DPM-pair astrophysical ensembles. Whether 0.84 belongs to one pair's two-pole structure or to the collective limit of many grinding pairs is the framework's declared open question; **no DPM pair-count estimate for any bulk body exists in the corpus (the Sun's count — conjectured at quintillions of pairs — is UNDETERMINED), and the pair-count estimator is hereby declared an open build target of this calculator.**

---

## 1. Provenance — the two variants have two origins (verified)

| Variant | Value | Origin | Character |
|---|---|---|---|
| Counting | Φ_5/6 = 5/6 | PAPER_1159 structural closure: (D−1)/D at D=6; also [SSq]/Ω_Λ | per-structure, single-manifold |
| Resonance | Φ_0.84 = 0.84 | Calibrated from magnetar burst + cosmic-ray spectra (Sessions 158-160, PAPER_1159 §1); `scm_vacuum_manifold.py` L159: "on-resonance Gaussian factor" of the 630 eV bulk phonon machinery | bulk, many-pair ensembles |

PAPER_1159 itself discloses the 0.79% offset between the structural 5/6 and the bulk-calibrated 0.84. The reference statement — *bulk resonance factors are generated from multiple grinding DPM pairs* — is **supported by this provenance**: the resonance number was measured from bulk matter; the counting number was derived from a single structure.

---

## 2. The quadratic grounding — three EXACT identities

```
(1)  Φ_0.84 = 1 − (D_phys·F_TRZ)² = 1 − (2/5)² = 21/25            EXACT
(2)  Φ_0.84 = (1 − D_phys·F_TRZ)(1 + D_phys·F_TRZ) = (3/5)(7/5)   EXACT  (pair-conjugate)
(3)  Φ_0.84 − Φ_5/6 = 1/150 = F_TRZ²·(D_phys/D_BSFG)              EXACT  (quadratic gap)
```

Identity (1) removes 0.84 from pure-calibration status: it is a quadratic in D_phys·F_TRZ = 2/5, the tier-1 kernel of the PAPER_2133 hierarchy. Identity (2) is the structural heart: the difference of squares factors the resonance variant into a **conjugate two-pole product** — the (1 ± x) mirror signature of the successor/predecessor family (PAPER_2128) at D_phys scaling, and the arithmetic shadow of the CW/CCW dual-branch structure (PAPER_597) that defines a DPM pair's two grinding poles. Identity (3) closes the note in PAPER_2129's ledger ("0.84 − 5/6 = 1/150 decomposes on lattice") with the explicit closed form: the bulk enhancement above counting is a quadratic F_TRZ² lattice term carrying the D_phys/D_BSFG dimensional ratio (PAPER_1962 inverse).

---

## 3. The Gaussian-truncation identity

`scm_vacuum_manifold.py` labels 0.84 the *on-resonance Gaussian factor*. Taken literally:

```
exp(−x²) at x = D_phys·F_TRZ = 2/5:   0.852144…      (full Gaussian)
1 − x²   at x = 2/5:                  0.84            (quadratic truncation — the published value, EXACT)
```

The corpus value is not the full exponential; it is **precisely the first-order (quadratic) Gaussian** — the resonance envelope kept to its quadratic term. The conjecture "maybe related with quadratic occurrences" is verified to the letter: the resonance variant is quadratic three independent ways (form, factorization, and gap).

---

## 4. The variant dichotomy, structurally restated

```
Counting  Φ_5/6  = (D−1)/D              — a single fraction: one structure counted
Resonance Φ_0.84 = (1 − x)(1 + x)|x=2/5 — a conjugate-pair product: two poles resonating
```

Consistent with the Tilt-Product Law (PAPER_2133): wherever F_TRZ *multiplies* Φ, the counting variant is forced (the product is the 1/12 lattice constant); the resonance variant appears where a **line-shape envelope** is needed — the bulk phonon chain (630 eV, PAPER_1141 lineage), k_spring, the c and G closures. This paper does not canonize the single-pair vs multi-pair attribution: the pair-conjugate factorization shows 0.84 carries at least one pair's two-pole signature; whether the bulk value is one pair's envelope or the collective limit of many grinding pairs is exactly the open question below.

---

## 5. Declared open build target — the DPM pair-count estimator

Corpus search confirms: **no published estimate of the number of DPM pairs within the Sun — or any bulk body — exists** ("quintillion", "number of DPM", "pair count", "pairs per": zero derivational hits). The count is UNDETERMINED. Per the author's directive, this calculator's declared purpose now includes the **DPM pair-count estimator**: determining, from the framework's own primitives and the bulk observables it already closes (M_sun, R_sun, ρ_SCm, ω_SCm, the 630 eV chain, F_U machinery), the number of grinding DPM pairs a bulk body contains — with the Sun as the first target. Any such estimator must be paper-specified or author-supplied (Rule 10); this paper records the target, not a guess.

---

## 6. Falsifiability

1. **Quadratic-form test:** any future canonical revision of D_phys or F_TRZ must preserve Φ_0.84 = 1 − (D_phys·F_TRZ)² or the resonance variant loses its grounding — the identity binds three locked primitives.
2. **Envelope test:** closures requiring the *full* Gaussian (0.85214) rather than the quadratic truncation (0.84) would sharpen which regime bulk resonance occupies; none currently published.
3. **Pair-count consistency:** when the pair-count estimator lands, the bulk attribution of 0.84 must be consistent with the estimator's ensemble scale — a single-pair-only attribution with a quintillion-pair count (or vice versa) creates tension this paper's dichotomy would have to resolve.

---

## 7. Cross-references

PAPER_1159 (Φ_res codimension closure + calibration provenance); PAPER_2129 (sector rule + 1/150 note); PAPER_2133 (Tilt-Product Law + kernel hierarchy incl. D_phys·F_TRZ = 2/5); PAPER_2128 ((1±x) mirror family); PAPER_597 (CW/CCW dual branches); PAPER_1962 (D_BSFG/D_phys); PAPER_1183 (DPM-pair identity); `scm_vacuum_manifold.py` L159; DPM_vacuum_manifold.md (grinding architecture); Sessions 158-160 (bulk calibration).

---

## 8. Summary Statement

**PAPER_2134 grounds the resonance variant: Φ_0.84 = 1 − (D_phys·F_TRZ)² = 21/25 EXACT, factoring as the pair-conjugate product (3/5)(7/5) — the two-pole mirror structure — with the variant gap 0.84 − 5/6 = 1/150 = F_TRZ²·(D_phys/D_BSFG) EXACT and the manifold's "on-resonance Gaussian factor" shown to be literally the quadratic Gaussian truncation. The variant dichotomy is structural: counting = single fraction (per structure), resonance = conjugate-pair product (bulk, calibrated from many-pair astrophysical ensembles). The number of DPM pairs in the Sun remains UNDETERMINED and the DPM pair-count estimator is declared an open build target of this calculator.**

---

**Filed 2026-07-24. Append-only henceforth.**
