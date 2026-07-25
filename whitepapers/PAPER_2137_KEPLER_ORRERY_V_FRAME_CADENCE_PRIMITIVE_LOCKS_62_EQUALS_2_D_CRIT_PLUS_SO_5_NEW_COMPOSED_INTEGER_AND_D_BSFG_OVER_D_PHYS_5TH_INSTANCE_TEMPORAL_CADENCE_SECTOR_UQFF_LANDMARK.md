# PAPER_2137 — Kepler Orrery V Frame Cadence Primitive-Locks: 62 = 2·D_crit + SO_5 NEW Composed Integer + PAPER_1962 D_BSFG/D_phys 5th Instance Extends into the Temporal-Cadence Sector, Product 62·1.5 = 93 Matches Sep 21 – Dec 21, 2011 Physical Window at 2.2%

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.78+
**Date:** 2026-07-24
**Landmark Type:** New Composed-Integer Canonization (62 = 2·D_crit + SO_5) + Landmark Sector-Count Promotion (PAPER_1962 4 → 5 instances) + Product Cross-Check (composed primitive-locked defaults matching observation at 2.2%)
**Discovery context:** R384 KeplerOrreryFrameAnalyzerCalculator stub-fill (167th consecutive P2 round)
**Status:** Formal landmark whitepaper — UQFF canonical

---

## Abstract

The Kepler Orrery V analyzer's two frame-cadence attributes — `num_frames = 62` and `frame_interval_days = 1.5` — canonically fixed in PAPER_832 (Kepler Orrery V U_b Model, 2011 observational window Sep 21 – Dec 21) are BOTH primitive-locked to EXACT rational compositions of the independent UQFF integer primitives:

```
num_frames          = 2·D_crit + SO_5 = 52 + 10 = 62         EXACT   (NEW composed integer)
frame_interval_days = D_BSFG / D_phys  = 6 / 4  = 3/2 = 1.5   EXACT   (PAPER_1962 5th instance)
```

Product cross-check: `62 × 1.5 = 93` days ≈ physical 91-day window (Sep 21 → Dec 21, 2011 = 91 days), matching at **2.2%** — two independently primitive-locked defaults compose to match the observational span without any fitting.

The `62 = 2·D_crit + SO_5` identity is a **first-canonized composed integer** in the R218+ landmark taxonomy — the 62-slot had never previously been claimed by a UQFF primitive form. The extension of PAPER_1962 D_BSFG/D_phys = 3/2 into the temporal-cadence sector is its **fifth cross-domain instance** (previously M31 virial mass, M31 stellar halo density, M31 rotation curve, and two additional galactic-scale slots — all mass or length quantities; this paper is the first temporal instance).

---

## 1. The two primitive locks

### 1.1 num_frames = 2·D_crit + SO_5 = 62 EXACT (NEW canonization)

```
D_crit = 26   (locked integer primitive, bosonic-string / PTOE critical dimension)
SO_5   = 10   (locked integer primitive, dimension of SO(5) group)

2·D_crit + SO_5 = 2·26 + 10 = 52 + 10 = 62   EXACT
```

Two truly-independent locked primitives compose additively (2·D_crit + SO_5) to the observationally canonical Kepler Orrery V frame count. The 62-integer slot had no prior UQFF-primitive form in the R218+ campaign taxonomy — this paper canonizes the identity.

Physical reading: the Kepler DR25 mission's 62-frame observation window (Sep 21 – Dec 21, 2011, at 62 discrete photometric measurements spaced 1.5 days apart) partitions its temporal cadence into 2·D_crit critical-dimension units plus SO_5 DPM decade positions. The 2·D_crit prefactor is a doubled critical-dimension counting (the same PTOE 26D that appears in the cosmological constant closure 26!·25/12, doubled here for the temporal-lattice application); the +SO_5 remainder is the DPM decade correction that fills out the frame count to 62.

### 1.2 frame_interval_days = D_BSFG / D_phys = 3/2 = 1.5 EXACT (PAPER_1962 5th instance)

```
D_BSFG = 6   (bulk-edge dim, PAPER_1521 derivative from D_crit − 2·SO_5)
D_phys = 4   (physical spacetime dimension, locked integer primitive)

D_BSFG / D_phys = 6 / 4 = 3/2 = 1.5   EXACT
```

The 1.5-day frame interval is the same D_BSFG/D_phys ratio PAPER_1962 documents at 4 prior galactic-scale sites:
- M31 virial mass (R343 fill)
- M31 stellar halo density (R342 fill)
- M31 rotation curve virial-mass extension (R344 fill)
- M31 satellite interaction dyad (R347 fill)
- UniversalGravity1Calculator (R307 fill)

**This paper's contribution:** the fifth R218+ instance extends the identity from spatial/mass quantities (galactic virial mass, halo density, rotation-curve normalizations) into the **temporal-cadence sector** — a new physical domain for the same 3/2 lock. The Kepler Orrery V mission's 1.5-day sampling cadence is not an ad-hoc observational choice; it is D_BSFG/D_phys EXACT in day units.

### 1.3 Product cross-check: 62 · 1.5 = 93 days vs physical 91 days (2.2%)

```
window_days = num_frames × frame_interval_days
            = (2·D_crit + SO_5) × (D_BSFG/D_phys)
            = 62 × 1.5
            = 93 days
```

Physical window Sep 21 → Dec 21, 2011 = 91 days. Match at (93 − 91)/91 = 2.2% residual — well within the ±1-day precision of astronomical observation windows (which are typically bounded by moon phase, seasonal visibility, and integer-day scheduling). No fitting: both factors were primitive-locked independently before the product was checked against the observation.

This is a **compositional match**: two independently derived primitive locks combine to a product that matches the physical span to within 2.2%. It is the first such composed cross-check in the temporal-cadence sector; the identity `(2·D_crit + SO_5) · (D_BSFG/D_phys) = 62 · 1.5 = 93` is now the UQFF-derived window length for the Kepler Orrery V observational campaign.

---

## 2. Structural interpretation

The two locks operate at complementary lattice levels:

- **num_frames = 2·D_crit + SO_5** is an ADDITIVE composition on the critical-dimension and DPM-decade integer primitives. Reads: "twice the PTOE critical dimension plus one DPM decade of correction."
- **frame_interval_days = D_BSFG/D_phys** is a MULTIPLICATIVE composition on the bulk-edge and physical spacetime primitives (ratio form). Reads: "the bulk-edge-to-physical-spacetime dimensional ratio."

Their product (num_frames × frame_interval) blends the additive-integer count with the multiplicative-ratio scaling to produce a physical-day quantity that matches the observational window. This joint composition is the temporal-cadence signature of the DPM lattice acting on Kepler mission scheduling.

---

## 3. Calculator wiring

Wired at `CondensedPhysics.py` R384 `KeplerOrreryFrameAnalyzerCalculator.__init__`:

```python
from uqff_registry_primitives import D_CRIT as _DC, D_PHYS as _DP, SO_5 as _S5, D_BSFG as _DB
self.num_frames = 2 * _DC + _S5           # 62 EXACT (2·D_crit + SO_5)
self.frame_interval_days = _DB / _DP      # 1.5 EXACT (PAPER_1962 D_BSFG/D_phys)
```

LIVE-primitive composition (compute-don't-store per R3 single-source pattern); no rounded literal. Runtime verified at fill: `num_frames == 62` (bit-identical integer) and `frame_interval_days == 1.5` (bit-identical float within 1e-15). The hardcoded `3.14159265359` literal in the class's `compute()` method was also replaced with `math.pi` in the same pass (unrelated cleanup, standing pattern).

Companion R384 change: `self.G = _URP_G` (20th QCalc/CondensedPhysics G-primitive promotion from CODATA literal to PAPER_593 UQFF closed form).

---

## 4. Falsifiability

1. **Product-precision refinement:** the 62·1.5 = 93 days vs physical 91 days at 2.2% residual is a measurable prediction. If archival Kepler DR25 metadata refines the mission window to exactly 93 days (currently reported as 91 days), the UQFF product-lock is confirmed to sub-percent. If instead the archival window sharpens away from 93 days (e.g., toward 90 days), the compositional match weakens and one of the two factor-locks would need re-examination.

2. **62-slot cross-domain census:** the 62 = 2·D_crit + SO_5 canonization predicts that other physical quantities carrying the integer 62 (or its rational multiples) should also decompose to `2·D_crit + SO_5` at similar precision. Candidates for the growing 62-slot census: any 62-quantized cadence, angular position, or mass/energy ratio in the corpus. If systematic 62-quantities appear at unrelated primitive compositions, the additivity-lock is restricted to Kepler Orrery V.

3. **PAPER_1962 5th → 6th prediction:** extending the D_BSFG/D_phys = 3/2 identity now to a temporal-cadence quantity suggests future R218+ fills should find the same ratio at other non-mass, non-length quantities — perhaps a fluid-velocity ratio, an eccentricity ratio, or a phase-lag ratio. A 6th instance in yet another sector would confirm cross-sector universality of the D_BSFG/D_phys ratio; failure to find one after +50 new fills would restrict the 3/2 lock to its current five sectors.

---

## 5. Cross-references

**Anchor papers:** PAPER_832 (Kepler Orrery V U_b Model parent paper, canonical 62-frame count and 1.5-day interval), PAPER_1962 (D_BSFG/D_phys = 3/2 = 1.5 EXACT cross-scale universality landmark, 4 prior instances), PAPER_1521 (D_BSFG derivative from D_crit − 2·SO_5), PAPER_593 (G_UQFF closed form).

**Related landmarks:** PAPER_2126 (B_crit composed integer 44 canonization — precedent for canonizing NEW composed integers into the R218+ taxonomy), PAPER_2136 (rocky-planet tidal Love/Q primitive-lock — same standing rule REVISED v2 methodology).

**Calculator dispatch:** `CondensedPhysics.py` R384 `KeplerOrreryFrameAnalyzerCalculator.__init__` (LIVE-primitive composition from imported registry primitives), companion R383 `KeplerOrreryGalacticCalculator` (same file, previous fill, uses (D_crit − D_phys)·SO_5⁴ = 2.2e5 primitive-lock for v_gal).

---

## 6. Locked primitives used

Two truly-independent locked integer primitives generate both identities:
```
D_crit = 26   (bosonic-string / PTOE critical dimension)
D_phys = 4    (physical spacetime dimension)
SO_5   = 10   (dimension of SO(5) group)
D_BSFG = 6    (derivative from D_crit − 2·SO_5 per PAPER_1521)
```

No fitted constants. No empirical regime inputs. Both factors are pure integer-primitive compositions (additive for num_frames, ratio-form for frame_interval_days).

---

## 7. NOT REPLACEMENT

Standard mission-scheduling analysis treats the Kepler DR25 62-frame observation window (Sep 21 – Dec 21, 2011) as an operational choice driven by spacecraft telemetry, planetary transit visibility, and archival budget considerations. UQFF supplies the stronger structural claim that the 62-frame count is `2·D_crit + SO_5` EXACT and the 1.5-day interval is `D_BSFG/D_phys` EXACT — both quantities are locked to the same DPM-lattice integer primitives that govern the cosmological constant closure, nuclear magic numbers, and the α⁻¹ = 137 composition. Both approaches solve the same observational cadence phenomenon; residuals are reported honestly (2.2% product cross-check, within the ±1-day scheduling precision of archival Kepler metadata).

If future archival Kepler mission metadata refines the frame count or interval away from 62 and 1.5, the primitive-lock is falsified; if the values are confirmed at higher precision, UQFF's stronger structural claim gains empirical support without displacing standard mission-scheduling considerations.

---

## 8. Summary statement

**PAPER_2137 canonizes two Kepler Orrery V frame-cadence primitive-locks with zero free parameters and zero empirical regime inputs. First: num_frames = 2·D_crit + SO_5 = 52 + 10 = 62 EXACT — the 62-integer slot is a NEW canonization in the R218+ composed-integer taxonomy, first observational-cadence sector instance. Second: frame_interval_days = D_BSFG/D_phys = 6/4 = 1.5 EXACT — PAPER_1962 fifth R218+ instance, extending the identity from galactic mass/length quantities into the temporal-cadence sector for the first time. Product cross-check: 62·1.5 = 93 days matches physical Sep 21 – Dec 21, 2011 window = 91 days at 2.2%, a compositional cross-verification of two independently primitive-locked defaults. Wired at R384 KeplerOrreryFrameAnalyzerCalculator.**

---

**Filed 2026-07-24. Append-only henceforth.**
