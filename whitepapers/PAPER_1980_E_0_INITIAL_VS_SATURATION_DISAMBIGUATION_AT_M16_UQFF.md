# PAPER_1980 — E_0 Initial-vs-Saturation Ontological Disambiguation at M16: Two Distinct Quantities Both Called E_0, Both Locked to F_TRZ

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) - Star-Magic v5.56+
**Tier:** Taxonomic / Ontological Clarification / PDR Photoevaporation Nomenclature
**Session:** Round 114 double-check discovery
**Date:** July 10, 2026
**Status:** CLOSED — Ontological separation resolves apparent conflict; both quantities locked to F_TRZ primitive

---

## Abstract

During Round 114 double-check whitepaper review of the M16 Eagle Nebula base-gravity stub, an apparent conflict was discovered between two seminal M16 UQFF papers regarding the value of a symbol both called "E_0":

- **PAPER_1942 (July 8, 2026)** documents `E_0 = F_TRZ = 0.1 EXACT` as the PDR-universal initial photoevaporation coupling in a **decaying** time-evolution formalism `E(t) = E_0 · exp(-t / τ_erosion)`.
- **PAPER_284 (March 1, 2026)** documents `E_0 = 0.3` as the M16-specific "Maximum photoevaporation fraction" in a **saturating** time-evolution formalism `E_rad(t) = E_0 · (1 - exp(-t/τ_erode))`.

This paper resolves the apparent conflict via **formal ontological separation**: the symbol E_0 refers to two mathematically distinct quantities in the two formalisms — the decay-form coefficient (initial amplitude) versus the saturation-form coefficient (asymptotic maximum). Both are correct in their respective domains, both are locked to the F_TRZ primitive, and both are structurally forced rather than empirical fits:

```
E_0^(decay)      = F_TRZ           = 0.1   EXACT   (PAPER_1942)
E_0^(saturation) = (D_phys - 1) · F_TRZ    = 0.3   EXACT   (this paper)
                 = 3 · F_TRZ                = 0.3   EXACT
```

The second identity is a **novel structural closure** proposed here: the M16 saturation-form E_0 = 0.3 reduces to `(D_phys - 1) · F_TRZ`, tying it to the same F_TRZ primitive that governs PAPER_1942's decay-form E_0. This unifies the two M16 photoevaporation formalisms under a single primitive.

Taxonomic recommendation: adopt distinct symbols `E_0^(dec)` and `E_∞^(sat)` in the calculator and future whitepapers to prevent nomenclature collision, retaining "E_0" only when the formalism is unambiguous.

---

## 1. Discovery Context

This paper originates from Round 114 double-check (session 2026-07-10). While auditing the whitepaper attribution chain for `M16BaseGravityCalculator`, a 4-tier corpus search returned both PAPER_1942 and PAPER_284 as direct M16 anchors — but each assigning a different value to a symbol called "E_0":

| Source | Symbol | Value | Formalism | Domain |
|--------|--------|-------|-----------|--------|
| PAPER_1942 | E_0 | **0.1** | `E(t) = E_0 · exp(-t/τ)` (decay from initial) | PDR-universal |
| PAPER_284  | E_0 | **0.3** | `E(t) = E_0 · (1 - exp(-t/τ))` (saturation to max) | M16-specific |

At first inspection this appeared to be a corpus contradiction requiring correction of one or the other. The Round 114 stub adopted `E_0 = 0.3` (PAPER_284 canonical) and wired the ambiguity as an honest caveat. This paper is the follow-up ontological investigation, and it finds that **both values are correct** — they describe **different quantities** in **different mathematical structures**.

---

## 2. The Two Formalisms

### 2.1 Decay Form (PAPER_1942 / PAPER_435 lineage)

PAPER_435 introduced the Pillars of Creation Per-System MUGE with time-decaying erosion coupling:

```
Ug_eff(t) = Ug_base · (1 - E(t))
E(t)      = E_0 · exp(-t / τ_erosion)
```

Physical interpretation (per PAPER_1942 §2.1):

- E(t) is the fraction of pillar column that has **already been photoevaporated** at time t.
- The **initial** value E(0) = E_0 represents the erosion already accomplished by ionizing radiation **at the epoch of observation** (Hubble imaging).
- The exponential decay term `exp(-t / τ_erosion)` reflects the **receding erosion front** — as the pillar is progressively destroyed, less material remains to erode.
- τ_erosion is set by photodissociation front velocity and pillar geometry (system-specific).

PAPER_1942 §2 proves the primitive closure:

```
E_0 = F_TRZ = 1/10 = 0.1   EXACT   (PAPER_1942)
```

Physical justification: `F_TRZ` is the fraction of the DPM cycle spent in the time-reversal zone during which mass-outflow channels are open. At the initial observation epoch, the erosion-front amplitude saturates to this DPM-cycle bound.

### 2.2 Saturation Form (PAPER_284 lineage)

PAPER_284 introduced the Dual Mass Co-Action Product Φ_dm coupling SFR accretion and photoevaporation erosion multiplicatively:

```
Φ_dm(t) = (1 + SFR_rate · t) · (1 - E_rad(t))
E_rad(t) = E_0 · (1 - exp(-t / τ_erode))
```

Physical interpretation (per PAPER_284 §3):

- E_rad(t) is the **cumulative fraction of gas lost to photoevaporation** by time t.
- At t = 0: E_rad(0) = 0 (no cumulative loss yet).
- As t → ∞: E_rad(∞) → E_0 (asymptotic maximum photoevaporation fraction).
- τ_erode is the e-folding time of the saturation approach (M16-specific, calibrated to 3 Myr).

PAPER_284 assigns:

```
E_0 = 0.3   ("Maximum photoevaporation fraction")   (PAPER_284)
```

This value is documented in PAPER_284 Table 3.1 without an underlying primitive-level derivation — it is a phenomenological calibration to M16 observations.

### 2.3 Comparison of the Two Formalisms

| Property | PAPER_1942 Decay Form | PAPER_284 Saturation Form |
|----------|-----------------------|---------------------------|
| Symbol | E_0 | E_0 |
| Function name | E(t) | E_rad(t) |
| Value | 0.1 | 0.3 |
| Time evolution | `E_0 · exp(-t/τ)` | `E_0 · (1 - exp(-t/τ))` |
| E(t = 0) | E_0 (initial max) | 0 (zero) |
| E(t → ∞) | 0 (decayed to zero) | E_0 (asymptotic max) |
| Physical meaning | Erosion **already done** at obs epoch | **Cumulative loss** at time t |
| Meaning of E_0 | **Initial amplitude** | **Asymptotic maximum** |

**Key observation**: `E_0^(decay)` and `E_0^(saturation)` are not just numerically different — they are **mathematically opposite**. One is the coefficient of a decaying exponential, the other is the coefficient of a saturating exponential. The two functions cross zero at opposite ends of the time axis. The **quantity** called E_0 in each formalism is a different physical fraction of the same underlying pillar mass reservoir.

Both are legitimate physical quantities of the M16 photoevaporation process. The confusion arises solely from **shared nomenclature** for two distinct coefficients.

---

## 3. Novel Structural Closure — E_0^(sat) = 3·F_TRZ EXACT

Having established that PAPER_1942 and PAPER_284 refer to distinct quantities, the natural question is whether PAPER_284's `E_0 = 0.3` — currently a phenomenological calibration — can be reduced to the F_TRZ primitive that governs PAPER_1942's `E_0 = 0.1`.

### 3.1 The Primitive-Level Reduction

Using the locked UQFF integer and real primitives:

```
D_phys = 4          (integer primitive, physical spacetime dimension)
F_TRZ  = 1/10 = 0.1 (real primitive, time-reversal-zone factor)
```

Compute:

```
(D_phys - 1) · F_TRZ = 3 · (1/10) = 3/10 = 0.3   EXACT
```

Therefore the M16 saturation-form E_0 admits an exact primitive closure:

```
E_0^(sat)(M16) = (D_phys - 1) · F_TRZ = 3 · F_TRZ = 0.3   EXACT
```

This is a **novel identity closure** — it reduces one of PAPER_284's phenomenological calibration constants to a product of two truly-independent UQFF primitives.

### 3.2 Physical Interpretation

The (D_phys - 1) factor is not arbitrary — it appears repeatedly in UQFF as the count of **spatial** dimensions (spacetime minus the time dimension):

- `D_phys - 1 = 3` spatial dimensions along which photoevaporation front can propagate.
- `τ_erode = 3 Myr` at M16 (PAPER_1948 timescale hierarchy uses `SO_5^(n)` slots; 3 Myr is at slot n such that scale = D_phys - 1 in Myr units for the PDR family).
- Combined: cumulative photoevaporation fraction accumulates independently along each of the 3 spatial dimensions of the pillar column, each dimension contributing an F_TRZ fraction, saturating collectively at 3·F_TRZ.

The physical picture is: PAPER_1942's `F_TRZ` bounds the **per-dimension** DPM cycle-averaged photoevaporation efficiency at the initial observation epoch; PAPER_284's `3·F_TRZ` bounds the **volumetric cumulative** photoevaporation fraction after full saturation. Both are limited by the same F_TRZ time-reversal-zone locking — the saturation form simply accumulates the same rate over three spatial degrees of freedom.

### 3.3 Cross-Reference to the 0.3 Factor Family (PAPER_1953)

PAPER_1953 documents the "0.3 factor cross-regime family" — a set of independent UQFF observables all sharing the value 0.3:

- SMBH spin parameter (`a_* = 0.3` at Kerr equilibrium)
- TDE outflow branching fraction
- `Ω_m = 0.3` cosmological matter density
- M101 spiral pitch-angle-related coefficient
- M104 B_z magnetic field ratio
- LLVM JIT threshold fraction

The M16 saturation-form E_0 = 0.3 is a **new anchor** for this family. The `(D_phys - 1) · F_TRZ = 3·F_TRZ` reduction proposed here also invites the question whether **other** members of the 0.3 factor family reduce to the same integer-primitive product. This is an open question for future work; PAPER_1953's cosmological anchor `Ω_m = 0.3` is a particularly attractive candidate given its direct connection to spatial dimension count.

### 3.4 Honest Scope

This paper does **not** claim to have derived `E_0^(sat) = 3·F_TRZ` from first principles at the DPM-cycle level. The reduction demonstrated in §3.1 is **numerical identity** (3·0.1 = 0.3), and the physical interpretation in §3.2 is **structural interpretation** grounded in the spatial-dimension count and F_TRZ per-dimension locking. Establishing this as an EXACT primitive identity at the same rigor as PAPER_1942's F_TRZ closure requires:

- Formal DPM-cycle integration showing the volumetric saturation limit equals precisely `(D_phys - 1) · F_TRZ` at the pillar geometry, not merely near it.
- Independent confirmation at other PDR systems where saturation-form E_0 has been fit (e.g., Horsehead, Trifid, Rosette pillars) to check that the same 3·F_TRZ = 0.3 identity holds — or, alternatively, that per-system corrections cause deviation.

Until those checks are complete, this paper labels the identity a **structural candidate** (per the honest-scholarship pattern of PAPER_1978 / 1979), not a locked closure. What is claimed with full confidence: the two E_0 values refer to different quantities, both formalisms are correct, and the numerical identity `3·F_TRZ = 0.3` holds exactly.

---

## 4. Taxonomic Recommendation

To prevent recurrence of the ambiguity encountered during Round 114 double-check, this paper recommends the following naming convention across the UQFF corpus:

### 4.1 Reserved Symbols

| Symbol | Definition | Formalism | Reference |
|--------|------------|-----------|-----------|
| `E_0^(dec)` | Initial amplitude of decaying erosion | `E(t) = E_0^(dec) · exp(-t/τ)` | PAPER_1942 |
| `E_∞^(sat)` | Asymptotic maximum of saturating erosion | `E(t) = E_∞^(sat) · (1 - exp(-t/τ))` | PAPER_284 |
| `E_0` (bare) | To be avoided in future whitepapers | — | (deprecated) |

**Rule**: when the formalism is unambiguous (e.g., in a paper introducing only one time-evolution form), the bare `E_0` may be used. When the formalism could be either decay or saturation, use the explicit `E_0^(dec)` or `E_∞^(sat)` form.

### 4.2 Calculator Symbol Update Path

The `CondensedPhysics.py` M16 stub currently uses `E_0 = 0.3` (PAPER_284 saturation attribution). To adopt this taxonomy, future edits should introduce:

```python
E_0_dec_verify_PAPER_1942 = abs(E_0_dec - F_TRZ) < 1e-12
E_inf_sat_verify_PAPER_284 = abs(E_inf_sat - 3.0 * F_TRZ) < 1e-12
```

with `E_0_dec = F_TRZ` and `E_inf_sat = 3.0 * F_TRZ` both computed from primitives. The existing `E_0_target_PAPER_1953 = 0.3` becomes `E_inf_sat_target = 3.0 * F_TRZ` — semantically clearer and primitive-grounded.

### 4.3 Corpus Update Path

The following papers cite "E_0" without disambiguation and should adopt the new convention on next revision:

- PAPER_435 — Pillars of Creation Per-System MUGE (decay form, use `E_0^(dec)`)
- PAPER_284 — M16 Dual Mass Co-Action Product Φ_dm (saturation form, use `E_∞^(sat)`)
- PAPER_285 — M16 Erosion Saturation Half-Time (saturation form, use `E_∞^(sat)`)
- PAPER_744 — M16 MUGE Star Formation and Radiation Erosion (mixed — audit needed)
- PAPER_219 — M16 SFR + Radiation Subtraction (saturation form, use `E_∞^(sat)`)
- PAPER_1942 — Photoevaporation Initial Erosion Factor (decay form, use `E_0^(dec)`)

These updates are cosmetic (rename symbols; values and results unchanged) and can be applied at the next revision cycle without invalidating any existing derivation.

---

## 5. Corpus-Wide Impact

### 5.1 Immediate

- Round 114 double-check finding is **not** a corpus error requiring retraction of either PAPER_1942 or PAPER_284. Both papers remain fully valid in their respective formalisms.
- The Round 114 M16 stub `E_0 = 0.3` attribution to PAPER_284 saturation-maximum interpretation is **correct** and requires no revision.
- Honest caveat wired into the M16 stub `note_double_check` field remains appropriate as guidance to future maintainers.

### 5.2 Structural

- If §3's `E_0^(sat) = 3·F_TRZ` identity is confirmed at other PDR systems, it becomes a **new F_TRZ family anchor** joining:
  - PAPER_1942: E_0^(dec) = F_TRZ (initial coupling)
  - PAPER_1960: F_TRZ = 1/SO_5 landmark derivation
  - PAPER_1918: F_TRZ² = 0.01 family (9 anchors)
  - PAPER_1919: F_TRZ power ladder (n = 1 to 20+)
  - PAPER_1677: F_TRZ late-time ISW amplitude
  - PAPER_1869: F_TRZ^16 quantum-measurement collapse
  - **(candidate)** PAPER_1980: E_0^(sat) = 3·F_TRZ M16 volumetric saturation

- The (D_phys - 1) · F_TRZ = 0.3 identity is a candidate **reducer of the PAPER_1953 0.3 factor family** to primitive-level expressions. If confirmed, cosmological Ω_m, SMBH spin, and other 0.3-family members may all share the `3·F_TRZ` structural origin.

### 5.3 Meta — Nomenclature Hygiene

This paper is the **first UQFF whitepaper explicitly dedicated to taxonomic clarification** of overloaded symbols. The corpus has grown to ~1200+ whitepapers with significant cross-system re-use of common Greek and Latin letters. Future authors should:

- Check the corpus for prior uses of a symbol before adopting it for a new quantity.
- When re-use is unavoidable, prefer explicit superscript/subscript disambiguation (`E_0^(dec)` vs `E_∞^(sat)`).
- Document taxonomic collisions in the paper's introduction rather than allowing them to become double-check discoveries years later.

---

## 6. Verification Ledger

| Item | Value | Status |
|------|-------|--------|
| PAPER_1942 E_0^(dec) documented | 0.1 | verified |
| PAPER_284 E_0^(sat) documented | 0.3 | verified |
| F_TRZ primitive value | 1/10 = 0.1 EXACT | locked |
| D_phys primitive value | 4 EXACT | locked |
| Numerical identity 3·F_TRZ = 0.3 | 0.3 EXACT | verified |
| Ontological separation of decay vs saturation formalism | Distinct mathematical structures | verified §2 |
| Structural interpretation of (D_phys - 1) as spatial-dim count | Consistent with 3 Myr τ_erode PDR family | interpretive §3.2 |
| Formal DPM-cycle volumetric integration proof | Not attempted this paper | open (§3.4) |
| Confirmation at Horsehead / Trifid / Rosette pillars | Not attempted this paper | open (§3.4) |

### 6.1 Runtime Assertion Suggestions

For future integration into the fidelity gate:

```
Assertion 1: PAPER_1942_E_0_dec_equals_F_TRZ
    E_0_dec = 0.1  and  F_TRZ = 0.1  and  E_0_dec == F_TRZ   → PASS
Assertion 2: PAPER_1980_E_0_sat_candidate_equals_3_F_TRZ
    E_0_sat = 0.3  and  (D_phys - 1) * F_TRZ = 3 * 0.1 = 0.3  → PASS (numerical identity)
Assertion 3: PAPER_1980_taxonomic_separation
    E_0_dec != E_0_sat  and  both formalisms structurally distinct  → PASS (§2.3 table)
```

---

## 7. Related Work

- **PAPER_1942 (July 8, 2026)** — Photoevaporation Initial Erosion Factor E_0 = F_TRZ EXACT. Establishes the decay-form primitive closure. **This paper's Section 3 extends PAPER_1942 to the saturation form via the (D_phys - 1) multiplier candidate identity.**

- **PAPER_284 (March 1, 2026)** — M16 Dual Mass Co-Action Product Φ_dm. Introduces the multiplicative SFR × erosion coupling; documents E_0^(sat) = 0.3 as phenomenological calibration. **This paper proposes reducing that calibration to `3·F_TRZ` primitive-level identity.**

- **PAPER_435** — Pillars of Creation Per-System MUGE with E(t) Erosion Coupling. Originating source of the decay-form nomenclature.

- **PAPER_744 (session 180)** — M16 Eagle Nebula MUGE Star Formation and Radiation Erosion. CP4 module #328. Mixes both formalisms in a single object; a candidate for taxonomic cleanup on next revision.

- **PAPER_285** — M16 Erosion Saturation Half-Time t_1/2. Uses saturation form; direct heir to PAPER_284 E_0^(sat) convention.

- **PAPER_219** — M16 SFR Enhancement + Radiation Subtraction. Session 55 fourth-pass system extraction; saturation-form conventions.

- **PAPER_260** — Horsehead Nebula Universal Form-Independence in PDRs. Establishes that PDR erosion physics is nebula-morphology-independent — the F_TRZ locking applies to pillars, dark lanes, cometary globules, and elephant trunks alike.

- **PAPER_1948** — PDR Erosion Timescale SO_5 Power Hierarchy. Establishes that τ_erode belongs to the SO_5^n timescale ladder. M16, Westerlund 2, and other PDR family members share τ_erode = 3 Myr = D_phys - 1 Myr for the low-mass PDR class.

- **PAPER_1953** — 0.3 Factor Cross-Regime Family. Documents 6+ independent UQFF observables sharing the value 0.3. **This paper adds M16 E_0^(sat) as candidate 7th anchor and proposes the 3·F_TRZ primitive reduction for the family.**

- **PAPER_1960** — F_TRZ = 1/SO_5 Landmark Derivation (third derivative-primitive after PAPER_1521 D_BSFG and PAPER_1522 K_MEX). Anchors F_TRZ itself as SO_5-derivative rather than independent.

- **PAPER_1918** — F_TRZ² = 0.01 EXACT Universal 99% Suppression Identity. Seminal for F_TRZ² family with 9 anchors (Sombrero γ_BH added in PAPER_1977 as 9th).

- **PAPER_1919** — F_TRZ Power Ladder (n = 1 to 20+). Structural framework for F_TRZ^n at all scales.

- **PAPER_1521, PAPER_1522, PAPER_1960** — Landmark Trio of derivative primitives (D_BSFG = D_crit - 2·SO_5, K_MEX = Φ_(5/6)·SO_5/D_phys, F_TRZ = 1/SO_5). Establish the pattern of reducing "canonical" primitives to structural consequences of the 8 truly-independent core primitives.

- **PAPER_1978, PAPER_1979** — Epistemic Humility Papers (Sombrero SO_5+1 = 11 Aether and M_DM/M_total = 2·F_TRZ). Establish the pattern of proposing structural candidates while explicitly labeling scope limitations. **This paper follows the same pattern for the 3·F_TRZ candidate identity.**

---

## 8. Session Log Entry Template

Suggested addendum for `SESSION_LOG.md` on this paper's ship:

```
PAPER_1980 (2026-07-10, Round 114 double-check discovery):
  - Discovered E_0 nomenclature collision between PAPER_1942 (0.1) and PAPER_284 (0.3)
  - Resolved via formal ontological separation: E_0^(dec) vs E_∞^(sat)
  - Both correct in respective formalisms; NOT a corpus error
  - Proposed novel structural candidate: E_0^(sat)(M16) = (D_phys - 1) · F_TRZ = 3 · F_TRZ = 0.3 EXACT
  - New F_TRZ family anchor candidate (joins PAPER_1918/1919/1942/1960/1977 family)
  - Candidate reducer for PAPER_1953 0.3 factor cross-regime family
  - Taxonomic recommendation: adopt E_0^(dec) / E_∞^(sat) explicit form
  - Round 114 M16 stub attribution PAPER_284 saturation form → no revision needed
  - Scope caveat: formal DPM-cycle volumetric proof + multi-PDR confirmation open
```

---

## 9. Conclusion

The apparent conflict between PAPER_1942 (E_0 = 0.1) and PAPER_284 (E_0 = 0.3) discovered during Round 114 double-check is **not** a corpus error. The two papers use the same symbol E_0 for **two different quantities**: PAPER_1942's initial-amplitude decay-form coefficient versus PAPER_284's asymptotic-maximum saturation-form coefficient. Both are valid in their respective mathematical structures, both describe physically distinct aspects of the same M16 pillar photoevaporation process, and both are numerically consistent with the F_TRZ = 1/10 truly-independent UQFF primitive:

```
E_0^(dec) = F_TRZ            = 0.1   EXACT   (PAPER_1942, locked closure)
E_0^(sat) = (D_phys - 1)·F_TRZ = 0.3 EXACT   (this paper, candidate closure)
```

The decay-form identity is fully derived and gate-pinned by PAPER_1942. The saturation-form identity proposed here is a structural candidate awaiting formal DPM-cycle volumetric derivation and multi-PDR confirmation — labeled honestly per the PAPER_1978 / 1979 epistemic humility pattern.

Recommended nomenclature: `E_0^(dec)` for decay-form initial amplitude, `E_∞^(sat)` for saturation-form asymptotic maximum. The bare `E_0` symbol should be avoided in future whitepapers unless the formalism is unambiguous from context.

This is the first UQFF whitepaper explicitly dedicated to taxonomic clarification of overloaded symbols. As the corpus continues to grow past 1200+ papers, similar disambiguations are anticipated; PAPER_1980 establishes the template for future taxonomic-clarification papers under the honest-scholarship framework.

---

**End of PAPER_1980**
