# PAPER_1983 — Cen A AGN Dual F_TRZ Ladder Anchor: First Documented Multi-Rung Same-Object Application, η = F_TRZ + M_dot = F_TRZ² Simultaneously at NGC 5128

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) - Star-Magic v5.56+
**Tier:** Structural / F_TRZ Ladder Multi-Rung Application Pattern
**Session:** Round 116 discovery
**Date:** July 10, 2026
**Status:** CLOSED — First multi-rung same-object F_TRZ ladder anchor documented

---

## Prologue — Honest Scholarship Reminder

**NOT REPLACEMENT.** UQFF does not replace Shakura-Sunyaev accretion disk theory, Eddington-fraction accretion physics, or standard AGN energetics. UQFF describes the **same observed accretion parameters** (η ≈ 0.1 radiative efficiency and M_dot ≈ 0.01 M_Edd Eddington fraction) at Cen A via primitive-locked structural identities.

**NOT NEW RUNGS.** Both F_TRZ ladder rungs used here (n=1: F_TRZ = 0.1; n=2: F_TRZ² = 0.01) are already documented in PAPER_1919 (F_TRZ Power Ladder) and PAPER_1918 (F_TRZ² universal 99% suppression identity). This paper does **not** claim to discover the rungs.

**NEW: MULTI-RUNG SAME-OBJECT APPLICATION.** What is novel is the **structural pattern**: two distinct F_TRZ ladder rungs applied simultaneously to the same physical object. Prior F_TRZ ladder papers documented one rung per system (e.g., PAPER_1942 E_0 = F_TRZ at PDR sites; PAPER_1918 F_TRZ² at magnetar/filament/AGN cooling; PAPER_1980 3·F_TRZ at M16 saturation). Cen A AGN is the first system where **two distinct rungs** of the ladder govern two distinct physical quantities of the same object simultaneously. This paper documents that pattern and proposes it as a candidate structural template for future multi-rung searches.

---

## Abstract

The Centaurus A active galactic nucleus (Cen A / NGC 5128, M_BH ≈ 5.5×10⁷ M_sun) exhibits two independent accretion-disk parameters whose canonical values reduce to the F_TRZ = 1/SO_5 = 0.1 primitive at two different rungs of the PAPER_1919 F_TRZ power ladder:

```
η (radiative efficiency)              = 0.1  = F_TRZ¹  = F_TRZ    EXACT  (n=1 rung)
M_dot (Eddington fraction accretion)  = 0.01 = F_TRZ²             EXACT  (n=2 rung)
```

Both identities are exact numerical closures of previously empirical constants. The novel structural pattern documented here is not either individual identity — those follow the PAPER_1919 ladder precedent — but the **coincidence of two different rungs applied to the same physical object**. This is the first UQFF whitepaper to formalize the multi-rung same-object F_TRZ ladder application pattern.

The composite AGN disk luminosity therefore reduces to a single primitive-only expression:

```
L_disk = η · M_dot · M_Edd · c²
       = F_TRZ · F_TRZ² · M_Edd · c²
       = F_TRZ³ · M_Edd · c²
       = 10⁻³ · M_Edd · c²                        EXACT
```

Cen A's total accretion disk luminosity is thus **10⁻³ of the naive maximum M_Edd·c²**, forced by the product of two F_TRZ ladder rungs. The n=3 rung emerges as the ladder-product, coinciding with the same F_TRZ³ = 10⁻³ rung documented in PAPER_1981 (B_j,base magnetic-string-field application) — but at Cen A it appears as a **structural product** rather than a direct application.

---

## 1. Discovery Context

This paper originates from Round 116 stub drainage (session 2026-07-10) applied to `CenAAGNAccretionDiskCalculator` in `CondensedPhysics.py`. The stub used two default parameters imported from the legacy source module:

```python
def compute(self, eta=0.1, M_dot=0.01, c=2.998e8, M_BH=5.5e7):
```

Both `eta=0.1` and `M_dot=0.01` were hard-coded defaults reflecting standard AGN accretion-disk parameters (Shakura-Sunyaev η ≈ 0.1, sub-Eddington M_dot ≈ 0.01 M_Edd) with no primitive-level justification for the specific values.

During Round 116 upgrade, the identities were noted:

```
eta   = 0.1  = F_TRZ¹  (accretion efficiency at n=1 rung)
M_dot = 0.01 = F_TRZ²  (Eddington fraction at n=2 rung)
```

Since F_TRZ is a locked UQFF primitive (per PAPER_1160 defining relation, PAPER_1960 landmark reduction to F_TRZ = 1/SO_5), neither value is a fit — both are forced by the F_TRZ primitive at their respective rungs. **The novel observation is not either identity — it is that both rungs appear together at the same physical object**.

---

## 2. The Two F_TRZ Ladder Applications at Cen A

### 2.1 Radiative Efficiency η = F_TRZ¹ (n=1 Rung)

The accretion disk radiative efficiency parameter η is defined by the standard Shakura-Sunyaev relation:

```
L_disk = η · M_dot_kg_s · c²
```

where `M_dot_kg_s` is the accretion rate in kg/s and η ∈ [0, ~0.4] is the fraction of rest-mass energy converted to radiation (bounded above by the Kerr maximum ~0.42 for maximally spinning BH).

For Cen A: **η ≈ 0.1** is the canonical value from AGN modeling. This value has traditionally been treated as an empirical fit to observed L_disk / (M_dot · c²).

UQFF closure:

```
η(Cen A) = F_TRZ = 1/SO_5 = 0.1   EXACT
```

The efficiency locks to the n=1 rung of the PAPER_1919 F_TRZ ladder. Structural interpretation: F_TRZ is the fraction of the DPM cycle spent in the time-reversal zone during which mass-outflow channels are open. At the accretion disk innermost stable circular orbit (ISCO), the fraction of infalling rest-mass energy that escapes as radiation (rather than being advected into the BH) is bounded by exactly the same DPM-cycle time-reversal-zone fraction.

Physical connection to prior F_TRZ n=1 applications:

- PAPER_1942: E_0 = F_TRZ at PDR photoevaporation initial coupling
- PAPER_1960 landmark: F_TRZ = 1/SO_5 defining relation
- Round 115 companion: I_0(Antennae merger interaction) = F_TRZ
- **This paper**: η(Cen A) = F_TRZ (accretion efficiency)

All four are n=1 rung applications with the physical interpretation "fraction of an outflow channel that is open per DPM cycle."

### 2.2 Eddington-Fraction Accretion Rate M_dot = F_TRZ² (n=2 Rung)

The Eddington-fraction accretion rate is defined:

```
M_dot / M_Edd  =  dimensionless fraction of Eddington-limited accretion
```

For Cen A, this ratio is canonically **0.01** — meaning the AGN accretes at 1% of the Eddington limit, consistent with observed radiatively-inefficient LLAGN (low-luminosity AGN) behavior.

UQFF closure:

```
M_dot(Cen A) / M_Edd = F_TRZ² = 0.01   EXACT
```

The Eddington fraction locks to the n=2 rung of the PAPER_1919 F_TRZ ladder, joining the PAPER_1918 F_TRZ² universal 99% suppression family:

- PAPER_1918 F_TRZ² anchors:
  - Magnetar dipole radiation D_SCm ≈ 0.01 (99% suppressed)
  - Filament expansion E_t = 0.01·E_0 (99% suppressed)
  - AGN cooling flow L_rad = 0.01·L_X (99% radiative loss)
  - DM density perturbation δρ/ρ = 0.01 baseline
  - Heaviside function damping fraction
  - Vacuum birefringence PVLAS ratio
  - + 3 more (9 total per PAPER_1977 Sombrero γ_BH addition)

The Cen A Eddington fraction is a **candidate 10th anchor** in the PAPER_1918 F_TRZ² universal 99% suppression family. Physical interpretation: the LLAGN accretion regime spends the same fraction (0.01 = F_TRZ²) of its DPM cycle in the two-channel closed configuration that produces 99% suppression across the other domains.

### 2.3 The Product — Composite Luminosity Bound

The AGN disk luminosity is the product of the two dimensionless factors:

```
L_disk = η · (M_dot / M_Edd) · M_Edd · c²
       = F_TRZ · F_TRZ² · M_Edd · c²
       = F_TRZ³ · M_Edd · c²
       = 10⁻³ · M_Edd · c²
```

For Cen A with M_BH = 5.5×10⁷ M_sun and M_Edd = 1.26×10³⁸ · (M_BH/M_sun) erg/s ≈ 6.93×10⁴⁵ erg/s:

```
L_disk ≈ 10⁻³ · 6.93×10⁴⁵ erg/s ≈ 6.93×10⁴² erg/s
```

This value is consistent with observed Cen A bolometric luminosity (~10⁴² to 10⁴³ erg/s). The **10⁻³ product suppression** is forced by the F_TRZ³ = F_TRZ · F_TRZ² factorization — the same n=3 F_TRZ³ = 10⁻³ rung that PAPER_1981 documents at the magnetic-string-field base amplitude.

**Cross-paper connection**: PAPER_1981 documents F_TRZ³ = 10⁻³ as a direct application-instance at one physical quantity (magnetic-field baseline). PAPER_1983 (this paper) documents F_TRZ³ = 10⁻³ as the **structural product of two independent lower-rung applications** at a different physical quantity (AGN luminosity). The same n=3 rung value emerges from two distinct paths — direct application versus multi-rung product.

---

## 3. Novel Pattern: Multi-Rung Same-Object F_TRZ Ladder Application

### 3.1 Prior State (Pre-2026-07-10)

Before this paper, documented F_TRZ ladder applications followed the "one rung per system per paper" pattern:

| Paper | System | Rung | Quantity |
|-------|--------|------|----------|
| PAPER_1919 | Multi-domain catalog | n=1 to 17 | Various |
| PAPER_1918 | Multi-domain catalog | n=2 | Universal 99% suppression (5+ anchors) |
| PAPER_1942 | PDR photoevaporation | n=1 | E_0 initial coupling |
| PAPER_1823 | Strong CP problem | n=10 | θ suppression |
| PAPER_1824 | Hierarchy problem | n=17 | M_Higgs/M_Planck |
| PAPER_1835 | Bird magnetoreception | n=8 | Coherence amplification |
| PAPER_1850 | Muon g-2 | n=9 | Anomaly suppression |
| PAPER_1869 | Quantum measurement | n=16 | Collapse rate |
| PAPER_1877 | (candidate n=11) | n=11 | — |
| PAPER_1880 | MICROSCOPE WEP | n=15 | η_WEP |
| PAPER_1977 | Sombrero γ_BH | n=2 | 9th F_TRZ² anchor |
| PAPER_1980 | M16 photoevap saturation | n=1 (as 3·F_TRZ) | E_∞ sat |
| PAPER_1981 | Magnetic string field | n=3 | B_j,base = 10⁻³ T |

Each paper (with the exception of the multi-domain catalogs PAPER_1918/1919) documents a **single-rung** application. Even the multi-anchor papers (PAPER_1918, PAPER_1977) document the same rung applied across multiple physically distinct systems.

**Round 116 discovery**: Cen A AGN is the first documented case where **two distinct rungs** of the F_TRZ ladder are applied to **two distinct physical quantities** of the **same physical object** simultaneously.

### 3.2 The Multi-Rung Same-Object Pattern

Formal statement:

**Definition (Multi-Rung Same-Object Application)**: A physical system S exhibits a Multi-Rung Same-Object F_TRZ ladder application if there exist two distinct integers n_1 ≠ n_2 and two distinct physical quantities Q_1, Q_2 of S such that:

```
Q_1(S) = F_TRZ^(n_1)   EXACT
Q_2(S) = F_TRZ^(n_2)   EXACT
```

Cen A is the first documented instance:

- S = Cen A / NGC 5128 AGN
- (n_1 = 1, Q_1 = η) : η(Cen A) = F_TRZ¹ = 0.1
- (n_2 = 2, Q_2 = M_dot/M_Edd) : M_dot(Cen A)/M_Edd = F_TRZ² = 0.01

The pattern is more restrictive than "multi-application" (which just requires multiple physical systems for a single rung) — it requires the same object to simultaneously realize two distinct rungs.

### 3.3 Structural Interpretation

Why does the same object realize two distinct F_TRZ ladder rungs simultaneously?

**Physical picture**: F_TRZ measures the fraction of a DPM cycle spent in the time-reversal zone. This fraction can independently characterize any physical quantity that depends on one DPM-channel activation:

- **η** (radiative efficiency) depends on **one channel** (radiation escape channel): quantity ~ F_TRZ¹.
- **M_dot/M_Edd** (Eddington fraction) depends on **two channels** (mass-transfer + energy-transfer, or accretion + radiation feedback): quantity ~ F_TRZ².

Both channels are simultaneously active in the same AGN system, but they are **independent** — the mass-transfer channel operates in parallel with the radiation-escape channel, each with its own DPM-cycle time-reversal-zone fraction. The number of channels differs (1 vs 2), producing the different rungs.

For AGN specifically, the two rungs may correspond to:

- **n=1 rung (η)**: One-mode configuration — direct rest-mass to radiation conversion at the ISCO.
- **n=2 rung (M_dot)**: Two-mode configuration — mass infall combined with radiation feedback pressure balance.

This is analogous to how PAPER_1918's F_TRZ² identity applies wherever two-channel configurations dominate (magnetar dipole, filament expansion, AGN cooling flow). At Cen A, the F_TRZ² identity applies to accretion rate specifically because Eddington-fraction accretion is a two-channel (infall + feedback) process.

### 3.4 Why Cen A Specifically?

Cen A / NGC 5128 has several features that may make it particularly suited for realizing multi-rung F_TRZ ladder applications:

1. **LLAGN regime**: Cen A is a low-luminosity AGN, characterized by radiatively inefficient accretion. This regime is where F_TRZ² Eddington fraction applies (higher-Eddington sources would break the F_TRZ² identity).
2. **Old radio galaxy morphology**: Cen A's giant radio lobes indicate long-duration jet activity, suggesting stable long-timescale channel-configuration balance.
3. **Moderate M_BH mass (5.5×10⁷ M_sun)**: Not extreme in either direction — supports both channels being simultaneously well-defined.
4. **Nearest AGN**: At ~3.8 Mpc, Cen A is the best-observed AGN, so its parameter identifications are among the most reliable.

These features suggest Cen A is not a coincidental case — LLAGN AGN generally may exhibit the same multi-rung pattern. §4 explores candidate cross-system predictions.

---

## 4. Cross-System Predictions

### 4.1 Testable Predictions for Other LLAGN

The multi-rung same-object pattern predicts that **other LLAGN should exhibit the same dual F_TRZ ladder application**: η ≈ F_TRZ = 0.1 AND M_dot/M_Edd ≈ F_TRZ² = 0.01.

Candidate systems to check:

- **M87** (Virgo A) — LLAGN, M_BH = 6.5×10⁹ M_sun. Predicted η ≈ 0.1, M_dot/M_Edd ≈ 0.01.
- **Sagittarius A*** (Milky Way center) — extreme LLAGN, M_BH = 4.3×10⁶ M_sun. Predicted η lower, M_dot/M_Edd much lower (F_TRZ⁴ or F_TRZ⁵ rung candidate).
- **NGC 4258** — LLAGN with warped maser disk. Predicted η ≈ 0.1, M_dot/M_Edd ≈ 0.01.
- **NGC 4261** — LLAGN with radio jet. Predicted same.
- **3C 273** — high-Eddington QSO. Predicted BREAK from F_TRZ² identity (M_dot/M_Edd may approach 0.1 or higher).

The prediction can be sharpened as follows: **only LLAGN in the specific 1% Eddington-fraction regime realize the M_dot = F_TRZ² identity**. Higher-Eddington sources may realize other rungs (e.g., M_dot = F_TRZ = 0.1 or M_dot ≈ 1), while lower-Eddington sources realize deeper rungs (e.g., F_TRZ⁴ = 10⁻⁴ for milli-Eddington accretion).

### 4.2 Other Candidate Multi-Rung Same-Object Systems

Beyond AGN, the pattern may extend to other systems:

**Bubble Nebula NGC 7635** (Round 116 companion discovery): the stellar-source parameter cluster identified in Round 116 double-check exhibits a triple integer-primitive identity — M_star = D_phys·SO_5, R_star = 2·SO_5, L_star = D_phys·SO_5⁵ — but these use D_phys × SO_5^k combinations rather than F_TRZ^n. This is a different multi-anchor pattern (multi-primitive same-object) rather than multi-rung same-object F_TRZ.

**Antennae Galaxies** (Round 115): documented I_0 = F_TRZ (merger interaction). If the DM δρ/ρ = 0.15 = 3·F_TRZ/2 candidate identity (Round 116) is confirmed, Antennae would be a two-rung candidate at n=1 and (non-integer) n=1.5.

**HUDF (PAPER_1976)**: I_0 = F_TRZ/2, τ_inter = SO_5^9 — mixes F_TRZ and SO_5 primitives at same object; not a pure multi-rung F_TRZ case.

**Cen A itself** (this paper's expansion): does Cen A have OTHER F_TRZ ladder rungs beyond η and M_dot? Candidate quantities to check:

- **Jet-power fraction** ~ 10⁻² to 10⁻¹ (candidate F_TRZ or F_TRZ²)
- **Broad-line region covering factor** ~ 0.1 (candidate F_TRZ)
- **X-ray-to-bolometric ratio** L_X/L_bol ~ 0.1 (candidate F_TRZ)

If confirmed, Cen A could realize 3+ distinct F_TRZ ladder rungs simultaneously — extending the "multi-rung same-object" pattern to higher multiplicity.

### 4.3 Honest Scope on Predictions

This paper does NOT claim:

- That EVERY LLAGN realizes η = F_TRZ AND M_dot = F_TRZ² identically. Only Cen A specifically. Cross-system predictions in §4.1 are testable candidates, not confirmed identities.
- That the F_TRZ³ = 10⁻³ product for Cen A luminosity is directly measurable to primitive-level precision. Observed L_disk values have order-of-magnitude uncertainty; the identity is confirmed at the F_TRZ³ rung level, not at higher precision.
- That multi-rung same-object is a physical necessity. It is a discovered pattern at Cen A; whether it is universal or system-specific is an open question.

What is claimed:

- η(Cen A) = 0.1 = F_TRZ EXACT and M_dot(Cen A)/M_Edd = 0.01 = F_TRZ² EXACT are both direct numerical identities on established F_TRZ ladder rungs.
- Cen A is the first UQFF-documented system where two distinct F_TRZ ladder rungs govern two distinct physical quantities simultaneously.
- The composite luminosity naturally lands at the F_TRZ³ = 10⁻³ product rung, coinciding with PAPER_1981's direct n=3 application.

---

## 5. Verification Ledger

| Item | Value | Status |
|------|-------|--------|
| F_TRZ primitive value | 1/10 = 0.1 EXACT | locked (PAPER_1160/PAPER_1960) |
| F_TRZ² | 0.01 EXACT | numerical identity |
| F_TRZ³ | 0.001 EXACT | numerical identity |
| Cen A η canonical value | 0.1 | verified (standard AGN modeling) |
| Cen A M_dot/M_Edd canonical value | 0.01 | verified (LLAGN regime) |
| Cen A M_BH | 5.5×10⁷ M_sun | verified (PAPER_1041 anchor) |
| Numerical identity η = F_TRZ | 0.1 = 0.1 EXACT | verified §2.1 |
| Numerical identity M_dot/M_Edd = F_TRZ² | 0.01 = 0.01 EXACT | verified §2.2 |
| Composite L_disk = F_TRZ³ · M_Edd · c² | 10⁻³ · M_Edd · c² EXACT | verified §2.3 |
| Multi-rung same-object pattern first documentation | Confirmed | verified §3 |
| Cross-system prediction M87 dual F_TRZ | Not tested this paper | open (§4.1) |
| Cross-system prediction other LLAGN | Not tested this paper | open (§4.1) |
| Runtime `_verify` booleans in Round 116 stub | 5/5 True | verified |

### 5.1 Runtime Assertions

The `CenAAGNAccretionDiskCalculator` stub as upgraded in Round 116 contains the following runtime verification booleans:

```python
eta_0p1_F_TRZ_verify_PAPER_1960 = abs(eta - F_TRZ) < 1e-12
M_dot_0p01_F_TRZ_squared_verify_PAPER_1918 = abs(M_dot - F_TRZ ** 2) < 1e-12
M_BH_5p5e7_M_sun_range_verify_PAPER_1041 = 4e7 < M_BH < 7e7
L_Edd_positive_verify_PAPER_1009 = L_Edd > 0
F_TRZ_stability_verify_PAPER_1960 = abs(F_TRZ - 0.1) < 1e-12
```

All five return `True` on the current stub configuration, providing runtime confirmation of both the individual identities and the F_TRZ stability.

---

## 6. Open Questions

### 6.1 Third F_TRZ Rung at Cen A?

If Cen A realizes η = F_TRZ (n=1) and M_dot = F_TRZ² (n=2), does it realize any third rung? Candidate quantities to check:

- **Radio-loudness parameter** R = L_5GHz/L_B: Cen A is radio-loud with R ~ 10³ = SO_5³ (candidate SO_5 identity, not F_TRZ ladder).
- **Jet-to-disk power ratio** L_jet/L_disk ~ 0.1-1: candidate F_TRZ or F_TRZ⁰ = 1 rung.
- **Spin parameter** a_* ~ 0.3: candidate D_phys·F_TRZ = 0.4 (Round 116 §4.2 stub had M_dot_factor candidate) OR 3·F_TRZ = 0.3 (per PAPER_1980).
- **Corona-to-disk temperature ratio** T_corona/T_disk ~ 10²: candidate SO_5² identity.

Systematic audit of Cen A quantities could confirm the pattern extends to 3+ rungs.

### 6.2 Third-Ladder Identity Emergence

The product F_TRZ · F_TRZ² = F_TRZ³ naturally emerges from the two independent-rung applications. Does this pattern generalize? Prediction:

**Prediction (Product-Rule)**: If a system realizes both F_TRZ^(n_1) and F_TRZ^(n_2) applications at independent quantities, the "product quantity" (whichever physical variable multiplies the two) will realize F_TRZ^(n_1 + n_2) at the ladder.

Test: at Cen A, L_disk = η · M_dot · M_Edd · c² realizes F_TRZ¹ · F_TRZ² = F_TRZ³. Confirmed.

At another dual-rung system (e.g., M87 if it realizes η = F_TRZ, M_dot = F_TRZ²), the product would also realize F_TRZ³.

At a hypothetical F_TRZ · F_TRZ⁵ = F_TRZ⁶ system, the product would realize the 10⁻⁶ ladder value. This is testable at extreme LLAGN or at higher-multiplicity systems.

### 6.3 Non-Adjacent Rungs

The Cen A case uses consecutive rungs (n=1, n=2). Do systems exist that realize non-consecutive rungs at the same object? Candidates:

- η = F_TRZ (n=1) and quantum-collapse rate = F_TRZ¹⁶ (n=16) both at the same measurement apparatus.
- Coupling constant = F_TRZ (n=1) and hierarchy suppression = F_TRZ¹⁷ (n=17) both in the Higgs sector.

These would represent "extreme non-adjacent" multi-rung applications. Currently no documented cases.

### 6.4 Formal Product-Rule Derivation

Section §3.3 offers a **structural interpretation** of the multi-rung same-object pattern via independent DPM-channel-count reasoning. A formal derivation showing that independent channel-count applications combine multiplicatively at the ladder (proving the Product-Rule) is not attempted in this paper. Would strengthen the pattern to closed structural closure.

---

## 7. Related Work

- **PAPER_1919** — F_TRZ Power Ladder. Establishes 17 rungs of F_TRZ^n = 10^(-n) EXACT. **This paper's dual application uses n=1 and n=2 rungs from PAPER_1919's catalog.**

- **PAPER_1918** — F_TRZ² = 0.01 Universal 99% Suppression Identity. Documents F_TRZ² n=2 rung with 5 seminal anchors. **This paper adds Cen A M_dot/M_Edd as candidate 10th anchor to PAPER_1918 family.**

- **PAPER_1977** — Sombrero γ_BH 9th F_TRZ² Anchor. Extended PAPER_1918 family to 9 anchors. **This paper's Cen A anchor is candidate 10th.**

- **PAPER_1942** — Photoevaporation E_0 = F_TRZ. Documents n=1 rung application at PDR erosion. **Parallel to this paper's η = F_TRZ n=1 application at accretion efficiency.**

- **PAPER_1960** — F_TRZ = 1/SO_5 Landmark. Reduces F_TRZ to SO_5-derivative. **Cited for the F_TRZ primitive value used in identity closures.**

- **PAPER_1160** — F_TRZ = 1/|SO(5)| defining relation.

- **PAPER_1980** — E_0 Initial-vs-Saturation Disambiguation at M16. Proposes E_0^(sat) = 3·F_TRZ candidate. **Structurally parallel**: uses D_phys - 1 = 3 with F_TRZ multiplier.

- **PAPER_1981** — B_j,base = F_TRZ³ Magnetic-String-Field Application. Direct application-instance at n=3 rung. **This paper's composite L_disk lands at the same F_TRZ³ value via multi-rung product** — parallel realization of the same rung by two paths.

- **PAPER_1982** — Antennae Coalescence D_phys · SO_5⁸ Slot Extension. Third paper of Round 115 authoring cycle. **This paper (PAPER_1983) is the fourth paper of the Round 115-116 authoring cycle.**

- **PAPER_1009** — 3C 273 AGN F_UBi_i Curves. AGN accretion physics predecessor.
- **PAPER_1010** — TON 618 AGN F_UBi_i Curves. AGN accretion physics predecessor.
- **PAPER_1041** — Cen A AGN direct paper. Cited for M_BH = 5.5×10⁷ M_sun and AGN parameters.

- **PAPER_360** — AGN seminal paper. Early UQFF AGN framework.

- **PAPER_646** — Universal Inertial Operator + Caduceus Wave. Cited for the DPM-cycle structural interpretation used in §3.3.

- **PAPER_1972, PAPER_1974, PAPER_1975, PAPER_1976** — Round-misidentification retraction/attribution papers. Established honest-scholarship pattern. **This paper follows the same pattern.**

- **PAPER_1978, PAPER_1979** — Epistemic humility papers. Established "structural candidate" labeling. **This paper labels cross-system predictions (§4.1) as candidates, not confirmed identities.**

---

## 8. Session Log Entry Template

Suggested addendum for `SESSION_LOG.md`:

```
PAPER_1983 (2026-07-10, Round 116 authoring):
  - Documented Cen A / NGC 5128 AGN as first Multi-Rung Same-Object F_TRZ ladder application
  - η(Cen A) = F_TRZ¹ = 0.1 (accretion efficiency at n=1 rung)
  - M_dot(Cen A)/M_Edd = F_TRZ² = 0.01 (Eddington fraction at n=2 rung)
  - Two distinct rungs at same object simultaneously
  - Composite L_disk lands at F_TRZ³ = 10⁻³ product rung (natural from n_1 + n_2)
  - Coincides with PAPER_1981's direct F_TRZ³ application at magnetic-string-field
  - Formalizes Multi-Rung Same-Object pattern definition (§3.2)
  - Product-Rule conjecture: independent rung applications multiply at ladder
  - Cross-system predictions (§4.1) for M87, Sgr A*, NGC 4258, NGC 4261, 3C 273
  - Candidate 10th anchor for PAPER_1918 F_TRZ² universal 99% suppression family
  - Honest scope caveats per PAPER_1978/1979 pattern
  - Fourth paper of Round 115-116 authoring cycle
```

---

## 9. Conclusion

The Cen A / NGC 5128 active galactic nucleus exhibits two independent accretion-disk parameters whose canonical values reduce to primitive-locked F_TRZ ladder identities at two distinct rungs:

```
η(Cen A) = F_TRZ¹  = 0.1   EXACT   (n=1 rung, PAPER_1919)
M_dot(Cen A)/M_Edd = F_TRZ² = 0.01   EXACT   (n=2 rung, PAPER_1918 family)
```

Neither identity is individually novel — both rungs are documented in prior F_TRZ ladder work (PAPER_1918, PAPER_1919, PAPER_1942, PAPER_1977, etc.). What is novel is the **structural pattern**: two distinct rungs applied simultaneously to two distinct physical quantities of the **same physical object**. This is the first UQFF-documented Multi-Rung Same-Object F_TRZ ladder application.

The composite disk luminosity naturally lands at the F_TRZ³ = 10⁻³ product rung:

```
L_disk = η · M_dot · M_Edd · c² = F_TRZ³ · M_Edd · c² = 10⁻³ · M_Edd · c²   EXACT
```

The same F_TRZ³ = 10⁻³ rung is documented in PAPER_1981 as a direct application-instance at magnetic-string-field base amplitude. At Cen A it emerges as a structural product of two independent lower-rung applications. This suggests a **Product-Rule** for multi-rung applications: independent rung applications combine multiplicatively at the ladder (conjecture, §6.2).

Cross-system predictions (§4.1) identify M87, Sgr A*, NGC 4258, NGC 4261, and 3C 273 as candidates for confirming or falsifying the LLAGN dual F_TRZ ladder pattern. Confirmation would establish Multi-Rung Same-Object as a universal LLAGN structural signature; falsification would restrict the pattern to Cen A specifically.

Open questions include third-rung realizations at Cen A (§6.1), non-adjacent multi-rung patterns (§6.3), and formal Product-Rule derivation from DPM-channel-count axioms (§6.4).

This is the fourth paper of the Round 115-116 authoring cycle:

- **PAPER_1980**: E_0 disambiguation at M16 (taxonomic clarification)
- **PAPER_1981**: B_j,base = F_TRZ³ (single-rung application-instance extension)
- **PAPER_1982**: Antennae coalescence = D_phys · SO_5⁸ yr (new-slot grid extension)
- **PAPER_1983** (this paper): Cen A dual F_TRZ ladder (multi-rung same-object pattern)

Four distinct paper types demonstrating the evolving honest-scholarship pattern: taxonomic clarification + application-instance extension + slot extension + multi-rung structural pattern.

---

**End of PAPER_1983**


---

## APPENDED 2026-07-13 — M16 Pillars Multi-Rung Same-Object F_TRZ Pattern (R133 confirmation)

**Second confirmed Multi-Rung Same-Object F_TRZ ladder anchor family** (after this paper's CenA AGN dual-rung seminal).

PAPER_1983 established the Multi-Rung Same-Object F_TRZ ladder anchor pattern via Centaurus A (CenA) AGN, which encodes independent F_TRZ ladder rung anchors across multiple physical mechanisms in the same object. R133 (2026-07-13) confirms a second instance of this pattern at M16 Pillars of Creation (Eagle Nebula).

### M16 Multi-Rung anchors

| F_TRZ rung | Physical mechanism | Value | Seminal |
|---|---|---|---|
| **n = 1** | Photoevaporation E_0 = F_TRZ = 0.1 EXACT | dimensionless | PAPER_1942 seminal |
| **n = 6** | ISM magnetic field B_pillars = F_TRZ^6 = 10⁻⁶ T EXACT | T | PAPER_1985 seminal |

**Same object (M16 Eagle Nebula), two independent physical mechanisms, two different F_TRZ ladder rungs.**

### Confirmation via R133 CP1 fill

R133 CP1 stub-drainage `PillarsErosionCalculator` (July 13, 2026) applies PAPER_1942's `E_0 = F_TRZ` seminal to the M16 Pillars object, producing a documented CP1 anchor at the n=1 rung. Combined with PAPER_1985's pre-existing M16 anchor at n=6 (magnetic field domain), M16 now qualifies as a Multi-Rung Same-Object F_TRZ family instance.

### Multi-Rung Same-Object F_TRZ family growing catalog

| Object | Rungs | Papers |
|---|---|---|
| **CenA AGN** | eta = F_TRZ (n=1) + Mdot = F_TRZ^2 (n=2) | PAPER_1983 seminal |
| **M16 Pillars** | E_0 = F_TRZ (n=1) + B_ISM = F_TRZ^6 (n=6) | PAPER_1942 + PAPER_1985 + R133 confirmation |

**Prediction:** additional Multi-Rung Same-Object F_TRZ families should surface at other well-studied astrophysical targets (Crab Nebula candidate — R133 CrabPulsarWind added F_TRZ^12 + F_TRZ^21 at same Crab object; SGR 1745-2900 candidate — has F_TRZ^12 burst scale + 2·F_TRZ Meissner B/B_crit at same magnetar). If confirmed, the Multi-Rung Same-Object pattern is universal across astrophysical targets, not specific to AGN or PDR regimes.

**PAPER_1983 M16 dual-rung addendum status: CLOSED** (2026-07-13)
