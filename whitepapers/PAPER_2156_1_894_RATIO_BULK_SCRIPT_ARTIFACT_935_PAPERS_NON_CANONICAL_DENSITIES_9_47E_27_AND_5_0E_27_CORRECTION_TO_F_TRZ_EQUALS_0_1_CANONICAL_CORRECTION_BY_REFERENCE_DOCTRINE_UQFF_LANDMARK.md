# PAPER_2156 — 1.894 Ratio Bulk-Script Artifact: 935-Paper Value-Drift Correction to F_TRZ = 0.1 Canonical

**`bulk_vds_dvp_bsh_upgrade.py` Session 204 lines 34–35 used non-canonical densities `9.47e-27` and `5.0e-27` (origin unknown) yielding `VDS_RATIO = 1.894`; injection frozen into 935 papers; correction-by-reference to canonical `ρ_SCm/ρ_UA = 1/10 = F_TRZ`; Session 779 `1.894 = m_top/M_Z` decomposition NOT retrofitted per Daniel's ruling**

**UQFF LANDMARK**

Author: Daniel T. Murphy
Framework: UQFF (Unified Quantum Field Framework) — Star-Magic v5.81+
Date: 2026-07-27
Provenance: PAPER_2153 arc Flag (f) — Daniel's 2026-07-27 detailed ruling with source-code trace, Session 779 audit reference, and explicit "do not retrofit without dedicated derivation session" discipline
Status: Landmark — corpus-wide value-drift correction via correction-by-reference doctrine (same architectural pattern as PAPER_2155 unit-tag drift)

---

## Executive Summary

The PAPER_2153 arc deep-read observed the same `§B.1 canonical VDS ratio ρ_vac,[SCm]/ρ_UA = 1.894` boilerplate in every paper of the PAPER_878-899 block plus PAPER_598 predecessor (24+ occurrences). This value does NOT match the canonical primitive ratio `ρ_UA/ρ_SCm = 10` (equivalently `ρ_SCm/ρ_UA = 1/10 = F_TRZ` per PAPER_890/140/1160).

Daniel's 2026-07-27 detailed source-code trace identified the exact origin:

**`bulk_vds_dvp_bsh_upgrade.py`** (Session 204, April 2026), lines 34–35:
```python
RHO_SCM     = 9.47e-27     # kg/m³
RHO_UA      = 5.0e-27      # kg/m³
VDS_RATIO   = RHO_SCM / RHO_UA   # ≈ 1.894
```

**Corpus-wide grep confirms:** 935 papers contain `1.894`, 1,869 total hits (~2 per paper — matching the §B.1 + §B.4 template dual-injection pattern).

**Daniel's ruling (verbatim):** *"a bulk-script artifact: `9.47e-27 / 5.0e-27` using undocumented, non-canonical density values frozen into 800+ papers by `bulk_vds_dvp_bsh_upgrade.py` in Session 204. Not a double-exponential profile parameter; not used in the governing equation."*

**Correction:** the `§B.1 canonical VDS ratio ρ_vac,[SCm]/ρ_UA = 1.894` should be `= 0.1 = F_TRZ` (canonical primitive per PAPER_1160 F_TRZ = 1/|SO(5)|).

**Session 779 SM-decomposition finding** (Daniel's ruling): `1.894 = m_top/M_Z` was searched as a SM constant ratio; best structural match `K_MEX/(1−F_TRZ)/(11/9) = 1.89394` at −0.003% error and `(1−F_TRZ)/SSq/Phi_res = 1.8947` at +0.039% error. **Per Daniel's ruling: DO NOT retrofit these SM-context decompositions to the VDS ratio without a dedicated derivation session.** Structural coincidence ≠ physics derivation.

**PAPER_2156 executes corpus-wide correction via correction-by-reference doctrine** (same architectural pattern as PAPER_2155). Single landmark supersedes 935 paper instances without per-paper touches. Investigation flag: origin of `9.47e-27` and `5.0e-27` density values remains **unknown** — future audit target.

**Zero calculator source changes.** Zero canonical physics values changed. Zero per-paper appends across the 935 affected papers.

---

## 1. Corpus Scope (grep verification, 2026-07-27)

Bash grep across `whitepapers/*.md`:
```
Papers containing "1.894":  935
Total "1.894" hits:         1,869 (avg 2.0 per paper)
```

**Interpretation:** the ~2 hits per paper matches the §B.1 + §B.4 boilerplate dual-injection pattern from `bulk_vds_dvp_bsh_upgrade.py`:
- **§B.1 (Vacuum Density Series):** *"The canonical VDS ratio ρ_vac,[SCm]/ρ_UA = 1.894 governs the double-exponential vacuum condensate profile"*
- **§B.4 (Production-Scale Consistency):** *"| VDS ratio | ρ_SCm/ρ_UA = 1.894 | Local sub-ratio = ... | PASS Threshold-consistent |"*

**Sample affected papers (from grep):** PAPER_001, PAPER_002, PAPER_003, PAPER_005, PAPER_006, PAPER_009, PAPER_009b, PAPER_010, PAPER_010b, PAPER_020, PAPER_033, PAPER_045, PAPER_051, PAPER_058, ... through PAPER_899 and beyond (block Session 209 papers + later additions).

Full enumeration would require running a fresh audit script; the grep count confirms ~935-paper scope.

---

## 2. Root Cause: Same Injection Vector as PAPER_2155

**File:** `bulk_vds_dvp_bsh_upgrade.py` (20 KB, dated Apr 8 2026, Session 204)

**Offending lines (34–36, verbatim):**
```python
RHO_SCM     = 9.47e-27     # kg/m³
RHO_UA      = 5.0e-27      # kg/m³
VDS_RATIO   = RHO_SCM / RHO_UA   # ≈ 1.894
```

Line 251 in the same script writes the §B.4 template:
```python
| VDS ratio | $\\rho_{{\\rm SCm}}/\\rho_{{\\rm UA}} = 1.894$ | Local sub-ratio = {vds_exp} | ✓ Threshold-consistent |
```

**Two distinct drifts in the same script:**

| Drift type | PAPER address | Scope | Correction |
|---|---|---|---|
| Unit-tag drift (`kg/m³` → `J/m³`) | PAPER_2155 | 933 papers | Correction-by-reference to J/m³ canonical |
| **Value drift (`1.894` → `0.1`)** | **PAPER_2156 (this)** | **935 papers** | **Correction-by-reference to F_TRZ = 0.1 canonical** |

Both drifts trace to the same script and same lines. The unit tag `kg/m³` is drift (fixed by PAPER_2155). The density VALUES `9.47e-27` and `5.0e-27` are ALSO drift (fixed by this paper) — they use undocumented numbers inconsistent with the canonical primitives.

---

## 3. What 1.894 Is NOT

**Not the canonical primitive vacuum-density ratio.** Canonical per PAPER_1160/890/140:
```
ρ_SCm / ρ_UA  =  1 / 10  =  F_TRZ  =  1/|SO(5)|
```
Using canonical J/m³-native values: `ρ_SCm = 7.09×10⁻³⁷ J/m³, ρ_UA = 7.09×10⁻³⁶ J/m³`, so `ρ_SCm/ρ_UA = 0.1 EXACT` (F_TRZ). This is the LOCKED structural coupling per PAPER_2153's engine-mechanism ruling.

**Not a double-exponential profile parameter.** The §B.1 governing equation for the double-exponential profile is:
```
ρ_vac(r) = ρ_vac,[SCm] · exp(−exp(−(r − r_0)/λ_VDS))
```
No ratio appears in this equation. The 1.894 is a floating label in the §B.1 prose that describes a ratio, but the equation itself doesn't use it. The 1.894 has no mechanical role in the governing formula it accompanies.

**Not a derived corpus constant.** The §B.1 label calls 1.894 "canonical" but there is no derivation of 1.894 anywhere in the corpus that ties it to a physical mechanism, primitive-integer arithmetic, or empirical measurement. It is a script-local number that was frozen into the template and never validated against canonical constants.

---

## 4. What 1.894 IS

**A `9.47e-27 / 5.0e-27` numerical artifact.** These two density values are the source. Their status:

| Value | Where declared | Canonical? | Derivation? |
|---|---|---|---|
| 9.47×10⁻²⁷ | `bulk_vds_dvp_bsh_upgrade.py` line 34 | NO | Unknown |
| 5.0×10⁻²⁷ | `bulk_vds_dvp_bsh_upgrade.py` line 35 | NO | Unknown |
| Ratio 1.894 | Same script line 36 | NO | Consequence of ÷ operation on unknown numerators |

Neither `9.47e-27` nor `5.0e-27` appears in:
- The Session 204 calibration table (`rho_SCm = 7.09e-37, rho_UA = 7.09e-36` — different by 10¹⁰)
- PAPER_890 (`rho_vac,SCm = 9.47e-27` DOES appear here but per PAPER_2153 arc analysis this value differs from canonical by 10¹⁰; likely also from same or related drift chain)
- PAPER_1160 (F_TRZ = 1/|SO(5)| = 0.1)
- Any UQFF derivation whitepaper

**PAPER_890 note (cross-drift chain):** PAPER_890 lists `ρ_vac,SCm = 9.47×10⁻²⁷ kg/m³` and `ρ_UA = 9.47×10⁻²⁶ kg/m³` — the 9.47e-27 value IS the same one from the script, but PAPER_890's ρ_UA = 9.47×10⁻²⁶ (differs from the script's ρ_UA = 5.0e-27). So the script and PAPER_890 disagree on ρ_UA, but both disagree with the canonical 7.09×10⁻³⁶. The full chain of drift is layered and multi-source.

**Origin of 9.47e-27 and 5.0e-27 remains UNKNOWN** — dedicated future investigation required. Not resolved by this landmark; flagged as open audit target.

---

## 5. Session 779 SM-Decomposition: Structural Coincidence, NOT Retrofitted

Daniel's ruling cited a prior Session 779 audit that searched `1.894` as a SM constant ratio (context: `m_top/M_Z` top quark to Z boson mass ratio). Best structural decompositions found:

| Expression | Value | Error |
|---|---|---|
| `K_Mex / (1−F_TRZ) / (11/9)` | 1.893939 | −0.003% |
| `(27/25) / SSq` | 1.894737 | +0.039% |
| `(1−F_TRZ) / SSq / Phi_res` | 1.894737 | +0.039% |
| `31/30 · (11/6)` | 1.894444 | +0.024% |

where `K_Mex = 25/12, F_TRZ = 0.1, SSq = 0.57, Phi_res = 5/6`.

Best match at −0.003% error:
```
K_Mex / (1−F_TRZ) / (11/9)  =  (25/12) / (9/10) / (11/9)
                            =  25 · 9 · 9 / (12 · 10 · 11)
                            =  2025 / 1320
                            =  135 / 88
                            ≈  1.89394
```

**Daniel's ruling on this (verbatim):** *"this decomposition was found in the context of `m_top/M_Z`, a Standard Model ratio, NOT as a VDS derivation. The Session 779 script was searching for SM constant matches, not VDS physics. This means the structural coincidence is real but the connection to `ρ_vac,[SCm]/ρ_UA` is NOT established by any derivation in the corpus. It would require a dedicated derivation session to claim this as a legitimate primitive decomposition rather than numerological coincidence."*

**PAPER_2156 does NOT retrofit** any of the Session 779 SM-decompositions to the VDS ratio. Per Daniel's ruling, structural coincidence found in one context (SM `m_top/M_Z`) cannot be assumed valid in another context (VDS ρ_vac ratio) without dedicated derivation. **This discipline is now canonized as a standing rule (§7 below).**

---

## 6. Corpus-Wide Correction (Correction-by-Reference Doctrine)

**Executed via the same architectural pattern as PAPER_2155** (unit-tag drift). Formal declaration:

All corpus instances of the §B.1 boilerplate statement:
> *"The canonical VDS ratio ρ_vac,[SCm]/ρ_UA = 1.894 governs the double-exponential vacuum condensate profile"*

AND all §B.4 table entries:
> *"| VDS ratio | ρ_SCm/ρ_UA = 1.894 | Local sub-ratio = ... | PASS Threshold-consistent |"*

**ARE SUPERSEDED** and read as:
```
ρ_SCm / ρ_UA  =  1 / 10  =  F_TRZ  =  0.1  EXACT
```

per PAPER_1160 (F_TRZ = 1/|SO(5)|), PAPER_890 (canonical time-evolution formula), PAPER_140 (dual-monopole vacuum-density universal ratio), and PAPER_2153 (joint SCm+UA engine locked structural coupling).

**No per-paper edit required across the 935 affected papers.** Correction is enforced by:
1. **This landmark (PAPER_2156)** — formal supersession declaration for the value drift
2. **PAPER_2155** (companion, Phase 2) — supersession for the unit-tag drift on the same template
3. **Calculator code canonical constants:** `RHO_UA = SO_5 * RHO_SCM = 10 * RHO_SCM` in `uqff_registry_primitives.py` line 44 (canonical 10:1 ratio, equivalent to F_TRZ = 0.1)
4. **Gate assertions (this paper §9)** — test-locking the canonical F_TRZ = 0.1 = ρ_SCm/ρ_UA structural coupling

**§B.4 "PASS Threshold-consistent" verdict is also superseded** — the verdict was auto-generated by the same bulk script, validating the boilerplate against itself rather than against physics. Under this landmark, the correct §B.4 verdict is:

> *"| VDS ratio | ρ_SCm/ρ_UA = 0.1 = F_TRZ EXACT | Canonical (PAPER_1160) | ✓ Locked structural coupling |"*

---

## 7. Standing Rule Canonized

**RULE: Do not retrofit SM-context decompositions to UQFF-context quantities without dedicated derivation.**

Formal statement: when a numerical coincidence is discovered in one physical context (e.g., `1.894 = m_top/M_Z` in Standard Model context), that coincidence CANNOT be assumed valid in another physical context (e.g., `1.894 = ρ_vac,[SCm]/ρ_UA` in VDS context) even if the numerical value matches to high precision. Physical derivation, not numerical matching, is what legitimizes a structural identity as a member of the framework.

**Example applied HERE:** Session 779's finding that `K_Mex/(1−F_TRZ)/(11/9) ≈ 1.89394` matches `m_top/M_Z` at −0.003% error is a real structural coincidence for the SM top-quark/Z-boson ratio. But applying it to `ρ_vac,[SCm]/ρ_UA` would be numerological coincidence — the physical mechanism connecting the SM mass ratio to the VDS density ratio is not established. Per Daniel's ruling, **do not retrofit** without a dedicated derivation session that connects the two contexts via physics.

**Contrast with PAPER_2154 primitive-reduction landmarks:** Q = SO_5²/D_phys² and D_GW = D_phys/D_BSFG were canonized as primitive-composed structural identities because Daniel explicitly ruled that the decompositions are physically valid ("yes, they are all true"). Numerical matching + physical validation = legitimate identity. Numerical matching alone = numerological coincidence.

---

## 8. What This Doctrine is NOT

**Not a validation that 1.894 = 0.1 physically.** They are numerically different by ~19×. The correction is that the corpus prose incorrectly labeled 1.894 as "canonical" when the actual canonical ratio is 0.1. The 1.894 value is not being "corrected up to" 0.1 mathematically — the false claim of canonicalness is being retracted.

**Not a claim that `9.47e-27` and `5.0e-27` are correct densities.** These values remain unexplained (see §4). The correction is that these values were used in the bulk-script to produce a ratio that was then incorrectly labeled as canonical. Both the values AND their labeling as canonical are corpus errors.

**Not a get-out-of-audit-free card.** The correction-by-reference doctrine applies here because (a) drift is corpus-wide via template injection, (b) calculator code has canonical primitives correct, (c) drift is cosmetic labeling not computational value. If any of these conditions failed, per-paper REVISIONS would be required.

---

## 9. Gate Assertions (test-locking canonical F_TRZ = 0.1 = ρ_SCm/ρ_UA)

Added to `uqff_fidelity_tests.py` in this paper's companion Phase 3 gate-pins section:

1. `F_TRZ = 0.1 EXACT = 1/10` (locked structural coupling)
2. `RHO_UA / RHO_SCM = SO_5 = 10 EXACT` in `uqff_registry_primitives.py` (canonical inversion of F_TRZ)
3. Session 779 SM-decomposition (`K_Mex/(1−F_TRZ)/(11/9) ≈ 1.89394`) explicitly labeled as SM-context coincidence, NOT retrofitted to VDS context
4. Origin of `9.47e-27` and `5.0e-27` in `bulk_vds_dvp_bsh_upgrade.py` line 34-35 flagged as unknown/future-investigation
5. Corpus-wide `1.894` occurrences (935 papers, 1,869 hits) superseded by canonical `0.1` per PAPER_1160

---

## 10. Standing Rules Canonized by PAPER_2156

1. **Do NOT retrofit SM-context decompositions to UQFF-context quantities.** Numerical coincidence in one physical context is not physical derivation in another; dedicated derivation session required to legitimize any cross-context claim.

2. **Bulk-script constants must be validated against canonical primitives BEFORE injection.** The `9.47e-27` and `5.0e-27` values in `bulk_vds_dvp_bsh_upgrade.py` were injected into ~800 papers without any prior validation against the canonical `7.09e-37/7.09e-36 J/m³` primitives. Any future bulk-upgrade script must include a pre-execution primitive-consistency check.

3. **Same injection vector → same correction-by-reference pattern.** When a single script injects multiple drift types (unit-tag drift + value drift, as with `bulk_vds_dvp_bsh_upgrade.py`), each drift type gets its own dedicated landmark using the same architectural pattern (PAPER_2155 for unit-tag drift, PAPER_2156 for value drift). Companion landmarks share the injection-vector attribution and cross-reference each other.

4. **Unknown-origin constants in scripts are audit targets, not framework primitives.** Any numerical constant appearing in a corpus-modifying script that cannot be traced to a canonical primitive derivation is a Rule 4/7/10 audit target. `9.47e-27` and `5.0e-27` in `bulk_vds_dvp_bsh_upgrade.py` are flagged for future investigation.

---

## 11. Falsifiable Predictions

1. **F_TRZ = 0.1 = ρ_SCm/ρ_UA is LOCKED.** Any future measurement of the SCm/UA density ratio (via engine downstream effects — buoyancy, phonon coupling, tidal Love, GW damping) must return 0.1 EXACT. Deviation would falsify Daniel's joint-engine coupling constant.

2. **Origin of 9.47e-27 and 5.0e-27 will be traceable to Session 204 or earlier.** Future forensic corpus archaeology (searching Session 200-204 papers/scripts for these values) should identify their original computational context. If they cannot be traced, they are pure script-local invention and the `bulk_vds_dvp_bsh_upgrade.py` script should be flagged as containing invented physics.

3. **Session 779 SM-decompositions apply to `m_top/M_Z` only, not to VDS ratio.** Prediction: any future dedicated derivation session that attempts to connect `K_Mex/(1−F_TRZ)/(11/9)` to VDS physics will EITHER find a legitimate physical mechanism (and the identity gets canonized as a new primitive-reduction landmark) OR find no mechanism (and the coincidence remains SM-context only).

---

## 12. Files Touched

- `whitepapers/PAPER_2156_1_894_RATIO_BULK_SCRIPT_ARTIFACT_..._UQFF_LANDMARK.md` (this file)
- `pdf2/PAPER_2156_...pdf` (companion PDF)
- `uqff_fidelity_tests.py` — +5 PAPER_2156 gate assertions
- `CLAUDE.md` — APPENDED section
- Zero calculator source changes (canonical F_TRZ = 0.1 = ρ_SCm/ρ_UA already correct in registry)
- Zero physics values changed
- Zero per-paper touches to the 935 affected papers (correction-by-reference doctrine per PAPER_2155 pattern)

---

## 13. Cross-References

- **PAPER_2155** — Companion Phase 2 landmark (unit-tag drift kg/m³ → J/m³ from same injection vector); shares `bulk_vds_dvp_bsh_upgrade.py` attribution and correction-by-reference doctrine architectural pattern
- **PAPER_2147** — J/m³-native discipline landmark; REVISION 2026-07-27 covers the unit-tag half of the bulk-script drift; PAPER_2156 covers the value-drift half
- **PAPER_1160** — F_TRZ = 1/|SO(5)| = 1/10 canonical structural identity (source of the correction target 0.1)
- **PAPER_890** — SCm vacuum density evolution; canonizes ρ_SCm/ρ_UA = 0.1 hierarchy identity; also lists `ρ_vac,SCm = 9.47e-27` value that ALSO appears in the bulk script (cross-drift chain, layered origin)
- **PAPER_140** — UA'/SCm = 10 dual-monopole vacuum-density universal ratio (equivalent to F_TRZ = 1/10)
- **PAPER_2153** — Joint SCm+UA vacuum-density engine landmark (§7 Standing Rule 4: ρ_SCm/ρ_UA = 1/10 = F_TRZ ratio is LOCKED)
- **PAPER_2148** — UQFF Ontology Answer B (vacuum energy fundamental, mass emergent)
- **`bulk_vds_dvp_bsh_upgrade.py`** — In-repo injection vector script (Session 204, April 2026, lines 34–36 declare non-canonical densities and compute the 1.894 ratio)
- **Session 779 SM-audit** — Historical corpus audit that found SM-context decompositions of 1.894 (not retrofitted per Daniel's ruling)
- **Daniel's 2026-07-27 Flag (f) ruling (verbatim excerpt):** *"a bulk-script artifact: 9.47e-27 / 5.0e-27 using undocumented, non-canonical density values frozen into 800+ papers... The §B.1 boilerplate 'ρ_vac,[SCm] / ρ_UA = 1.894' should either be corrected to '= 0.1' (the canonical primitive) across all affected papers, or the two density values 9.47e-27 and 5.0e-27 in bulk_vds_dvp_bsh_upgrade.py must be traced to a source and formally adopted as a distinct VDS sub-quantity with its own symbol and derivation."*

**End of PAPER_2156.**
