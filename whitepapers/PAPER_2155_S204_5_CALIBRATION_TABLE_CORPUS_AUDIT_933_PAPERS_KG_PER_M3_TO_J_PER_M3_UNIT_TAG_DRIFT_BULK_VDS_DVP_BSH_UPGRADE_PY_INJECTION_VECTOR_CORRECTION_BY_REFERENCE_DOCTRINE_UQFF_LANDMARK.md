# PAPER_2155 — S204.5 Calibration Table Corpus Audit: 933-Paper kg/m³ → J/m³ Unit-Tag Drift Correction

**`bulk_vds_dvp_bsh_upgrade.py` (Session 204, April 2026) identified as the template-injection vector; correction-by-reference doctrine canonized; 918 of 957 unit-labeled instances corrected without per-paper touches; code layer already correct**

**UQFF LANDMARK**

Author: Daniel T. Murphy
Framework: UQFF (Unified Quantum Field Framework) — Star-Magic v5.81+
Date: 2026-07-27
Provenance: PAPER_2153 arc Flag (g) — Daniel's 2026-07-27 ruling on massless-SCm unit-tag drift + `_audit_rho_scm_units.json` (2026-05-15, 61 KB, in-repo)
Status: Landmark — corpus-wide paperwork correction executed via correction-by-reference doctrine

---

## Executive Summary

Daniel's 2026-07-26 ruling canonized that UA and SCm are massless and cannot be expressed in kg — only in J (energy density). Subsequent PAPER_2153 arc deep-read (24 papers, PAPER_878-901) noted the Session 204 calibration table appearing verbatim in 11+ papers with `kg/m³` tags.

The pre-existing `_audit_rho_scm_units.json` file (already in the repo, run 2026-05-15) revealed the actual scope: **933 papers, 918 kg/m³ hits, 39 J/m³ hits, 14 papers with internal conflict**. The drift is a corpus-wide template-injection issue from `bulk_vds_dvp_bsh_upgrade.py` (Session 204).

**PAPER_2155 executes the correction via correction-by-reference doctrine** — a single landmark supersedes all 933 affected paper instances without requiring 933 individual REVISION appendices. Daniel's ruling: *"You do NOT need to touch any of the 933 individual papers. The unit tag in the prose is cosmetic documentation. The computations are correct — the numbers 7.09e-37 and 7.09e-36 are dimensionally consistent as J/m³ in all code. The papers are wrong only in their label string, not in their arithmetic. A corpus-wide supersession by landmark is the cleanest resolution consistent with Rule 9 (append-only) and the established PAPER_2147 pattern."*

**Zero calculator source changes.** Zero physics values changed. Zero per-paper appends across the 933 affected papers. Complete supersession via this landmark + PAPER_2147 REVISION 2026-07-27 companion + calculator code's already-correct J/m³ constants + gate assertions test-locking the code.

---

## 1. Audit Findings (from `_audit_rho_scm_units.json`, in-repo)

The audit was run 2026-05-15 (61 KB JSON, machine-generated). Corpus-wide grep for `7.09e-37` value with nearby unit-label detection.

**Full findings:**

| Metric | Count | % of unit-labeled |
|---|---|---|
| Papers containing `7.09e-37` value | 933 | — |
| Total hits with explicit unit labels | 957 | 100.0% |
| `kg/m³` labels (DRIFT) | **918** | **95.9%** |
| `J/m³` labels (correct) | **39** | 4.1% |
| Papers with internal conflict (BOTH units) | 14 | — |
| Papers `kg/m³` only | 898 | — |
| Papers `J/m³` only | 20 | — |
| Hits with NO unit label nearby | 29 | — |

**Interpretation:**
- The vast majority (95.9%) of unit-labeled instances carry the drift `kg/m³` tag
- Only 4.1% carry the correct `J/m³` tag (papers authored after PAPER_2147's canonization or by hand rather than template)
- 14 papers contain both units (internal inconsistency — presumably where hand-authored J/m³ content was juxtaposed with template-injected kg/m³ content)

**Sample affected papers (from audit):**
- `PAPER_001_GW170817_UQFF_Damping_Analysis.md` (kg/m³ only)
- `PAPER_002_GW190425_Mass_Gap_Interpretation.md` (kg/m³ only)
- `PAPER_036_FUBii_Buoyancy_Variant1_Archimedes_UQFF.md` (both units — conflict)
- `PAPER_145_UQFF_MUGE_Compression_Cycle3_Unified_Framework_12Term_Resonance.md` (both units — conflict)
- `PAPER_1066_UQFF_Lagrangian_Derivation.md` (J/m³ only — correct)
- `PAPER_1109_26Level_Vacuum_Density_Ladder_Ramanujan.md` (J/m³ only — correct)

Full enumeration of all 933 papers remains in `_audit_rho_scm_units.json`.

---

## 2. Root Cause: `bulk_vds_dvp_bsh_upgrade.py` Session 204 Injection Vector

**File:** `bulk_vds_dvp_bsh_upgrade.py` (20 KB, dated Apr 8 2026, Session 204)

**Offending lines (34–35):**
```python
RHO_SCM     = 9.47e-27     # kg/m³
RHO_UA      = 5.0e-27      # kg/m³
```

**Impact:** The script bulk-injected the §B "VDS/DVP/BSH Deep Synthesis" appendix (including the §B.4 "Production-Scale Consistency" table listing `rho_SCm = 7.09 × 10⁻³⁷ kg/m³`) into PAPER_001 through PAPER_877 in a single pass. The `kg/m³` label was baked into the template as a fixed string. Every paper that received the bulk-upgrade received the drifted tag.

**Note:** the specific density VALUES `9.47e-27` and `5.0e-27` in the script are ALSO non-canonical (Daniel's Flag (f) ruling: bulk-script artifact using undocumented density values inconsistent with canonical constants). That value-drift issue is addressed separately in PAPER_2156 (Phase 3 in-flight — 1.894 ratio audit). PAPER_2155 addresses ONLY the unit-tag drift on the canonical values.

---

## 3. Framework State — Code Layer Already Correct

Daniel's ruling documented that the CALCULATOR CODE already uses J/m³ correctly (drift is exclusively in whitepaper prose). Verified locations:

| File | Content | Status |
|---|---|---|
| `QCalc_Performance.py` | ρ_SCm declaration | J/m³ ✓ |
| `CondensedPhysics_InputData.py` | ρ_SCm declaration | J/m³ ✓ |
| `ScmVacuumDensityModule.cpp` | ρ_SCm declaration | J/m³ ✓ |
| `_smoke_qv_bundle9.py` | Test values | J/m³ ✓ |
| `provenance_recorder.py` | Metadata | J/m³ ✓ |
| `scm_vacuum_manifold.py` | Proof function | Prints `[MASSLESS]` explicitly ✓ |

**The framework code layer is CORRECT. The drift is EXCLUSIVELY in whitepaper prose templates.** Since the numbers `7.09e-37` and `7.09e-36` are dimensionally consistent in the code as J/m³ (energy densities of the massless SCm and UA), and the calculators use them correctly, the drift is COSMETIC — a label-string error in documentation, not a computational error.

---

## 4. Correction-by-Reference Doctrine (canonical ruling)

**Daniel's ruling (verbatim):** *"You do NOT need to touch any of the 933 individual papers. The unit tag in the prose is cosmetic documentation. The computations are correct — the numbers 7.09e-37 and 7.09e-36 are dimensionally consistent as J/m³ in all code. The papers are wrong only in their label string, not in their arithmetic. A corpus-wide supersession by landmark is the cleanest resolution consistent with Rule 9 (append-only) and the established PAPER_2147 pattern."*

**Formal doctrine:** when a corpus-wide drift affects `N` papers due to a single template-injection vector (e.g., a bulk-upgrade script), and the underlying calculator code is correct, a SINGLE landmark declaration + companion REVISION append to the standing rule paper (PAPER_2147 in this case) is SUFFICIENT to correct all N papers without touching them individually.

**Enforced by three overlapping mechanisms:**
1. **This landmark (PAPER_2155)** — formal supersession declaration
2. **PAPER_2147 REVISION 2026-07-27** — extends the standing rule paper's authority to cover the S204.5 calibration table specifically
3. **Calculator code J/m³ constants** — the code IS the ground truth; the papers only carry documentation strings
4. **Gate assertions (this paper §6)** — test-locking the J/m³ units in the calculator code prevents future drift

**Explicit supersession:** all S204.5 calibration table instances anywhere in the whitepaper corpus that list:
- `rho_SCm = 7.09 × 10⁻³⁷ kg/m³`
- `rho_UA = 7.09 × 10⁻³⁶ kg/m³`

**ARE SUPERSEDED** and read as:
- `ρ_SCm = 7.09 × 10⁻³⁷ J/m³` (massless SCm vacuum energy density)
- `ρ_UA = 7.09 × 10⁻³⁶ J/m³` (massless UA vacuum energy density)

**No per-paper edit required.** No REVISION append required for any of the 933 affected papers. The correction is complete via this landmark's declaration.

---

## 5. What This Doctrine is NOT

**Not a get-out-of-audit-free card.** The correction-by-reference doctrine applies ONLY when:
- The drift affects a corpus-wide template (i.e., a bulk-script injection, not per-paper hand-authoring)
- The calculator code layer is already correct (numbers correct, only strings drifted)
- The drift is cosmetic/label-only (no computational impact)
- A companion REVISION to an existing standing-rule landmark (e.g., PAPER_2147) formally extends the standing rule's authority

**When per-paper REVISIONS ARE required:** if the drift affects computational values (not just labels), or if the calculator code is also wrong, or if there is no existing standing-rule landmark to extend, individual per-paper REVISIONS remain necessary.

**Not silent drift correction.** Every future reader of the 933 affected papers must be able to trace back to this landmark and PAPER_2147 REVISION 2026-07-27 to see the corrected reading. That's why the correction is DECLARED explicitly here (not hidden), and why the audit JSON (`_audit_rho_scm_units.json`) remains in the repo as the machine-generated record of scope.

---

## 6. Gate Assertions (test-locking the code layer)

To ensure the code layer remains correct in perpetuity (preventing future re-drift), the following gate assertions are added to `uqff_fidelity_tests.py` in this paper's companion Phase 2 gate-pins section:

```
1. RHO_SCM constant in canonical calculator modules is 7.09e-37 J/m³ (not kg/m³)
2. RHO_UA constant in canonical calculator modules is 7.09e-36 J/m³
3. RHO_UA / RHO_SCm = 10.0 EXACT (F_TRZ = 1/10 canonical coupling per PAPER_890/140/1160)
4. Calculator modules explicitly declare J/m³ unit in their constant declarations
5. `scm_vacuum_manifold.py` proof function returns [MASSLESS] verdict for SCm/UA
```

These assertions test-lock the code's J/m³-native discipline. Any future accidental reversion of a canonical calculator constant to `kg/m³` label would fail the gate.

---

## 7. Standing Rules Canonized by PAPER_2155

1. **Template-injection audit is standing procedure:** any bulk-upgrade script that injects prose into multiple papers (like `bulk_vds_dvp_bsh_upgrade.py`) must be audited BEFORE execution to ensure the injected template is drift-free. Retrospective audits (like `_audit_rho_scm_units.json`) should be run when template-injection drift is suspected.

2. **Correction-by-reference is preferred over mass-append:** when a corpus-wide drift affects a large number of papers due to template injection, a single landmark supersession is preferred over N individual per-paper appends. Rule 9 (append-only) is satisfied by ONE authoritative supersession, not N minor appends.

3. **Code layer is ground truth for computational values:** when the calculator code and whitepaper prose disagree on a computational value (constant, unit, formula), the code layer takes precedence. Whitepapers are documentation — code is the framework's operational truth.

4. **In-repo audit files are canonical evidence:** files like `_audit_rho_scm_units.json` produced by systematic corpus scans are canonical evidence for framework-state claims. Future audits should preserve their audit-file outputs rather than re-running scans that produce different snapshots.

---

## 8. Falsifiable Predictions

1. **`_audit_rho_scm_units.json` re-run after PAPER_2155 canonization should show the same 933-paper scope.** If a future re-run shows a different scope, either (a) more papers have been added to the corpus (expected), or (b) some existing papers have been silently modified (unexpected — would be a Rule 9 violation flag).

2. **Framework code layer will remain J/m³ correct.** Gate assertions in §6 test-lock this. Any future gate failure on the J/m³ assertions would indicate someone reverted a canonical constant.

3. **No per-paper appends across the 933 affected papers.** Any future session that attempts to individually revise the 933 papers is going against the correction-by-reference doctrine — a Rule 9 architectural violation.

---

## 9. Files Touched

- `whitepapers/PAPER_2155_S204_5_CALIBRATION_TABLE_CORPUS_AUDIT_..._UQFF_LANDMARK.md` (this file)
- `pdf2/PAPER_2155_...pdf` (companion PDF)
- `whitepapers/PAPER_2147_J_PER_M3_NATIVE_..._UQFF_LANDMARK.md` — REVISION 2026-07-27 append (companion, extends standing rule authority)
- `uqff_fidelity_tests.py` — +5 PAPER_2155 gate assertions (RHO_SCM/RHO_UA J/m³ test-lock)
- `CLAUDE.md` — APPENDED section
- Zero calculator source changes (code layer already correct)
- Zero physics values changed
- Zero per-paper touches to the 933 affected papers (correction-by-reference doctrine)

---

## 10. Cross-References

- **PAPER_2147** — Original J/m³-native discipline landmark (2026-07-25); REVISION 2026-07-27 companion extends its authority
- **PAPER_2148** — UQFF Ontology Answer B (vacuum energy fundamental, mass emergent — SCm massless per this ontology)
- **PAPER_2153** — Joint SCm+UA vacuum-density engine (§7 Standing Rule 2 flags SCm direct-detection claims — massless SCm consistent)
- **`_audit_rho_scm_units.json`** — In-repo audit record (2026-05-15, 61 KB, machine-generated), canonical evidence for the 933-paper scope
- **`bulk_vds_dvp_bsh_upgrade.py`** — In-repo injection vector script (Session 204, April 2026, lines 34–35 non-canonical densities with drift kg/m³ label)
- **PAPER_2156** (Phase 3 in-flight) — 1.894 ratio bulk-script artifact correction; separate from unit-tag drift but same injection vector
- **PAPER_890/140/1160** — Canonical ρ_SCm/ρ_UA = 1/10 = F_TRZ structural coupling
- **Calculator canonical J/m³ constants:** `QCalc_Performance.py`, `CondensedPhysics_InputData.py`, `ScmVacuumDensityModule.cpp`, `provenance_recorder.py`, `scm_vacuum_manifold.py`
- **Daniel's 2026-07-26 massless-SCm ruling:** *"UA and SCm are massless, so it can not be expressed in terms of (kg), only J energy; this was AI drift."*
- **Daniel's 2026-07-27 Flag (g) ruling:** correction-by-reference doctrine per PAPER_2147 REVISION + new landmark pattern

**End of PAPER_2155.**
