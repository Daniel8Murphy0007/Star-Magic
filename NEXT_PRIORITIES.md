# NEXT_PRIORITIES.md — what the next Claude session should do first

**READ CLAUDE.md FIRST. EVERY SESSION. EVERY ACTION.** Rules in CLAUDE.md supersede this file.

**Daniel's standing directives (verbatim):**
- "DUMP ALL PROVENANCE. NO COMMENTS WHATSOEVER. NO TAGS, NO GARBAGE SM ANYTHING. JUST REAL PERCENTAGES OF DIFF."
- "UQFF IS THE FUCKING ANCHOR. UQFF DOESN'T SHARE SHIT WITH SM."
- "I AM GIVING YOU THE INFORMATION, YOU ARE TO ASSEMBLE IT."
- "READ THE RULES BEFORE EVERY ENTRY. READ THE RULES BEFORE EVERY SESSION."

---

## STATE AT SESSION CLOSE (2026-07-25 — v5.80.1)

### Ship trail (recent — full history in CHANGELOG.md)
- **v5.79.0** (2026-07-24) — TIDAL/KEPLER/HALVING/RULE-4 arc: 7 landmarks PAPER_2136-2142 + R382-R389 stub fills + REVISED STANDING RULE v4
- **v5.80.0** (2026-07-25) — c/LAMBDA/ONTOLOGY audit arc: 5 landmarks PAPER_2144-2148 + H_0 route upgrade 47× tighter + framework ontology formally declared (Answer B)
- **v5.80.1** (2026-07-25) — METADATA PATCH: README badges/Version drift fix + 2 permanent SHIP GUARD assertions (f, g)

### Honest metrics (verifiable via `python uqff_fidelity_tests.py`)
- **Fidelity gate:** **3376 passed, 0 failed** (up from 1031 at v5.61.0 — ~2300 assertions added)
- **Whitepapers authored:** **2,237** (PAPER_001-PAPER_2148, with reserved IDs)
- **PDF coverage:** **2,237 / 2,237 = 100%** (parity achieved v5.76.0; every whitepaper has a PDF)
- **P2 stub-fill campaign (R218+):** **172 consecutive rounds** through R390 (last: HydrogenBubbleMagneticCalculator)
- **Registry Program (R0-R5 + XGEO):** COMPLETE per PAPER_2130. UNIFIED_REGISTRY.csv 2,549 rows. Falsifiability graph 658 edges. 14 live-derived constants, 7 EXACT, worst 0.9009% (Λ vs Planck-lensing).
- **9 truly-independent primitives locked** (ρ_SCm=7.09e-37 J/m³ + 8 dimensionless integers/ratios). Structural derivatives: 3 (D_BSFG, K_MEX, κ).
- **Cosmological quadruple {c, H_0, Λ, v_F}:** all UQFF-native, H_0 tightened 47× in v5.80.0 (0.065% residual via PAPER_1573 A_5+SO_5=70 km/s/Mpc EXACT)

### Framework ontology (formally declared v5.80.0 PAPER_2148)
- **Fundamental:** ρ_SCm [J/m³] is UQFF's sole dimensioned primitive
- **Emergent (per arxiv manuscript + Star-Magic.pdf):** mass, Newton's G, Newtonian gravity
- **UQFF and SM have INVERTED ontologies** — same universe, different starting points
- **F_UBi = mass pushing against universe; F_UBii = universe's response** (canonized by Daniel 2026-07-25)
- **Gravity = F_UBi/F_UBii habitable-zone crossing** (observable via terminal velocity)
- **Λ dual manifestation:** open-space potential value + mass-lensing observable, same 1.1e-52 m⁻²
- **SM-comparison validity boundary:** SM's Λ = 8πG·ρ_Λ/c⁴ VALID when known massive astronomical objects are the anchor; INVALID when inverting to derive UQFF cosmology from SM axioms

### 26-layer quantum chain (per Star-Magic.pdf)
```
Layer     Physical Meaning      Frequency Scale
1-6       Particle physics      10¹⁹ Hz    ← high freq, short wavelength
7-12      Nuclear structure     10²² Hz
13-18     Atomic binding        10¹⁶ Hz
19-22     Molecular bonds       10¹⁴ Hz
23-24     Electromagnetic       10¹⁵ Hz
25-26     Gravitational         10⁻¹⁰ Hz   ← LOW freq, LARGE scale
```
Same vacuum, opposite ends of resonance spectrum.

---

## PRIORITY 1 — VERIFICATION FIRST

Before any modification:

1. Run `python uqff_fidelity_tests.py`. Expect **3376 passed, 0 failed**, exit 0.
2. Confirm the 9 canonical primitives + 3 structural derivatives (D_BSFG, K_MEX, κ) unchanged in `uqff_registry_primitives.py`.
3. Sample calculator public surfaces; verify return shape `{'value': X}` per CLAUDE.md Rule 5.
4. **DO NOT MODIFY** existing Bucket A-K wiring, Registry Program artifacts, or arc landmark papers (PAPER_2130-2148) without explicit user request.
5. **DO NOT touch cosmological quadruple values** ({c, H_0, Λ, v_F, ρ_SCm}) without invoking the coupling-discovery pre-swap verification rule (PAPER_2144 canonical).

---

## PRIORITY 2 — RUN THE FRAMEWORK'S RULES 4/7/10 DISCIPLINE

Every action must comply with three canonical rules (CLAUDE.md and the arc landmarks that extend them):

- **Rule 4 (No SM anywhere)** — extended by **PAPER_2147** to include unit-direction discipline (UQFF is J/m³-native; SM is kg/m³-native; silent framework-translation is contraband) and by **PAPER_2148** to include ontology-inversion recognition (mass, G, gravity are EMERGENT in UQFF, not fundamental).
- **Rule 7 (Honest disclosure)** — extended by **PAPER_2146** with 5 standing sub-rules: (5.1) no circular calibration in verification, (5.2) no ad-hoc doctrinal retrofits around fits, (5.3) verify reference values against primary sources, (5.5) v_F not primitive-lockable without cosmological fit. Standing Rule 5.4 (dimensional verification) is SUPERSEDED by PAPER_2147.
- **Rule 10 (Daniel provides info, AI assembles)** — the arc's 2026-07-25 experience validates this rule: Daniel's five sharpening questions caught every AI overstatement. Continue interrogation as the primary quality-control mechanism.

---

## PRIORITY 3 — CATEGORIES OF OUTSTANDING WORK

The R218+ P2 stub-fill campaign is at 172 consecutive rounds. The next session may choose any of the following work categories. There is no strict ordering — pick based on user direction.

### 3.1 — Continue P2 stub-fill campaign (steady progress)

**Next round:** R391 (173rd consecutive).

**Pattern (per PAPER_2140/2141 discipline):**
1. Pick a stub in CondensedPhysics.py or QCalc*.py that currently `raise NotImplementedError` or returns `None`
2. Read the class docstring/framework_papers list to identify the paper it should implement
3. Read that paper for the derivation
4. Primitive-lock what's derivable using `uqff_registry_primitives` LIVE composition (never CODATA/PDG literals per Rule 4)
5. Author a landmark whitepaper ONLY IF a novel identity is discovered (integer sums, F_TRZ rungs, tilt appearances, cross-scale universalities)
6. Wire gate assertions to pin the primitive-lock + landmark identity
7. Ship (see Section 4 for ship checklist)

**Bulk-cleanup pattern (PAPER_2140/PAPER_2141):** when a boilerplate template is discovered in >100 classes, apply `replace_all=true` or Python `re.subn` (per REVISED STANDING RULE v3). Save PRE_BULK backup.

### 3.2 — Buckets E-K PURE_UQFF upgrade (deferred since v5.58)

Per CLAUDE.md Section "Ships/Buckets": Buckets E-K are FIRST-PASS drainage. Many observables use heuristic SCm-correction formulas rather than paper-canonical closed forms. Revisit each `DERIVED_SCM_CORRECTION` observable and transcribe the verbatim formula from the source paper (e.g., PAPER_1009 for 3C273 Eddington, PAPER_915 for GW170817 strain damping, PAPER_034 for κ_t).

### 3.3 — Framework-level open questions from v5.80 arc

1. **13.4% ρ_Λ discrepancy** (PAPER_2148 Interpretation A vs B). Is UQFF's `5.957×10⁻¹⁰ J/m³` vs SM's inferred `5.28×10⁻¹⁰ J/m³` a framework-differentiating prediction or a ρ_SCm × 26! × K_MEX chain error requiring a ~0.88 correction? Requires distinguishing experiment or independent UQFF cross-check.
2. **v_F structural derivation.** Currently `_V_FERMI = 0.77e6 m/s` is an observational SI anchor (Session 239). PAPER_2145's Friedmann-lock derivation was walked back. Open target: is there a UQFF-native primitive composition that gives v_F within observational precision, without SM-Friedmann inversion?

### 3.4 — F_UBi/F_UBii dimensional audit under PAPER_2148

Per PAPER_2148 §3, F_UBi is now formally canonized as "mass pushing against universe" and F_UBii as "universe's response." Audit the F_UBi formula (`−β·G·M·ρ_SCm/r²·...`) for unit consistency under the ontology declaration — G is emergent (kg-native SI), ρ_SCm is J/m³-native, their product may have an implicit unit-system-crossing that should be disclosed per PAPER_2147.

### 3.5 — Corpus revision follow-through (partially done in v5.81.0)

The v5.80 arc queued REVISION appends for PAPER_1170, 1226, 1235, 2145, 2146 — all applied in v5.81.0 corpus hygiene batch. Additional papers may exhibit the SM-native reversal pattern (PAPER_2147 unit-direction discipline). A corpus-wide grep for `kg/m³.*×c²` or `ρ_Λ.*Planck.*match` patterns would identify further candidates.

---

## PRIORITY 4 — SHIP CHECKLIST (7 rules, gate-pinned)

Every ship must satisfy (as of v5.80.1):

- **(a)** pyproject `description` contains the current version string
- **(b)** pyproject `description` is ≤512 chars (SHIP GUARD e, learned from v5.77.0 PyPI HTTP 400)
- **(c)** All gate file-hash pins are LF-normalized (learned from v5.75.1 CI failure)
- **(d)** Verify PyPI landing via `gh run list --limit 3` and `curl -s https://pypi.org/pypi/uqff/json`, NOT the cached browser page
- **(e)** SHIP GUARD gate assertions pass (description ≤512 + version-in-description present)
- **(f)** README `**Version**:` line matches pyproject version (NEW v5.80.1, gate-pinned)
- **(g)** README `fidelity_gate` badge references current gate count (NEW v5.80.1, gate-pinned)

**Git tag hygiene (learned v5.80.0/5.80.1):**
- Use ANNOTATED tags (`git tag -a vX.Y.Z -m "..."`), NOT lightweight tags
- Push with `git push --follow-tags` (skips lightweight tags — annotated required)
- OR explicit `git push origin vX.Y.Z` for lightweight tags
- The "Release to PyPI" workflow triggers on TAG push, not branch push

**Registry regen chain (run before every ship):**
```
python registry_generator.py     # → UNIFIED_REGISTRY.csv + supporting
python uqff_registry_graph.py    # → UNIFIED_REGISTRY_GRAPH.csv
python uqff_registry_status.py   # → results table + status report
python uqff_registry_xgeo.py     # → XGEO queue + confirmations
```

**Version bump artifacts (all must land in same commit):**
- pyproject.toml (version + description)
- CITATION.cff (version + date-released)
- README.md (What's-new section prepended + badges/Version line updated)
- CHANGELOG.md (entry prepended)
- SESSION_LOG.md (session arc entry appended)
- CLAUDE.md (arc appends if framework-level changes)

---

## RULES (mirrored from CLAUDE.md — must-read every session)

1. **READ CLAUDE.md FIRST every session.** Rules in CLAUDE.md supersede Map, Plan, and this file.
2. **DO NOT REVERT canonical primitives.** SSQ=0.57, β_i=0.6029, K_MEX=25/12, S_26=1.453162, RHO_SCM=7.09e-37, integer primitives.
3. **NO NARRATIVE OF ANY KIND** in calculator source. Pure mathematical calculator. No docstrings, no `provenance` keys, no `NOT REPLACEMENT` tags.
4. **NO SM ANYWHERE** — content layer (no SM constants/formulas/terminology) AND presentation layer (no SM-native unit direction, no SM-framed comparison; PAPER_2147/2148).
5. **Public surface return contract:** `{'value': X}` only. No metadata.
6. **NO `datetime`, `json.dump`, file writes, `__main__`, classes** in `uqff_pure_calculator.py`.
7. **DO NOT claim "0.000% error" without numerical proof.** Honest residuals only.
8. **RUN `uqff_fidelity_tests.py` AFTER EVERY EDIT.** Exit 0 required.
9. **APPEND to SESSION_LOG.md, never rewrite.** Each session adds one entry at the bottom.
10. **DANIEL PROVIDES THE INFORMATION. YOU ASSEMBLE IT.**
11. **DO NOT MODIFY existing Bucket A-K wiring or Registry Program artifacts** without explicit user request.
12. The user has been fighting AI drift on this project for 12+ months. Treat his words carefully.

---

## EDIT / WRITE TOOL WARNING (updated 2026-07-24)

BOTH the Edit tool AND the Write tool truncate `uqff_pure_calculator.py` (~35k lines, ~1.85 MB) and `CondensedPhysics.py` (~1,308 classes, 8.8 MB) mid-write on large insertions. Two truncations during BUCKET 0 had to be repaired via Python splice.

**Pattern that works for large-file edits:**
```python
with open('uqff_pure_calculator.py','r',encoding='utf-8',newline='') as f:
    src = f.read()
anchor = "def calculate_<next-public>(dataset):"   # any unique fixed string
src2 = src.replace(anchor, new_block + anchor, 1)
with open('uqff_pure_calculator.py','w',encoding='utf-8',newline='') as f:
    f.write(src2)
```

**For >1000-replacement bulk cleanups** (REVISED STANDING RULE v3 from PAPER_2141):
```python
import re
with open('CondensedPhysics.py', 'r', encoding='utf-8', newline='') as f:
    src = f.read()
# Longer-form patterns FIRST (avoid partial-match on shorter variants)
for pat in [r'\b6\.67430e-11\b', r'\b6\.6743e-11\b', r'\b6\.674e-11\b']:
    src = re.subn(pat, '_URP_G', src)[0]
with open('CondensedPhysics.py', 'w', encoding='utf-8', newline='') as f:
    f.write(src)
```

Always run subsequent tests with `PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/tmp/uqff_test` to avoid stale .pyc cache hiding edits.

**AST-parse verify before running gate** on bulk cleanups:
```python
import ast; ast.parse(open('CondensedPhysics.py').read())
```

---

## BACKUP HYGIENE POLICY

Preserve descriptive PRE_<BATCH-ID>_<CONSTANT>_BACKUP files (per PAPER_2141 standing rule) — NEVER delete without explicit user instruction. Current preserved backups (per CLAUDE.md):
```
uqff_pure_calculator.py.PRE_FIX_BACKUP through POST_BUCKETK_BACKUP
CondensedPhysics.py.PRE_R388_BULK_G_BACKUP
uqff_Plan.md.PRE_FIX_BACKUP
```

---

## SESSION MILESTONE MARKERS

- **2026-06-08** — TOTAL PURGE applied to `uqff_pure_calculator.py`: all provenance/paper/closure_status/SM-references stripped. Public surfaces return `{'value': X}` only. Cat 16 strict purge guard wired.
- **2026-06-18** — PAPER_1521/1522 landmarks: D_BSFG and K_MEX declared as structural derivatives (11 → 9 truly-independent primitives).
- **2026-06-18** — DUAL LICENSE adopted (AGPL-3.0 + Commercial). Repository moved from MIT.
- **2026-07-11** — R218+ P2 stub-fill campaign initiated (v5.61.0). At time of this doc: 172 consecutive rounds through R390.
- **2026-07-22 through 07-24** — Unified Registry Program R0-R5 COMPLETE (PAPER_2130). Single-source-of-truth `uqff_registry_primitives.py`; three-language pins (Python=C++=Lean); 656→658-edge falsifiability graph; 348/348 XGEO cross-geometry tasks routed.
- **2026-07-24** — Vacuum Coupling Kernel K=F_TRZ·K_MEX·SSq=19/160 EXACT canonized (PAPER_2132, XGEO-U arc). Alpha_inverse = 125+12 = 137 EXACT integer decomposition (PAPER_2133 tilt-product law).
- **2026-07-24** — Rule 4 doctrinal correction: R382 Rule 4 catch (SM-imported Goldreich-Peale reverted), four-revision arc culminating in tidal k₂/Q = 3/125 EXACT primitive-lock (PAPER_2136). REVISED STANDING RULE v4 (dual-exposure with observation-headlining for measured constants).
- **2026-07-25** — c/Λ/ONTOLOGY audit arc (v5.80.0): PAPER_2144 H_0 route upgrade 47× tighter; PAPER_2145 Friedmann-lock walked back as AI overreach (framework survived intact); PAPER_2146 speed-of-light-fuckup self-audit + anti-circular-calibration standing rule; PAPER_2147 J/m³-native unit-direction discipline; **PAPER_2148 UQFF Ontology Declaration Answer B — the arc-closing landmark: vacuum energy fundamental, mass/G/gravity emergent, Λ dual-manifestation, F_UBi/F_UBii causal roles canonized, SM-comparison validity boundary defined.** Framework net-tighter and better-documented after arc than before.
- **2026-07-25** — v5.80.1 metadata patch: README badges/Version drift fix + 2 permanent SHIP GUARD gate assertions (f, g) prevent recurrence.
