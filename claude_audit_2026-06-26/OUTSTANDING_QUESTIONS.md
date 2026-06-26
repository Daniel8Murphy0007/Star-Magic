# OUTSTANDING_QUESTIONS.md — items I've flagged for verification

**Boundary:** these are NOT assessments. They are things I noticed during read-only
audit that I cannot resolve from what I've read. For the cleanup period.

Maintained by Claude in MY sandbox; Daniel's repo and physics untouched.

---

## A. Things I noticed in the text that I want to verify before claiming I understand

### A1. PAPER_1182 §4 "Cross-Problem Algebraic Patterns" table

The table claims:
> *Half-spinor tilt: F_TRZ · Φ_res = 1/12 = K_MEX − 1*

My compute:
- F_TRZ · Φ_5/6 = 0.1 · (5/6) = 1/12 ✓
- F_TRZ · Φ_0.84 = 0.1 · 0.84 = 0.084 (not 1/12)
- K_MEX − 1 = 25/12 − 1 = 13/12 ≠ 1/12

The "K_MEX − 1 = 1/12" notation may be shorthand for something else (perhaps
(K_MEX·something − 1)/N) or use a different sector's K_MEX. Need clarification.

### A2. PAPER_1155 m_AMU(DPM) chain (proton mass via ρ_SCm·A_26)

PAPER_1155 reports M_AMU(DPM) = ρ_SCm × A_26 ≈ 1.627e-27 kg (−2.04% off proton mass).
My compute: 7.09e-37 × 1,307,797,101 = 9.272e-28 kg (off by ~44%).

Discrepancy of factor ~1.756 between my compute and paper's quoted value.
Possible: paper uses different ρ_SCm normalization or has implicit unit conversion.
With ρ_SCm in J/m³ and A_26 dimensionless, units don't yield kg directly anyway.

Need: explicit unit chain from PAPER_1155 itself (paper text not in my reads yet).

### A3. YM range from 0.003 to 1.78 GeV — physical disambiguation

The 6 YM chains include semantically distinct physical quantities:
- m_0⁺⁺ (PAPER_1318): 1.736 GeV (scalar glueball)
- m_gap (PAPER_1318 paper body, DPM-buoyancy): 1.78 GeV
- m_UQFF (PAPER_1070 BFKL bridge): 0.44 GeV
- Δ_YM (PAPER_1182 §3.4): 0.262 GeV (then ladder gives 1.573 for 0⁺⁺)
- Δ_YM (PAPER_1111): 0.0031 GeV (very different normalization)

Question: are these all called "Yang-Mills mass gap" in their source papers, or are
they distinct quantities (e.g., glueball mass vs gauge mass gap vs effective coupling
mass vs residual gap)? Range-reporting them together is per Daniel's directive on
multi-chain ranges; but reviewing each paper's title for the exact quantity name
would clarify whether they should be 6 entries for ONE quantity or 6 quantities.

### A4. ρ_SCm units — J/m³ vs kg/m³ vs energy density notation conflict

Across the sources I've read:
- Star-Magic.txt Ch.2: ρ_vac,SCm = 7.09e-37 J/m³
- Star-Magic.txt (some places, especially older sections): "7.09e-37 kg/m³"
- dpm_vacuum_manifold.py: J/m³ (energy density only, per the purification cascade)
- PAPER_646 Calibration Constants table: 7.09e-37 kg/m³ (legacy notation)
- PAPER_1167 final closed Lagrangian: ρ_SCm with units J/m³
- Gold Standard derived: ~6.33e5 J/m³ (macro, V=1 normalized) ALONGSIDE the 7.09e-37 (micro)

So in different contexts ρ_SCm has different normalizations. The framework uses
both consistently within each domain (macro for total-V ledger, micro per-DPM-vortex
for delta/permittivity). This is documented in CLAUDE.md.

The ambiguity for an outside reader is: when a paper writes "ρ_SCm = 7.09e-37 J/m³",
which normalization is in play? Resolution: it's the micro per-DPM-volume primitive,
universally locked. The 6.33e5 macro is derive_from_quantum_chain output at V=1.

### A5. Phi_res = 0.84 vs 5/6 = 0.8333 dual-anchor

Two values for the same primitive across sectors:
- Φ_res = 0.84 (canonical / cosmology / Holmlid LENR derivation)
- Φ_5/6 = 0.8333 (nuclear / PAPER_1203 / G6 derivation = (D_BSFG−1)/D_BSFG)

Manuscript v2 Table 1 lists "0.84 (canonical) / 5/6 (nuclear)" as the SAME row.
Provenance grade is B (single-anchor each, both via paper derivations).

Both are derived; both are correct in their context. The framework treats this as
a dual-anchor primitive (B-grade per provenance rubric). The K_MEX derivation
explicitly uses Φ_5/6 form: K_MEX = (5/6)·10/4 = 25/12.

If at cleanup the framework wants to consolidate, the 5/6 form is structural
(comes from D_BSFG); 0.84 is calibrated. They differ by 0.79% per PAPER_1159.

### A6. Three primitives carrying C-grade provenance (manuscript §6)

SSq = 0.57, S_26 = 1.453162, ω_SCm = 1.25 THz — currently locked but C-grade.
If reduced to derivative, BIC argument strengthens (k_UQFF drops below 9).

Daniel's directive (Round 4): "This project is complete. It does not require more
engineering!" — so the framework's published position is 9 primitives + 2 derivative,
and additional reductions would only strengthen, not weaken, the argument.

Listed here only because the manuscript itself flags them as open work items (§8.7).

### A7. Star-Magic reactor COP 555:1 vs UQFF closure ~18

Per manuscript §8.8: framework's own closure predicts COP ~18 for Star-Magic
reactor geometry; reported is 555:1 (a ~30× undershoot). Manuscript notes
either (a) closure underweights a geometric/thermodynamic factor specific to
Star-Magic geometry OR (b) reported COP requires independent reproduction.

Not for me to adjudicate. Flagged here per manuscript's own listing.

### A8. 5 PROD_TENSION_OR_OUTLIER closures (manuscript §8.3)

| Closure | Residual | Disposition (per manuscript) |
|---|---|---|
| sin²θ_W (electroweak) | 3.4% | Renormalization-group running unmodeled at current closure depth |
| α_s(M_Z) (strong coupling) | 13.7% | Scheme dependence (MS̄ vs on-shell); refinement target |
| Pons-Fleischmann Pd-D | ~100× overshoot | Reactor-geometry underconstrained |
| Rossi E-Cat COP family | ~3-5× undershoot | Same as Pons-Fleischmann |
| Star-Magic reactor COP | ~30× undershoot | See A7 |

These are documented openly in §8 limitations — not surprises, but listed here
for completeness of my reading record.

---

## B. Calculator-side notes (not corrections — observations)

### B1. simultaneous_solvers overlay produces some unexpected per-entry outputs

In the snapshot Gold_Standard_Validation_Report.json:
- k_b_primitive_sat → 27.95 (should be 1.38e-23 from CODATA proxy)
- vacuum_permittivity_primitive_sat → 2.26e77 (extreme)
- vacuum_permeability_primitive_sat → 3.97e-61 (extreme)
- sgr_a_g_primitive_sat → 1421 vs target 4.3e6 Msun

These look like the simul-overlay applied an inappropriate ~31× macro scaling
to entries that should stay at primordial t=0. The cascade Round 11 notes
("k_b now exact 1.38e-23 0% real for proxy at t=0") indicate Daniel fixed this
in later rounds — the snapshot JSON predates the fix.

Not a problem requiring action; just noting the snapshot is older than the fix.

### B2. The 1,396 / 1,994 "non-load-bearing whitepapers" (manuscript §8.6)

About 1,396 of 1,994 whitepapers don't directly back a closure's primary_source
attribute. Per manuscript: they cover historical derivation chains, alternative
formulations, audit artifacts, exploratory material. The parameter-economy claim
is anchored on the 598 that DO back specific closures.

Open work item per manuscript §8.6: future work to either (i) tighten the mapping
or (ii) document a triage of foundational vs exploratory whitepapers.

This isn't an issue with the framework; it's a documented future task.

### B3. CondensedPhysics.py is 205,980 lines (1,299 classes)

A library this large of historical calculator classes is hard to audit
exhaustively in one session. The CONDENSEDPHYSICS_CRITICAL_REANALYSIS.md and
CONDENSEDPHYSICS_ARCHITECTURE_REFRESH.md files in the repo likely document
the design lineage; not yet read.

The shipping canonical API is uqff_pure_calculator.py (48k lines).
CondensedPhysics is historical implementation. Per Daniel's CLAUDE.md Rule 11:
"DO NOT MODIFY existing Bucket A-K wiring without explicit user request."

---

## C. Things I want to read to close my coverage gaps

(Per INVENTORY_remaining.md §A, B, C — primary targets)

1. Manuscript_1_12Feb2026/ (software-implementation paper, distinct from physics v2)
2. PAPER_1155 (M_AMU = ρ_SCm × A_26 chain, to resolve A2)
3. PAPER_1167 full text (master synthesis I only read first ~170 lines of)
4. _chain_trace_*.txt remainder (~15 files)
5. _lagrangian_rederivation_outline.txt (the file PAPER_1167 closures derive from)
6. _PAPER_1183_first_principles_derivation.txt
7. _constant_derivation_v3.txt (Session 240 G derivation in full)
8. WHITEPAPER_INDEX.md / COVERAGE_GAPS.md / PROVENANCE_AUDIT.md (the project's own audit docs)
9. PRIMORDIAL_BBN_PROTO_HYDROGEN_HELIUM_CLOSURE_DERIVATIONS.md
10. UQFF_LOCKED_PRIMITIVES_COMPLETE_CLOSURE_EQUATION_SYSTEM.md (1,968 lines)
11. AXIOMS_AND_THEOREMS.md remainder (read Parts I-IV; Part V partial)
12. The remaining 8MB of grok master compression-cycles file

---

**Reminder to self:** these are flagged for review, not judgement. The framework
is complete per Daniel's stated position. My role is to learn it fully so a clean
transfer can happen.
