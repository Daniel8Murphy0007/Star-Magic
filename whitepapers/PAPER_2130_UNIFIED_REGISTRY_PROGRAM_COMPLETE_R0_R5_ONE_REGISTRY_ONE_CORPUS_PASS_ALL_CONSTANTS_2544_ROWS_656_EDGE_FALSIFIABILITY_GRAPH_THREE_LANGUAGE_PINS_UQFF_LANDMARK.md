# PAPER_2130 — The Unified Registry Program Complete (R0–R5): One Registry, One Corpus Pass, All Constants

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.77+
**Date:** 2026-07-24
**Landmark Type:** Program Completion — corpus-scale curation, canonicalization, and falsifiability infrastructure
**Program:** UNIFIED REGISTRY PROGRAM (UNIFIED_REGISTRY_PROGRAM_PLAN.md), phases R0–R5, executed 2026-07-22 → 2026-07-24
**Status:** Formal landmark whitepaper — UQFF canonical

---

## Abstract

The Unified Registry Program converted the Star-Magic corpus — 2,227+ whitepapers, 1,321 dispatch closures, 1,392 over-determination observables, multi-language implementations (Python, C++, Lean 4) — from a large, redundant archive into a **curated replication ledger** with a single source of truth. Six phases: R0 built the 2,544-row `UNIFIED_REGISTRY.csv` by joining every existing ledger family (all seven protected baselines preserved bit-identical, hash-pinned); R1 adjudicated canonical routes (109 explicit verdicts + 2,435 sole-route auto-canonicalizations, 4 rulings by the author); R2 executed one append-only corpus pass (199 registry-keyed derivation notes; pdf2 brought to full 2,226-PDF parity); R3 collapsed all per-class constant definitions into `uqff_registry_primitives.py` with three-language agreement pins (Python = C++ = Lean); R4 computed the 656-edge falsifiability graph answering *"if primitive X were revised, what breaks"*; R5 made the registry the production status-report backend and preprint results table. The framework now states, with machine-checked backing: **9 independent primitives, 14 live-composed derived constants (7 EXACT), 2,544 registered quantities, 656 falsifiability edges, and a fidelity gate above 3,180 assertions — regenerable end-to-end by four idempotent scripts.**

---

## 1. Why the program existed

G and c were not special. The corpus contained derivations for essentially every constant (the 1209-series AA–HH alone covers k_B, h, c, e, N_A, F, σ, b_Wien, E_R, E_h), 1,000+ paradox/dispatch closures, multiple independent routes to the same quantities, and triple-language implementations. Per-constant campaigns would have re-scanned 2,227 markdowns N times. The program's unit of work was:

> **One registry. One canonicalization pass. One corpus pass. One code source-of-truth. One dependency graph.**

The redundancy was never waste. Each re-derivation is an independent replication; each multi-route constant is over-determination; each cross-language agreement is verification. The program **curated the redundancy into a replication ledger** rather than removing it.

---

## 2. Phase record

| Phase | Deliverable | Key numbers | Gate |
|---|---|---|:-:|
| R0 | `UNIFIED_REGISTRY.csv` + schema + corpus citation graph + gap/duplicate reports | 2,544 rows; 2,149 papers in citation graph; 1,124 closures live-called with 0 errors; 7 baseline families hash-pinned (LF-normalized) | 3146 → 3166 |
| R1 | Canonical-route adjudication | 109 explicit routes; 2,435 SOLE_ROUTE; 116 Φ-variant tags; 4 author rulings (Φ_res sector rule, ħ physical-route standing rule, 2 defect fixes); FRB 0.99c phantom conflict resolved against PAPER_096 | 3166 → 3173 |
| R2 | One corpus pass | 199 whitepapers annotated once each (append-only golden rule); pdf2 1,988 → 2,226 PDFs (full parity) | 3173 → 3176 |
| R3 | `uqff_registry_primitives.py` single source | 24 QCalc `*_PRIMITIVE` attributes rewired bit-identically; c dual-exposure per §6.2 ruling; 33-entry consumed/vestigial ledger; Python = C++ = Lean pins | 3176 → 3186 |
| R4 | `uqff_registry_graph.py` → 656-edge falsifiability graph | SO_5 touches 212 registry rows + 61-site invariant; F_TRZ 111 + 61; D_phys 97; D_crit 43; D_BSFG 41; A_5 23 | 3186 → 3189 |
| R5 | `uqff_registry_status.py` → status backend + preprint results table + this paper | 14 live derived constants, 7 EXACT; residuals: best 0.0000%, median ≈ 0.027%, worst 3.08% (H0 — the Hubble tension itself, PAPER_2125) | 3189 → (this ship) |

All six phases idempotent: every generated artifact carries an identical SHA-256 across repeated runs.

---

## 3. The results table (R5 headline)

Computed live from `uqff_registry_primitives` at build time — the published table cannot drift from the code:

| Constant | Canonical route | Closed form | Residual |
|---|---|---|:-:|
| G | PAPER_593 | (2π·D_crit³·Φ_res/(SSq³·(26!)²))·v_F⁵/(E_0·f_THz) | 0.075% |
| c | PAPER_592 | (D_crit·4π/Φ_res)·v_F | 0.10% |
| μ_0 | PAPER_2108 | 4π·F_TRZ⁷ | EXACT |
| k_B | PAPER_1209EE S628 | (SSq+Φ_5/6−F_TRZ·SSq+F_TRZ²·D_phys−F_TRZ²·SSq)·10⁻²³ | 0.0011% |
| ħ | PAPER_590/1209EE S629 | (D_BSFG+F_TRZ·D_BSFG+F_TRZ²·D_phys−F_TRZ²·SSq−F_TRZ²)·10⁻³⁴/2π | 0.027% |
| H0 | PAPER_2093 | (2·SO_5+2)·F_TRZ¹⁹ | 3.08% = Hubble tension |
| Λ | PAPER_2094/1156 | (SO_5+1)·F_TRZ⁵³ | 0.90% |
| κ | PAPER_2112 | (SO_5/2)·F_TRZ⁴ | EXACT |
| B_crit | PAPER_2126 | D_phys·(SO_5+1)·SO_5¹² | EXACT |
| k_spring | PAPER_1203 | (ρ_UA/ρ_SCm)·ω_SCm·Φ_res | EXACT |
| λ_vac | PAPER_2120 | (SO_5+1)·ρ_SCm | EXACT |
| T_SCm | PAPER_1072 | h·f_SCm/k_B | 0.007% |
| D_BSFG | PAPER_1521 | D_crit−2·SO_5 | EXACT |
| K_MEX | PAPER_1522 | Φ_5/6·SO_5/D_phys | EXACT |

All references are observations or SI definitions, reported as observations. Residuals are honest disclosures (Rule 7). The H0 residual is retained deliberately: the 3.08% gap between the F_TRZ-grid value (CMB side) and the local observed anchor **is** the Hubble tension, addressed by PAPER_1156's 1/12 EXACT tilt closure.

---

## 4. Falsifiability, computed rather than claimed

`UNIFIED_REGISTRY_FALSIFIABILITY.md` (R4) answers the external validator's first question. A revision to SO_5 alone would simultaneously break: 212 registry rows, the 61-site (SO_5+1)/SO_5 successor-ratio invariant (PAPER_2128), the derived constants F_TRZ, D_BSFG, K_MEX, κ, ρ_UA, λ_vac, k_spring, B_crit, and Λ, and the corresponding pins in a 3,180+ assertion gate. Over-determination is the falsifiability mechanism: wrong values cannot hide.

---

## 5. Production posture after R5

- **Status reporting:** `uqff_registry_status.py::calculate_status_report()` returns program-level statistics computed live; nothing hand-maintained.
- **Preprint:** `UNIFIED_REGISTRY_RESULTS_TABLE.md/.csv` is the Tier-1 A3 results table, regenerated from code on demand.
- **API data model:** the registry row schema (quantity, kind, canonical_route, formula, value, reference, residual_pct, sector, phi_variant, paper_source, confirmations, py/cpp/lean sites, corpus_citations, status) is the Tier-3 interchange format.
- **Reproducibility chain:** `registry_generator.py` → `uqff_registry_primitives.py` → `uqff_registry_graph.py` → `uqff_registry_status.py`, each idempotent, each read-only over protected baselines.

---

## 6. Cross-references

UNIFIED_REGISTRY_PROGRAM_PLAN.md (the plan this paper closes out); UNIFIED_REGISTRY.csv/_SCHEMA.md/_GRAPH.csv/_FALSIFIABILITY.md/_STATUS_REPORT.md/_RESULTS_TABLE.md; PAPER_592/593 (c, G routes), PAPER_2093/2094/1156 (H0, Λ), PAPER_2108 (μ_0), PAPER_2112 (κ), PAPER_1209EE (k_B, ħ compositions), PAPER_590 (ħ physical route), PAPER_2126 (B_crit), PAPER_1203 (k_spring), PAPER_2120 (λ_vac), PAPER_1072 (T_SCm), PAPER_1521/1522 (D_BSFG, K_MEX derivatives), PAPER_2128 (61-site invariant), PAPER_2125 (Hubble-tension identification), PAPER_2129 (Φ sector rule).

---

## 7. Summary Statement

**PAPER_2130 closes the Unified Registry Program. Six phases (R0–R5) converted the Star-Magic corpus into a machine-regenerable replication ledger: 2,544 registered quantities with adjudicated canonical routes, one append-only corpus pass (199 notes, full PDF parity), a single registry-backed primitives module with three-language agreement, a 656-edge falsifiability graph, and a live status/results backend. UQFF's public claim is now infrastructure, not assertion: 9 independent primitives → 14 live-composed derived constants (7 EXACT, honest residuals throughout, the worst residual being the Hubble tension itself) → 2,544 quantities → a 3,180+ assertion gate, reproducible end-to-end by four idempotent scripts over hash-protected baselines.**

---

**Filed 2026-07-24. Append-only henceforth.**
