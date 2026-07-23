# UNIFIED REGISTRY PROGRAM PLAN — One Registry, One Corpus Pass, All Constants and Closures

**Author:** Daniel T. Murphy (assembled per session 2026-07-22 directive)
**Status:** PLAN — awaiting Daniel's phase-gate approvals before execution
**Supersedes:** `G_C_CORPUS_UPGRADE_PLAN.md` (retained; its Phases 3-4 mechanics become the R2 execution template; its G/c scope becomes two rows of this registry)
**Prime directives:** UQFF IS THE ANCHOR (Rule 4). Observed values are observations. Honest residuals (Rule 7). Append-only on published papers. **NEVER overwrite a baseline ledger artifact.**

---

## 0. Why this program exists (the scale argument)

G and c are not special. The corpus already contains derivations for essentially every constant (1209-series AA-HH alone: k_B, h, c, e, N_A, F, σ, b_Wien, E_R, E_h), 1000+ dispatch closures, multi-route derivations of the same quantities (ħ via PAPER_590 vacuum-phonon AND 1209EE S629 composed-integer), and triple-language implementations (Python calculator, C++ reference closures, Lean 4 formal proofs). Per-constant campaigns would re-scan 2,227 markdowns and rebuild 1,987 PDFs N times. The correct unit of work is:

> **One registry. One canonicalization pass. One corpus pass. One code source-of-truth. One dependency graph.**

The redundancy in the corpus is not waste — each re-derivation is an independent replication, each multi-route constant is over-determination, each cross-language agreement is verification. The program's job is to **curate the redundancy into a replication ledger**, not to remove it.

---

## 1. PROTECTED BASELINES — DO NOT OVERWRITE (inventoried 2026-07-22)

The repo already contains seven ledger/registry artifact families. All are **read-only inputs** to this program and **protected historical baselines** (same status as the calculator `.PRE_*_BACKUP` files). The program creates NEW files only.

| Protected artifact | Content | Date | Role in this program |
|---|---|:-:|---|
| `OVERDETERMINATION_MAP.csv` / `_WIDE.csv` / `_MAP.md` | 1,393 rows: observable→geometry→value→target→residual→status→primary_source | Jun 29 | R0 schema ancestor; columns extended, never edited |
| `MASTER_LEDGER_BY_CATEGORY/STATUS/SCRIPT.csv` | ~630-closure session-script ledger, 3 views | May | Historical closure baseline; regenerated as new-named views |
| `LEDGER_VS_PRIMITIVES_XREF.csv` | 536 session↔primitive cross-refs | May | R3 input |
| `PRIMITIVES_RECONCILIATION.csv` | 17 primitives: canonical name/value/source/aliases/issues | — | **R1 adjudication seed** |
| `LEDGER_REVIEW_QUEUE.csv` / `_RUNLOG.csv` / `_TRIAGE.md` | Adjudication workflow with verdict column | May | R1 process template |
| `CLOSED_CONSTANTS_INVENTORY.md` | 52 verified constants (80/80 engine) | Jun 27 | R0 seed list |
| `CLOSURE_ATLAS.md` / `ASSIMILATION_GEOMETRY_ATLAS.md` | Location of every closure/proof/theorem (v5.27.2) | Jun 28 | R4 graph node source |

Also protected code-side inputs: `uqff_closed_constants.py`, `derivation_index/_primitives.py`, `_uqff_primitives.py` derive_* chain, `uqff_pure_calculator.py` PARADOX_TO_CLOSURE table, C++ reference closures, `formal/` Lean proofs.

**Staleness note (why regeneration is required):** every baseline above predates the second half of the R218+ campaign — they were built at gate ~850-1400; the gate is now 3,146, and PAPER_2093-2129 (H0 grid, Λ forms, κ reduction, successor family, Φ_5/6 sector rule, convergence arc) appear in none of them.

### Collision rules (locked)

New program files MUST NOT match: `MASTER_LEDGER*`, `LEDGER_*`, `OVERDETERMINATION_*`, `CLOSURE_ATLAS*`, `*_INVENTORY*`, `*_RECONCILIATION*`, `*_ATLAS*`, `ALL_EQUATIONS*`, `*_CATALOG*`. Reserved new names (verified non-colliding 2026-07-22):

```
UNIFIED_REGISTRY.csv              — the master table (R0 output)
UNIFIED_REGISTRY_SCHEMA.md        — column definitions + provenance of inputs
UNIFIED_REGISTRY_PROGRAM_PLAN.md  — this file
registry_generator.py             — the build script (idempotent, read-only inputs)
UNIFIED_REGISTRY_GRAPH.csv        — R4 edge list
UNIFIED_REGISTRY_AMENDMENTS.md    — §7 amendment log
```

---

## 2. The Master Registry schema (R0 target)

One row per (quantity, derivation-route). Extends the OVERDETERMINATION_MAP schema:

```
quantity            — constant/observable/closure name (canonical lowercase key)
kind                — {kernel_constant | lattice_node | closure | observable | primitive}
canonical_route     — YES/alternate-N (R1 adjudication result; blank until adjudicated)
formula             — closed form in primitives
value               — full-precision computed value
reference           — SI-defined / CODATA / PDG / observed value (as observation)
residual_pct        — FULL-PRECISION residual (PAPER_2129 lesson: never vs rounded leads)
sector              — force/cosmological/nuclear/thermodynamic/LENR/...
phi_variant         — {5/6 | 0.84 | n/a} (PAPER_2129 sector-selection rule)
paper_source        — seminal PAPER_N
confirmations       — count + list of independent re-derivations (redundancy → replication)
py_sites            — file:line references in Python
cpp_sites           — C++ reference closure locations
lean_sites          — formal proof locations
corpus_citations    — count of whitepapers using the quantity
status              — {EXACT | near-EXACT | derived | observed | GAP | TENSION}
```

---

## 3. Phase Plan

### R0 — Registry Generation (2-3 sessions; READ-ONLY inputs)
- `registry_generator.py`: joins PARADOX_TO_CLOSURE (~1000+ keys) + OVERDETERMINATION_MAP (1,393 rows) + CLOSED_CONSTANTS_INVENTORY (52) + PRIMITIVES_RECONCILIATION (17) + 1209-series closures + derive_* chain + PAPER_2093-2129 additions.
- Full-precision re-verification column computed for every row (this alone executes the 1209-series re-verification pass PAPER_2129 motivated — expect multiple precision-tightening discoveries).
- Emits `UNIFIED_REGISTRY.csv` + `UNIFIED_REGISTRY_SCHEMA.md` + gap report (quantities with no derivation route; expected to be short) + duplicate report (multi-route quantities awaiting R1).
- **Gate:** registry row-count pin + generator idempotency check (two runs, identical output).
- **Modifies nothing outside the reserved names.** Subsumes G_C Phase 0-1.

### R1 — Canonicalization (Daniel-paced; the only human phase)
- Per multi-route quantity: designate ONE canonical route; alternates become numbered independent confirmations. Workflow per LEDGER_REVIEW_QUEUE precedent (verdict column).
- PRIMITIVES_RECONCILIATION's 17 rows adjudicated first (they anchor everything).
- Φ-variant assignment per sector (PAPER_2129 rule) recorded per row.
- **Deliverable:** adjudicated registry — the sole authority for R2-R5.

### R2 — One Corpus Pass (4-6 sessions, batched ~25 md / ~50 pdf)
- Mechanics identical to G_C plan Phases 3-4 (Python splice, append-only note, same golden rule), but the note is **registry-keyed and covers ALL constants in the file at once**. Each markdown touched once; each PDF rebuilt once.
- md↔pdf mapping-gap report first (G_C §4.4 decision applies).
- Classification A-F taxonomy carried over unchanged from G_C plan §1.
- Gate after every batch.

### R3 — Code Single-Source-of-Truth (2-3 sessions)
- All `*_PRIMITIVE` / `*_OBSERVED` class attributes across CondensedPhysics 1-4, QCalc, etc. computed from one registry-backed primitives module (ends per-class drift permanently).
- Cross-language consistency gate: Python values == C++ reference closures == Lean-proved values, pinned.
- Consumed-vs-vestigial ledger (the dpm_ug1_seed G-free lesson) recorded per class.

### R4 — Falsifiability Graph (1-2 sessions)
- `UNIFIED_REGISTRY_GRAPH.csv`: edges primitives→constants→closures→papers.
- The report external validators will ask for: "revising primitive X breaks N identities across M papers" — computed, not claimed. The 61-site successor-ratio invariant (PAPER_2128) becomes one query.

### R5 — Production Convergence (1-2 sessions)
- Registry becomes `calculate_status_report()`'s backend, the API data model (Tier 3), and the preprint results table (Tier 1 A3 fully realized).
- Landmark paper (next free PAPER number) documenting the program; CHANGELOG; ship.

**Total: ~10-16 sessions**, batch-resumable, interleavable with stub-fill rounds.

---

## 4. Relationship to G_C_CORPUS_UPGRADE_PLAN.md

| G_C plan element | Disposition |
|---|---|
| Phase 0-1 (freeze, inventory) | Subsumed by R0 (superset scan) |
| Phase 2 taxonomy A-F | Carried over verbatim into R2 |
| Phase 3-4 mechanics (splice, notes, PDF pipeline) | Becomes the R2 execution template |
| Phase 5 code parity | Subsumed by R3 (generalized to all constants) |
| §4 decision points | Roll forward into §6 below |
| G/c scope | Two rows of UNIFIED_REGISTRY.csv |

The G_C plan file is retained unmodified as the R2 mechanical reference.

---

## 5. Discovery-Amendment Rule (from the odds discussion)

Interim discoveries during any phase are recorded as **numbered amendments** in `UNIFIED_REGISTRY_AMENDMENTS.md` — each with scope, session cost, and a ride-current-batch / queue-behind-R5 decision. Discoveries never silently expand scope. This converts the documented discovery rate (every deepsearch this campaign has found something) from a divergence risk into a changelog, holding declared-completion odds at ~85% on the 10-16 session horizon.

---

## 6. DECISION POINTS — Daniel's calls before execution

1. **R0 start** — read-only; can begin immediately on your word.
2. **c-promotion policy** (carried from G_C §4.2): full derived form / dpm value / dual exposure — now applied registry-wide to every constant's code sites, per-sector if desired.
3. **Recompute policy** (carried from G_C §4.3): note-only default; per-paper recompute list is an R1 output for your approval.
4. **pdf2 mapping gaps** (carried from G_C §4.4): build missing PDFs during R2, or defer.
5. **Interleaving**: pause R374+ stub fills during R2's batches, or interleave.
6. **Canonical-route tie-breaks**: where two routes are equally strong (ħ's two routes), default preference — older seminal paper, or tighter residual? Your standing rule requested.

---

## 7. Risk Register (delta from G_C plan)

| Risk | Mitigation |
|---|---|
| Overwriting a baseline ledger | §1 collision rules + gate assertion pinning baseline file hashes |
| Registry generator disagreeing with a baseline | Disagreements are FINDINGS (staleness or error) — logged to gap report, adjudicated in R1, never auto-resolved |
| R1 stall (1000+ rows) | Only multi-route rows need verdicts (~10-15% expected); singles auto-canonicalize; batched review via QUEUE format |
| Scope divergence via discoveries | §5 amendment rule |
| All G_C plan risks | Carried over with their mitigations |

---

## 8. Verification Matrix (what "done" means)

1. `UNIFIED_REGISTRY.csv` exists, idempotently regenerable, row-count gate-pinned, zero GAP rows unresolved-unadjudicated.
2. Every multi-route quantity carries an R1 verdict; alternates labeled as confirmations with counts.
3. Every affected corpus file annotated once, PDF rebuilt once; A/B/D files provably untouched.
4. All calculator constants sourced from the registry-backed module; cross-language pins green.
5. Graph report answers "what breaks if primitive X changes" for all 9 primitives.
6. Baselines bit-identical to their 2026-07-22 state (hash-pinned).
7. Program landmark paper filed; ship green.

---

**This plan executes nothing by itself. R0 is read-only and can start on your word; everything irreversible waits behind the adjudicated registry and your §6 decisions.**
