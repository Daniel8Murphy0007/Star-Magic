# G/c CORPUS-WIDE UPGRADE PLAN — UQFF Gravity + Speed-of-Light Derivations

**Author:** Daniel T. Murphy (assembled per session 2026-07-22 directive)
**Status:** PLAN — awaiting Daniel's phase-gate approvals before execution
**Scope:** whitepapers/ (2,227 .md) audit → pdf2/ (1,987 .pdf) upgrade → code parity → gate expansion
**Prime directive:** UQFF IS THE ANCHOR (Rule 4). Observed values are observations, not SM. Honest residuals only (Rule 7). Append-only on published papers (Rule 9 discipline extended to whitepapers).

---

## 0. The Canonical Anchors (fixed reference for every phase)

Both derivations are production-status, parameter-free, audited (Sessions 237-240), and wired:

```
PAPER_592 — Speed of Light from Pre-Mass Triad Equilibrium
  c_UQFF = (26·4π/Φ_res)·v_F = 2.995e8 m/s          (0.13% vs observed)
  third SI anchor: v_F = 0.77e6 m/s Fermi-velocity proxy (INDEPENDENT of c;
  dpm_vacuum_manifold.py lines 3701/4896/5224 — circularity closed Session 239)

PAPER_593 — Gravitational Constant from Void Coupling
  G_UQFF = (2π·26³·Φ_res/(SSq³·(26!)²))·v_F⁵/(E_0·f_THz)
         = 6.668991909557279e-11                      (0.08% vs observed)
  SI-clean anchors: E_0 = 1e-20 J, f_THz = 1.25e12 Hz, v_F, + 26! barrier

Wired at:  uqff_pure_calculator.py lines 20-24 (module-level C_LIGHT, G_NEWTON)
           _uqff_primitives derive_G_newton() / derive chain
           CondensedPhysics.py QCalc classes R364-R373 (G promoted 2026-07-22)
Gate pins: BUCKET 0 (c within 0.5%, G within 0.5%) + G-PROMOTION block (bit-identity)
```

Corpus recon (2026-07-22):

| Metric | Count |
|---|:-:|
| whitepapers/*.md | 2,227 |
| pdf2/*.pdf | 1,987 |
| .md files containing G literals (6.674e-11 / 6.6743e-11 / 6.67e-11) | 168 |
| .md files containing c literals (3e8 / 2.998e8 / 2.995e8 / 299792458) | 61 |
| .md files already citing PAPER_592/593 | 36 |

Estimated adjudication load: ≤ 229 files (union of G+c sets, overlap expected), ~10% of corpus.

---

## 1. Classification Taxonomy (applied to every literal occurrence)

Every G/c literal occurrence in every file gets exactly one class:

| Class | Meaning | Action |
|:-:|---|---|
| **A** | Canonical derivation source (PAPER_592, PAPER_593, their audit notes, PAPER_590/591 siblings) | **NO CHANGE** — these are the anchors |
| **B** | Cites the derivation correctly (references PAPER_592/593, uses derived value or notes residual) | **NO CHANGE** — already compliant |
| **C** | Uses CODATA literal as **input** to a UQFF formula/derivation chain | **UPGRADE** — annotate with UQFF closed form + residual; flag if downstream numbers shift |
| **D** | Uses literal as **observational anchor** for residual comparison (measured G/c as the observation UQFF is tested against) | **NO CHANGE** — observations are observations, not SM (Rule 4 clarification); optionally tag "anchor: observed" |
| **E** | Documents calculator internals that carried the literal (pre-promotion code listings, class references) | **UPGRADE** — append derivation note matching the 2026-07-22 code promotion |
| **F** | Ambiguous / mixed usage | **ESCALATE to Daniel** — no automated action |

**Golden rule for C-class:** published downstream numbers are NEVER silently recomputed. The upgrade appends a standardized derivation note stating the UQFF-derived value, the residual, and (where material) the % shift the derived value would produce. Recomputation of a published result happens only on Daniel's explicit per-paper instruction.

**Standardized note block (appended to C/E-class files):**

```
---
## G/c DERIVATION NOTE (appended 2026-07-XX, corpus-wide upgrade)
This paper uses [G = 6.674e-11 | c = 3e8 | ...] at [locations]. Per PAPER_593
[/PAPER_592], the UQFF parameter-free derivation gives
[G_UQFF = 6.66899e-11 (0.08%) | c_UQFF = 2.995e8 (0.13%)]. The value above is
retained as [the observational anchor | the original published input]; the UQFF
derivation is canonical. [Downstream shift if adopted: X.XX% on quantities Q.]
---
```

---

## 2. Phase Plan

### PHASE 0 — Freeze + Backup (½ session)
- `git status` clean checkpoint; tag or note current HEAD in SESSION_LOG.
- Backup manifests: `whitepapers_G_C_PREUPGRADE_FILELIST.txt` (name+size+mtime of all 2,227 md and 1,987 pdf). No file copies needed — git is the backup; the manifest is the tamper-evidence.
- Confirm gate baseline (3138/0) recorded.
- **Deliverable:** frozen baseline note in SESSION_LOG.

### PHASE 1 — Automated Inventory (1 session)
- Script `audit_g_c_corpus.py` (repo tools/, or outputs scratch) scans whitepapers/*.md for the literal families:
  - G: `6.674e-11`, `6.6743e-11`, `6.67e-11`, `6.67430`, `6.669e-11`, `6.66899e-11`
  - c: `3e8`, `3.0e8`, `2.998e8`, `2.995e8`, `299792458`, `299,792,458`, `c^2 =`, `c² =`
  - context window ±2 lines per hit
- Auto-classifies obvious cases: A (filename match), B (PAPER_592/593 cited within 10 lines), D (context contains "anchor", "observed", "residual", "CODATA", "vs measured").
- Emits `G_C_AUDIT_MANIFEST.csv`: `file, line, literal, class_auto, context, action_proposed`.
- Emits summary: counts per class, list of all F-class (ambiguous) rows.
- **Deliverable:** manifest CSV + summary. NO file modifications in this phase.
- **GATE CHECK:** manifest row count pinned in a new gate assertion (audit reproducibility).

### PHASE 2 — Adjudication (Daniel, 1-2 sessions)
- Daniel reviews the manifest — primarily the C and F rows (estimated ≤ 80 files).
- Per-row verdict recorded in a `class_final, recompute_yn` column.
- **Decision points requiring Daniel's explicit call (see §4).**
- **Deliverable:** adjudicated manifest. This is the sole authority for Phase 3.

### PHASE 3 — Markdown Upgrades (2-3 sessions, batched)
- Batches of ~25 files. Python splice scripts only (Edit-tool truncation risk on large files per CLAUDE.md warning).
- C-class: append standardized derivation note; inline `[G/c → see Derivation Note]` markers only where Daniel approved.
- E-class: append note referencing the 2026-07-22 code promotion (compute-don't-store).
- D-class (if tagging approved): add one-line `anchor: observed` tag; else untouched.
- After each batch: gate run (must stay green), batch log appended to SESSION_LOG.
- **Deliverable:** upgraded whitepapers, per-batch logs, zero gate regressions.

### PHASE 4 — pdf2 Regeneration (2-4 sessions, batched)
- For every markdown upgraded in Phase 3 that has a pdf2 counterpart: rebuild the PDF **in place, same filename**, now carrying the derivation note.
- Mapping pass first: `md ↔ pdf2` filename join (expect imperfect mapping — 2,227 md vs 1,987 pdf; emit `PDF2_MAPPING_GAPS.txt` for papers with md but no pdf and vice versa — Daniel decides whether gaps get new PDFs).
- Rebuild pipeline: reportlab, standard house style (Helvetica/Courier, color tables, subscript/superscript XML tags — never Unicode sub/superscripts), `_build_log.txt` appended per batch.
- Batch size ~50 PDFs; verify each batch: file count parity, spot-open 3 PDFs per batch, size sanity (no zero-byte).
- **Deliverable:** upgraded pdf2 folder, mapping-gap report, build logs.

### PHASE 5 — Code Parity Sweep (1-2 sessions)
- Extend the R364-R373 G-promotion to ALL earlier R218+ fills still carrying `G_PRIMITIVE = 6.674e-11` or `self.G = 6.674e-11` (R329-R347 era: CompressionUg1/Ug3, M51 suite, M31 suite, UFEUgGravityMode, ~12-15 classes) — same live closed form.
- c-promotion decision (see §4.2) applied uniformly: all `C_PRIMITIVE = 3e8` / `self.c = 3e8` sites.
- Consumed-vs-vestigial audit per class (dpm_ug1_seed G-free finding): tag which classes actually consume G/c, so instance counts in future papers count consumption, not decoration.
- **Deliverable:** corpus-uniform code; consumed/vestigial ledger.

### PHASE 6 — Gate Expansion (½ session)
- New assertions: (a) zero remaining raw `= 6.674e-11` assignments in CondensedPhysics calculator classes (regression guard, allowlist for D-class anchor dicts); (b) same for `= 3e8` per the §4.2 decision; (c) manifest reproducibility pin; (d) derived-G/derived-c bit-identity pins against uqff_pure_calculator module values.
- **Deliverable:** gate ~3150+, exit 0.

### PHASE 7 — Landmark Documentation + Ship (1 session)
- PAPER_2129 (or next free number): "Corpus-Wide G/c Derivation Upgrade — PAPER_592/593 Propagation Across 2,227 Whitepapers + 1,987 PDFs" — full statistics, class-by-class counts, consumed/vestigial ledger, honest-residual tables.
- SESSION_LOG entry; CHANGELOG entry; version bump (candidate v5.74.0); ship staging per the established PowerShell master-branch procedure (defensive index.lock removal, `git push origin master`).
- **Deliverable:** shipped release with the upgrade as its headline.

---

## 3. Sequencing & Effort Summary

| Phase | Effort | Blocking dependency |
|:-:|---|---|
| 0 Freeze | ½ session | none |
| 1 Inventory | 1 session | Phase 0 |
| 2 Adjudication | 1-2 sessions (Daniel) | Phase 1 manifest |
| 3 Markdown | 2-3 sessions | Phase 2 verdicts |
| 4 pdf2 | 2-4 sessions | Phase 3 per-file completion |
| 5 Code parity | 1-2 sessions | §4 decisions (can parallel Phase 3/4) |
| 6 Gate | ½ session | Phases 3-5 |
| 7 Paper + ship | 1 session | all |

Total: ~8-13 sessions. Phases 3-5 are batch-resumable — safe to interleave with continued stub-fill rounds if desired.

---

## 4. DECISION POINTS — Daniel's calls before execution

**4.1 — G value policy for D-class (observational anchors).** Keep 6.674e-11 as "observed anchor" untouched (recommended; observations are observations), or tag every instance with `anchor: observed`?

**4.2 — c promotion policy in code.** Three options:
  (a) **Full promotion:** `C_PRIMITIVE = (26·4π/0.84)·0.77e6 = 2.9947e8` everywhere — pure UQFF, shifts c-consuming outputs ~0.18%;
  (b) **dpm precedent:** source from `dpm_vacuum_manifold._C_LIGHT = 2.998e8` (R226 MasterBuoyant precedent, 0.07% from prior literal);
  (c) **Dual exposure:** keep 3e8 as `C_OBSERVED` + add `C_DERIVED_PRIMITIVE` live form, closures consume the derived one only on your instruction.
  Recommendation: (a) for vestigial/unconsumed sites immediately, your call per-class where consumed.

**4.3 — Recompute policy.** Which (if any) C-class papers get downstream numbers recomputed with derived G/c, vs. note-only? Default: note-only, append-only.

**4.4 — pdf2 mapping gaps.** Papers with md but no PDF (and orphan PDFs): build missing PDFs during Phase 4, or defer?

**4.5 — Interleaving.** Pause the R374+ stub-fill campaign during the upgrade, or interleave (batches between rounds)?

---

## 5. Risk Register

| Risk | Mitigation |
|---|---|
| Edit-tool truncation on large md/py files | Python splice scripts exclusively for batch edits (CLAUDE.md protocol) |
| Silent alteration of published numbers | Golden rule: append-only notes; recompute only on per-paper instruction (§4.3) |
| Gate regression mid-batch | Gate run after every batch; batch size small enough to bisect |
| md↔pdf mapping mismatches | Phase 4 mapping pass with gap report before any rebuild |
| Instance-count contamination (vestigial vs consumed) | Phase 5 ledger; future papers count consumption |
| Fastly/PyPI cache confusion at ship | Verify via `gh run view` job status, not the PyPI page (v5.72.0 lesson) |
| Session interruption mid-phase | Every phase leaves resumable artifacts (manifest, batch logs); SESSION_LOG entry per batch |

---

## 6. Verification Matrix (what "done" means)

1. Every one of the ≤229 affected md files carries exactly one class; C/E files carry the derivation note; A/B/D files provably untouched (manifest diff).
2. Every upgraded md's pdf2 counterpart rebuilt; mapping-gap report resolved per §4.4.
3. Zero raw CODATA G assignments in calculator classes outside the D-class allowlist; c per §4.2.
4. Gate exit 0 with new regression guards; bit-identity pins on derived G/c.
5. PAPER_2129 filed with full statistics; SESSION_LOG + CHANGELOG entries; ship green.

---

**This plan executes nothing by itself. Phase 1 (read-only inventory) can start on your word; Phases 2+ wait on the manifest and your §4 decisions.**
