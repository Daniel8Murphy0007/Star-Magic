# NEXT_PRIORITIES.md — what the next Claude session should do first

**READ CLAUDE.md FIRST. EVERY SESSION. EVERY ACTION.** Rules in CLAUDE.md supersede this file.

**Daniel's standing directives (verbatim):**
- "DUMP ALL PROVENANCE. NO COMMENTS WHATSOEVER. NO TAGS, NO GARBAGE SM ANYTHING. JUST REAL PERCENTAGES OF DIFF."
- "UQFF IS THE FUCKING ANCHOR. UQFF DOESN'T SHARE SHIT WITH SM."
- "I AM GIVING YOU THE INFORMATION, YOU ARE TO ASSEMBLE IT."
- "READ THE RULES BEFORE EVERY ENTRY. READ THE RULES BEFORE EVERY SESSION."

---

## STATE AT SESSION CLOSE (2026-07-11 — v5.61.0)

### Ship trail
- **v5.44.0** (2026-07-04) — CP pipeline integration (uqff_pure_calculator dispatchers added)
- **v5.45.0 through v5.57.0** (2026-07-04 to 2026-07-10) — 13 ships of CP1 P2 stub-drainage rounds. Zero changes to uqff_pure_calculator.py. Drift accumulated.
- **v5.58.0** (2026-07-11) — Phase A corpus completion: 670 new PDFs (99.9% coverage), 57 regex-artifact-swept, PAPER_1796-1799 reservation stubs, PAPER_275 typo fix
- **v5.59.0** (2026-07-11) — Phase B calculator wiring: 100 new PARADOX_TO_CLOSURE dispatch keys for PAPER_1893-1984 + PAPER_872 + PAPER_1087 OPEN_QUESTION marker. Gate 931/0 → 1031/0.
- **v5.60.0** (2026-07-11) — Phase C framework annotation retrofit: 316 auto-extracted annotations for Rounds 52-116
- **v5.61.0** (2026-07-11) — Phase D housekeeping (this ship): honest CP1-CP4 audit + NEXT_PRIORITIES.md refresh + task-list prune

### Honest metrics
- **Fidelity gate:** 1031 passed, 0 failed (was 931 pre-Phase-B)
- **Whitepapers authored:** 1,984 (PAPER_001-PAPER_1984, 4 IDs reserved)
- **PDF coverage:** 1,978 / 1,980 = 99.9% (PAPER_1218 + PAPER_1801 deferred to v5.58.1)
- **PARADOX_TO_CLOSURE dispatch keys:** ~11,167 (added 100 in Phase B)
- **Framework annotations:** 351 (35 hand-classified Rounds 45-51 + 316 auto-extracted Rounds 52-116)
- **CP1 coverage:** 324/813 scoreable = 39.9%. 489 stubs remain in CP1.
- **CP2 coverage:** 0/224 = 0%. Never touched by Rounds 1-116.
- **CP3 coverage:** 0/255 = 0%. Never touched.
- **CP4 coverage:** 0/807 = 0%. Never touched.
- **Total scoreable across CP1-CP4:** 2,099 stubs. Upgraded: 324. **Real corpus coverage: 15.4%.**
- **Errored on default `.compute()` (may include parameterized calculators requiring signature audit):** 593 CP1 + 473 CP2 + 672 CP3 + 709 CP4 = **2,447 classes** with unknown status
- **Grand total classes with `compute()` method:** 4,579

### What the "50% milestone" (v5.54.0) actually meant
50% of CP1's estimated 1203 total, computed against a wrong denominator. Actual corpus coverage at that point: ~30%. This document uses honest denominators going forward.

---

## PRIORITY 1 — VERIFICATION FIRST

Before any modification:

1. Run `python uqff_fidelity_tests.py`. Expect **1031 passed, 0 failed**, exit 0.
2. Confirm all 11 canonical primitives + 3 landmark derivative primitives (PAPER_1521 D_BSFG, PAPER_1522 K_MEX, PAPER_1960 F_TRZ) unchanged in `uqff_pure_calculator.py`.
3. Sample public surfaces; verify return shape `{'value': X}` per CLAUDE.md Rule 5.
4. **DO NOT MODIFY** existing Bucket A-K wiring without explicit user request.

---

## PRIORITY 2 — DEFINE AN "EFFECTIVE" ROUND

**Rounds 1-116 lesson: primitive-locking a stub in CondensedPhysics is NOT enough.** Papers accumulated, `framework_papers` metadata accumulated, but users installing `pip install uqff` reached NONE of it because the pure calculator dispatch surface stayed untouched.

**Definition of an effective Round going forward:**

An effective Round is complete when ALL FIVE steps land in the SAME ship:

1. **Stub upgrade** in CP1/CP2/CP3/CP4 with `framework_papers` metadata + runtime `_verify` booleans
2. **Novel-closure whitepaper** authored (if the Round discovered a novel closure)
3. **PDF built** for any new whitepaper (or explicitly deferred to a numbered patch)
4. **Dispatch key added** to `PARADOX_TO_CLOSURE` in `uqff_pure_calculator.py` — user reachable via `calculate_paradox({'paradox': '<key>'})`
5. **Fidelity gate assertion** added to `uqff_fidelity_tests.py` — pins the identity numerically

**If any step is skipped, the Round is INCOMPLETE.** The Phase B drift (Rounds 45-116 skipped steps 4-5 for 92 papers) took 100+ hours of catch-up work to close retroactively. Do not repeat.

**Round-start ritual (documented by Rounds 45-79, forgotten by Rounds 100-116):**
- Region-safety pre-check (verify ORB_OLBERS_PARAMS + SOURCE57 + SOURCE71 present)
- Regex-pattern verified (Round 42 misfire lesson — use `class X\b[^\n]*?:` + explicit next-class-name lookahead)
- Framework annotations in stub return dict (`backbone / method / shells_used / CPCH / spine / time_frame` — Round 45-79 convention, dropped Rounds 100-116)

---

## PRIORITY 3 — CATCH UP THE REMAINING CP CORPUS

**Scope disclosure:** 1,775 unequivocal stubs across CP1-CP4 remain. Plus 2,447 errored classes needing signature audit. Total unresolved: ~4,222 classes.

At the historical 5-stubs-per-Round pace: **355+ rounds of stub drainage**, or roughly 8-12 months at the historical cadence.

### Suggested ordering

**Round 117 onward — CP1 completion first** (489 stubs remaining):
- Complete CP1 P2 stub drainage before starting CP2/CP3/CP4
- Apply the 5-step Effective Round definition every time
- Ship every 5 rounds (matches v5.47-v5.56 cadence)

**Round ~215 onward — CP2 first pass** (224 stubs):
- CP2 = `CondensedPhysics2.py` — a separate module Rounds 1-116 never touched
- Same 5-step discipline

**Round ~260 onward — CP3 first pass** (255 stubs)

**Round ~310 onward — CP4 first pass** (807 stubs — largest)

**Round ~470 onward — 2,447 errored classes audit**:
- Which are legitimate parameterized calculators (like `AetherMetricQCalcCalculator(M=1e30, r=1e6)`) — they need signature-aware Round work
- Which are truly broken stubs — bucket into drainage queue

### Whitepaper authoring in parallel

Each Round finding a novel closure produces a whitepaper. Historical pace: ~1 whitepaper per 3-5 Rounds. Extrapolating: **~120 more whitepapers** during CP1-CP4 completion, taking corpus from PAPER_1984 to roughly PAPER_2100.

---

## PRIORITY 4 — OPEN ITEMS FROM v5.58-v5.60

### PAPER-level
- **PAPER_1218** — PDF not built. LaTeX error: `Paragraph ended before \text@ was complete` from Unicode superscripts in Higgs branching-ratio table. Fix: wrap formula cells in inline math. → **v5.58.1 patch**
- **PAPER_1801** — PDF not built. LaTeX error: `Missing $ inserted` before `\quad`. Fix: trace unclosed math-mode delimiter earlier in file. → **v5.58.1 patch**
- **PAPER_1087** — Unit erratum marked OPEN_QUESTION in v5.59.0. Resolution requires clarifying κ units in §3 table vs abstract. → Future work; not blocking

### Framework-annotation depth
- 316 Phase C auto-extracted annotations have only `framework_papers` field populated. Full 10-field classification (`backbone / method / shells_used / CPCH / spine / F_U_zero_shell / time_frame / candidate_closures_flagged`) requires per-stub physics review. → Future work

### Calculator-surface expansion
- `calculate_cp_call_UQFF` currently routes to CP1-4 via generic dispatch. Wiring CP1 primitive-locked identities as dedicated `calculate_*` surfaces (instead of only via `PARADOX_TO_CLOSURE`) would provide a second discovery path. → Future work

### Task-list housekeeping
- Task list has 400+ items, most historical. Phase D pruning was performed for actives. Future sessions should prune completed rounds' fine-grained tasks aggressively.

---

## RULES (mirrored from CLAUDE.md)

1. READ CLAUDE.md FIRST EVERY SESSION.
2. DO NOT REVERT canonical primitives.
3. NO NARRATIVE in calculator.
4. NO SM ANYWHERE.
5. PUBLIC SURFACE: `{'value': X}` only.
6. NO datetime/json/file-writes/`__main__`/classes in calculator.
7. NO "0.000% error" without numerical proof.
8. RUN FIDELITY GATE AFTER EVERY EDIT.
9. APPEND to SESSION_LOG, never rewrite.
10. DANIEL PROVIDES INFORMATION. YOU ASSEMBLE IT.
11. DO NOT MODIFY existing Bucket A-K wiring without explicit user request.
12. Daniel has been fighting AI drift for 10+ months.
13. **NEW (v5.61.0): DEFINE EFFECTIVE ROUND per Priority 2 above. 5 steps land in same ship or the round is incomplete.**
14. **NEW (v5.61.0): CP1-CP4 status honestly disclosed. Do not celebrate "50% coverage" against a wrong denominator.**

---

## EDIT / WRITE TOOL WARNING

- Edit tool truncates files > ~1 MB silently. Confirmed multiple times.
- Write tool truncates the same way.
- For `uqff_pure_calculator.py` (3.2 MB, CRLF line endings), `CondensedPhysics.py` (8.2 MB), and similarly large files: **use Python heredoc + `replace()`** instead:

```python
with open('uqff_pure_calculator.py','r',encoding='utf-8',newline='') as f:
    src = f.read()
# ... mutate src via .replace() ...
# ALWAYS assert tail integrity before write:
assert src.rstrip().endswith(']')  # or appropriate sentinel
with open('uqff_pure_calculator.py','w',encoding='utf-8',newline='') as f:
    f.write(src)
```

**Also: preserve CRLF explicitly** — pass `newline=''` to both open() calls or CRLF gets normalized to LF and every line of the file counts as modified.

---

## BACKUP HYGIENE POLICY (v5.61.0 formalization)

Fourteen `.PRE_*` backups accumulated on `uqff_pure_calculator.py` between 2026-06-08 and 2026-06-16. Policy:

**Keep:**
- `PRE_PHASE2_BACKUP`, `POST_BUCKET0_BACKUP`, `POST_BUCKETA_BACKUP` through `POST_BUCKETK_BACKUP` — reference points for the Bucket A-K wiring
- `PRE_PURIFY_BACKUP` — pre-narrative-purge state
- `PRE_RESTORE_BACKUP` — post-purge, pre-canonical-restoration state

**May prune** (after 6-month cold storage):
- Intermediate `PRE_FIX_BACKUP` variants that predate the current TOTAL PURGE state
- Any backup where git already has the same commit-tagged state as a ship version (v5.44.0, v5.47.0, etc.)

**Never prune:**
- Any backup marked "canonical" or "landmark" in its filename
- Any backup Daniel has explicitly named in a session

Prune only after fresh gate run confirms current state is stable AND Daniel explicitly authorizes each pruned file.

---

## SESSION MILESTONE MARKERS

Historically the project has celebrated milestones. Honest markers going forward:

- **50% CP1** — 407/813 = 407 stubs. Currently at 324. **~83 stubs (~17 Rounds) away.**
- **75% CP1** — 610/813. **~286 stubs (~57 Rounds) away.**
- **100% CP1** — 813/813. **~489 stubs (~98 Rounds) away.**
- **All-CP 50%** — 1,050/2,099. **~726 stubs (~145 Rounds) away.**
- **All-CP 100%** — 2,099/2,099. **~1,775 stubs (~355 Rounds) away.**

At the sustainable 6 Rounds/day pace observed in the best week (mid-June 2026), 100% CP corpus is ~60 sessions away.

---

## END OF NEXT_PRIORITIES.md v5.61.0


---

## APPENDED 2026-07-11 — Session Rounds 117-120 candidate log

### Deferred novelty candidates worth future consolidation

**F_TRZ ladder anchors accumulated this session:**

| n | Rung | Current anchors | Candidate additions | Session |
|---|---|---|---|---|
| 6 | F_TRZ^6 = 10^-6 | Pillars ISM (PAPER_1985) | — | R117 discovery |
| 8 | F_TRZ^8 = 10^-8 | Birds (PAPER_1835), solar wind (PAPER_588), Crab outer (PAPER_1986), + 5 more via N-regime | rho_SCm,vac base (PAPER_043/049), QEC phonon floor (PAPER_1056), CMB-S4 mu upper (PAPER_1180), Helix omega_0 (PAPER_070), Uranus/Neptune (PAPER_1078/1079) | R118 + deeper double-check |
| 10 | F_TRZ^10 = 10^-10 | Strong CP (PAPER_1823), MOND a_0 (PAPER_1855) | **Crab quantum uncertainty Delta_x = 1e-10 m (R120 candidate; illustrative default, not measured)** | R120 |

**SO_5-power ladder anchors accumulated this session:**

| Slot | Current anchors | Candidate additions | Session |
|---|---|---|---|
| SO_5^3 = 1000 | v_superwind M82 (PAPER_784) | — | pre-session |
| 2*SO_5^3 = 2000 | v_wind Antennae (PAPER_1972), Westerlund 2 + NGC 3603 (PAPER_1911 seminal) | Rings v_wind (R120 numerical coincidence at lensing scale) | R119 + R120 |
| SO_5^7 = 10^7 | — | **NGC 2525 M_BH = (N_CH/D_phys)*SO_5^7 (PAPER_1985)** | R117 discovery |
| SO_5^10 = 10^10 | Hubble time = 10 Gyr (PAPER_1952, PAPER_1955) | **CompressedMode r = 10^10 m (R120 candidate; length not time)** | R120 |
| SO_5^30 = 10^30 | NS mass ~2.8x10^30 kg (PAPER_148, PAPER_1944) | **CompressedMode M = 10^30 kg (R120 candidate; illustrative default)** | R120 |

### 1.683 prefactor investigation (deferred from R117 deeper double-check)

- PAPER_462 documents rho_SCm/rho_UA = 1.683e-97 (three-leg proofset Leg 2)
- PAPER_463 documents Bohr ground-state E_0 = 1.683 x 10^-37 J
- **1.683 appears in TWO independent contexts** — worth investigating whether 1.683 itself is a structural coefficient of primitives

### Meta-catalog papers spawned this session

- **PAPER_1985** (R117 dual discovery) — Pillars F_TRZ^6 + NGC 2525 mass
- **PAPER_1986 Draft 2** (R118 + deeper double-check) — F_TRZ^8 three-regime -> N-regime (8 anchors)
- **PAPER_1987** (R118 deeper double-check spawn) — 2/3 EXACT supercomposite catalog (8 domains)
- **PAPER_1988** (R119 double-check spawn) — bipartite Sum_Ug closure with D_phys anchor

### Actionable next-session items

1. Consider a **PAPER_1919 revision** to comprehensively catalog n=8 (now 8+ anchors) and add the R120 candidate at n=10 (Crab Delta_x)
2. Consider a **PAPER_1955 revision** to promote SO_5^10 length application (previously time-only)
3. **1.683 prefactor deep-dive** (deferred from R117)
4. **Continue Round 121+** stub drainage at existing rhythm
