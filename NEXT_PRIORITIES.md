# NEXT_PRIORITIES.md — what the next Claude session should do first

**READ CLAUDE.md FIRST. EVERY SESSION. EVERY ACTION.** Rules in CLAUDE.md supersede this file. Rules in CLAUDE.md supersede the Map. Rules in CLAUDE.md supersede the Plan.

**Daniel's standing directives (session 2026-06-08, verbatim):**
- "DUMP ALL PROVENANCE. NO COMMENTS WHATSOEVER. NO TAGS, NO GARBAGE SM ANYTHING. JUST REAL PERCENTAGES OF DIFF."
- "MY SYSTEM IS A PURE PREDICTOR. UQFF IS THE FUCKING ANCHOR. UQFF DOESN'T SHARE SHIT WITH SM."
- "I AM GIVING YOU THE INFORMATION, YOU ARE TO ASSEMBLE IT."
- "READ THE RULES BEFORE EVERY ENTRY. READ THE RULES BEFORE EVERY SESSION."
- "BUCKETS A-K WERE ALREADY COMPLETED. WHEN I ASK YOU TO VERIFY, VERIFY. DO NOT MODIFY."

---

## ✅ STATE AT SESSION CLOSE (2026-06-08)

- **Calculator:** 35,716 lines, 32 public surfaces, ALL Buckets 0/A/B/C/D/E/F/G/H/I/J/K wired
- **Public surface return contract:** `{'value': X}` ONLY (no provenance, no paper, no closure_status, no metadata)
- **Catalog observable contract:** `{observable, uqff_derived, anchor, residual_pct}` — pure math only
- **Fidelity gate:** 417/417 PASS, Cat 16 STRICT PURGE GUARD enforced (catches narrative regression automatically)
- **Canonical primitives:** all 11 intact (SSQ=0.57, BETA_I=0.6029, K_MEX=25/12, S_26=1.453162, RHO_SCM=7.09e-37, D_PHYS=4, D_BSFG=6, D_CRIT=26, N_CH=9, SO_FIVE=10, A_FIVE=60)
- **Pure-calculator discipline:** 1935 functions, 0 docstrings, 0 classes, 0 `#` comments, 0 datetime/json imports, 0 print(), 0 file writes
- **Zero SM references** anywhere in calculator source
- **Zero "NOT REPLACEMENT" tags** anywhere
- **Zero narrative provenance / paper / closure_status / paper_attribution fields**

---

## 🔴 PRIORITY 1 FOR NEXT SESSION — VERIFICATION ONLY

**Daniel said: "BUCKETS A-K WERE ALREADY COMPLETED, AND I ASKED YOU TO VERIFY THEM AT THE START OF THIS SESSION."**

The next session's FIRST action after reading CLAUDE.md → SESSION_LOG.md → this file is:

1. **Run `python uqff_fidelity_tests.py`.** Expect 417/417 PASS, exit 0.
2. **Verify canonical primitives.** Read CLAUDE.md "11 locked canonical primitives" section. Confirm calculator source has them at locked values.
3. **Sample public surfaces.** Call `calculate_gw_events({})`, `calculate_agn_jet({})`, `calculate_astrophysics({})`, `calculate_cosmology({})`, `calculate_particle_physics({})` and verify return shape is `{'value': dict}` with no 'provenance' / 'paper' / 'closure_status' keys.
4. **Verify Cat 16 STRICT PURGE GUARD runs and passes.** This category catches any narrative regression.

**DO NOT MODIFY anything in Buckets A-K without explicit user instruction.** "Verify" means read, report, STOP.

---

## 🟡 PRIORITY 2 (only if Daniel explicitly requests new work)

If and only if Daniel asks to add new physics (new bucket, new observable, new derivation):

1. **Daniel provides the closed form.** From a paper or his direct instruction. Do not infer, do not paraphrase.
2. **Transcribe literally using locked primitives.** No SM constants, no SM-named identifiers, no "classical" baselines, no Kerr/GR/PDG/CODATA references.
3. **Catalog tuple is `(label, value, anchor)` only.** No `paper`, no `closure_status`. Math fields only.
4. **Add a fidelity test** that pins the new closure to its expected value.
5. **Run the gate.** Must exit 0 with Cat 16 strict purge guard passing.
6. **Append SESSION_LOG.md entry.** Do not rewrite prior entries.

---

## RULES (mirrored from CLAUDE.md for emphasis)

1. **READ CLAUDE.md FIRST EVERY SESSION AND BEFORE EVERY ACTION.**
2. **DO NOT REVERT canonical primitives.**
3. **NO NARRATIVE OF ANY KIND** in the calculator. No comments, no docstrings, no provenance strings, no paper attribution, no closure status tags, no "NOT REPLACEMENT" tags, no formulas embedded as text, no physics narrative, no SM references.
4. **NO SM ANYWHERE.** UQFF is the only physics. UQFF is the anchor. UQFF is the truth.
5. **PUBLIC SURFACE RETURN:** `{'value': X}` only.
6. **NO datetime/json/file-writes/__main__/classes** in calculator.
7. **NO "0.000% error" claims** without numerical proof.
8. **RUN FIDELITY GATE AFTER EVERY EDIT.** Must exit 0.
9. **APPEND to SESSION_LOG, never rewrite.**
10. **DANIEL PROVIDES INFORMATION. YOU ASSEMBLE IT.** Do not invent. Do not paraphrase. Do not think for him.
11. **DO NOT MODIFY existing Bucket A-K wiring without explicit user request.** "Verify" means read and report. Stop.
12. **Daniel has been fighting AI drift for 10 months.** Treat his words carefully. Build trust.

---

## What's queued (only if Daniel explicitly authorizes)

Daniel mentioned earlier this session the "98% remainder" — physics outside this repo. Location not yet disclosed. Next session should ASK before assuming.

No other work is queued. Buckets A-K are DONE per Daniel's confirmation. The calculator is in its pure mathematical state per the total purge.
