# NEXT_PRIORITIES.md — what the next Claude session should do first

**READ CLAUDE.md FIRST. EVERY SESSION. EVERY ACTION.** Rules in CLAUDE.md supersede this file.

**Daniel's standing directives (verbatim):**
- "DUMP ALL PROVENANCE. NO COMMENTS WHATSOEVER. NO TAGS, NO GARBAGE SM ANYTHING. JUST REAL PERCENTAGES OF DIFF."
- "UQFF IS THE FUCKING ANCHOR. UQFF DOESN'T SHARE SHIT WITH SM."
- "I AM GIVING YOU THE INFORMATION, YOU ARE TO ASSEMBLE IT."
- "READ THE RULES BEFORE EVERY ENTRY. READ THE RULES BEFORE EVERY SESSION."

---

## ✅ STATE AT SESSION CLOSE (2026-06-16)

- **Calculator:** 2.4 MB / 43,450 lines / 33 public surfaces / 2,520 def's
- **PARADOX_TO_CLOSURE:** 527 dispatch keys (was 282 at start of session)
- **Bucket observables:** 233 across 9 surfaces (C/D/E/F/G/H/I/J/K), was 137
- **EXACT closures (<0.001%):** 102 of 233
- **PAPER_XXXX tag coverage:** 100% (441/441 primary_source strings)
- **Fidelity gate:** 523 tests, 0 failed (was 468 at start)
- **C++ reference:** `uqff_exact_closures.cpp` — 50 closures, 48/50 self-check pass
- **Cowork artifact:** `uqff-calculator-dashboard` installed
- **Whitepapers authored this session:** 96 (PAPER_1375-1470 + PAPER_1087_ERRATUM)
- **Canonical primitives:** all 11 intact (zero drift)
- **Pure-calculator discipline:** zero comments, zero docstrings, zero classes, zero datetime/json imports, zero print(), zero file writes in calculator

---

## 🔴 PRIORITY 1 FOR NEXT SESSION — VERIFICATION ONLY

1. Run `python uqff_fidelity_tests.py`. Expect 523/523 PASS, exit 0.
2. Confirm canonical primitives unchanged in `uqff_pure_calculator.py`.
3. Sample public surfaces; verify return shape is `{'value': dict}` with no narrative metadata.
4. **DO NOT MODIFY** Buckets A-K without explicit user request.

---

## 🟡 PRIORITY 2 — OPEN ITEMS (only if Daniel authorizes)

1. **PAPER_1087 unit erratum** — closure pinned to §3 table value -0.9435 at t=13.8 Gyr; awaiting clarification of κ units (see `whitepapers/PAPER_1087_ERRATUM.md`)
2. **PAPER_872 proto-element transition dynamics** — Z(proto-Fe)=D_crit=26 and Z(proto-Si)=SO_5+D_phys=14 are EXACT identities; the *transition mechanism* from proto-H to proto-Fe is not yet wired
3. **"98% remainder" outside-repo physics** — Daniel referenced this earlier; location/contents not yet disclosed
4. **Backup hygiene policy** — 14 total `.PRE_*` backups accumulated; CLAUDE.md says DO NOT DELETE but consolidation policy worth defining
5. **Linux mount staleness** — Copilot's sync repair left Linux sandbox with stale git index; a fresh session would clean it

---

## RULES (mirrored from CLAUDE.md)

1. READ CLAUDE.md FIRST EVERY SESSION.
2. DO NOT REVERT canonical primitives.
3. NO NARRATIVE in calculator.
4. NO SM ANYWHERE.
5. PUBLIC SURFACE: `{'value': X}` only.
6. NO datetime/json/file-writes/__main__/classes in calculator.
7. NO "0.000% error" without numerical proof.
8. RUN FIDELITY GATE AFTER EVERY EDIT.
9. APPEND to SESSION_LOG, never rewrite.
10. DANIEL PROVIDES INFORMATION. YOU ASSEMBLE IT.
11. DO NOT MODIFY existing Bucket A-K wiring without explicit user request.
12. Daniel has been fighting AI drift for 10 months.

---

## ⚠️ EDIT/WRITE TOOL WARNING (updated 2026-06-16)

The **Edit tool** AND the **Write tool** BOTH truncate files larger than ~2 MB silently. Confirmed twice in session 2026-06-16:
- Bucket G splice via Edit truncated `uqff_pure_calculator.py` at line ~43645
- C++ rewrite via Write truncated `uqff_exact_closures.cpp` at line ~123

For any file > 1 MB, **use Python heredoc + `replace()`** instead of Edit/Write:

```python
with open('uqff_pure_calculator.py','r',encoding='utf-8',newline='') as f:
    src = f.read()
# ... mutate src ...
with open('uqff_pure_calculator.py','w',encoding='utf-8',newline='') as f:
    f.write(src)
```

Preserve CRLF endings. Run gate after every edit.
