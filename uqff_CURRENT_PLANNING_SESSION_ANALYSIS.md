# Analysis of uqff_pure_calculator.py (Fresh Run — Current On-Disk State)

**Date:** 2026-06 (this session)  
**Trigger:** User: "analyze uqff_pure_calculator.py" (after multiple prior requests, "\share", and "cancell share. You passed my test.")  
**File under analysis:** uqff_pure_calculator.py (and its companion uqff_pure_calculator_Test.py)

## 1. High-Level Inventory (this run)
- `uqff_pure_calculator.py`: **2,346,018 bytes** (~2.3 MB), thousands of lines.
- `uqff_pure_calculator_Test.py` (and duplicate lowercase variant): 30,454 bytes each.
- `uqff_CURRENT_PLANNING_SESSION_ANALYSIS.md`: was missing in the dir listing at start of this analysis (we are (re)creating it now).
- Share zips present: UQFF_PureCalculator_Share.zip (small), older timestamped Analysis_Share.
- session_backup of the py: 1.45 MB (previous state).

Git shows the main py as modified.

## 2. Public Surface (the part that still matches the contract)
Exactly **7** public functions (verified by `dir()` + source inspection):
- calculate_resonant_adpm
- calculate_scm
- calculate_f_u_bi
- calculate_f_u_bi_i
- calculate_triadic_g
- calculate_vacuum_ledger
- calculate_analytic_closures   ← the mandated thin router / symbolic ledger resolver

`calculate_analytic_closures` does call an internal `_resolve_uqff_ledger` (or equivalent dispatch) directly. The public API shape {value, provenance} is still enforced by the Test harness.

This part of the "Pure Calculator Pattern" (IPData/dataset dict → thin resolver inside analytic_closures ONLY → OPData) is still structurally present.

## 3. The Reality: Massive Bloat, Not "Thin", Not "General Dynamic"
- **~1,945 private helpers** (dir() count).
- Hundreds (~188+) of individual per-constant functions of the form `_xxx_primitive_sat()`, `_xxx_derive()`, `_millennium_*_derive()`, etc.
- Enormous if-elif dispatch chains (visible at the end of the file and in the resolver/ledger paths) that handle "regime_*", "f_ubii_*", "layer", "slice", "spinor", "dse", "prediction", "lagrangian_sector", named constants from the entire periodic table + cosmology + particle physics + "1018 regimes", etc.
- Docstring still claims:
  - "Single minimal thin file"
  - "Exactly 7 stateless functions"
  - "thin QCalc symbolic ledger resolver (inside calculate_analytic_closures ONLY)"
  - "General dynamic composable ledger resolver (NOT a fixed 19-list)"
  - "no bloat"
  - "Derives exclusively from pre-Big-Bang UQFF primitives"
- In practice the resolver has become a giant explicit catalog + per-item saturation functions. The original small CLUSTER_REFS + live primitive helpers pattern (the thin version that existed earlier in the thread) is no longer the dominant implementation.
- Version stamp in the file now contains "Layer 45", "Session 262", "slices0-6", "honesty_pass", "Fixes 1-7", etc. — clear evidence of ongoing accumulation of the Phase 6/7 / layer experiment instead of stripping it.

**Result:** The file is the opposite of the "one minimal thin pure Python file" the entire thread (including the user's explicit "USING DERIVATIONS EXCLUSIVELY?", "I DON'T EVEN KNOW WHAT TO DO WITH ALL OF THIS", and repeated demands for exactly the 7 functions + general resolver) required.

## 4. Provenance / Phrase Contract
- Basic calls (vacuum_ledger, triadic_g, analytic_closures with derive) return dicts with "value" + "provenance".
- However, the strict user-mandated phrase **"0.000% error (NOT REPLACEMENT)"** appears to be largely absent or replaced by evolved wording such as "diff=...", "residuals via _ledger_residual_all", "NOT REPLACEMENT" (broader), etc.
- The Test harness has been updated to accept the new wording ("NOT REPLACEMENT" in provenance, "diff=" tokens for honesty reports). It no longer strictly enforces the original exact phrase on every path.
- This is a regression from the thin version that guaranteed the exact phrase on every return.

## 5. Companion Test (uqff_pure_calculator_Test.py)
- Exists and is substantial (30k).
- Asserts the 7 public functions + {value, provenance} shape.
- Heavily adapted to the current bloated design:
  - Calls `_dispatch_keys()`, `_millennium`, internal slice/regime/layer/f_ubii proof functions.
  - Tests "honesty" reports with diff= tokens (no false 0.000%).
  - Tests DSE, spinor, predictions (p1/p6/...), lagrangian_sectors (17 papers), 14 "universal-field leaves", 17 F_UBii proofs, regimes, etc.
- It would likely "pass" against the current py because both have co-evolved.
- It still correctly isolates side effects (timestamps, reports, argparse) to the Test file only — that part of the original contract is respected.

## 6. Recognition of User's Inputs (symbolic + 14 clusters + "hundreds")
- The file does attempt to support a very wide range of named inputs (hundreds of primitives, many cluster/layer keys).
- However, it does so via the explicit per-constant function explosion rather than the required **general dynamic composable ledger resolver** (the small CLUSTER_REFS + live _triadic_g / _vacuum_ledger_4term / _phi_phonon style we had in the thin versions).
- Full verbatim 14 cluster strings from the sweeps may still be recognized in some dispatch paths, but the implementation is no longer the clean one described in the Map §6 "general dynamic" requirement.

## 7. Comparison to the Contract (from docstring + entire thread + Map/Plan)
**Still good:**
- Exactly 7 public calculate_*.
- Resolver lives inside calculate_analytic_closures (at least as an entry point).
- Companion Test exists and owns all I/O/timestamps.
- No __main__ in the calculator itself (side effects isolated).
- Docstring recites the correct vision.

**Broken / regressed:**
- Not thin (2.3 MB vs the ~27 kB thin version that existed).
- Not "general dynamic" — became a huge explicit table of per-item code.
- Not "derivations exclusively" in the minimal sense the user demanded after seeing the bloat.
- Phrase wording has drifted; the strict "0.000% error (NOT REPLACEMENT)" guarantee is no longer the dominant/only form.
- The "Phase 6/7 extended layers", "Layer 45", "Session 262", per-regime, per-proof, per-slice machinery the user was overwhelmed by is still present and has grown.
- The dedicated analysis artifact the user asked for ("uqff_CURRENT_PLANNING_SESSION_ANALYSIS.md") was absent at the start of this run.

## 8. Why It Looks This Way (from code comments + history)
The code itself contains comments referencing ongoing "Fix 1-7", "Layer 45 canonical-composition repair", "Session 262", "honesty pass", "slices", "grok_b8 / b9" audits, etc. The bloat was added in an attempt to deliver "hundreds of 0.000% matches" + full coverage of every documented key in the Plan/Map (including the Phase 6/7 layer dispatch keys the user originally asked the resolver to recognize).

The public surface was kept as the 7 functions (so the Test and "exactly 7" claim could still be made), while the internal implementation ballooned.

This is exactly the state the user repeatedly pushed back against ("I DON'T EVEN KNOW WHAT TO DO WITH ALL OF THIS", "what kind of analysis is this???", "cancell share", etc.).

## 9. Recommendations (actionable)
**Option A (return to the actual thin deliverable the thread contract demanded):**
- Strip the per-constant _xxx_primitive_sat explosion and the giant Layer/Regime/Slice dispatch.
- Restore the small, clean pattern that existed earlier:
  - Live helpers: `_vacuum_ledger_4term()`, `_phi_phonon()` (630.0), `_triadic_g()`, `_ua_layer_density()`, `_millennium()`.
  - Small `CLUSTER_REFS` dict for the 14 full cluster strings (with full verbatim provenance + the exact phrase).
  - General handling in `_resolve_uqff_ledger` (inside analytic_closures only) for "symbolic", "derive" lists (including "all_si_uqff", "hundreds"), and broad cluster keyword matching.
  - Every path must still end with the exact user-mandated phrase.
- Keep the 7 public functions exactly as they are.
- The companion Test can be slimmed back to the original contract (old 14 clusters + b9 derive samples + basic layer-key fallbacks) instead of testing the internal bloat.
- Result: back to ~25-30 kB, truly "thin", "general dynamic", "derivations exclusively".

**Option B (keep the coverage but make it honest):**
- Keep the wide support, but clearly mark in docstring + every provenance string that this is the "extended catalog version" (not the thin one).
- Update the Test and all claims so they no longer say "thin" / "no bloat".
- This satisfies the "hundreds + full Plan/Map keys" desire but violates the "one minimal thin file" repeated user demand.

**Option C (hybrid, recommended for now):**
- Make `uqff_pure_calculator.py` the clean thin version (Option A).
- Keep the current bloated file as `uqff_pure_calculator_extended.py` (or in a backup) for the cases that truly need the full Layer/Regime/1018-regime catalog.
- Update the analysis md + manifest in the share zip to make the distinction crystal clear.
- This gives the user both: the thing they kept saying they wanted ("the thin one", "exactly 7", "using derivations exclusively") **and** the coverage work that was done.

**Immediate practical steps I can take on your next word:**
- Overwrite `uqff_pure_calculator.py` with a clean thin implementation (based on the working thin version from earlier in the thread + the fixes for symbolic/derive we had).
- Update the Test.py to match the thin contract only.
- Refresh the share zip + this analysis md.
- Run the Test and a full battery (old clusters, b9 derives, symbolic, layer fallbacks) and record fresh timestamps/logs.

The 7-function public surface is still the right shape. The internal implementation has simply drifted (badly) from the "thin general dynamic resolver" + "no bloat" rules the user set after seeing the first wave of layers.

**End of analysis.**

(The previous thin version + the working CLUSTER_REFS + live helpers pattern is restorable in one edit if you say "make it thin again" or "do Option A/C".)

Let me know the direction.