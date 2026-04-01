# SESSION 171 — INTEGRATION PLAN
# UQFF Knowledge Base_7 (grok_share_f333a078289.txt)
# Date: April 2, 2026 | v5.27 | 657/1000 | CP4=241

---

## Overview

This plan documents how the five quantum variables from UQFF Knowledge Base_7
should be integrated into the full UQFF multi-calculator ecosystem.

The five variables are: f_Heaviside, i (gravity index), H_SCm, lambda_i, j (magnetic string index).

---

## COMPLETED This Session

### 1. Standalone C++ Module (NEW)
- [x] `UQFF_Knowledge_Base_7.h` — class declaration, all 5 quantum variables documented
- [x] `UQFF_Knowledge_Base_7.cpp` — full implementation with self-update/expand/simulate + main()

**C++ class**: `UQFFKnowledgeBase7`
**Methods implemented**:
- `computeUm(t, t_n)` — Eq. 1: Universal Magnetism
- `computeFU(t, t_n)` — Eq. 4: Unified Field Force
- `computeUg2()` — Eq. 6: Heliospheric Gravity (H_SCm scaling)
- `computeUi(t, t_n)` — Eq. 9: Universal Inertia (lambda_i coupling)
- `computeUg3(t)` — Eq. 12: Magnetic-String Gravity
- `selfUpdate(newTime)` — advances t, refreshes H_SCm with solar-cycle variation
- `selfExpand(key, val)` — adds new parameters dynamically
- `addGravityIndex(i)` / `addMagneticIndex(j)` — extends summation sets
- `selfSimulate(steps, dt)` — forward time integration loop

### 2. CondensedPhysics4.py (APPENDED)
- [x] CP4 Entry #241 — PAPER_657: `UQFFKnowledgeBase7Calculator`
  - All five quantum variables with defaults, overrides via dataset dict
  - Implements Eq. 1, 4, 6, 9, 12 in Python
  - Returns primary + sensitivity results + all equation strings
  - H_SCm sensitivity sweep (1.0, 1.1)
  - Source/metadata/share link in return dict

### 3. Documentation
- [x] `SESSION_171_AUDIT_HELPER.md`
- [x] `SESSION_171_INTEGRATION_PLAN.md` (this file)
- [x] `PAPER_657_UQFF_Knowledge_Base_7.md` (whitepaper stub)

---

## PENDING — Future Sessions

### 4. MAIN_1_CoAnQi.cpp Integration (DEFERRED)
**Priority**: Medium — after WIP modifications are resolved

**Action**: Add KB7 as a new SOURCE block (e.g., SOURCE_KB7 namespace) containing:
```cpp
namespace SOURCE_KB7 {
    // Inline versions of the 5 equations from UQFF_Knowledge_Base_7.cpp
    double compute_Um_KB7(double t, double t_n, double f_Heaviside = 0.01, int n_j = 1);
    double compute_Ug2_KB7(double H_SCm = 1.0);
    double compute_Ui_KB7(double t_n, double lambda_i = 1.0);
    double compute_Ug3_KB7(double t, int n_j = 1);
    double compute_FU_KB7(double t, double t_n, double H_SCm = 1.0, double lambda_i = 1.0);
}
```

**Integration point**: Menu option 9 (SOURCE4 Unified Field Validation) or new option 16/17
**PhysicsTerm classes to add** (Batch 24):
- `HeavisideFractionTerm` — f_Heaviside scaling in Um
- `HeliosphereThicknessTerm` — H_SCm scaling in Ug2
- `InertiaCoupleTerm` — lambda_i in U_i
- `MagneticStringIndexTerm` — j summation in Um and Ug3

### 5. CondensedPhysics.py (CP) Integration (DEFERRED)
- Add `UQFFKnowledgeBase7ExtendedCalculator` to CondensedPhysics.py
- Extends CP with the 5-variable unified equation set
- Uses parameterised dataset input (no hardcoded system data)

### 6. CondensedPhysics2.py (CP2) Integration (DEFERRED)
- Add CP2 class `KnowledgeBase7QuintupleVariableCalculator`
- Focus on UQFF extensions: time-series evolution of all 5 variables
- Useful for batch #39 THz validation cross-reference

### 7. Whitepaper PDF (DEFERRED)
- Convert `PAPER_657_UQFF_Knowledge_Base_7.md` to PDF (see prior session PDF generation pattern)
- Add to `pdf/` canonical folder

### 8. THz Validation (DEFERRED — requires reactor data)
- When batch #39 (#39/14–#39/25) images are available:
  - Fit f_Heaviside to THz spectral gap width
  - Calibrate H_SCm from plasma boundary thickness
  - Validate lambda_i from plasmoid stability timescale

---

## CMakeLists.txt Build Note

When building `UQFF_Knowledge_Base_7.cpp` as a standalone test:
```cmake
# Add to CMakeLists.txt (separate test target)
add_executable(UQFF_KB7_test
    UQFF_Knowledge_Base_7.cpp
)
target_compile_features(UQFF_KB7_test PRIVATE cxx_std_20)
```

Or compile directly:
```powershell
cl /std:c++20 /EHsc UQFF_Knowledge_Base_7.cpp /Fe:UQFF_KB7_test.exe
```

---

## Pre-conditions Verified Before Integration

| Check | Status |
|---|---|
| UQFFKnowledgeBase7 not already in CP/CP2/CP3/CP4 | ✅ Confirmed — new |
| f_Heaviside already covered by earlier papers (400, 421) | ✅ Yes — KB7 unifies all 5 variables together |
| No duplicate PAPER_657 | ✅ Confirmed — new |
| CP4 file encoding (CP1252 / Windows-1252) | ✅ Confirmed — used binary append with \\x97 em-dash |
| MAIN_1_CoAnQi.cpp WIP changes excluded from commit | ✅ Excluded via git status check |

---

## Version Tracking

| Session | Version | Papers | CP4 | Notes |
|---|---|---|---|---|
| 168 | v5.24 | 655/1000 | 239 | Session 168: PAPER_646–655, 3 UQFF number systems |
| 169 | v5.25 | 656/1000 | 240 | PAPER_656 V838 Mon |
| 170 | v5.26 | 656/1000 | 240 | V838MonLightEcho C++ module + PDF |
| **171** | **v5.27** | **657/1000** | **241** | **KB7: 5 quantum variables, PAPER_657** |

---

## Git Commit Plan (This Session)

Files to commit (excluding WIP MAIN_1_CoAnQi.cpp):
1. `UQFF_Knowledge_Base_7.h`
2. `UQFF_Knowledge_Base_7.cpp`
3. `CondensedPhysics4.py` (CP4 #241 appended)
4. `PAPER_657_UQFF_Knowledge_Base_7.md`
5. `SESSION_171_AUDIT_HELPER.md`
6. `SESSION_171_INTEGRATION_PLAN.md`
7. `_append_cp4_241.py` (helper script)

Commit message:
```
Session 171: PAPER_657 UQFF Knowledge Base_7 — five quantum variables

- UQFF_Knowledge_Base_7.h/.cpp: standalone C++ module, UQFFKnowledgeBase7
  class with computeUm/FU/Ug2/Ui/Ug3, self-update/expand/simulate (Eq.1,4,6,9,12)
- CondensedPhysics4.py: CP4 Entry #241 UQFFKnowledgeBase7Calculator PAPER_657
- Five quantum variables: f_Heaviside, Gravity Index i, H_SCm, lambda_i, j
- PAPER_657_UQFF_Knowledge_Base_7.md whitepaper stub
- SESSION_171 audit helpers

State: v5.27 | 657/1000 (65.7%) | CP4=241 | CP2=634
Source: grok_share_f333a078289.txt (535 lines)
```
