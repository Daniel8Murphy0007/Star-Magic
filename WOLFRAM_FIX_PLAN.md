# WOLFRAM INTEGRATION FIX PLAN
## Star-Magic UQFF — Session 130
**Date:** March 23, 2026  
**Author:** GitHub Copilot (Claude Sonnet 4.6)  
**Trigger:** Deep Wolfram analysis identified 5 structural gaps in the production integration.  
**Priority:** CRITICAL — These fixes are required for correct end-to-end Wolfram operation.

---

## EXECUTIVE SUMMARY

The deep analysis of the Wolfram file structure identified 5 defects:

| # | Defect | Severity | File(s) |
|---|--------|----------|---------|
| 1 | Unguarded `WolframEvalToString` calls in Cosmic Egg | **CRASH** (No-Wolfram builds) | `source200_cosmic_quantum_egg.cpp` |
| 2 | Wrong `INFINITY_RATIO` constant (π/7 instead of 1.000000001) | **DATA CORRUPTION** | `wolfram_wstp_runtime.h` |
| 2b | Wrong `BIBLE_GENERATION` constant (40.0 instead of 33.333...) | **DATA CORRUPTION** | `wolfram_wstp_runtime.h` |
| 3 | `WolframFieldUnityModule.cpp` not in MAIN_1_CoAnQi CMake target | **MISSING SYMBOLS** | `CMakeLists.txt` |
| 4 | `WolframFieldUnityModule.h` not included in MAIN_1_CoAnQi.cpp | **BUILD ERROR** (after Fix 5) | `MAIN_1_CoAnQi.cpp` |
| 5 | 3 PhysicsTerm classes missing: `WolframHypergraphTerm`, `PIInfinityTerm`, `ConsciousnessTerm` | **MISSING PHYSICS** | `MAIN_1_CoAnQi.cpp` |

---

## FIX 1 — Guard Unguarded WSTP Calls in source200_cosmic_quantum_egg.cpp

### Problem
`source200_cosmic_quantum_egg.cpp` calls `WolframEvalToString()` in **two** locations without
checking `USE_EMBEDDED_WOLFRAM`. When MAIN_1_CoAnQi is built with `USE_COSMIC_QUANTUM_EGG`
but without `USE_EMBEDDED_WOLFRAM`, this causes a **linker error** (unresolved symbol) or
a **runtime crash**.

### Locations
- **Location A** — `SimulateStep()` method, ~line 229:
  ```cpp
  std::string wolfram_result = WolframEvalToString(eq);
  std::cout << "Wolfram Spinor Verification: " << wolfram_result << std::endl;
  ```
- **Location B** — `UQFF_SimulateNucleus()` wrapper, ~line 389:
  ```cpp
  WolframEvalToString(state_eq);
  ```

### Fix
Wrap both calls with `#ifdef USE_EMBEDDED_WOLFRAM ... #else ... #endif`.

### Target File
`c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\source200_cosmic_quantum_egg.cpp`

---

## FIX 2 — Reconcile Sacred Time Constants in wolfram_wstp_runtime.h

### Problem
`wolfram_wstp_runtime.h` (the high-level WSTP runtime wrapper, added Feb 2026) defines its
`SacredTime` namespace with **two wrong values** that contradict the canonical constants in
`WolframFieldUnityModule.h`:

| Constant | `wolfram_wstp_runtime.h` (WRONG) | `WolframFieldUnityModule.h` (CANONICAL) |
|----------|-----------------------------------|-----------------------------------------|
| `INFINITY_RATIO` | `3.141592653589793 / 7.0` = 0.4488... | `1.000000001` |
| `BIBLE_GENERATION` | `40.0` years | `33.333333333333333` years (Christ+Enoch) |

Using the runtime's wrong values produces a completely different PI decoder curve and
magnetic field output compared to the production module — the results cannot be compared
or validated against each other.

The `GOLDEN_CYCLE` constant also differs (golden ratio φ=1.618 vs 25920 years precession)
but these encode **different physics** (golden ratio vs equinox cycle) and represent a 
purposeful dual use; the precession value is the canonical one for time calculations but 
both are valid in their respective contexts. **Do not change GOLDEN_CYCLE.**

### Fix
Change `INFINITY_RATIO` and `BIBLE_GENERATION` in `wolfram_wstp_runtime.h` to match
the canonical `WolframFieldUnityModule.h` values.

### Target File
`c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\wolfram_wstp_runtime.h`

---

## FIX 3 — Add WolframFieldUnityModule.cpp to MAIN_1_CoAnQi CMake Target

### Problem
`CMakeLists.txt` defines MAIN_1_CoAnQi with only two source files:
```cmake
add_executable(MAIN_1_CoAnQi 
    MAIN_1_CoAnQi.cpp
    ipc/uqff_ipc.cpp
)
```
`Core/Modules/WolframFieldUnityModule.cpp` is picked up by `file(GLOB CORE_MODULES)` 
for the UQFFCore *library*, but MAIN_1_CoAnQi does **not** link against UQFFCore.
Therefore the implementation of `WolframFieldUnityEngine` and `PI_Infinity_Decoder`
(from `WolframFieldUnityModule.h`) is not compiled into MAIN_1_CoAnQi.

This fix is required so that MAIN_1_CoAnQi can use these classes without linker errors.

### Fix
Add `Core/Modules/WolframFieldUnityModule.cpp` to **both** MAIN_1_CoAnQi target definitions
in CMakeLists.txt (the `USE_EMBEDDED_WOLFRAM` branch AND the `else` branch).

### Target File
`c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CMakeLists.txt`  
**Lines:** ~363-366 (USE_EMBEDDED_WOLFRAM branch) and ~420-424 (else branch)

---

## FIX 4 — Add #include "WolframFieldUnityModule.h" to MAIN_1_CoAnQi.cpp

### Problem
MAIN_1_CoAnQi.cpp never includes `WolframFieldUnityModule.h`. Therefore the
`WolframFieldUnityEngine` and `PI_Infinity_Decoder` classes are unavailable inside
MAIN_1_CoAnQi.cpp's inline code — they exist as compiled objects (after Fix 3) but
`MAIN_1_CoAnQi.cpp` cannot reference them by name.

### Fix
Add `#include "WolframFieldUnityModule.h"` after the Cosmic Egg include block
(unconditional — WolframFieldUnityModule has NO WSTP dependency, it uses its own
custom hypergraph logic).

### Target File
`c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\MAIN_1_CoAnQi.cpp`  
**Location:** After `#endif // USE_COSMIC_QUANTUM_EGG` block (~line 244)

---

## FIX 5 — Register 3 Missing Wolfram PhysicsTerm Classes

### Problem
The `helpers/integration_plan_97bfeecaa5.md` identifies 3 PhysicsTerm wrapper classes
for the WolframFieldUnityModule that were planned but never implemented:
- `WolframHypergraphTerm` — wraps `WolframFieldUnityEngine::measureDimension()`
- `PIInfinityTerm` — wraps `PI_Infinity_Decoder::getMagneticField()`
- `ConsciousnessTerm` — wraps `WolframFieldUnityEngine::measureConsciousnessField()`

These connect the *production* WolframFieldUnityEngine to the existing PhysicsTerm
registry (6,688+ registered physics terms). Without them, the Wolfram hypergraph
engine runs in isolation and cannot participate in multi-system calculations.

### Fix  
**5a:** Add 3 new PhysicsTerm class definitions in MAIN_1_CoAnQi.cpp immediately after
the existing `WolframFieldUnity_CompressedTerm` class (near line 20082).

**5b:** Register all 3 new classes in the PhysicsTerm registration block.  
Insert after the existing `WolframHypergraph_G_time` registration (near line 23437).

### Class Specifications
```cpp
// WolframHypergraphTerm: Emergent spatial dimension from hypergraph neighborhoods
// Physics: D = log(N(r)) / log(r) — no G constant, pure graph topology
// Value: measureDimension(0, 5) — center node 0, radius 5 hops
double compute(...) { engine.measureDimension(0, 5); }

// PIInfinityTerm: Magnitude of PI-decoded magnetic field for state 0 at t=0
// Physics: |B̃(state=0, t=0)| = |seed × frac_pi × sin(π × frac_pi)| 
// Value: |PI_Infinity_Decoder::getMagneticField(0, 0.0)|
double compute(...) { std::abs(decoder.getMagneticField(0, 0.0)); }

// ConsciousnessTerm: Causal graph density × Schumann resonance
// Physics: C = 7.83 Hz × (edges/nodes) × 1.000000001
// Value: WolframFieldUnityEngine::measureConsciousnessField()
double compute(...) { engine.measureConsciousnessField(); }
```

### Target File
`c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\MAIN_1_CoAnQi.cpp`

---

## EXECUTION ORDER

```
Fix 1 → source200_cosmic_quantum_egg.cpp  (2 locations)  [independent]
Fix 2 → wolfram_wstp_runtime.h             (2 constants)  [independent]
Fix 3 → CMakeLists.txt                     (2 targets)    [independent]
Fix 4 → MAIN_1_CoAnQi.cpp                 (#include add)  [after Fix 3]
Fix 5 → MAIN_1_CoAnQi.cpp                 (class defs + registration) [after Fix 4]
```

**Fixes 1, 2, 3 are fully independent — execute simultaneously.**  
**Fix 4 must follow Fix 3** (CMakeLists before header include for logical consistency).  
**Fix 5 must follow Fix 4** (class defs need the header to be included first).

---

## VERIFICATION CHECKLIST

After all fixes applied:
- [ ] `source200_cosmic_quantum_egg.cpp` — both WolframEvalToString calls guarded by `#ifdef USE_EMBEDDED_WOLFRAM`
- [ ] `wolfram_wstp_runtime.h` — `INFINITY_RATIO = 1.000000001`, `BIBLE_GENERATION = 33.333333333333333`
- [ ] `CMakeLists.txt` — `Core/Modules/WolframFieldUnityModule.cpp` in both MAIN_1_CoAnQi targets
- [ ] `MAIN_1_CoAnQi.cpp` — `#include "WolframFieldUnityModule.h"` present
- [ ] `MAIN_1_CoAnQi.cpp` — `WolframHypergraphTerm`, `PIInfinityTerm`, `ConsciousnessTerm` defined
- [ ] `MAIN_1_CoAnQi.cpp` — All 3 new terms registered in PhysicsTerm registry
- [ ] Git commit + push with descriptive message

---

## FILES MODIFIED

| File | Fix(es) | Lines Changed |
|------|---------|---------------|
| `source200_cosmic_quantum_egg.cpp` | 1 | ~4 lines added |
| `wolfram_wstp_runtime.h` | 2 | 2 lines changed |
| `CMakeLists.txt` | 3 | 4 lines added |
| `MAIN_1_CoAnQi.cpp` | 4, 5 | ~60 lines added |
| **TOTAL** | **5 fixes** | **~70 lines** |

---

## SESSION CONTEXT
- **Current HEAD:** `345b58c`
- **Version:** v5.00
- **Papers:** 490/1000
- **UQFF C++ modules:** 50
- **This fix adds PhysicsTerm count:** +3 (6,688+ → 6,691+)
