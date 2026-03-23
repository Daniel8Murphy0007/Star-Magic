# Integration Plan — grok_share_97bfeecaa5.txt Modules
# Session 128 | Star-Magic v5.00
# 7 new UQFF modules → Codebase integration roadmap
# Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

## Summary
This plan describes how the 7 new UQFF modules (from grok_share_97bfeecaa5.txt) integrate
into the Star-Magic codebase across all 6 tiers:
- Tier 2 Calculators (MAIN_1_CoAnQi.cpp + CondensedPhysics2.py)
- Tier 1 GUI (source2.cpp query field → results display)
- Tier 5 Storage (CondensedPhysics_OutputData.py recall)
- Tier 6 Docs (PAPER_484–490 whitepapers)

---

## Phase 1: C++ Integration into MAIN_1_CoAnQi.cpp

### 1.1 Header Inclusion Block
Add to MAIN_1_CoAnQi.cpp near line 100 (alongside other module includes):

```cpp
// v5.00 — Session 128 UQFF Modules (97bfeecaa5)
#include "UQFFCalculationsModule.h"
#include "UQFFBuoyancySNRModule.h"
#include "UQFFCassiniBuoyancyModule.h"
#include "UQFFMultiAstroSystemsModule.h"
#include "UQFFEightAstroSystemsModule.h"
#include "UQFFNineteenAstroSystemsModule.h"
#include "WolframFieldUnityModule.h"
```

### 1.2 New Menu Option (Option 16 / after SOURCE4)
Add new menu entry after the existing Option 15 (SOURCE4 Unified Field Validation):

```
16. Session128 Multi-System UQFF Validation [NEW]
    - Runs all 7 new modules: Calculations(5), BuoyancySNR(5), Cassini(3),
      MultiAstro(11), EightAstro(8), NineteenAstro(19), WolframUnity(hypergraph)
    - Outputs compact table: system name + F_compressed + g_poly magnitude
17. Exit (was 16)
```

### 1.3 Wolfram Integration
WolframFieldUnityModule is independent of Wolfram WSTP (it uses custom hypergraph logic).
Register it as an alternative "Wolfram-style" computation that runs without WSTP:
- Add to the "No Wolfram Build" conditional block
- WolframFieldUnityEngine.evolveMultiway() uses OpenMP → add to CMakeLists.txt OpenMP section

### 1.4 PhysicsTerm Registration
Register new physics terms from each module into PhysicsTerm registry
(MAIN_1_CoAnQi.cpp line 411+, existing PhysicsTermRegistry):

| Module | New Terms |
|--------|-----------|
| UQFFCalculations | UQFFDPMTerm, UQFFMagnetismTerm, UQFFEtaTerm |
| UQFFBuoyancySNR | UQFFBuoyancyMasterTerm, UQFFLENRTerm, UQFFFRelTerm |
| UQFFCassini | UQFFCassiniBuoyTerm, UQFFCassiniTHzTerm |
| UQFFMultiAstro | UQFFMultiCompressedTerm, UQFFMultiResonanceTerm |
| UQFFEightAstro | UQFFEightCompressedTerm (with proof comments) |
| UQFFNineteen | UQFF26DPolyTerm, UQFF26DGravTerm |
| WolframUnity | WolframHypergraphTerm, PIInfinityTerm, ConsciousnessTerm |

### 1.5 CMakeLists.txt Updates
```cmake
# Add new source files (already in Core/Modules/)
target_sources(MAIN_1_CoAnQi PRIVATE
    Core/Modules/UQFFCalculationsModule.cpp
    Core/Modules/UQFFBuoyancySNRModule.cpp
    Core/Modules/UQFFCassiniBuoyancyModule.cpp
    Core/Modules/UQFFMultiAstroSystemsModule.cpp
    Core/Modules/UQFFEightAstroSystemsModule.cpp
    Core/Modules/UQFFNineteenAstroSystemsModule.cpp
    Core/Modules/WolframFieldUnityModule.cpp
)
# OpenMP for WolframFieldUnityModule parallel evolution
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MAIN_1_CoAnQi OpenMP::OpenMP_CXX)
endif()
```

---

## Phase 2: Python Calculator Integration (CondensedPhysics2.py)

### 2.1 New Calculator Classes (add to CP2, v4.3.9+)

Add the following 7 calculator classes. Each follows strict CondensedPhysics.py rules:
- No hardcoded data inside the calculator
- Receives datasets from source2.cpp via APIFetch.py
- Outputs equation sets to CondensedPhysics_OutputData.py

#### Class 1: UQFFCalculationsCalculator
```python
class UQFFCalculationsCalculator:
    """UQFF 5-system calculations: U_g1, U_g3, U_m, E, eta for astrophysical objects"""
    def compute(self, dataset: dict) -> dict:
        # dataset keys: r, B, z, sfr, t_age, Z (atomic number)
        # Returns: F_U_g1, F_U_g3, U_m, E_field, eta (neutron rate)
```

#### Class 2: UQFFBuoyancySNRCalculator
```python
class UQFFBuoyancySNRCalculator:
    """F_U_Bi_i master buoyancy equation: F_LENR + F_act + F_DE + F_res + F_neutron + F_rel"""
    def compute(self, dataset: dict) -> dict:
        # dataset keys: r, B, z, v_exp, eta_base, Q_wave, t_age
```

#### Class 3: UQFFCassiniBuoyancyCalculator
```python
class UQFFCassiniBuoyancyCalculator:
    """Complex-valued Cassini ring buoyancy: U_Bi, U_Ii, U_Mi, THz_hole, q_scope"""
    def compute(self, dataset: dict) -> dict:
        # dataset keys: orbital_r, ring_r, B, rotation_period, wind_vel
```

#### Class 4: UQFFMultiAstroCalculator
```python
class UQFFMultiAstroCalculator:
    """Simultaneous Compressed/Resonance/Buoyancy solutions + DPM creation rate"""
    def compute(self, dataset: dict) -> dict:
        # dataset keys: r, sfr, B, z, t_age
```

#### Class 5: UQFFEightAstroProofCalculator
```python
class UQFFEightAstroProofCalculator:
    """UQFF 7-step equation proofs for star-forming regions: F_compressed + F_resonance"""
    def compute(self, dataset: dict) -> dict:
        # Returns: F_compressed, F_resonance, proof_compressed, proof_resonance
```

#### Class 6: UQFFNineteen26DCalculator
```python
class UQFFNineteen26DCalculator:
    """26D polynomial framework: g(r) = Sum(i=1..26) E_DPM,i/r_i^2 x f_TRZ_i x f_Um_i"""
    def compute(self, dataset: dict) -> dict:
        # dataset keys: M_0, r, sfr, B, z, t_age
        # Returns: g_poly, F_compressed, F_resonance, taylor_coeffs[26]
```

#### Class 7: WolframFieldUnityCalculator
```python
class WolframFieldUnityCalculator:
    """Wolfram hypergraph + PI infinity decoder + sacred time constants"""
    def compute(self, dataset: dict) -> dict:
        # dataset keys: r, sfr, evolution_depth
        # Returns: emergent_dimension, buoyant_gravity, consciousness_field
```

### 2.2 CP2 Version Bump
- Current: v4.3.8 (548 classes)
- Target: v4.3.9 (555 classes after adding 7 new calculators)

---

## Phase 3: source2.cpp GUI Integration

### 3.1 New Query Support
Tab 1 (UQFF Calculator) — add new calculation type picker entries:
- "UQFF 26D Polynomial" → routes to UQFFNineteen26DCalculator
- "Wolfram Field Unity" → routes to WolframFieldUnityCalculator
- "SNR Buoyancy (UQFF)" → routes to UQFFBuoyancySNRCalculator

### 3.2 Results Tab
Tab 7 (Compute Results) — new display rows for:
- g_poly (complex, 26D sum)
- F_U_Bi_i (master buoyancy)
- Emergent Dimension (Wolfram)
- Consciousness Field
- PI Infinity Decoder output (26 states)

---

## Phase 4: CondensedPhysics_OutputData.py Recall Storage

Add recall entries for Session 128 computations:
```python
SESSION_128_MODULES = [
    "UQFFCalculations_5systems",
    "UQFFBuoyancySNR_5systems",
    "UQFFCassiniBuoyancy",
    "UQFFMultiAstro_11systems",
    "UQFFEightAstro_8systems",
    "UQFFNineteen26D_19systems",
    "WolframFieldUnity"
]
```

---

## Phase 5: Whitepaper Documentation (PAPER_484–490)

| Paper | Title | Key Content |
|-------|-------|-------------|
| PAPER_484 | UQFFCalculations 5-System Module | U_g1, U_g3, U_m, E, η for M82/IC418/CanisMajor/NGC6302/NGC7027 |
| PAPER_485 | UQFF SNR Buoyancy Master Equation | F_LENR+F_act+F_DE+F_res+F_neutron+F_rel for SN1006/EtaCar/etc. |
| PAPER_486 | UQFFCassini Complex Ring Buoyancy | THz hole, q-scope, U_Bi/U_Ii/U_Mi for Saturn rings (complex<double>) |
| PAPER_487 | UQFFMultiAstro 11-System Triad | Simultaneous Compressed/Resonance/Buoyancy + DPM creation |
| PAPER_488 | UQFFEightAstro Star-Forming Proof | 7-step equation proofs for AFGL5180/NGC346/LMC×5/NGC2174 |
| PAPER_489 | UQFF 19-System 26D Polynomial | Full 26D gravity sum: g=Σ E_DPM,i/r_i² framework breakthrough |
| PAPER_490 | Wolfram Field Unity + Sacred Time | Hypergraph rules, PI decoder (312), Mayan/Biblical/Schumann constants |

---

## Phase 6: Integration Tracker Update

Update `INTEGRATION_TRACKER.csv` with 7 new entries:

| SOURCE_FILE | STATUS | UNIQUE_PHYSICS | NOTES |
|-------------|--------|----------------|-------|
| UQFFCalculationsModule | INTEGRATED | 5 | U_g1,U_g3,U_m,E,eta — 5 systems |
| UQFFBuoyancySNRModule | INTEGRATED | 7 | F_U_Bi_i master, 6 components — 5 SNR systems |
| UQFFCassiniBuoyancyModule | INTEGRATED | 5 | Complex U_Bi/Ii/Mi/THz/qscope — Saturn |
| UQFFMultiAstroSystemsModule | INTEGRATED | 4 | 3-mode simultaneous + DPM — 11 systems |
| UQFFEightAstroSystemsModule | INTEGRATED | 3 | Star-forming proof system — 8 systems |
| UQFFNineteenAstroSystemsModule | INTEGRATED | 8 | 26D polynomial framework — 19 systems |
| WolframFieldUnityModule | INTEGRATED | 7 | Hypergraph+PI decoder+sacred time — 7 rules |

**Total new physics terms: ~39 (across 7 modules)**
**New systems tracked: 52**
**Build status:** Files created, CMake integration pending (Phase 1)

---

## Priority Order

1. ✅ **COMPLETE** — .h and .cpp files created (14 files)
2. ✅ **COMPLETE** — Helper docs created
3. ✅ **COMPLETE** — Whitepapers PAPER_484–490
4. 🔵 **NEXT SESSION** — MAIN_1_CoAnQi.cpp: add #includes + menu option 16
5. 🔵 **NEXT SESSION** — CMakeLists.txt: add 7 new source files
6. 🔵 **NEXT SESSION** — CondensedPhysics2.py v4.3.9: add 7 new calculator classes
7. 🔵 **NEXT SESSION** — source2.cpp: new Tab 1 query type entries + Tab 7 result rows
8. 🔵 **NEXT SESSION** — INTEGRATION_TRACKER.csv update
