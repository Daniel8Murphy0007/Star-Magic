# Wolfram WSTP Phase Roadmap
**Version:** 1.0 (Session 132, March 2026)  
**Purpose:** Integration plan and workflow guide for Wolfram WSTP symbolic verification

---

## Current State Summary

| Component | Status | File | Notes |
|-----------|--------|------|-------|
| WSTP Bridge | ✅ WORKING | `source174_wolfram_bridge_embedded.cpp` | WSOpenArgv + wolfram.exe -mathlink |
| UQFF Export Prototype | ✅ WORKING | `source175_uqff_wolfram_export.cpp` | ExportFullUQFFPrototype() |
| Auto Full Export | ✅ WORKING | `source176_auto_full_uqff.cpp` | AutoExportFullUQFF() + ostringstream |
| Wolfram Field Unity | ✅ WORKING | `source177_wolfram_field_unity.cpp` | Hypergraph rules |
| Grok API | ✅ WORKING | `source178_grok_api.cpp` | xAI integration |
| Wolfram Engine 14.3 | ✅ ACTIVATED | `wolfram.exe` | Free edition, activated |
| Build Flags | ✅ APPLIED | `CMakeLists.txt` | AVX2, /GL, LTCG, /Gy, C++20 |
| Physics Terms | ✅ 5,034+ | `wolfram_physics_terms_5034.txt` | Exceeds original 3,000 goal |
| WOLFRAM_TERM macros | ✅ INJECTED | All source*.cpp files | Via `add_wolfram_terms.py` |

---

## WSTP Launch Configuration (CANONICAL)

This is the **only working launch method** for free Wolfram Engine 14.3 on Windows.  
Do NOT change this — it was validated through extensive debugging in Nov 2025.

```cpp
// In source174_wolfram_bridge_embedded.cpp InitializeWolframKernel()
const char* argv[] = {
    "CoAnQi",
    "-linkmode", "launch",
    "-linkname",
    "\"C:\\Program Files\\Wolfram Research\\Wolfram Engine\\14.3\\wolfram.exe\" -mathlink",
    nullptr
};
int argc = 5;
int err = 0;
ws_link = WSOpenArgv(ws_env, argc, const_cast<char**>(argv), &err);
```

**Key Rule:** Use `wolfram.exe`, NOT `WolframKernel.exe`.  
The free Engine requires the wrapper executable for WSTP authentication.

---

## Build and Run Procedure

### 1. Build (Release, x64)
```powershell
# Configure (Visual Studio 2022, x64)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Build with AVX2 + LTCG + C++20 (requires MSVC v14.44+)
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi
```

### 2. Run Auto-Export (Menu Option)
```powershell
# Set library path for WSTP runtime
$env:PATH = "C:\Program Files\Wolfram Research\Wolfram Engine\14.3\SystemFiles\Libraries\Windows-x86-64;" + $env:PATH

# Run calculator
.\build_msvc\Release\MAIN_1_CoAnQi.exe
# Select: Option 10 (Auto-export full UQFF to Wolfram)
# Or:     Option 9  (WSTP kernel interface)
```

### 3. Expected Output Sequence
```
[WSTP] Initializing environment...
[WSTP] Launching kernel...
[DEBUG] Launch command: "C:\Program Files\Wolfram Research\..." -mathlink
Kernel launched and connected.
=== AUTO FULL UQFF EXPORT (bulletproof) ===
Exporting all 175+ source files with WOLFRAM_TERM definitions...
Built 5034 terms (XX.X MB) — sending to kernel...
Simplifying...
=== WOLFRAM VERDICT ===
[result]
```

---

## Remaining Work (Prioritized)

### Priority 1: Full FullSimplify Run with 5,034 Terms
The 5,034-term Wolfram export has not yet produced a stable verdict (FullSimplify
ran for 2+ hours in Nov 2025 before the session ended). To get the verdict:

1. Build Release x64 (verified above)
2. Run option 10 (Auto-export)
3. Leave running — **expected 4–48 hours** on Ryzen 5 5600G / 128 GB RAM
4. Watch Task Manager: WolframKernel.exe memory should climb to 60–110 GB
5. Copy verdict line when it appears

**Tip:** Add a heartbeat print to AutoExportFullUQFF() to confirm progress:
```cpp
std::cout << "FullSimplify in progress... " << elapsed_seconds << "s elapsed\r";
std::cout.flush();
```

### Priority 2: Dataset Validation Suite
Connect UQFF predictions to observational data:

```python
# Proposed: validate_uqff_observations.py
# Compare UQFF g(r,t) values against:
datasets = [
    ("Planck 2018", "CMB power spectrum χ² residuals"),
    ("DESI 2024",   "BAO scale distance moduli"),
    ("JWST CEERS",  "High-z galaxy SED fits"),
    ("Gaia DR4",    "Parallax-corrected distance moduli"),
]
# Method: compute F_U_Bi_i for each dataset system, compare to observed
# Output: χ² per dataset, residual plot, confidence intervals
```

**Estimated new terms from validation:** +400 observational constraint terms

### Priority 3: Hypergraph Multiway Export (source177)
Export UQFF field equations as Wolfram Language hypergraph evolution rules:

```mathematica
(* Framework for source177_wolfram_field_unity.cpp equivalent in WL *)
rules = {
  {Ug1[r_], Ug2[r_], Ug3[r_,t_]} :> {FU[r,t]},
  {FU[r_,t_], Ubi[r_]}           :> {F_U_Bi_i[r,t]},
  (* ... 26D extension rules ... *)
};
WolframModelEvolutionObject[rules, init, steps]
```

**Estimated new terms:** +312 hypergraph evolution terms

### Priority 4: Extend add_wolfram_terms.py Coverage
Current injector only processes `source*.cpp`. Extend to:
- `CondensedPhysics2.py` classes → extract WOLFRAM_TERM-equivalent docstrings
- `QCalc.py` master equations → Mathematica-formatted lambda expressions
- `index.js` physics functions → `WOLFRAM_TERM = "(* JS system *)"` export

---

## Menu Reference (Current Build)

| Option | Action | WSTP? |
|--------|--------|-------|
| 9 | WSTP kernel interface (manual Mathematica input) | Yes |
| 10 | Auto-export full UQFF to Wolfram (5,034+ terms → FullSimplify) | Yes |
| 11 | Run Wolfram Field Unity Simulation (source177 rules) | Yes |
| 14 | Test Grok AI Integration (source178) | No (HTTP) |
| 15 | SOURCE4 Unified Field Validation (UQFF + MUGE Compressed + Resonance) | No |

---

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Kernel launch failed (error 11)` | Not using wolfram.exe | Change to `wolfram.exe -mathlink` in argv[4] |
| `Waiting for WSActivate` (frozen) | Listener mode triggered | Use direct launch mode (not `-linkmode listen`) |
| `Connect to link:` dialog appears | Interactive mode | Add `-nogui` flag to launch string |
| Memory stuck at ~6 GB | 32-bit build | Ensure x64 platform in VS2022 |
| Slow string build (hours) | Old `string +=` loop | Use `std::ostringstream` (already applied) |
| Kernel exits immediately | Engine not activated | Run `wolframscript -activate` once |
| `error C2065 'argc'` | Variable declared after use | Declare `int argc = 5;` before argv array |
| `WOLFRAM_TERM macro redefinition` | Duplicate in included files | Warning only — not an error, safe to ignore |

---

## File Inventory (Wolfram Integration)

```
source174_wolfram_bridge_embedded.cpp  (6,459 B)  — WSTP bridge (WSOpenArgv)
source175_uqff_wolfram_export.cpp      (2,823 B)  — ExportFullUQFFPrototype()
source176_auto_full_uqff.cpp           (9,900 B)  — AutoExportFullUQFF() + ostringstream
source177_wolfram_field_unity.cpp     (10,896 B)  — WolframFieldUnityModule (hypergraph)
source178_grok_api.cpp                 (9,105 B)  — Grok xAI API integration
wolfram_wstp_runtime.h                (22,945 B)  — WSTP runtime header
wolfram_wstp_runtime.cpp              (22,890 B)  — WSTP runtime implementation
wolfram_physics_classes.cpp             (5.1 MB)  — Full wolfram physics registry
wolfram_sources_bridge.cpp             (783  KB)  — Included in MAIN_1_CoAnQi.cpp
wolfram_physics_terms_5034.txt         (3.6  KB)  — 5034-term reference list
wolfram_physics_terms_FULL.txt          (364 KB)  — Full terms with equations
add_wolfram_terms.py                  (1,414 B)  — WOLFRAM_TERM macro injector
WOLFRAM_INTEGRATION_STATUS.md        (11,400 B)  — MSVC build guide
WOLFRAM_INTEGRATION_COMPLETE.md      (12,000 B)  — Integration milestone docs
WOLFRAM_LAST_EVALUATION_REPORT.md    (20,000 B)  — Last kernel run report
WOLFRAM_FIX_PLAN.md                   (9,500 B)  — Fix tracking (resolved)
QCalc_Wolfram_Extensions.py           (277 KB)  — Python Wolfram extensions
```

---

## Key Commits (Wolfram Integration History)

| Commit | Date | Description |
|--------|------|-------------|
| f1abfd8 | Nov 20, 2025 | Add physics term inventories and extraction tools |
| dfceccc | Nov 2025 | Document complete file integration and Wolfram WSTP status |
| aab8760 | Nov 2025 | Activate Wolfram WSTP integration (C++20, MSVC) |
| 93301ac | Nov 2025 | Workspace update: VS Code restart files, restore points |
| 146a3c4 | Nov 2025 | MSVC C++20 compatibility fixes + AVX2/LTCG/Gy in CMakeLists.txt |

---

*See `GROK_SHARE_84A767D3_ANALYSIS.md` for the full analysis of the Nov 2025 debugging session.*  
*See `PROGRESS_TO_3000.md` for current physics term count status.*
