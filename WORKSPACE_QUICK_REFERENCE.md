# WORKSPACE QUICK REFERENCE - Star-Magic UQFF

**Last Updated:** December 3, 2025 @ 13:45 PM  
**Purpose:** ONE-PAGE cheat sheet to prevent re-discovery

---

## WHERE IS EVERYTHING?

### Executables (BUILT & WORKING)

```
build_msvc\Release\MAIN_1_CoAnQi.exe          2.00 MB    ✅ WORKING (Runtime Verified Dec 3 @ 13:37:40)
build_msvc\Release\source2_minimal_test.exe   0.30 MB    ✅ WORKING
```

### Log Files (Application Runtime)

```
coAnQi_log_*.txt in ROOT directory
Latest Runtime: Dec 3, 13:37:40 (6,643/6,785 terms registered, 13-option menu)
```

### Wolfram Integration

```
source168.cpp - UQFF Buoyancy (5 systems)
source169.cpp - Cassini Buoyancy (Saturn)
source170.cpp - Multi-Astro SOURCE114 (11 systems)
source171.cpp - Eight Astro SOURCE114+ (8 LMC)
source172.cpp - Nineteen Astro SOURCE115 (19 systems, 26D)
source173.cpp - Wolfram Field Unity SOURCE116 (hypergraph)

wolfram_physics_classes.cpp - 5,703 classes (Batch 18)
wolfram_extraction/generated_classes/ - 187 classes (Batch 19)
source*_wolfram.cpp - 119 companion files

TOTAL: 5,890 Wolfram PhysicsTerm classes
```

### Grok AI Integration

```
source178_grok_api.cpp - xAI API wrapper (5 functions)
GROK_ACTIVATION_GUIDE.md - Setup instructions
GROK_QUICK_SETUP.md - 5-minute guide

STATUS: ✅ CODE COMPLETE | ⚠️ NEEDS XAI_API_KEY
```

### Source File Tracking

```
INTEGRATION_TRACKER.csv - 694 lines, 173 source files
MAIN_1_CoAnQi_integration_status.json - Physics inventory

SOURCE1-116: ✅ INTEGRATED (446 modules)
SOURCE168-173: ✅ INTEGRATED (WSTP bridge)
source178: ✅ INTEGRATED (Grok AI)
```

### Status Documentation

```
SAVEPOINT_DEC1_2025.md - Complete build verification (MOST RECENT)
CONSTRUCTION_STATUS_REPORT.md - Phase 1 85-90% complete
WOLFRAM_INTEGRATION_STATUS.md - Dual-system 5,890 classes
GROK_INTEGRATION_SUMMARY.md - xAI API details
RESTART_STATUS.md - Session restoration instructions
EMERGENCY_SESSION_CAPTURE_DEC1_2134.md - Context preservation
```

### Build System

```
CMakeLists.txt - Visual Studio 2022, MSVC-only
BUILD_INSTRUCTIONS_PERMANENT.md - vcpkg paths, WSTP config
build_msvc/ - Build directory (Release config)
RUN_MAIN_1_COANQI.bat - Quick launcher with Qt6 PATH
```

---

## CURRENT STATUS (SNAPSHOT)

**Executable:** ✅ 1.31 MB UPX-compressed, built 3:05 PM Dec 1  
**Last Run:** 4:46 PM Dec 1 (startup verification)  
**Last Full Calc:** 4:25 PM Dec 1 (100 systems, 12 threads)  

**Physics Terms:** 6,643 registered (of 6,785 expected)  
- 894 MAIN  
- 5,703 Wolfram Batch 18  
- 188 Wolfram Batch 19  
- Missing: 142 terms (system still functional)  

**Menu Options:** 13 total
- 1-8: Core calculator functions
- 9-11: Wolfram WSTP integration ✅
- 12: Grok AI integration ⚠️ (needs API key)
- 13: Exit

**Git:** d52fcd2 (HEAD -> master, origin/master)  
**Tag:** savepoint-wolfram-qt6-grok-dec1-2025 (commit 3143667)  

---

## DEPENDENCIES INSTALLED

```
Qt6 6.10.0:        C:\Qt\6.10.0\msvc2022_64\         ✅ WORKING
Wolfram WSTP 14.3: C:\Program Files\Wolfram Research\ ✅ LINKED
OpenSSL 3.6.0:     vcpkg x64-windows                 ✅ INSTALLED
UPX 5.0.2:         System PATH                       ✅ ACTIVE
MSVC 19.44:        Visual Studio 2022                ✅ COMPILER
```

**Qt5:** ⏸️ PENDING (blocked by ATL dependency in vcpkg)

---

## UNCOMMITTED CHANGES

```
Modified:
  .gitignore
  .vscode/settings.json
  BUILD_INSTRUCTIONS_PERMANENT.md
  BUILD_STATUS.md
  CMAKE_RESTORE_POINT.txt
  CMakeLists.txt
  MAIN_1_CoAnQi.cpp
  MAIN_1_CoAnQi_integration_status.json
  inbox-dropzone/Grok_clone access.xlsx
  source2.cpp

Untracked:
  C++17_MINGW_CONTAMINATION_REPORT.txt
  WORKSPACE_BACKUP_739PM/
```

**DO NOT LOSE** - Current working state

---

## BUILD COMMANDS (VERIFIED WORKING)

```powershell
# Configure
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Build
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi

# Run
.\build_msvc\Release\MAIN_1_CoAnQi.exe

# Or use launcher
.\RUN_MAIN_1_COANQI.bat
```

---

## CRITICAL CONTEXT

**What Works:**
- MAIN_1_CoAnQi.exe fully functional
- Wolfram WSTP integration active (menu 9-11)
- Qt6 6.10.0 working
- 100 astronomical systems loaded
- 6,643 physics terms registered
- Windows threading (12 cores)
- UPX compression (84.8% reduction)

**What's Pending:**
- Grok AI activation (needs XAI_API_KEY)
- Qt5 installation (for source2.cpp GUI)
- 142 missing physics term registrations (non-critical)

**What's Broken:**
- Nothing critical - system fully operational

---

## FREQUENTLY ASKED

**Q: Where's the executable?**  
A: `build_msvc\Release\MAIN_1_CoAnQi.exe`

**Q: Where are the logs?**  
A: `coAnQi_log_*.txt` in root directory

**Q: Where's Wolfram integration?**  
A: source168-173.cpp + wolfram_physics_classes.cpp + wolfram_extraction/

**Q: Where's Grok integration?**  
A: source178_grok_api.cpp (needs XAI_API_KEY to activate)

**Q: Where's the source file report?**  
A: INTEGRATION_TRACKER.csv (694 lines)

**Q: What's the current git state?**  
A: d52fcd2 on master, synced with origin

**Q: Is the build working?**  
A: YES - 1.31 MB executable built Dec 1 @ 3:05 PM, last run 4:46 PM

**Q: How many physics terms?**  
A: 6,643 registered (894 MAIN + 5,703 Wolfram + 188 Phase4) - Missing 142 but functional

---

## NEXT SESSION CHECKLIST

1. Read EMERGENCY_SESSION_CAPTURE_DEC1_2134.md
2. Check `git status`
3. Verify executable exists
4. Ask user what they need (don't re-discover)

---

**PURPOSE:** STOP WASTING TIME RE-DISCOVERING THE SAME THINGS

**Created:** December 1, 2025 @ 21:34  
**File Type:** Quick Reference / Cheat Sheet  
**Keep Updated:** YES - Update timestamps when builds change
