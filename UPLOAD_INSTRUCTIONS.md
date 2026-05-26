# Star-Magic Upload & Setup Instructions (Windows / PowerShell)

**Repository**: https://github.com/Daniel8Murphy0007/Star-Magic  
**Primary Application**: `MAIN_1_CoAnQi.exe` (C++ console with full UQFF + Grok API + Verification Framework)  
**Date**: May 2026 (post Phase 3 verification)

---

## 1. Clone the Repository

```powershell
git clone https://github.com/Daniel8Murphy0007/Star-Magic.git
cd Star-Magic
```

---

## 2. Python Dependencies (`python -m pip install`)

**Core project requirements** (CondensedPhysics pipeline, Grok API, QCalcGeom, etc.):

```powershell
python -m pip install -r requirements.txt
```

**Verification Framework specific** (26D geometric matching, clone snapshots, VDS_CVP_BH26, QCalcGeom.py helpers):

```powershell
python -m pip install -r clones_archive/requirements.txt
```

This installs:
- numpy, scipy, requests (mandatory for geometric solvers, 26D/VDS_CVP_BH26 analysis, GrokAPI.py calls from source178)

Optional (recommended for full QCalcGeom / plotting / symbolic work):
```powershell
python -m pip install matplotlib sympy
```

> **IMPORTANT**: Run both `pip install` commands above. The verification framework (Phase 3) and Grok integration will not function correctly without them.

---

## 3. Build the C++ Application (CMake + MSVC)

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

The resulting `MAIN_1_CoAnQi.exe` (and supporting DLLs) will be in `build/Release/` (or your configured output dir).

**CMake will automatically copy** the required Python helpers next to the .exe:
- GrokAPI.py
- APIKeyManager.py
- (and any other .py files declared in CMakeLists.txt POST_BUILD rules)

---

## 4. Python File Co-Location / PATH Requirements

Several features rely on launching Python subprocesses from the C++ console:

- Grok API (source178_grok_api.cpp → `python GrokAPI.py` / `APIKeyManager.py`)
- Unified Geometric Verification Framework (VerificationOrchestrator → `python clones_archive/clone_snapshot_writer.py`, QCalcGeom.py, Wolfram bridges)
- 26D geometric checkpoint generation

**Rule**: The following .py files **must** be either:
- Co-located in the same directory as `MAIN_1_CoAnQi.exe`, **or**
- On your system `PATH` when the .exe runs.

Files that must resolve:
- `GrokAPI.py`
- `APIKeyManager.py`
- `QCalcGeom.py`
- `clones_archive/clone_snapshot_writer.py`
- Any other verification / Grok Thread import scripts you invoke from the menu

If Python cannot find them, the "Configure Grok API Key", "Test Grok AI Integration", and "Run 26D Geometric Verification Batch (Grok Threads)" menu items will report errors.

**Quick test from the built .exe directory**:
```powershell
python --version
python -c "import requests, numpy, scipy; print('OK')"
python APIKeyManager.py --status
```

---

## 5. Grok API Key Activation (xAI)

After building and running `MAIN_1_CoAnQi.exe`:

1. From the interactive menu, choose **"Configure Grok API Key"**.
2. Paste your key from https://console.grok.x.ai (or xAI platform).
3. The key is persisted in `grok_api_config.json` (next to the .exe) via `APIKeyManager.py`.
4. Use **"Test Grok AI Integration"** to verify a live call.

Fallback: You can also run directly:
```powershell
python APIKeyManager.py set YOUR_KEY_HERE
python APIKeyManager.py --status
```

See `GROK_ACTIVATION_GUIDE.md` for full troubleshooting (PATH issues, config file location, Windows launcher variants `py`/`python`/`py -3`).

---

## 6. Running the Unified Geometric Verification Framework (Phase 3)

After the pip installs and build:

- In `MAIN_1_CoAnQi.exe`, select the menu option for **"Run 26D Geometric Verification Batch (Grok Threads)"** (or invoke directly via the orchestrator once wired).
- Or run the snapshot writer POC:
  ```powershell
  python clones_archive/clone_snapshot_writer.py --example
  ```
  This produces timestamped archives under `clones_archive/YYYYMMDD_HHMMSS/` containing full 26D DPMVars, QCalcGeom checkpoints, Wolfram events, VDS_CVP_BH26 branches, and the 9 mandatory test results (Primordial, Gold Standard, First Primitives, Derivations, Variational Stationarity, First Axiom, Audit Outputs, Closure Whitepapers, Closure Test).

See `VERIFICATION_CONTRACT.md` (v0.2) for the authoritative spec, including:
- Grok Threads reframed as the primary 1000+ clone history source
- Chandra https://chandra.harvard.edu/photo/category/misc.html as external data
- QCalcGeom.py + Wolfram as the geometric checkpoint engines
- VDS_CVP_BH26 26D variational models

---

## 7. IMPORTANT: Preservation of All TUI Session Work

**All additions and changes made to all files since the start of this TUI thread MUST be preserved across any future commits, pushes, resets, or uploads.**

This includes (but is not limited to):

- **Grok API Activation** (b5b735bf): `source178_grok_api.cpp` (printGrokAPIStatus, configureGrokAPIKey, rewritten callGrokAPI), `MAIN_1_CoAnQi.cpp` (menu cleanups + status call), `CMakeLists.txt` (copy rules), `GROK_ACTIVATION_GUIDE.md` (full rewrite removing Qt Network myth).
- **Unified Geometric Verification Framework (Phase 3)** (cfc4d4a3 + follow-ups): Full `VERIFICATION_CONTRACT.md` update resolving all 5 Phase 2 open questions with exact language (Chandra misc, QCalcGeom+Wolfram, VDS_CVP_BH26, the complete 9-item named test suite, "The 1000+ clones *are* the Grok Threads... primary authoritative historical source"), new `VerificationOrchestrator.h` + `VerificationOrchestrator.cpp` (DPMVars26D with VDS, GeometricCheckpoint, CloneSnapshot with is_grok_thread + 9-test results map, MasterProof, runVerificationBatch, verifyGeometricMatch, persistArchive, requestGrokExplanation hook), `clones_archive/clone_snapshot_writer.py` (full featured, --example producing NGC2264 26D snapshot with all sections + provenance), `clones_archive/requirements.txt`, menu integration comment/hook in `MAIN_1_CoAnQi.cpp` next to SelfModifier clone case, all supporting 26D / Grok Thread artifacts.
- **Python pip install documentation** (this file + the inline header added to `clone_snapshot_writer.py` + `clones_archive/requirements.txt`).
- Any other files touched, created, or modified during the entire TUI session (including but not limited to PLAN.md updates, source172 26D work, observational_systems_config.h references, etc.).

**Never perform a hard reset, clean, or selective staging that would discard the verification framework, Grok activation wiring, or these instructions.** The entire post-activation verification implementation and documentation are part of the canonical state of the repository.

---

## 8. Quick Sanity Checklist After Clone + Install

- [ ] `python -m pip install -r requirements.txt` succeeded
- [ ] `python -m pip install -r clones_archive/requirements.txt` succeeded
- [ ] `python -c "import numpy, scipy, requests"` works
- [ ] CMake build completed without errors
- [ ] `MAIN_1_CoAnQi.exe` launches and shows Grok status at startup
- [ ] "Configure Grok API Key" + "Test Grok AI Integration" work
- [ ] `python clones_archive/clone_snapshot_writer.py --example` produces a valid timestamped 26D JSON archive with all 9 test results and provenance fields
- [ ] All files listed in Section 7 are present and unmodified from their TUI-session state

---

If you encounter issues with Python subprocess launches from the .exe, double-check co-location / PATH and re-run the two pip install commands above.

This document fulfills the explicit request: "git commit, push origin master. Then add python -m pip install. Keep all additions/changes made to all files since the start of this TUI thread."
