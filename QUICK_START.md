# QUICK START - Star-Magic UQFF Calculator

**Last Updated: December 4, 2025 @ 18:58 PM**  
**Current Commit:** `5a6346f` (Phase 32: Qt Networking + Virgo Cluster Integration - Dec 4 @ 18:58)  
**Runtime Verified:** MAIN_1_CoAnQi.exe @ 17:40:40 (6,643+ terms registered)

## ✅ Current State: FULLY OPERATIONAL

- **Wolfram WSTP:** ✅ ENABLED (menu options 9-11)
- **Qt6 6.10.0:** ✅ ENABLED (QCoreApplication + QNetworkAccessManager)
- **Grok AI:** ✅ ACTIVE (xAI grok-2-1212)
- **Tracing:** ✅ ACTIVE (uqff_trace.log)
- **Build:** ✅ SUCCESS (1.35 MB executable, built Dec 4 @ 18:58)
- **MSVC:** 19.44.35207, C++20
- **Virgo Cluster:** ✅ INTEGRATED (source82_wolfram.cpp)

---

## 🚀 Run the Calculator (Right Now)

```powershell
.\build_msvc\Release\MAIN_1_CoAnQi.exe
```

**Menu Options Available (12 total):**
- Options 1-8: Core UQFF calculations (single, parallel, clone, custom, dynamic, simulations, stats, optimize)
- Options 9-11: Wolfram WSTP integration (Export UQFF, Auto-Export, Evaluate)
- Option 12: Exit

---

## 🔧 Rebuild from Scratch

```powershell
# Clean
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue

# Configure (with Wolfram)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 -DUSE_EMBEDDED_WOLFRAM=ON

# Build
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi -j 8

# Run
.\build_msvc\Release\MAIN_1_CoAnQi.exe
```

**Build Time:** ~30-60 seconds (depending on system)

---

## 🔄 Restore This Save Point (if needed)

```powershell
# Option 1: Using tag (recommended)
git checkout savepoint-wolfram-qt6-grok-dec1-2025

# Option 2: Using commit hash
git checkout 33fbd3c1240a55f2fa2baf943e1461c8f360c7ff

# Option 3: Create new branch from save point
git checkout -b working-dec1-2025 savepoint-wolfram-qt6-grok-dec1-2025

# Then rebuild
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 -DUSE_EMBEDDED_WOLFRAM=ON
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi -j 8
```

---

## 🧪 Test Wolfram Integration

```powershell
.\build_msvc\Release\MAIN_1_CoAnQi.exe
```

Then select:
- **Option 9:** WSTP kernel interface
- **Option 10:** Auto-export full UQFF to Wolfram
- **Option 11:** Run Wolfram Field Unity Simulation

Wolfram Engine will launch automatically via WSTP.

---

## 🤖 Test Grok AI Integration

```powershell
# Verify API key is set
$env:XAI_API_KEY

# If not set, configure your xAI API key (obtain from https://console.x.ai)
[System.Environment]::SetEnvironmentVariable('XAI_API_KEY', 'your-xai-api-key-here', 'User')

# Run calculator
.\build_msvc\Release\MAIN_1_CoAnQi.exe

# Select Option 12: Test Grok AI Integration
```

Expected output: Grok response confirming C++ physics simulation assistance.

---

## 📊 Verify Build Configuration

```powershell
# Check Wolfram is enabled
Get-Content build_msvc\CMakeCache.txt | Select-String "USE_EMBEDDED_WOLFRAM"
# Should show: USE_EMBEDDED_WOLFRAM:BOOL=ON

# Check Qt6 is found
Get-Content build_msvc\CMakeCache.txt | Select-String "Qt6_DIR"
# Should show: C:/Qt/6.10.0/msvc2022_64/lib/cmake/Qt6

# Check executable exists
Test-Path build_msvc\Release\MAIN_1_CoAnQi.exe
# Should show: True

# Check executable size
(Get-Item build_msvc\Release\MAIN_1_CoAnQi.exe).Length / 1MB
# Should show: ~1.31
```

---

## 📝 View Tracing Output

```powershell
# Run calculator (tracing starts automatically)
.\build_msvc\Release\MAIN_1_CoAnQi.exe

# After running, view trace log
Get-Content uqff_trace.log | Select-Object -Last 50
```

Trace file includes:
- System initialization events
- Physics term registration (6,809 terms)
- Calculation spans with duration
- Performance metrics

---

## 🆘 Troubleshooting

### Menu only shows options 1-9 (Wolfram missing)

```powershell
# Wolfram is disabled, rebuild with:
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 -DUSE_EMBEDDED_WOLFRAM=ON
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi
```

### Grok API not responding

```powershell
# Check API key
$env:XAI_API_KEY

# Set if missing
[System.Environment]::SetEnvironmentVariable('XAI_API_KEY', 'your-key-here', 'User')

# Restart PowerShell to load new environment variable
```

### Build fails with "FATAL" errors

```powershell
# This is fixed in save point - TraceLevel enum uses TRACE_ prefix
# If you see this error, you're not on the save point commit
git checkout savepoint-wolfram-qt6-grok-dec1-2025
```

### Build fails with duplicate symbol errors (source178)

```powershell
# This is fixed in save point - CMakeLists.txt doesn't compile source178 separately
# If you see this error, you're not on the save point commit
git checkout savepoint-wolfram-qt6-grok-dec1-2025
```

---

## 📚 Documentation Files

- **SAVEPOINT_DEC1_2025.md** - Complete save point documentation
- **BUILD_INSTRUCTIONS_PERMANENT.md** - Detailed build guide (updated Dec 1)
- **TRACING_GUIDE.md** - Tracing system usage
- **README.md** - Project overview
- **ENHANCEMENT_GUIDE.md** - Self-expanding framework guide

---

## 🎯 What This Build Includes

### Physics Components
- **6,809 Physics Terms** (894 MAIN + 5,703 Wolfram + 212 new)
- **145+ Astronomical Systems** (ESO137, NGC1365, Vela, Carina, etc.)
- **26D Compressed Gravity** (26-layer quantum states)
- **19-System Master Equations** (SOURCE115/source172)
- **Wolfram Hypergraph** (SOURCE116/source173)

### Integration Features
- **Wolfram WSTP** - Direct kernel communication
- **Qt6 Networking** - Grok API calls
- **xAI Grok-2-1212** - AI-assisted diagnostics
- **AI Toolkit Tracing** - Performance monitoring
- **Windows Threading** - Parallel calculations

### Build Optimizations
- **UPX Compression** - 84.8% size reduction (9.07 MB → 1.31 MB)
- **MSVC /Os** - Favor small code
- **Whole Program Optimization** - /GL + /LTCG
- **Global Data Optimization** - /Gw
- **Function-level Linking** - /Gy

---

**Status:** ✅ READY TO USE  
**Date:** December 1, 2025 3:05 PM  
**Commit:** `33fbd3c` (tagged as `savepoint-wolfram-qt6-grok-dec1-2025`)
