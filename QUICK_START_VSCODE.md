# VS Code C++ Quick Reference

**Last Updated**: December 4, 2025 @ 18:58 PM (Phase 32: Qt Networking + Virgo Cluster)  
**Latest Commit**: 5a6346f - Phase 32 workspace update  
**Executable**: 1.35 MB (MAIN_1_CoAnQi.exe, runtime verified 17:40:40)

## 🚀 Quick Start Commands

### Build

- `Ctrl+Shift+B` → Select build task
  - "C++: Build CoAnQi (6,477 classes)"
  - "C++: Build with UPX compression"
  - "C++: Clean rebuild"
- **Visual Studio 2022 Build:**

  ```powershell
  cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64
  cmake --build build_msvc --config Release --target MAIN_1_CoAnQi
  ```

- **Manual Compression:** `compress_upx.bat` (UPX 5.0.2 --best --lzma)

### Debug

- `F5` → Start debugging
- `F9` → Toggle breakpoint
- `F10` → Step over
- `F11` → Step into
- `Shift+F11` → Step out
- `Shift+F5` → Stop debugging

### Code Navigation

- `Ctrl+Click` → Go to definition
- `Alt+F12` → Peek definition
- `Shift+F12` → Find all references
- `F2` → Rename symbol
- `Ctrl+Space` → Trigger IntelliSense

### Formatting

- `Ctrl+S` → Save and auto-format
- `Shift+Alt+F` → Format document
- `Ctrl+K Ctrl+F` → Format selection

## 📊 Current Framework

**6,643 Physics Terms Active** (MAIN_1_CoAnQi.cpp - Runtime Verified Dec 3 @ 13:37:40)

### Core Components

- **MAIN_1_CoAnQi.cpp**: 102,672 lines
- **Physics Terms**: 6,643 registered (6,785 expected, 97.9% rate)
- **Modules Integrated**: SOURCE1-116 (446 unique terms)
- **Executable**: 2.00 MB (built Dec 3 @ 13:21:45)
- **Build System**: MSVC 19.44.35207, C++20, CMake 3.31.0
- **Runtime Status**: ✅ OPERATIONAL (13-option menu, Wolfram + Grok AI active)
- **Build status**: 0 compiler warnings (all suppressions active)
- **9-option menu**: Calculate, Analyze, Optimize, Simulate, Self-Modify, Export, Import, Wolfram WSTP, Exit

### Module System

- Source163: Multi-System UQFF (8 terms)
- Source164: Nebula UQFF (7 terms)
- Source165: Buoyancy UQFF (9 terms)
- Source166: Quantum 26-State (12 terms)
- Source167: UQFF Core June 2025 (8 terms)
- **163 total source modules** (Source1-167 with gaps)

### Compiler & Build

- **Compiler**: MSVC 14.44.35207 (Visual Studio 2022 Professional)
- **C++ Standard**: C++20
- **Optimization**: Release-MaxCompress (/Os /GL /Gw /Gy /Zc:inline /GF /Oi)
- **Warning Level**: /W3 with 6 suppressions (/wd4100 /wd4244 /wd4267 /wd4305 /wd4996 /wd4005)
- **CMake**: Updated Nov 22, 2025 (Phase 27, UPX 5.0.2 integration)
- **Latest Executable**: MAIN_1_CoAnQi.exe (1.17 MB compressed, 11/22/2025 22:06:54 + UPX)
- **UPX Compression**: 7.95 MB → 1.17 MB (85.3% reduction, --best --lzma)

### Validation Status

- ✅ Quantum sum=351
- ✅ U_m=7.97e-22 T
- ✅ eta=9.48e12 n/s
- ✅ All 273 PhysicsTerm classes integrated

## 🔧 Configurations

All located in `.vscode/`:

- `c_cpp_properties.json` - IntelliSense
- `tasks.json` - Build tasks (Node.js auto-run removed)
- `launch.json` - Debug configs
- `settings.json` - C++ settings (auto-reset disabled)
- `extensions.json` - Recommended extensions

### Recent Changes (Nov 22, 2025 - Phase 27)

- ✅ All compiler warnings fixed (0 warnings, 6 suppressions added)
- ✅ UPX 5.0.2 integration (--best --lzma compression)
- ✅ 85.3% executable compression achieved (7.95 MB → 1.17 MB)
- ✅ Grok AI integration (Qt6::Network, source178_grok_api.cpp)
- ✅ Wolfram WSTP fully functional (6,477 physics classes accessible)
- ✅ Git commit 144f9b8 pushed to origin/master
- ✅ compress_upx.bat script created for manual compression

### Backup Files (Nov 13, 2025 @ 9:30 PM)

- `restore_point_13nov2025_930pm.cpp`
- `MAIN_1_CoAnQi_backup_13nov2025_930pm.cpp`
- `CoAnQi_enhancements.cpp`
- `13nov2025_backup_930pm/` directory
- `cmake_backup.md`
- `CoAnQi_log_13nov2025_930pm.txt`

**Status**: ✅ Complete - Ready for development | ⚠️ Git push pending (branch divergence)
