# PAPER_505: MSVC Release-MaxCompress Build Profile

**Session:** 137 | **Source:** grok_share_84a767d3.txt (lines 1500–1900, 2700–2800)
**Date:** November 2025 — verified commit 146a3c4
**Machine:** Trident01 — AMD Ryzen 5 5600G, 128 GB RAM, Windows 11

---

## §1.1 Abstract

The Release-MaxCompress build profile is the canonical optimized production configuration for MAIN_1_CoAnQi.cpp under MSVC (Visual Studio 2022). It combines whole-program optimization, link-time code generation, size-favoring code generation, and post-build UPX compression to produce the smallest possible fully-functional x64 executable. This paper documents each flag, its effect, and the achieved result: **1.75 MB** for a 108,000-line C++20 monolith.

---

## §1.2 Complete Flag Table

| Category | Flag | Effect | Size/Speed Gain |
|----------|------|--------|----------------|
| Standards | `/std:c++20` | Full C++20 conformance | — |
| Standards | `/permissive-` | Strict standards mode | — |
| Standards | `/Zc:__cplusplus` | Correct `__cplusplus` macro value | — |
| Architecture | `/arch:AVX2` | 256-bit AVX2 SIMD (Ryzen 5 5600G) | +10–20% speed |
| Optimization | `/GL` | Whole Program Optimization | WPO enabled |
| Optimization | `/Os` | Favor small code size | −15–25% size |
| Optimization | `/Gw` | Global data optimization | −10–20% size |
| Optimization | `/GF` | String pooling | minor |
| Optimization | `/Gy` | Function-level linking (COMDAT) | −5–15% size |
| Optimization | `/Oi` | Enable intrinsic functions | +5–15% speed |
| Linker | `/LTCG` | Link-Time Code Generation | WPO finalization |
| Linker | `/OPT:REF` | Remove unreferenced code/data | −20–30% size |
| Linker | `/OPT:ICF` | Identical COMDAT folding | −5–10% size |
| Post-build | UPX `--best --ultra-brute` | Runtime decompressor | −70% size |

---

## §1.3 CMake Implementation

```cmake
if(MSVC)
    # Base C++20 setup
    add_compile_options(/permissive-)
    add_compile_options(/Zc:__cplusplus)

    # Release-MaxCompress flags (only in Release config)
    add_compile_options("$<$<CONFIG:Release>:/arch:AVX2>")
    add_compile_options("$<$<CONFIG:Release>:/GL>")
    add_compile_options("$<$<CONFIG:Release>:/Os>")
    add_compile_options("$<$<CONFIG:Release>:/Gw>")
    add_compile_options("$<$<CONFIG:Release>:/GF>")
    add_compile_options("$<$<CONFIG:Release>:/Gy>")
    add_compile_options("$<$<CONFIG:Release>:/Oi>")
    add_link_options("$<$<CONFIG:Release>:/LTCG>")
    add_link_options("$<$<CONFIG:Release>:/OPT:REF>")
    add_link_options("$<$<CONFIG:Release>:/OPT:ICF>")
endif()

# Optional post-build UPX compression
find_program(UPX_EXECUTABLE upx)
if(UPX_EXECUTABLE)
    add_custom_command(TARGET MAIN_1_CoAnQi POST_BUILD
        COMMAND ${UPX_EXECUTABLE} --best --ultra-brute
            "$<TARGET_FILE:MAIN_1_CoAnQi>"
        COMMENT "UPX compressing executable"
    )
endif()
```

---

## §1.4 Compression Chain Analysis

Starting from uncompressed debug build:

```
Debug x64 (no opt)           →  ~50 MB
+ /GL /LTCG /OPT:REF /OPT:ICF  → ~17.5 MB  (−65%)
+ /Os                           → ~13.3 MB  (−24%)
+ /Gw /GF /Gy                   → ~11.3 MB  (−15%)
+ /arch:AVX2 /Oi                → ~11.0 MB  (−2%, mainly speed)
= Pre-UPX Release binary        → ~11.0 MB
+ UPX --best --ultra-brute      →  ~1.75 MB  (−84%)
= FINAL MAIN_1_CoAnQi.exe       →   1.75 MB  ✅
```

---

## §1.5 Build Commands

```powershell
# Configure (Visual Studio 17 2022, x64)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Build Release target
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi

# Clean rebuild
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64
cmake --build build_msvc --config Release
```

---

## §1.6 Verification

After build, verify:
```powershell
(Get-Item "build_msvc\Release\MAIN_1_CoAnQi.exe").Length / 1MB  # Should be ~1.75
cmake --build build_msvc --config Release 2>&1 | Select-String "error C|Finished"
```

---

## §1.7 Equations

```
Compression ratio      C_ratio = size_final / size_debug
LTCG gain              G_LTCG  = (size_no_LTCG - size_with_LTCG) / size_no_LTCG
AVX throughput factor  F_AVX2  = throughput_AVX2 / throughput_SSE2   ≈ 2× (256-bit width)
UPX expansion ratio    R_UPX   = time_decompression / time_cold_start  (< 50 ms typical)
```

---

## §1.8 Citation

Source: grok_share_84a767d3.txt, lines 1500–1900, 2700–2800
Confirms: 1.75 MB executable at commit 146a3c4 (Nov 22, 2025)
Machine: AMD Ryzen 5 5600G, 128 GB RAM, Windows 11, MSVC v14.44.35219
Paper number: PAPER_505
