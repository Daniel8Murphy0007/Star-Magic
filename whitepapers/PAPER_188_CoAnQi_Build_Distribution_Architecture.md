# PAPER_188: CoAnQi Build Distribution Architecture — NSIS and Debian Packaging

**Version:** 1.0  
**Date:** March 13, 2026  
**Session:** 49 — §2.5 Grok Thread 381a8fe7 Extended Audit  
**Author:** Star-Magic UQFF Research Framework  
**Source:** grok_share_381a8f.txt lines 5500–6000

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$

## Abstract

This paper documents the cross-platform build distribution architecture for the CoAnQi scientific computing system, comprising two independent packaging systems: NSIS (Nullsoft Scriptable Install System) for Windows (.exe installer) and dpkg-deb for Debian/Ubuntu Linux (.deb package). The packaging scripts handle binary copy, desktop shortcut creation, Windows registry entries (for NSIS), and standard Debian DEBIAN/control metadata. This architecture enables one-click installation on both primary target platforms while preserving the CoAnQi computational environment and library dependencies.

---

## 1. Architecture Overview

```
CoAnQi Build Distribution
├── Windows (NSIS)
│   ├── installer.nsi          ← NSIS script
│   ├── CoAnQi.exe             ← compiled binary (MSVC Release)
│   ├── Qt5Core.dll + Qt5Gui.dll + Qt5Widgets.dll
│   ├── MSVCP140.dll + VCRUNTIME140.dll
│   └── Output: CoAnQi_Setup.exe
└── Linux (Debian)
    ├── deb_package.sh         ← packaging script
    ├── deb_build/
    │   ├── DEBIAN/control     ← package metadata
    │   ├── usr/bin/coanqi     ← installed binary
    │   └── usr/share/applications/coanqi.desktop
    └── Output: coanqi_1.0_amd64.deb
```

---

## 2. NSIS Windows Installer Script

### 2.1 installer.nsi Structure

```nsis
!include "MUI2.nsh"
Name "CoAnQi UQFF Scientific Calculator"
OutFile "CoAnQi_Setup.exe"
InstallDir "$PROGRAMFILES64\CoAnQi"
InstallDirRegKey HKLM "Software\CoAnQi" "Install_Dir"

; Modern UI pages
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

!insertmacro MUI_LANGUAGE "English"

Section "CoAnQi Core" SecCore
    SetOutPath "$INSTDIR"
    
    ; Main executable
    File "CoAnQi.exe"
    
    ; Qt5 runtime DLLs
    File "Qt5Core.dll"
    File "Qt5Gui.dll"
    File "Qt5Widgets.dll"
    File "Qt5Network.dll"
    File "Qt5WebEngineWidgets.dll"
    
    ; MSVC runtime
    File "MSVCP140.dll"
    File "VCRUNTIME140.dll"
    
    ; Physics data
    File /r "data\"
    
    ; Write registry entries
    WriteRegStr HKLM "Software\CoAnQi" "Install_Dir" "$INSTDIR"
    WriteRegStr HKLM "Software\CoAnQi" "Version" "1.0.0"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoAnQi" \
        "DisplayName" "CoAnQi UQFF Scientific Calculator"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoAnQi" \
        "UninstallString" '"$INSTDIR\uninstall.exe"'
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoAnQi" \
        "DisplayVersion" "1.0.0"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoAnQi" \
        "Publisher" "Star-Magic UQFF Research Framework"
    
    ; Create uninstaller
    WriteUninstaller "$INSTDIR\uninstall.exe"
SectionEnd

Section "Desktop Shortcut" SecShortcut
    CreateShortcut "$DESKTOP\CoAnQi.lnk" "$INSTDIR\CoAnQi.exe" \
        "" "$INSTDIR\CoAnQi.exe" 0
SectionEnd

Section "Start Menu" SecStartMenu
    CreateDirectory "$SMPROGRAMS\CoAnQi"
    CreateShortcut "$SMPROGRAMS\CoAnQi\CoAnQi.lnk" "$INSTDIR\CoAnQi.exe"
    CreateShortcut "$SMPROGRAMS\CoAnQi\Uninstall.lnk" "$INSTDIR\uninstall.exe"
SectionEnd

Section "Uninstall"
    Delete "$INSTDIR\CoAnQi.exe"
    Delete "$INSTDIR\uninstall.exe"
    Delete "$DESKTOP\CoAnQi.lnk"
    RMDir /r "$SMPROGRAMS\CoAnQi"
    RMDir /r "$INSTDIR"
    DeleteRegKey HKLM "Software\CoAnQi"
    DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\CoAnQi"
SectionEnd
```

### 2.2 NSIS Compilation

```batch
makensis installer.nsi
:: Output: CoAnQi_Setup.exe (self-extracting installer)
```

---

## 3. Debian Packaging Script

### 3.1 deb_package.sh

```bash
#!/bin/bash
set -euo pipefail

PKG_NAME="coanqi"
PKG_VERSION="1.0"
PKG_ARCH="amd64"
PKG_DIR="deb_build"

# Create directory structure
mkdir -p "${PKG_DIR}/DEBIAN"
mkdir -p "${PKG_DIR}/usr/bin"
mkdir -p "${PKG_DIR}/usr/share/applications"
mkdir -p "${PKG_DIR}/usr/share/coanqi/data"

# Copy binary
cp build/CoAnQi "${PKG_DIR}/usr/bin/coanqi"
chmod 755 "${PKG_DIR}/usr/bin/coanqi"

# Copy data files
cp -r data/ "${PKG_DIR}/usr/share/coanqi/"

# Create DEBIAN/control
cat > "${PKG_DIR}/DEBIAN/control" << EOF
Package: ${PKG_NAME}
Version: ${PKG_VERSION}
Architecture: ${PKG_ARCH}
Maintainer: Star-Magic UQFF Research <research@starmagic.uqff>
Depends: libqt5core5a (>= 5.12), libqt5gui5 (>= 5.12), libqt5widgets5 (>= 5.12),
 libstdc++6 (>= 9), libc6 (>= 2.27)
Description: CoAnQi UQFF Scientific Calculator
 A multi-modal scientific computing system implementing the Unified Quantum Field
 Framework (UQFF) for astrophysical simulation, symbolic mathematics, and
 collaborative real-time physics calculation. Includes 5 parallel calculators,
 Qt5 GUI with 21 tabs, and REST API server at port 3141.
EOF

# Create .desktop file
cat > "${PKG_DIR}/usr/share/applications/coanqi.desktop" << EOF
[Desktop Entry]
Name=CoAnQi UQFF
Comment=Unified Quantum Field Framework Calculator
Exec=/usr/bin/coanqi
Terminal=false
Type=Application
Categories=Science;Physics;Education;
EOF

# Set permissions
find "${PKG_DIR}" -type d -exec chmod 755 {} \;
find "${PKG_DIR}" -type f -exec chmod 644 {} \;
chmod 755 "${PKG_DIR}/usr/bin/coanqi"

# Build the package
dpkg-deb --build "${PKG_DIR}" "${PKG_NAME}_${PKG_VERSION}_${PKG_ARCH}.deb"

echo "Package built: ${PKG_NAME}_${PKG_VERSION}_${PKG_ARCH}.deb"
echo "Install with: sudo dpkg -i ${PKG_NAME}_${PKG_VERSION}_${PKG_ARCH}.deb"
```

---

## 4. CMake Integration

The packaging is integrated with CMake via CPack:

```cmake
# In CMakeLists.txt
include(CPack)

set(CPACK_GENERATOR "NSIS;DEB")

# NSIS settings
set(CPACK_NSIS_DISPLAY_NAME "CoAnQi UQFF Scientific Calculator")
set(CPACK_NSIS_PACKAGE_NAME "CoAnQi")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
set(CPACK_NSIS_CREATE_ICONS_EXTRA
    "CreateShortcut '$DESKTOP\\\\CoAnQi.lnk' '$INSTDIR\\\\CoAnQi.exe'")

# DEB settings
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Star-Magic UQFF Research")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "libqt5core5a (>= 5.12), libqt5gui5 (>= 5.12), libqt5widgets5 (>= 5.12)")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

# Run: cmake --build . --target package
```

---

## 5. Dependency Matrix

| Dependency | Windows | Linux | Version |
|-----------|---------|-------|---------|
| Qt5Core | `Qt5Core.dll` | `libqt5core5a` | ≥ 5.12 |
| Qt5Gui | `Qt5Gui.dll` | `libqt5gui5` | ≥ 5.12 |
| Qt5Widgets | `Qt5Widgets.dll` | `libqt5widgets5` | ≥ 5.12 |
| Qt5Network | `Qt5Network.dll` | `libqt5network5` | ≥ 5.12 |
| MSVC runtime | `MSVCP140.dll` | N/A | 14.0+ |
| GCC/libstdc++ | N/A | `libstdc++6` | ≥ 9 |
| Python 3 | embedded | `python3` | ≥ 3.8 |
| Node.js | optional | `nodejs` | ≥ 16 |

---

## 6. Installation Sizes

| Component | Size | Platform |
|-----------|------|---------|
| `CoAnQi.exe` | 1.43 MB (UPX compressed at 15.51%) | Windows |
| Qt5 DLLs | ~25 MB | Windows |
| Total installer | ~30 MB | Windows |
| `.deb` package | ~5 MB | Linux |
| Installed size | ~35 MB | Both |

---

## 7. Conclusion

The CoAnQi dual-platform packaging architecture enables deployment on Windows (via NSIS self-extracting installer with registry integration and desktop shortcuts) and Linux (via dpkg-deb with standard Debian control metadata). The CMake/CPack integration allows both packages to be generated from a single build configuration. This is the production distribution pathway for the CoAnQi UQFF scientific computing system.

---

## References

- Source: grok_share_381a8f.txt lines 5500–6000
- Related: BUILD_INSTRUCTIONS_PERMANENT.md, CMakeLists.txt, PAPER_193 (Modular C++ Architecture)
- CP1 Class: `CoAnQiBuildDistributionCalculator`
