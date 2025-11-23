@echo off
REM UPX Compression Script for MAIN_1_CoAnQi.exe
REM Created: November 22, 2025
echo.
echo === UPX COMPRESSION SCRIPT ===
echo Target: build_msvc\Release\MAIN_1_CoAnQi.exe
echo.

cd /d "%~dp0"

if not exist "build_msvc\Release\MAIN_1_CoAnQi.exe" (
    echo ERROR: MAIN_1_CoAnQi.exe not found!
    echo Build the executable first: cmake --build build_msvc --config Release
    pause
    exit /b 1
)

echo BEFORE COMPRESSION:
dir "build_msvc\Release\MAIN_1_CoAnQi.exe" | findstr "MAIN_1"
echo.

echo Running UPX compression (--best --lzma)...
C:\Tools\upx\upx.exe --best --lzma --verbose "build_msvc\Release\MAIN_1_CoAnQi.exe"

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo WARNING: UPX compression failed! Trying fallback (--best only)...
    C:\Tools\upx\upx.exe --best "build_msvc\Release\MAIN_1_CoAnQi.exe"
)

echo.
echo AFTER COMPRESSION:
dir "build_msvc\Release\MAIN_1_CoAnQi.exe" | findstr "MAIN_1"
echo.
echo === COMPRESSION COMPLETE ===
pause
