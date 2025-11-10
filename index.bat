@echo off
REM Star-Magic UQFF Computational Engine v2.1 Launcher
REM Enhanced Edition with 106 Systems and Sources 4-6 Integration
echo Starting UQFF Computational Engine v2.1...
node index.js %*
if %ERRORLEVEL% EQU 0 (
    echo.
    echo Computation complete.
) else (
    echo.
    echo Error during execution (Exit code: %ERRORLEVEL%)
)
