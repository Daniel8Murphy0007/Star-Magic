@echo off
setlocal enabledelayedexpansion
set "VS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat"
if exist "%VS_PATH%" (
    call "%VS_PATH%"
) else (
    echo VS2022 not found at expected path
    exit /b 1
)
cd /d "C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\.venv\Scripts\python.exe build_uqff_core.py
