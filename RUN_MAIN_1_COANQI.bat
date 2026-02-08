@echo off
REM Quick launcher for MAIN_1_CoAnQi.exe with Qt6 PATH
REM Created: December 1, 2025
REM Adds Qt6 bin directory to PATH before launching

SET PATH=C:\Qt\6.10.0\msvc2022_64\bin;%PATH%
cd /d "%~dp0build_msvc\Release"
MAIN_1_CoAnQi.exe
pause
