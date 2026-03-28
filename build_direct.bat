@echo off
setlocal
set "MSVC=C:\Program Files\Microsoft Visual Studio\18\Professional\VC\Tools\MSVC\14.50.35717"
set "SDK=C:\Program Files (x86)\Windows Kits\10"
set "SDKVER=10.0.26100.0"
set "PATH=%MSVC%\bin\Hostx64\x64;%SystemRoot%\System32"
set "INCLUDE=%MSVC%\include;%SDK%\Include\%SDKVER%\ucrt;%SDK%\Include\%SDKVER%\shared;%SDK%\Include\%SDKVER%\um"
set "LIB=%MSVC%\lib\x64;%SDK%\Lib\%SDKVER%\ucrt\x64;%SDK%\Lib\%SDKVER%\um\x64"
cd /d c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic
cl.exe /std:c++20 /Zc:__cplusplus /O2 /EHsc /D_USE_MATH_DEFINES /DQCALCGEOM_STANDALONE /I. /Fe:qcalcgeom_phaseE.exe QCalcGeom.cpp _qcg_main.cpp
echo EXIT_CODE=%ERRORLEVEL%
