@echo off
call "C:\Program Files\Microsoft Visual Studio\18\Professional\VC\Auxiliary\Build\vcvars64.bat" >nul
cl /nologo /std:c++17 /EHsc /O2 /MD /DNOMINMAX _qcg_main.cpp QCalcGeom.cpp /Fe:_qcg_test.exe
