call "C:\Program Files\Microsoft Visual Studio\18\Professional\VC\Auxiliary\Build\vcvars64.bat"
cd /d "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
cl /std:c++17 /O2 /EHsc /D_USE_MATH_DEFINES /DQCALCGEOM_STANDALONE /I. /Fe:qcalcgeom_phaseE.exe QCalcGeom.cpp _qcg_main.cpp
