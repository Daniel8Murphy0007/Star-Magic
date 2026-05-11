// Standalone test harness for QCalcGeom::runQCalcGeomTests()
// Session 254/255 — verifies T71-T73 (G1 V(UA) closure tests).
// Build: g++ -std=c++17 -I. _qcalcgeom_test_runner.cpp QCalcGeom.cpp -o _qcalcgeom_tests.exe
// Run:   ./_qcalcgeom_tests.exe
#include "QCalcGeom.h"
int main() {
    QCALCGEOM::runQCalcGeomTests();
    return 0;
}
