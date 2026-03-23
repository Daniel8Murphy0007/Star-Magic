// Ug3DiskVectorModule.h
// Disk-geometry Ug3 string-rotation field: Ug3_disk = G M_disk / r² × (h/r) scalar correction.
// For flat disk (scale height h ≪ r): Ug3_disk ≈ Ug3_sphere × (h/r). Milky Way: h/r ≈ 0.07.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L8265

#ifndef UG3_DISK_VECTOR_MODULE_H
#define UG3_DISK_VECTOR_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class Ug3DiskVectorModule {
private:
    std::map<std::string, double> variables;

public:
    Ug3DiskVectorModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeDiskUg3();                 // Ug3 = G M_disk/r² × (h/r) [m/s²]
    double computeScaleHeightRatio();        // h/r disk flattening factor
    double computeDiskSurfaceDensity();      // Σ = M_disk / (π r²) [kg/m²]
    double computeToomreQ();                 // Q = c_s κ / (π G Σ) stability
    std::string getEquationText();
    void printVariables();
};

#endif // UG3_DISK_VECTOR_MODULE_H
