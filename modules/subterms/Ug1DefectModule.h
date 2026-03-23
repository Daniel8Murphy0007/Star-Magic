// Ug1DefectModule.h
// Computes Ug1 defect correction when dipole field has structural imperfection:
// Ug1_corr = Ug1 × (1 - δ_def) where δ_def = fractional field defect ∈ [0,1].
// Magnetar crust fracture: δ_def ≈ 0.05–0.15. Feeds UQFF calibration pipeline.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L8093

#ifndef UG1_DEFECT_MODULE_H
#define UG1_DEFECT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class Ug1DefectModule {
private:
    std::map<std::string, double> variables;

public:
    Ug1DefectModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeDefectFraction();          // δ_def from Ug1 nominal vs. observed
    double computeCorrectedUg1(double Ug1);  // Ug1_corr = Ug1(1 - δ_def) [m/s²]
    double computeDefectEnergy();            // E_def = mass × g × δ_def × r [J]
    double computeRepairTimescale();         // τ_rep = r / v_Alfven [s]
    std::string getEquationText();
    void printVariables();
};

#endif // UG1_DEFECT_MODULE_H
