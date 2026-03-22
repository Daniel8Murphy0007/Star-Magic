// UFEOrbModule.h
// Unified Field Equation (UFE) for Red Dwarf Reactor Plasma Orb Experiment.
// 496 frames at 33.3 fps; 6 batch types; 26 quantum levels.
// Key physics: t^- = -t_n * exp(pi - t_n) time transformation.
// Computes UP(t) and FU for plasmoid dynamics with SCm/UA vacuum energies.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef UFE_ORB_MODULE_H
#define UFE_ORB_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

enum class BatchType {
    BATCH_31,       // Frame 301, t=9.03 s
    BATCH_39,       // Frame 451, t=13.53 s
    EARLY_SEQUENCE, // Photo #9, t=0.24 s
    MID_SEQUENCE,   // Batch 30 end, t=8.73 s
    LATE_SEQUENCE,  // Batch 39/6, t=13.68 s
    GENERIC
};

class UFEOrbModule {
private:
    std::map<std::string, double> variables;
    BatchType current_batch;

    double computeTminus(double t_n);
    double computeUgSum(double t, double r);
    double computeUmSum(double t, double r);
    double computeMetricTerm();
    double computeUbTerm(double t_minus);
    double computeFUExtension();
    double computeVacEnergy(const std::string& type);
    double computePlasmoidCount(double timestamp);

public:
    UFEOrbModule(BatchType batch = BatchType::GENERIC);

    void setBatch(BatchType batch);

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeUP(double t);  // Full UP(t) = Ug_sum + Um_sum + metric + Ub + vac*spin
    double computeFU(double t);  // FU = UP + extension

    std::string getEquationText();
    std::string getSolutions(double t);

    void printVariables();
};

#endif // UFE_ORB_MODULE_H
