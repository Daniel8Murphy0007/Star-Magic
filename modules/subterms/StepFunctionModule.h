// StepFunctionModule.h
// Heaviside step function θ(x): returns 0 for x < 0, 1 for x ≥ 0.
// Used in piecewise UQFF regime switching, f_TRZ boundary transitions, and Ug window functions.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L7212

#ifndef STEP_FUNCTION_MODULE_H
#define STEP_FUNCTION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StepFunctionModule {
private:
    std::map<std::string, double> variables;

public:
    StepFunctionModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeHeaviside(double x);       // θ(x): 0 (x<0), 1 (x≥0)
    double computeSigmoid(double x, double alpha); // Smooth θ: 1/(1+exp(-α x))
    double computeWindow(double x, double a, double b); // θ(x-a) - θ(x-b) boxcar
    double computeSignFunction(double x);    // sgn(x) = 2θ(x) - 1
    std::string getEquationText();
    void printVariables();
};

#endif // STEP_FUNCTION_MODULE_H
