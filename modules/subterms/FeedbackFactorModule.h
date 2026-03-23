// FeedbackFactorModule.h
// Computes environmental feedback factor F_env = f_AGN + f_SN + f_SF in UQFF MUGE.
// F_env modulates base gravity: g_base × (1 + F_env(t))
// Typical: F_env ≈ 0.0 for most systems; enhanced for AGN environments.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L3094

#ifndef FEEDBACK_FACTOR_MODULE_H
#define FEEDBACK_FACTOR_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class FeedbackFactorModule {
private:
    std::map<std::string, double> variables;

public:
    FeedbackFactorModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeF_AGN();       // AGN jet feedback contribution
    double computeF_SN();        // Supernova feedback contribution  
    double computeF_SF();        // Star formation feedback
    double computeF_env();       // Total: f_AGN + f_SN + f_SF
    std::string getEquationText();
    void printVariables();
};

#endif // FEEDBACK_FACTOR_MODULE_H
