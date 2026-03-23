// ScmReactivityDecayModule.h
// Rate equation for [SCm] reactivity loss: d[SCm]/dt = -k_r [SCm]
// Solution: [SCm](t) = [SCm]_0 exp(-k_r t). k_r ≈ 1e-18 s⁻¹; [SCm]_0 = 0.57.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L6096

#ifndef SCM_REACTIVITY_DECAY_MODULE_H
#define SCM_REACTIVITY_DECAY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ScmReactivityDecayModule {
private:
    std::map<std::string, double> variables;

public:
    ScmReactivityDecayModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeSCmAtTime(double t);       // [SCm](t) = [SCm]_0 exp(-k_r t)
    double computeDecayRate();               // k_r = -1/[SCm] × d[SCm]/dt [s⁻¹]
    double computeHalfLife();               // t_½ = ln(2)/k_r [s]
    double computeCurrentReactivity();       // [SSq] = [SCm](t_now) [unitless]
    std::string getEquationText();
    void printVariables();
};

#endif // SCM_REACTIVITY_DECAY_MODULE_H
