// UgIndexModule.h
// Tracks gravitational sub-field indices i ∈ [1,4] across n ∈ [1,26] dimensional layers.
// Ug_{i,n} = G M / r² × Q_i × [UA]_i × [SCm]_i — indexer for 26D UQFF gravity array.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L4049

#ifndef UG_INDEX_MODULE_H
#define UG_INDEX_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <array>

class UgIndexModule {
private:
    std::map<std::string, double> variables;

public:
    UgIndexModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeUg(int i, int n);              // Returns Ug_{i,n} [m/s²]
    double computeUgSum(int n);                  // Sum over i=1..4 for layer n
    double computeFullSum();                     // Total Σ_n Σ_i Ug_{i,n}
    std::string getEquationText();
    void printVariables();
};

#endif // UG_INDEX_MODULE_H
