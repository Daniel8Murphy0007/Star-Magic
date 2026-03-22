// SMBHBinaryUQFFModule.h
// Modular C++ implementation of UQFF frequency-domain gravity for SMBH Binary (LISA target).
// Novel physics: frequency-domain MUGE a = f_total * lambda_Planck / (2*pi).
// Params: M1=4e6 Msun, M2=2e6 Msun, t_coal=1.555e7 s, r_init=0.1 ly, z=0.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SMBHBINARY_UQFF_MODULE_H
#define SMBHBINARY_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class SMBHBinaryUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeRsep(double t);      // orbital separation as function of time
    double computeFsuper();
    double computeFfluid();
    double computeFquantum();
    double computeFAether();
    double computeFUg1(double r);
    double computeFUg2();
    double computeFUg3(double r);
    double computeFUg4(double t);
    double computeFtotal(double t, double r);

public:
    SMBHBinaryUQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    // Returns acceleration a = f_total * lambda_Planck / (2*pi)
    double computeG(double t, double r);
    double computeCoalescenceTime();
    std::string getEquationText();
    void printVariables();
};

#endif // SMBHBINARY_UQFF_MODULE_H
