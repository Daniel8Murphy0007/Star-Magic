// OuterFieldBubbleModule.h
// Models expanding outer vacuum field bubble: r_bubble = r_0 × exp(H t) where H is the
// local Hubble constant. Describes [UA] outer boundary growth. At t=10 Gyr, r ≈ 1.28 r_0.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L5530

#ifndef OUTER_FIELD_BUBBLE_MODULE_H
#define OUTER_FIELD_BUBBLE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class OuterFieldBubbleModule {
private:
    std::map<std::string, double> variables;

public:
    OuterFieldBubbleModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeBubbleRadius(double t);    // r = r_0 exp(H t) [m]
    double computeExpansionRate(double t);   // dr/dt = H r [m/s]
    double computeUAEnergy(double t);        // E_UA decays as r expands [J]
    double computeBubbleVolume(double t);    // V = (4/3)π r³ [m³]
    std::string getEquationText();
    void printVariables();
};

#endif // OUTER_FIELD_BUBBLE_MODULE_H
