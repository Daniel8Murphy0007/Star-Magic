"""
gen_fubiicalc_header.py
Returns the C++ file header, includes, BuoyancyEquation struct, and
printSection() utility for F_U_Bi_i_QCalc.cpp.
"""


def get_header() -> str:
    return r"""// ======================================================================
// F_U_Bi_i_QCalc.cpp  -  Universal Buoyancy Equation Catalogue
// Author  : Daniel T. Murphy
// Project : Star-Magic / UQFF (Unified Quantum Field Framework)
// Version : v4.83   |   Date: 2026-03-21
// ----------------------------------------------------------------------
// SOURCES:
//   grok_share_c020496d9e.txt  (6,955 lines)
//   BB_C_Equations_04Sept2025.pdf
//   UQFF+Framework_Progress_Completion_Calibration_22Sept2025.pdf
//   UQFF_99_9_Complete_14Sept2025.pdf
// ----------------------------------------------------------------------
// PHYSICS:
//   F_U_Bi_i = Universal Buoyancy (outside->inside)
//              Governs ALL mass/non-mass actions at every scale.
//   F_U_Bi   = Universal Buoyancy (inside->outside) [negligibility factor]
//   Together -> suspension and perpetual motion.
//
//   MASTER INTEGRAL (grok_share_c020496d9e.txt, lines 1-60):
//   F_U_Bi_i = Int_0^{x2} [
//     -F0
//     + (me*c^2/r^2)*DPM_momentum*cos(theta)
//     + (GM/r^2)*DPM_gravity
//     + rho_vac[UA]*DPM_stability
//     + k_LENR*(omega_LENR/omega0)^2
//     + k_act*cos(omega_act*t)
//     + k_DE*L_X
//     + 2qB0V*sin(theta)*DPM_resonance*P_pol
//     + k_neutron*sigma_n
//     + k_rel*(Ecm_eff/Ecm)^2
//     + k_UV*L_UV  +  k_mm*L_mm*f_mm
//   ] dx
//
//   KEY NUMERICAL RESULTS:
//     F_U_Bi_i (LENR dominant)  ~= +2.11e208 N  (x2 ~= -1.35e172 m)
//     F_U_Bi_i (F_rel dominant) ~= -8.31e211 N
//     F_rel  =  4.31e33 N   (2024 LEP validation)
//     F_LENR =  1.56e36 N
//
// SECTIONS:
//   A  29 per-system g_UQFF equations (astrophysical systems)
//   B   6 Compressed UQFF backbone + Triadic Master equations
//   C  10 Sub-equations  (Um, [SSq], t_n, f_Ub, Ug2, vacuum series)
//   D  12 F_U_Bi_i master integral component force equations
//   E  79 unique F_UBii variant types  (BB_C_Equations catalogue)
//   F  68 unique Um   variant types
//   G  25 numerical solutions and calibration constants
//   H   7 Lambda-CDM / MOND comparison and validation notes
// ======================================================================

// Requires: C++20 or later
// Build:  cl /EHsc /std:c++20 F_U_Bi_i_QCalc.cpp /Fe:F_U_Bi_i_QCalc.exe
//         g++ -std=c++20 -O2 F_U_Bi_i_QCalc.cpp -o F_U_Bi_i_QCalc

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <numbers>   // C++20: std::numbers::pi etc.
#include <ranges>    // C++20: range views

// -------------------------------------------------------------------
// BuoyancyEquation: core data struct for every catalogued equation
// -------------------------------------------------------------------
struct BuoyancyEquation {
    int         id;               // Catalogue number
    std::string name;             // Short identifier
    std::string expression;       // Mathematical expression (ASCII notation)
    std::string context;          // Physical context / domain
    std::string validation_note;  // Key validation comment or cross-reference
    double      numerical_result; // Computed value (0.0 = not numerically resolved)
    std::string units;            // SI or natural units
    std::string system;           // Astrophysical system or domain
    std::string section;          // Section letter A-H
};

// -------------------------------------------------------------------
// printSection: display all equations in a named section
// -------------------------------------------------------------------
static void printSection(const std::string& title,
                         const std::vector<BuoyancyEquation>& eqs)
{
    const int W = 90;
    std::cout << "\n" << std::string(W, '=') << "\n";
    std::cout << "  SECTION: " << title
              << "  [" << eqs.size() << " equations]\n";
    std::cout << std::string(W, '-') << "\n";
    for (const auto& eq : eqs) {
        std::cout << "\n  [" << std::setw(5) << eq.id << "]  "
                  << eq.name << "  (" << eq.section << ")\n";
        std::cout << "    EXPR : " << eq.expression << "\n";
        std::cout << "    CTX  : " << eq.context    << "\n";
        if (!eq.validation_note.empty())
            std::cout << "    NOTE : " << eq.validation_note << "\n";
        if (eq.numerical_result != 0.0)
            std::cout << "    NUM  : "
                      << std::scientific << std::setprecision(4)
                      << eq.numerical_result
                      << "  " << eq.units << "\n";
    }
    std::cout << std::string(W, '-') << "\n";
}

// Forward declarations
std::vector<BuoyancyEquation> buildSectionA();
std::vector<BuoyancyEquation> buildSectionB();
std::vector<BuoyancyEquation> buildSectionC();
std::vector<BuoyancyEquation> buildSectionD();
std::vector<BuoyancyEquation> buildSectionE();
std::vector<BuoyancyEquation> buildSectionF();
std::vector<BuoyancyEquation> buildSectionG();
std::vector<BuoyancyEquation> buildSectionH();

"""
